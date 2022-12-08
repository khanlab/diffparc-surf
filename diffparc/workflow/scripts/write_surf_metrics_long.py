import pandas as pd
import nibabel as nib
import numpy as np

# this loads the cifti dlabel and dscalar files and creates a long table

# axis 0 has the dense data (ie scalars or labels),
#  and axis 1 the different surfaces (left striatum, right striatum;
#   coded as CORTEX_LEFT, CORTEX_RIGHT, so they can be visualized in wb_view)

df = pd.DataFrame()

struc_to_hemi_dict = {
    "CIFTI_STRUCTURE_CORTEX_LEFT": "L",
    "CIFTI_STRUCTURE_CORTEX_RIGHT": "R",
}

reduce_func_dict = {"inout": np.mean, "surfarea": np.sum}


dlabel_nib = nib.load(snakemake.input.dlabel)
dlabel_data = dlabel_nib.get_fdata().T

label_dict = dlabel_nib.header.get_axis(0).label[0]
del label_dict[0]  # remove background label

# loop over metrics:
for dscalar_nii, metric in zip(snakemake.input.dscalars, snakemake.params.metrics):

    dscalar_nib = nib.load(dscalar_nii)
    dscalar_data = (
        dscalar_nib.get_fdata().T
    )  # we have to transpose in order for later indices to work..

    # loop over structures (ie CORTEX_LEFT, CORTEX_RIGHT)
    for (brain_struct, data_indices, brain_model) in dscalar_nib.header.get_axis(
        1
    ).iter_structures():

        hemi = struc_to_hemi_dict[brain_struct]
        scalar = dscalar_data[data_indices]
        label = dlabel_data[data_indices]

        # loop over labels (parcels)
        for labelnum, (labelname, rgba) in label_dict.items():

            if np.any(label == labelnum):
                value = reduce_func_dict[metric](
                    scalar[label == labelnum]
                )  # mean for inout, sum for surfarea
            else:
                print(
                    f"warning: {labelname} is empty, subject={snakemake.wildcards.subject}, setting to 0"
                )
                value = 0

            row = dict()

            row["subject"] = [snakemake.wildcards.subject]

            if "session" in snakemake.wildcards._names:
                row["session"] = [snakemake.wildcards.session]

            row["label"] = [labelname]
            row["metric"] = [metric]
            row["hemi"] = [hemi]
            row["value"] = [value]

            df_row = pd.DataFrame.from_dict(row)
            df = pd.concat((df, pd.DataFrame.from_dict(row)), ignore_index=True)

# write to output file
df.to_csv(snakemake.output.csv, index=False)
