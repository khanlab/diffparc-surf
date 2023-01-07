import pandas as pd
import nibabel as nib
import numpy as np

# this loads the maxprob dseg files for each hemi, and creates a wide table
#  1 row per subject

df = pd.DataFrame()

# dict for creating the dataframe
row = dict()

# sets the index column as sub-{subject} or sub-{subject}_ses-{session}
row[snakemake.params.index_col_name] = snakemake.params.index_col_value

dseg_nib_dict = {
    "Left": nib.load(snakemake.input.maxprob_L),
    "Right": nib.load(snakemake.input.maxprob_R),
}

# loop over hemis
for hemi, dseg_nib in dseg_nib_dict.items():

    dseg_data = dseg_nib.get_fdata()
    voxel_vol = np.product(dseg_nib.header.get_zooms())

    print(f"voxel_vol: {voxel_vol}")

    # loop over parc labels
    for i, labelname in enumerate(snakemake.params.parc_list):
        labelnum = i + 1  # i is zero-based
        parc_vol = np.sum(dseg_data == labelnum) * voxel_vol
        row[f"{hemi}_{labelname}"] = [parc_vol]

df = pd.DataFrame.from_dict(row)

# write to output file
df.to_csv(snakemake.output.csv, index=False)
