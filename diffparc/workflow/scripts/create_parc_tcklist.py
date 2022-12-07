import nibabel as nib
import numpy as np

# get the label data from the gifti
label_array = nib.load(snakemake.input.label_gii).agg_data()

# get the indices that are equal to label_num (offset by 1 since .tck files start at 1 instead of 0)
indices = np.ix_(label_array == snakemake.params.label_num) + np.array(1)


with open(snakemake.output.tcklist, "w") as f:

    for i in indices.flat:
        f.write(snakemake.params.tck_filename.format(index=i))
        f.write("\n")
