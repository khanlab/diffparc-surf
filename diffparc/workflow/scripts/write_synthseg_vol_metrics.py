import pandas as pd
import nibabel as nib
import numpy as np

# this loads the synthseg dseg file and creates a wide table
#  1 row per subject

df = pd.DataFrame()

# dict for creating the dataframe
row = dict()

# sets the index column as sub-{subject} or sub-{subject}_ses-{session}
row[snakemake.params.index_col_name] = snakemake.params.index_col_value

dseg_nib = nib.load(snakemake.input.dseg)
labels_df = pd.read_table(snakemake.input.labels_tsv)

# print(labels_df)

dseg_data = dseg_nib.get_fdata()
voxel_vol = np.product(dseg_nib.header.get_zooms())

# print(f"voxel_vol: {voxel_vol}")

# loop over parc labels
for labelnum, labelname in zip(labels_df.label, labels_df.name):
    parc_vol = np.sum(dseg_data == labelnum) * voxel_vol
    #    print(f'labelnum {labelnum} is {labelname}, vol = {parc_vol}')
    row[f"{labelname}"] = [parc_vol]

df = pd.DataFrame.from_dict(row)

# write to output file
df.to_csv(snakemake.output.csv, index=False)
