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
metric_nib = nib.load(snakemake.input.metric)
labels_df = pd.read_table(snakemake.input.labels_tsv)

# print(labels_df)

dseg_data = dseg_nib.get_fdata()
metric_data = metric_nib.get_fdata()

# print(f"metric shape: {metric_data.shape}")
# print(f"dseg shape: {dseg_data.shape}")

# loop over parc labels
for labelnum, labelname in zip(labels_df.label, labels_df.name):
    if np.sum(dseg_data == labelnum) == 0:
        print(f"{labelname} is missing")
    parc_val = np.mean(metric_data[dseg_data == labelnum])
    row[f"{labelname}"] = [parc_val]

df = pd.DataFrame.from_dict(row)

# write to output file
df.to_csv(snakemake.output.csv, index=False)
