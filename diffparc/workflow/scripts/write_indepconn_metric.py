import pandas as pd
import numpy as np

conn = dict()
df = pd.DataFrame()
row = dict()

# sets the index column as sub-{subject} or sub-{subject}_ses-{session}
row[snakemake.params.index_col_name] = snakemake.params.index_col_value

conn["Left"] = np.loadtxt(snakemake.input.csv_left, delimiter=",", skiprows=1)
conn["Right"] = np.loadtxt(snakemake.input.csv_right, delimiter=",", skiprows=1)

for hemi in ["Left", "Right"]:
    for i, labelname in enumerate(snakemake.params.target_labels):
        arr = conn[hemi][:, i]
        value = arr[np.nonzero(arr)].mean()  # get mean from nonzero values
        row[f"{hemi}_{labelname}"] = [value]


df = pd.DataFrame.from_dict(row)

# write to output file
df.to_csv(snakemake.output.csv, index=False)
