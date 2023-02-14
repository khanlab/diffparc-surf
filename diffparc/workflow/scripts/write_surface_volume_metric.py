import pandas as pd
import numpy as np

df = pd.DataFrame()
row = dict()

# sets the index column as sub-{subject} or sub-{subject}_ses-{session}
row[snakemake.params.index_col_name] = snakemake.params.index_col_value
for vol_txt, col_name in zip(snakemake.input.vol_txts, snakemake.params.col_names):
    value = np.loadtxt(vol_txt)
    row[col_name] = [value]


df = pd.DataFrame.from_dict(row)

# write to output file
df.to_csv(snakemake.output.csv, index=False)
