import pandas as pd
import pyvista as pd

df = pd.DataFrame()
row = dict()

# sets the index column as sub-{subject} or sub-{subject}_ses-{session}
row[snakemake.params.index_col_name] = snakemake.params.index_col_value

for surf, col_name in zip(snakemake.input.surfs, snakemake.params.col_names):
    value = pv.load(surf).volume
    row[col_name] = [value]


df = pd.DataFrame.from_dict(row)

# write to output file
df.to_csv(snakemake.output.csv, index=False)
