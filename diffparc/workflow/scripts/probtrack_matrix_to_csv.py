import pandas as pd

df = pd.read_table(snakemake.input[0],sep='\s+',names=snakemake.params.col_headers)
df.to_csv(snakemake.output[0],index=False)

