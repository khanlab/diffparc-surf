import pandas as pd

pd.concat([pd.read_csv(csv) for csv in snakemake.input.csvs]).to_csv(
    snakemake.output.csv, index=False
)
