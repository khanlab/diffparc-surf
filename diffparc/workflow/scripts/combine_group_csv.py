import pandas as pd

from glob import iglob
from pathlib import Path

csv_fpaths = str(Path(snakemake.output.combined_csv).parent / "*")

for idx, csv_fpath in enumerate(iglob(csv_fpaths)):
    col_suffix = csv_fpath.split("/")[-1].split("_")  # Split fname by comp.
    col_suffix[-1] = col_suffix[-1].split(".")[0]  # Remove extension

    # Filter through list to keep only label, method, and metric
    col_suffix = [
        comp
        for comp in col_suffix
        if not (
            comp.startswith("sub")
            or comp.startswith("seedspervertex")
            or comp.startswith("desc")
        )
    ]
    col_suffix = "_".join(comp for comp in col_suffix)

    tmp_df = pd.read_csv(csv_fpath)
    tmp_df.rename(
        columns=lambda col: f"{col}_{col_suffix}"
        if col.lower() not in ["subj"]
        else col,
        inplace=True,
    )

    if idx == 0:
        combined_df = tmp_df.copy()
    else:
        combined_df = pd.merge(left=combined_df, right=tmp_df, how="outer", on="subj")

    combined_df.to_csv(snakemake.output.combined_csv, index=False)
