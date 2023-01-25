#!/usr/bin/env python 
import matplotlib
import matplotlib.pyplot as plt
from nilearn import plotting


def plot_synthseg(vol_nii, synthseg_dseg, out_png, wildcards):
    """Plot synthseg dseg overlaid on subject T1w"""
    # Use non-gui backend
    matplotlib.use("Agg")

    # Make output directory if it doesn't exist
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)

    # Setup figure
    fig, ax = plt.subplots(1, 1, figsize=(9, 3.5))

    # Plot overlay
    plotting.plot_roi(
        roi_img=synthseg_dseg,
        bg_img=vol_nii,
        draw_cross=False,
        axes=ax,
    )

    # Finalize and save fig
    fig.suptitle(str(Path(out_png).name.strip(".png")[:-2]))
    fig.savefig(out_png, dpi=200)


if __name__ == "__main__":
    plot_synthseg:
        vol_nii=snakemake.input.vol_nii,
        synthseg_dseg=snakemake.input.synthseg_dseg,
        out_png=snakemake.output.png,
    )