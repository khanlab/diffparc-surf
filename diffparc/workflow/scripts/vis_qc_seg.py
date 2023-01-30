#!/usr/bin/env python
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from nilearn import plotting, surface
from pathlib import Path


def get_surf_bounds(surf_mesh):
    """Get coordinate limits to maximize display of surface"""
    coords, _ = surface.load_surf_mesh(surf_mesh)

    mins = np.floor(np.min(coords, axis=0))
    maxs = np.ceil(np.max(coords, axis=0))

    return mins, maxs


def plot_surface(fig, fig_gs, surf_mesh, surf_roi):
    """Plot surfaces - lateral + medial (striatum) OR dorsal + ventral (vtasnc)"""
    # Get bounds
    mins, maxs = get_surf_bounds(surf_mesh)

    # Lateral / Dorsal
    ax_surf_left = fig.add_subplot(fig_gs[2], projection="3d")

    # Set up views
    if "hemi-L" in surf_mesh:
        hemi = "left"
        view = "lateral" if "striatum" in surf_mesh else "dorsal"
        title = "Lateral" if "striatum" in surf_mesh else "Dorsal"
    elif "hemi-R" in surf_mesh:
        hemi = "right"
        view = "medial" if "striatum" in surf_mesh else "ventral"
        title = "Medial" if "striatum" in surf_mesh else "Ventral"
    else:
        raise ValueError("Unknown hemisphere.")

    plotting.plot_surf_roi(
        surf_mesh=surf_mesh,
        roi_map=surf_roi,
        view=view,
        title=title,
        hemi=hemi,
        axes=ax_surf_left,
    )
    ax_surf_left.set_xlim((mins[0], maxs[0]))
    ax_surf_left.set_ylim((mins[1], maxs[1]))
    ax_surf_left.set_zlim((mins[2], maxs[2]))

    # Medial / Ventral
    ax_surf_right = fig.add_subplot(fig_gs[3], projection="3d")

    # Set up views
    if "hemi-R" in surf_mesh:
        hemi = "right"
        view = "lateral" if "striatum" in surf_mesh else "dorsal"
        title = "Lateral" if "striatum" in surf_mesh else "Dorsal"
    elif "hemi-L" in surf_mesh:
        hemi = "left"
        view = "medial" if "striatum" in surf_mesh else "ventral"
        title = "Medial" if "striatum" in surf_mesh else "Ventral"
    else:
        raise ValueError("Unknown hemisphere.")

    plotting.plot_surf_roi(
        surf_mesh=surf_mesh,
        roi_map=surf_roi,
        view=view,
        title=title,
        hemi=hemi,
        axes=ax_surf_right,
    )
    ax_surf_right.set_xlim((mins[0], maxs[0]))
    ax_surf_right.set_ylim((mins[1], maxs[1]))
    ax_surf_right.set_zlim((mins[2], maxs[2]))


def plot_vol(fig, fig_gs, vol_roi, vol_nii, ax_title):
    """Plot volume with dseg overlay"""
    ax_nii = fig.add_subplot(fig_gs[0:2])
    plotting.plot_roi(
        roi_img=vol_roi,
        bg_img=vol_nii,
        draw_cross=False,
        axes=ax_nii,
        title=ax_title,
    )


def plot_qc(surf_mesh, surf_roi, vol_nii, vol_roi, out_png):
    """Plot QC figure with volume and surface overlays"""
    # Use non-gui backend
    matplotlib.use("Agg")

    # Make output directory if it doesn't exist
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)

    # Setup figure and main grid
    fig = plt.figure(figsize=(18, 7))
    gs_parent = gridspec.GridSpec(1, 2, figure=fig)

    # Plot left hemi
    gs_lh = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_parent[0])
    plot_vol(fig, gs_lh, vol_roi[0], vol_nii, "Left hemi")
    plot_surface(fig, gs_lh, surf_mesh[0], surf_roi[0])

    # Plot right hemi
    gs_rh = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_parent[1])
    plot_vol(fig, gs_rh, vol_roi[1], vol_nii, "Right hemi")
    plot_surface(fig, gs_rh, surf_mesh[1], surf_roi[1])

    # Finalize and save figure
    fig.suptitle(str(Path(out_png).name.strip(".png")[:-2]))
    fig.savefig(out_png, dpi=200)


if __name__ == "__main__":
    plot_qc(
        surf_mesh=snakemake.input.surf_mesh,
        surf_roi=snakemake.input.surf_roi,
        vol_nii=snakemake.input.vol_nii,
        vol_roi=snakemake.input.vol_roi,
        out_png=snakemake.output.png,
    )
