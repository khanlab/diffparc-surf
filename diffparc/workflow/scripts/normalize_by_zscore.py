import nibabel as nib
import numpy as np

nib_metric = nib.load(snakemake.input.metric)
metric_vol = nib_metric.get_fdata()

mask_vol = nib.load(snakemake.input.mask).get_fdata()


mean = np.mean(metric_vol[mask_vol == 1])
std = np.std(metric_vol[mask_vol == 1])

metric_vol = (metric_vol - mean) / std

out_nib = nib.Nifti1Image(
    metric_vol, affine=nib_metric.affine, header=nib_metric.header
)
out_nib.to_filename(snakemake.output.metric)
