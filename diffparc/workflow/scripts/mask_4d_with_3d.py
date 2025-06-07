import nibabel as nib
import numpy as np

# Load 3D mask image
mask_img = nib.load(snakemake.input.mask)
mask = mask_img.get_fdata()

# Load 4D data image
data_img = nib.load(snakemake.input.data)
data = data_img.get_fdata()

# Check that shapes are compatible
if data.shape[:3] != mask.shape:
    raise ValueError(f"Shape mismatch: mask {mask.shape} vs data {data.shape[:3]}")

# Multiply each volume with the mask
# Broadcasting the 3D mask across the 4th dimension
result = data * mask[..., np.newaxis]

# Save the result as a 4D NIfTI file
result_img = nib.Nifti1Image(result, affine=data_img.affine, header=data_img.header)
nib.save(result_img, snakemake.output.data)
