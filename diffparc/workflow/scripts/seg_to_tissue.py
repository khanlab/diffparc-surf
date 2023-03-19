import nibabel as nib
import pandas as pd
import numpy as np

df = pd.read_table(snakemake.input.label_tsv)

tissue_lut = snakemake.params.tissue_lut

dseg_nib = nib.load(snakemake.input.dseg)
dseg_vol = dseg_nib.get_fdata()
tissue_vol = np.zeros(dseg_vol.shape)

#extract tissue label
for tissue_label_name, tissue_label_num in tissue_lut.items():

    labels = df.loc[df['tissue'] == tissue_label_name].label 

    for label in labels:
        tissue_vol[dseg_vol==label] = tissue_label_num

#save the tissue dseg file
out_nib = nib.Nifti1Image(tissue_vol,affine=dseg_nib.affine,header=dseg_nib.header)
out_nib.to_filename(snakemake.output.dseg)


