import nibabel as nib
import numpy as np
from os.path import join
from glob import glob

stat = snakemake.wildcards.statsacross

sample_txts = sorted(glob(join(snakemake.input.sampledti_dir, "*.txt")))
darray = np.zeros((len(sample_txts), 1))

for i, txt in enumerate(sample_txts):
    temp_arr = np.loadtxt(txt, delimiter=" ", skiprows=1)

    if stat == "mean":
        darray[i] = temp_arr.mean()
    elif stat == "max":
        darray[i] = temp_arr.max()
    elif stat == "min":
        darray[i] = temp_arr.min()
    elif stat == "median":
        darray[i] = np.median(temp_arr)


sample_darray = nib.gifti.GiftiDataArray(
    data=darray, datatype="NIFTI_TYPE_FLOAT32"  # default intent is NIFTI_INTENT_NONE
)
gifti = nib.GiftiImage()
gifti.add_gifti_data_array(sample_darray)

gifti.to_filename(snakemake.output.gii_metric)
