import nibabel as nib
import numpy as np

stat = snakemake.wildcards.statsacross

# load label gii
label_array = nib.load(snakemake.input.label_gii).agg_data()

darray = np.zeros(label_array.shape)

for i, parc in enumerate(snakemake.params.parcs):
    lbl = i + 1
    in_txt = snakemake.params.sample_txt_file.format(parc=parc, **snakemake.wildcards)
    temp_arr = np.loadtxt(in_txt, delimiter=" ", skiprows=1)

    if stat == "mean":
        darray[label_array == lbl] = np.nanmean(temp_arr)
    elif stat == "max":
        darray[label_array == lbl] = np.nanmax(temp_arr)
    elif stat == "min":
        darray[label_array == lbl] = np.nanmin(temp_arr)
    elif stat == "median":
        darray[label_array == lbl] = np.nanmedian(temp_arr)

sample_darray = nib.gifti.GiftiDataArray(
    data=darray, datatype="NIFTI_TYPE_FLOAT32"  # default intent is NIFTI_INTENT_NONE
)

gifti = nib.GiftiImage()
gifti.add_gifti_data_array(sample_darray)

gifti.to_filename(snakemake.output.gii_metric)
