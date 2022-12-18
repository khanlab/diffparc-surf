import nibabel as nib
import numpy as np

conn = np.loadtxt(snakemake.input.csv, delimiter=",", skiprows=1)

#pre-normalize connectivity if enabled:
if snakemake.params.normalize_percentile:
    for j in range(conn.shape[1]):
        
        norm_by = np.percentile(conn[:,j],snakemake.params.normalize_percentile)
        seeds_per_vertex = snakemake.params.seeds_per_vertex #to rescale back after normalizing
        conn[:,j] = conn[:,j] / norm_by * seeds_per_vertex


conn_darray = nib.gifti.GiftiDataArray(
    data=conn, intent="NIFTI_INTENT_TIMESERIES", datatype="NIFTI_TYPE_INT32"
)

gifti = nib.GiftiImage()
gifti.add_gifti_data_array(conn_darray)
gifti.to_filename(snakemake.output.gii_metric)
