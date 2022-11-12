import pandas as pd
import nibabel as nib




for dscalar_nii in snakemake.input.dscalars:
    dscalar_nib = nib.load(dscalar_nii)
    print(dir(dscalar_nib))
    print(dscalar_nib.get_fdata().shape)

    dscalar_data = dscalar_nib.get_fdata().T #we have to transpose in order for later indices to work..

    axes = dscalar_nib.header.get_axis(1) 
    for (brain_struct,data_indices,brain_model) in axes.iter_structures(): #loops over structures (ie CORTEX_LEFT, CORTEX_RIGHT)
        print(brain_struct)
        hemi_data = dscalar_data[data_indices]
    

