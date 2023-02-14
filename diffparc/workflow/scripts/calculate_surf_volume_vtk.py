import pyvista as pv
import numpy as np

vol = pv.read(snakemake.input.surf_vtk).volume
print(vol)
volarr = np.zeros(1)
volarr[0] = vol
np.savetxt(snakemake.output.txt, volarr)
