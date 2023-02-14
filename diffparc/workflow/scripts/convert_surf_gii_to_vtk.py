import nibabel as nib
import numpy as np
import pyvista as pv

# read gifti data
nib_gii = nib.load(snakemake.input.surf_gii)

print(dir(nib_gii))
points_garray = nib_gii.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0]
tri_garray = nib_gii.get_arrays_from_intent("NIFTI_INTENT_TRIANGLE")[0]


# gifti has Nx3 for triangles
tri_data = tri_garray.data
points_data = points_garray.data

print(f"shape of tri_data: {tri_data.shape}")
print(f"shape of points_data: {points_data.shape}")

print("tri_data")
print(tri_data[:5, :])
print("points_data")
print(points_data[:5, :])

# we want to have 4xN (where first row is all 3's, indicating the number of points in the triangle)
# append the 3's, then transpose:
faces = np.hstack((3 * np.ones((tri_data.shape[0], 1)), tri_data)).astype("int16")

# points data we just have to transpose:
vertices = points_data

# print(faces[:,:5])
print(faces[:5, :])
print(vertices[:5, :])

print(f"shape of faces: {faces.shape}")
print(f"shape of vertices: {vertices.shape}")

surf = pv.PolyData(vertices, faces)

surf.save(snakemake.output.surf_vtk)
