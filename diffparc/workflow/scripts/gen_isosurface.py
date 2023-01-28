import pyvista as pv
import nibabel as nib
import numpy as np
from vtk import vtkNIFTIImageReader


class NIFTIReader(pv.BaseReader):
    _class_reader = vtkNIFTIImageReader


reader = NIFTIReader(snakemake.input.nii)
vol = reader.read()

affine = nib.load(snakemake.input.nii).affine

# keep the translation portion only, i.e. replace top-left 3x3 with identity
affine[:3, :3] = np.eye(3)


# gifti_coordsys = nib.gifti.GiftiCoordSystem(dataspace='NIFTI_XFORM_UNKNOWN',
#                            xformspace='NIFTI_XFORM_MNI_152',
#                            xform=affine)

# the contour function produces the isosurface
surface = vol.contour([snakemake.params.threshold])

surface = surface.decimate(float(snakemake.params.decimate_percent) / 100.0)

# faces from pyvista surface are formatted with number of verts each row
# reshape and remove the first col to get Nx3
faces = surface.faces
faces = faces.reshape((int(faces.shape[0] / 4), 4))[:, 1:4]

points = surface.points

# transform points with affine so they are in world coordinates
points_h = np.vstack((points.T, np.ones((1, points.shape[0]))))
xform_points_h = affine @ points_h
xform_points = xform_points_h[:3, :].T


tri_darray = nib.gifti.GiftiDataArray(
    data=faces, intent="NIFTI_INTENT_TRIANGLE", datatype="NIFTI_TYPE_INT32"
)

points_darray = nib.gifti.GiftiDataArray(
    data=xform_points, intent="NIFTI_INTENT_POINTSET", datatype="NIFTI_TYPE_FLOAT32"
)

gifti = nib.GiftiImage()
gifti.add_gifti_data_array(points_darray)
gifti.add_gifti_data_array(tri_darray)

gifti.to_filename(snakemake.output.surf_gii)
