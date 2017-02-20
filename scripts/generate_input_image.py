import numpy as np
import vtk

radius = 100.0 # mm
spacing = 2.0 # mm
extents = [120, 120, 120]

grid = vtk.vtkImageData()
grid.SetSpacing(spacing, spacing, spacing)
grid.SetDimensions(extents)

prolif_data = vtk.vtkDoubleArray()
prolif_data.SetNumberOfTuples(extents[0]*extents[1]*extents[2])
prolif_data.SetName("proliferating")

tumour_data = vtk.vtkDoubleArray()
tumour_data.SetNumberOfTuples(extents[0]*extents[1]*extents[2])
tumour_data.SetName("tumour")

centre = np.array((extents[0]*spacing/2.0, extents[1]*spacing/2.0, extents[2]*spacing/2.0))
for kdx in range(extents[2]):
    for jdx in range(extents[1]):    
        for idx in range(extents[0]):
            index =idx + extents[0]*jdx + extents[0]*extents[1]*kdx
            loc = np.array((idx*spacing, jdx*spacing, kdx*spacing))
            if np.linalg.norm(loc-centre)<=radius:
                tumour_data.SetTuple1(index, 1.0)
                prolif_data.SetTuple1(index, 1.0e9)
            else:
                tumour_data.SetTuple1(index, 0.0)
                prolif_data.SetTuple1(index, 0.0) 

grid.GetPointData().AddArray(tumour_data)   
grid.GetPointData().AddArray(prolif_data)  

writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName("/home/grogan/rad_100.vti")
writer.SetInputData(grid)
writer.Write()               