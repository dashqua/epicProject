import pyvista as vtki
import numpy as np
import os
# source https://pswpswpsw.github.io/posts/2018/09/blog-post-modify-vtk-openfoam/
foldername = os.path.basename(os.getcwd())
#print(foldername)
grid = vtki.UnstructuredGrid("./VTK/"+foldername+"_0.vtk")
#print(grid.points)
#print(grid.cell_arrays)
#print(grid.point_arrays)
p_cell = grid.cell_arrays['p']
p2_cell = p_cell**2
grid._add_cell_array(p2_cell, 'p2')
grid.save("./VTK/"+foldername+"_alt_0.vtk")
