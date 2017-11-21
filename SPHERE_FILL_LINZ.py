#!/usr/bin/python
"""
Thi script is used for creating sphere clumps to fill the STL mesh in LIGGGHTS
simulations. The method used is based on the following paper:
* Amberger S, Friedl M, Goniva C, Pirker C, Approximation of Objects By Spheres
for multisphere simulations in DEM

https://www.researchgate.net/publication/259693012

@author: Dr Wei Chen

"""

import numpy as np
import trimesh
from scipy.spatial import Delaunay
import math
import time
import sys

class STL_MESH:
    def __init__(self, file_name):
        self.MESH_OBJ = trimesh.load_mesh(file_name)
        self.FACES = self.MESH_OBJ.faces
        self.VERTICES = self.MESH_OBJ.vertices
        self.BOUNDS = self.MESH_OBJ.bounds
        self.CNVX_HULL_FACES = self.MESH_OBJ.convex_hull.faces
        self.BOUND_BOX = self.MESH_OBJ.bounding_box
        
def GRID_MESH(mesh_bounds):
    '''
    This function is designed to create the uniform grid based on the input
    boundary limits of the mesh. Resolution is kept at 10 X 10 X 10
    '''
    
    return
    
def GRID_DATA(grid, mesh_obj):
    '''
    This function is designed to compute that if the grid point is within the
    mesh. If so, the distance to all faces of the mesh
    '''
    
    return

def COMPUTE_DIST(grid_point, FACES):
    '''
    This function is designed to compute the distance between the grid point to
    all surface on a mesh if it is within the mesh.
    '''
    dist = np.zeros(len(FACES))
    
    # compute vector
    np
    
    return


def SPH_FILL(mesh_obj):
    
    # construct the gird for the mesh first
    grid = GRID_MESH(mesh_obj.BOUNDS)
    # compute grid information
    grid_data = GRID_DATA(grid, mesh_obj)
    
    return


def main():
    # Print command line arguments
    for arg in sys.argv[1:]:
        print('Currently processing' + arg)
        tmp_mesh = STL_MESH(arg)
        SPH_FILL(tmp_mesh)

if __name__ == "__main__":
    main()
    
    