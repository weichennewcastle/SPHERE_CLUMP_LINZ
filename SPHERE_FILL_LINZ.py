#!/usr/bin/python
"""
Thi script is used for creating sphere clumps to fill the STL mesh in LIGGGHTS
simulations. The method used is based on the following paper:
* Amberger S, Friedl M, Goniva C, Pirker C, Approximation of Objects By Spheres
for multisphere simulations in DEM

https://www.researchgate.net/publication/259693012

@author: Dr Wei Chen

"""


import trimesh
import math
import sys
import itertools

import numpy as np

from scipy.spatial import Delaunay


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
    boundary limits of the mesh. Resolution is kept at a fixed resolution
    '''
    resolution = 20
    # Getting x,y,z lower bounds
    x_min = mesh_bounds[0][0]
    y_min = mesh_bounds[0][1]
    z_min = mesh_bounds[0][2]
    # Getting x,y,z upper bounds
    x_max = mesh_bounds[1][0]
    y_max = mesh_bounds[1][1]
    z_max = mesh_bounds[1][2]
    # Construct range
    x_range = np.linspace(x_min, x_max, resolution + 2)
    y_range = np.linspace(y_min, y_max, resolution + 2)
    z_range = np.linspace(z_min, z_max, resolution + 2)
    # Construct grid point numpy array
    grid_points = np.zeros((pow(resolution, 3), 3))
    cnt = 0
    for gp in itertools.product(x_range[1:-1], y_range[1:-1], z_range[1:-1]):
        grid_points[cnt][0] = gp[0]
        grid_points[cnt][1] = gp[1]
        grid_points[cnt][2] = gp[2]
        cnt += 1
    return grid_points


def GRID_IN_HULL(grid_points, mesh_obj_vertices):
    """
    Test if points in `p` are in `hull`
    `grid_point` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    # Construct a boolean array
    grid_in_mesh = np.zeros(len(grid_points), dtype=bool)
    # Construct a mesh based on the vertices on the mesh
    if not isinstance(mesh_obj_vertices, Delaunay):
        hull = Delaunay(mesh_obj_vertices)
    # Loop through points to get the status of in/out of the mesh
    cnt = 0
    for p in grid_points:
        if hull.find_simplex(p) >= 0:
            grid_in_mesh[cnt] = True
        cnt += 1
    return grid_in_mesh


def POINT_DIST_TRI(p_coord, FACES):

    return


def COMPUTE_DIST(grid_points, grid_in_mesh, FACES):
    '''
    This function is designed to compute the distance between the grid point to
    all surface on a mesh if it is within the mesh.

    Check if the grid point in with in the mesh.
    if yes -> Calculate the distance
    if no -> skip
    '''
    # Construct a numpy array for grid point distance to triangle distance
    grid_face_dist = np.zeros((len(grid_points), len(FACES)))

    cnt = 0
    for in_out_flag, gp in zip(grid_in_mesh, grid_points):
        if in_out_flag = True:
            grid_face_dist[cnt] = POINT_DIST_TRI(gp, FACES)
        else:
            continue
    return grid_face_dist


def SPH_FILL(mesh_obj):
    '''
    This is actually the main part of the function
    '''
    # construct the gird for the mesh first
    grid_points = GRID_MESH(mesh_obj.BOUNDS)
    # compute if grid points are within the mesh, output a boolean array
    grid_in_mesh = GRID_IN_HULL(grid_points, mesh_obj.VERTICES)
    # Compute the distance if a point is within the mesh
    grid_distance = COMPUTE_DIST(grid_points, grid_in_mesh, mesh_obj.FACES)

    return


def main():
    # Print command line arguments
    for arg in sys.argv[1:]:
        print('Currently processing' + arg)
        tmp_mesh = STL_MESH(arg)
        SPH_FILL(tmp_mesh)


if __name__ == "__main__":
    main()
