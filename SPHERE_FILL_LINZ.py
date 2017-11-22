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
from numpy.core.umath_tests import inner1d
from scipy.spatial.distance import cdist

class STL_MESH:
    def __init__(self, file_name):
        self.MESH_OBJ = trimesh.load_mesh(file_name)
        self.FACES = self.MESH_OBJ.faces
        self.VERTICES = self.MESH_OBJ.vertices
        self.BOUNDS = self.MESH_OBJ.bounds
        self.CNVX_HULL_FACES = self.MESH_OBJ.convex_hull.faces
        self.BOUND_BOX = self.MESH_OBJ.bounding_box
        self.VOLUME = self.MESH_OBJ.volume


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

def pointsToTriangles(points, triangles):

    with np.errstate(all='ignore'):
        # Unpack triangle points
        p0,p1,p2 = np.asarray(triangles).swapaxes(0, 1)
        # Calculate triangle edges
        e0 = p1-p0
        e1 = p2-p0
        a = inner1d(e0, e0)
        b = inner1d(e0, e1)
        c = inner1d(e1, e1)
        # Calculate determinant and denominator
        det = a*c - b*b
        invDet = 1. / det
        denom = a-2*b+c
        # Project to the edges
        p  = p0-points[:,np.newaxis]
        d = inner1d(e0,p)
        e = inner1d(e1,p)
        u = b*e - c*d
        v = b*d - a*e
        # Calculate numerators
        bd = b+d
        ce = c+e
        numer0 = (ce - bd) / denom
        numer1 = (c+e-b-d) / denom
        da = -d/a
        ec = -e/c
        # Vectorize test conditions
        m0 = u + v < det
        m1 = u < 0
        m2 = v < 0
        m3 = d < 0
        m4 = (a+d > b+e)
        m5 = ce > bd
        #
        t0 =  m0 &  m1 &  m2 &  m3
        t1 =  m0 &  m1 &  m2 & ~m3
        t2 =  m0 &  m1 & ~m2
        t3 =  m0 & ~m1 &  m2
        t4 =  m0 & ~m1 & ~m2
        t5 = ~m0 &  m1 &  m5
        t6 = ~m0 &  m1 & ~m5
        t7 = ~m0 &  m2 &  m4
        t8 = ~m0 &  m2 & ~m4
        t9 = ~m0 & ~m1 & ~m2
        #
        u = np.where(t0, np.clip(da, 0, 1), u)
        v = np.where(t0, 0, v)
        u = np.where(t1, 0, u)
        v = np.where(t1, 0, v)
        u = np.where(t2, 0, u)
        v = np.where(t2, np.clip(ec, 0, 1), v)
        u = np.where(t3, np.clip(da, 0, 1), u)
        v = np.where(t3, 0, v)
        u *= np.where(t4, invDet, 1)
        v *= np.where(t4, invDet, 1)
        u = np.where(t5, np.clip(numer0, 0, 1), u)
        v = np.where(t5, 1 - u, v)
        u = np.where(t6, 0, u)
        v = np.where(t6, 1, v)
        u = np.where(t7, np.clip(numer1, 0, 1), u)
        v = np.where(t7, 1-u, v)
        u = np.where(t8, 1, u)
        v = np.where(t8, 0, v)
        u = np.where(t9, np.clip(numer1, 0, 1), u)
        v = np.where(t9, 1-u, v)
        # Return closest points
        return (p0.T + u[:, np.newaxis]*e0.T + v[:, np.newaxis]*e1.T).swapaxes(2,1)

def POINT_DIST_TRI(p_coord, FACES):
    '''
    This function is designed to calculate the distance between a grid point to
    the input stl mesh
    If the projected normal is with the triangle, use the normal distance
    If the projected normal is outside the triangle, use nearest edge of tri
    '''
    for tri in FACES:

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

def CONSTRUCT_TRI(vertices, faces):
    '''
    This function takes the vertices and faces of a triangular mesh to construct
    a numpy array with mesh coordinates
    '''
    triangles = np.zeros((len(faces), 3, 3,))
    cnt = 0
    for fc in faces:
        p0 = [vertices[fc[0]][0], vertices[fc[0]][1], vertices[fc[0]][2]]
        p1 = [vertices[fc[1]][0], vertices[fc[1]][1], vertices[fc[1]][2]]
        p2 = [vertices[fc[2]][0], vertices[fc[2]][1], vertices[fc[2]][2]]
        triangles[cnt] = [p0, p1, p2]
        cnt += 1
    return triangles

def FILTER_GRID_POINTS(grid_in_mesh, grid_points):
    '''
    This function is designed to filter the grid points within the mesh and
    return the coordinates of points only in the mesh
    '''
    f_grid_points = np.zeros((sum(grid_in_mesh), 3))
    cnt = 0
    for in_flag, gp in zip(grid_in_mesh, grid_points):
        if in_flag:
            f_grid_points[cnt] = gp
            cnt += 1
    return f_grid_points

def GRID_MESH_DIST(f_grid_points, grid_tri_nearest, faces):
    '''
    This function is designed to calculate each grid point to the surfaces on
    the trimesh
    '''
    grid_distance = np.zeros((len(f_grid_points), len(faces)))
    cnt = 0
    for gp, gtn in zip(f_grid_points, grid_tri_nearest):
        tmp = np.zeros((len(faces), 3))
        tmp_gp = tmp + gp
        # Calculat the euclidean distance between two sets of matrix
        grid_distance[cnt] = cdist(tmp_gp, gtn, 'euclidean')[0]
        cnt += 1
    return grid_distance

def grid_update(eligible_grid, overlap):
    '''
    This function dynamically update the grid point array to be always having
    largest distance on top and descending order
    '''
    current_gp = eligible_grid[0]
    # Construct a tmp numpy array
    first_grid = np.zeros((len(eligible_grid), 4)) + current_gp
    # find the distance between first grid and other grid
    tmp_dist = cdist(first_grid[:, 0:3], eligible_grid[:, 0:3], 'euclidean')[0]
    # Construct a numpy arrray for the new eligible_grid
    tmp_dist_flag = tmp_dist <= current_gp[3]*overlap/100
    # ------------Below we update the eligible grid array
    eligible_grid = np.delete(eligible_grid, (0), axis=0)
    new_eligible_grid = np.zeros((len(eligible_grid)-sum(tmp_dist_flag), 4))
    cnt = 0
    for tdf, egd in zip(tmp_dist_flag, eligible_grid):
        if ~tdf:
            new_eligible_grid[cnt] = egd
            cnt += 1
    return new_eligible_grid

def FIND_SPHERE(grid_distance, r_thresh, f_grid_points, sphere_num, overlap):
    '''
    This function is designed to find all eligible spheres nominated by the user
    '''
    # find the largest amount the min distnace from grid to trimesh
    min_dist = np.array([grid_distance.min(1)])
    # Construc a new array
    dist_grid = np.concatenate((f_grid_points, min_dist.T), axis = 1)
    # Sort he dist_grid array to start finding the spheres from the largest to
    # increase filled volume
    dist_grid_sorted = np.flipud(dist_grid[dist_grid[:,3].argsort()])
    # from r_thresh to gradually reduce the distance to get Spheres
    # Once it is under the thresh, -> fill the sphere; -> delete grids
    # -> substract total sphere number
    # define the sphere array [x, y, z, radius]
    dist_flag = dist_grid_sorted[:, 3] <= r_thresh
    eligible_grid = np.zeros((sum(dist_flag), 4))
    cnt = 0
    for gp_f, gp in zip(dist_flag, dist_grid_sorted):
        if gp_f:
            eligible_grid[cnt] = gp
            cnt += 1
    # Define sphere array
    spheres = np.zeros((sphere_num, 4))
    sphere_cnt = 0
    while sphere_cnt < sphere_num:
        spheres[sphere_cnt] = eligible_grid[0]
        sphere_cnt += 1
        eligible_grid = grid_update(eligible_grid, overlap)
    return spheres

def SPH_FILL(mesh_obj, sphere_num, overlap, r_thresh):
    '''
    This is actually the main part of the function
    '''
    # construct the gird for the mesh first
    grid_points = GRID_MESH(mesh_obj.BOUNDS)
    # compute if grid points are within the mesh, output a boolean array
    grid_in_mesh = GRID_IN_HULL(grid_points, mesh_obj.VERTICES)
    # Filter the grid points so that the only within a mesh is retained
    f_grid_points = FILTER_GRID_POINTS(grid_in_mesh, grid_points)
    # Compute the distance if a point is within the mesh
    # Firstly, contruct the trianglular mesh
    triangles_coord = CONSTRUCT_TRI(mesh_obj.VERTICES, mesh_obj.FACES)
    # Calculate the nearest points between grid points and all triangles
    grid_tri_nearest = pointsToTriangles(f_grid_points, triangles_coord)
    # Calculate the distances between grid points and triangles
    # Return is a matrix with grid point to all mesh distances
    grid_distance = GRID_MESH_DIST(f_grid_points, grid_tri_nearest, mesh_obj.FACES)
    # Also added a threshold to have 1/10 of the equivalent radius of particle
    radius_thresh = pow((3*mesh_obj.VOLUME)/4/np.pi, 1/3)*r_thresh
    sphere_pop = FIND_SPHERE(grid_distance, radius_thresh, f_grid_points, sphere_num, overlap)
    return sphere_pop


def main():
    # Print command line arguments
    # stl_f        : Particle stl name
    # sphere_number: Total number spheres allowed
    # overlap      : Overlap percentage, 100 being no overlap
    stl_f = sys.argv[1]                 # 'test.stl'
    sphere_number = sys.argv[2]         # 5 20 40 ...
    overlap = sys.argv[3]               # 110 100 90 85 ...
    r_thresh = sys.argv[4]              # 0.5 0.4 0.3 ...
    # run the following functions
    print('Currently processing' + stl_f)
    tmp_mesh = STL_MESH(stl_f)
    SPH_FILL(tmp_mesh, sphere_number, overlap, r_thresh)

    '''
    # Below is the code for testing the sphere results in Paraview python

    '''

if __name__ == "__main__":
    main()
