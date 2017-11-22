#!/usr/bin/python
"""
Thi script is used for creating sphere clumps to fill the STL mesh in LIGGGHTS
simulations. The method used is based on the following paper:
* Amberger S, Friedl M, Goniva C, Pirker C, Approximation of Objects By Spheres
for multisphere simulations in DEM

https://www.researchgate.net/publication/259693012

@author: Dr Wei Chen

"""

from paraview.simple import *

fop = open('C:\Users\wchen\Desktop\Sphere_Clump_Gen-master\sphere.txt', 'r')
lines = fop.readlines()
for line in lines:
    tmp = line.strip('\n')
    tmp_str = tmp.split(',')
    Sphere(Center=[float(tmp_str[0]), float(tmp_str[1]), float(tmp_str[2])], Radius = float(tmp_str[3]))

fop.close()
