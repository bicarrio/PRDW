# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 14:25:15 2014

@author: crozas
"""

import netCDF4
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import os
import prdw_tsunami_tools as ptt

#Folder name where the sww files are
folder_sww_files = 'Conchali'

#Main directory
main_dir = os.getcwd()  

#Read netCDF file
infile = netCDF4.Dataset(main_dir + '\\' + folder_sww_files + '\\' + 'velocities.sww')

x, y, z, triangles = ptt.extract_sww_grid(infile)

#Read variables
time = infile.variables['time'][:]
stage = infile.variables['stage'][:]
u = infile.variables['u'][:]
v = infile.variables['v'][:]
vel = infile.variables['velocity'][:]
timestep = infile.timestep

#Define timestep to plot
time_to_plot = 25*60 #time to plot in seconds
t_plot = time_to_plot/timestep

#Define variables to plot
surf = ptt.surface(stage[t_plot], z)
triang_wet = tri.Triangulation(x, y, triangles)
ptt.mask_dry_triangles(surf, triang_wet, triangles)

#Inicialize figure
fig = plt.figure()
ax = plt.gca()
ax.set_aspect('equal')

#Plot topobathymetry
levels=np.linspace(min(z),max(z),101)
topo = ax.tricontourf(x, y, triangles, z, levels, cmap='terrain')

#Plot stage
wave = ax.tricontourf(triang_wet, surf, [-8, 0, 8], cmap='bwr')
