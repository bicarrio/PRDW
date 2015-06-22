# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 11:10:01 2015

@author: crozas
"""

import netCDF4
import matplotlib.pyplot as plt
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

#Plot topobathymetry
fig = plt.figure()
ax = plt.gca()
ax.set_aspect('equal')

#Plot all cells (dry)
levels=np.linspace(min(z),max(z),101)
topo = ax.tricontourf(x, y, triangles, z, levels, cmap='terrain')
tcbar = fig.colorbar(topo, pad = 0.05)
topo_contours = ax.tricontour(x, y, triangles, z, colors='k')

##z=0
#level0=np.linspace(-0.2,0.2,3)
#z0 = ax.tricontour(x, y, tri_mesh_pointers, z, [0],cmap='binary',lw=6)

def set_axis_limits(lmts,dx=500):
    x1 = dx*(np.trunc(lmts[0]/dx))
    x2 = dx*(np.trunc(lmts[1]/dx)+1)
    
    return [x1,x2]

def resadjust(ax, xres=None, yres=None):
    """
    Send in an axis and I fix the resolution as desired.
    """

    if xres:
        start, stop = ax.get_xlim()
        ticks = np.arange(start, stop + xres, xres)
        ax.set_xticks(ticks)

    if yres:
        start, stop = ax.get_ylim()
        ticks = np.arange(start, stop + yres, yres)
        ax.set_yticks(ticks)

res = 2000
xl = set_axis_limits([min(x), max(x)], res)
yl = set_axis_limits([min(y), max(y)], res)
ax.axis(np.append(xl,yl))

#resadjust(ax,xres=2000,yres=2000)

ax.set_xticks(np.arange(xl[0],xl[1]+res,res))
ax.set_yticks(np.arange(yl[0],yl[1]+res,res))

ax.set_xticklabels(np.arange(xl[0],xl[1]+res,res), rotation = 90)
ax.set_yticklabels(np.arange(yl[0],yl[1]+res,res))

plt.xlabel('Easting [m]')
plt.ylabel('Northing [m]')
plt.title('Topobatimetria')
plt.tight_layout()
plt.show()

plt.savefig('Topobathymetry.png', format='png', dpi=300)
#plt.close("all")