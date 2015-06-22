# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 12:24:07 2015

Package for post-processing and plotting tsunami model results.

@author: crozas
"""

import numpy as np

#-----------------------------------------------------BOF-------------------------------------------------
def extract_sww_grid(swwfile):
    """   
    Read grid data from sww file.
    swwfile: netCDF4 object
    """
    x = swwfile.variables['x'][:] +  swwfile.xllcorner
    y = swwfile.variables['y'][:] +  swwfile.yllcorner
    z = swwfile.variables['elevation'][:]
    triangles= swwfile.variables['volumes'][:]   
        
    return x,y,z,triangles
#-----------------------------------------------------EOF------------------------------------------------- 

#-----------------------------------------------------BOF------------------------------------------------- 
def surface(eta, z, drytol=-1e-2):
   """
   Identify wet nodes
   
   Return a masked array containing the surface elevation only in wet cells.
   Surface is eta = h+topo
   """
   h = eta-z
   water = np.ma.masked_where(h<=drytol, eta)
   return water
#-----------------------------------------------------EOF------------------------------------------------- 

#-----------------------------------------------------BOF-------------------------------------------------    
def mask_dry_triangles(eta, triang_wet, triangles):
    """
    Mask off undesired dry triangles
    
    eta: masked array [:] surface elevation
    triang_wet: triangulation, triangles vertices to be masked in wet triangles
    triangles: array[:,3] triangles vertices
    """
    eta_mask_indices = np.where(eta.mask)
    tri_mask = np.any(np.in1d(triangles, eta_mask_indices).reshape(-1,3), axis=1)
    triang_wet.set_mask(tri_mask)
#-----------------------------------------------------EOF------------------------------------------------- 

#-----------------------------------------------------BOF-------------------------------------------------    
def mask_wet_triangles(eta, triang_dry, triangles):
    """
    Mask off undesired wet triangles
    
    eta: masked array [:] surface elevation
    triang_dry: trinagulation, triangles vertices to be masked in dry triangles
    triangles: array[:,3] triangles vertices
    """
    eta_mask_inverted = [not i for i in eta.mask]
    eta_mask_indices = np.where(eta_mask_inverted)
    tri_mask = np.any(np.in1d(triangles, eta_mask_indices).reshape(-1,3), axis=1)
    triang_dry.set_mask(tri_mask)
#-----------------------------------------------------EOF------------------------------------------------- 
