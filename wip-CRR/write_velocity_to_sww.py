# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:12:16 2015

@author: crozas
"""

import netCDF4
import numpy as np

#Read netCDF file
infile = netCDF4.Dataset('flood_fixed.sww')

#Define timesteps to plot
#Mejora (CRR): crear una funcion que lea los dt del archivo original y calcule el dt para los plots.
t = infile.variables['time']
dt = t[1]-t[0]
#Plot every 1minute
timestep_plot = int(60/dt) #Datos cada 10s en archivo original, se guardan cada 1 minuto para video y resultados.
t_id_plot = np.arange(0,len(t),timestep_plot)

#Read variables
x = infile.variables['x']
y = infile.variables['y']
z = infile.variables['elevation']
triangles = infile.variables['volumes']

ti = t[t_id_plot]
st = infile.variables['stage'][t_id_plot,:]
mx = infile.variables['xmomentum'][t_id_plot,:]
my = infile.variables['ymomentum'][t_id_plot,:]

#Define new variables h, u, v, vel
n_t_plot = len(ti)
n_points = len(infile.dimensions['number_of_points'])

h = np.zeros((n_t_plot,n_points))
u = np.zeros((n_t_plot,n_points))
v = np.zeros((n_t_plot,n_points))
vel = np.zeros((n_t_plot,n_points))

#Calculate h, u, v, vel
for i in range(n_t_plot):
    h_aux = st[i,:]-z[:]
    u_aux = mx[i,:]/h_aux
    v_aux = my[i,:]/h_aux
    u_aux[np.isnan(u_aux)] = 0
    v_aux[np.isnan(v_aux)] = 0
    u_aux[np.isinf(u_aux)] = 0
    v_aux[np.isinf(v_aux)] = 0
    
    h[i,:] = h_aux
    u[i,:] = u_aux
    v[i,:] = v_aux    
    vel[i,:] = np.sqrt(u[i,:]*u[i,:]+v[i,:]*v[i,:])

#Write new netCDF file with calculated values
outfile = netCDF4.Dataset('velocities.sww',mode='w',format='NETCDF4')

#Write sww file attributes
outfile.institution = 'PRDW'
outfile.author = 'CRR'
outfile.starttime = infile.getncattr('starttime')
outfile.xllcorner = infile.getncattr('xllcorner')
outfile.yllcorner = infile.getncattr('yllcorner')
outfile.zone = infile.getncattr('zone')
outfile.false_easting = infile.getncattr('false_easting')
outfile.false_northing = infile.getncattr('false_northing')
outfile.datum = infile.getncattr('datum')
outfile.projection = infile.getncattr('projection')
outfile.units = infile.getncattr('units')

#Define dimensions
number_of_volumes = outfile.createDimension('number_of_volumes', len(infile.dimensions['number_of_volumes']))
number_of_triangle_vertices = outfile.createDimension('number_of_triangle_vertices', len(infile.dimensions['number_of_triangle_vertices']))
number_of_vertices = outfile.createDimension('number_of_vertices', len(infile.dimensions['number_of_vertices']))
number_of_points = outfile.createDimension('number_of_points', n_points)
number_of_timesteps = outfile.createDimension('number_of_timesteps', None)

#Copy unchanged variables in sww file
def copy_netcdf_variable(varin, var_name,nc_out=outfile):
    outVar = nc_out.createVariable(var_name, varin.dtype, varin.dimensions)
    #Copy variable attributes    
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    outVar[:]=varin[:]

x = copy_netcdf_variable(x,'x',outfile)
y = copy_netcdf_variable(y,'y',outfile)
z = copy_netcdf_variable(z,'elevation',outfile)
triangles = copy_netcdf_variable(triangles,'volumes',outfile)

#Write new variables to sww file
time = outfile.createVariable('time', np.float32, ('number_of_timesteps',))
time[:] = ti[:]

stage = outfile.createVariable('stage', np.float32, ('number_of_timesteps','number_of_points'))
stage[:] = st[:]

u_out = outfile.createVariable('u', np.float32, ('number_of_timesteps','number_of_points'))
u_out[:] = u[:]
u_out.units='m/s'

v_out = outfile.createVariable('v', np.float32, ('number_of_timesteps','number_of_points'))
v_out[:] = v[:]
v_out.units='m/s'

vel_out = outfile.createVariable('velocity', np.float32, ('number_of_timesteps','number_of_points'))
vel_out[:] = vel[:]
vel_out.units='m/s'

#Close files
infile.close(); outfile.close(); print('Datasets are closed!')