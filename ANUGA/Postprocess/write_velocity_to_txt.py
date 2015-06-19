# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:12:16 2015

@author: crozas
"""

# tried to include a bit of exception handling--------------------
try:
    import netCDF4
except:
    import Scientific.IO.NetCDF as netCDF4
# -----------------------------------------------------------------
import numpy as np

#Read netCDF file
infile = netCDF4.Dataset('flood_fixed.sww')
#print(nc)

#Define timesteps to plot
#Mejora (CRR): crear una funcion que lea los dt del archivo original y calcule el dt para los plots.
t = infile.variables['time']
dt = t[1]-t[0]
#Plot every 1minute
timestep_plot = int(60/dt) #Datos cada 10s en archivo original, se guardan cada 1 minuto para video y resultados.
t_plot = np.arange(0,len(t),timestep_plot)

#Read variables
x = infile.variables['x']
y = infile.variables['y']
z = infile.variables['elevation']
triangles = infile.variables['volumes']

times=[144,309,324,418]
stage = infile.variables['stage'][times,:]
mx = infile.variables['xmomentum'][times,:]
my = infile.variables['ymomentum'][times,:]

#stage = infile.variables['stage'][t_plot,:]
#mx = infile.variables['xmomentum'][t_plot,:]
#my = infile.variables['ymomentum'][t_plot,:]

#Define new variables h, u, v, vel
n_t_plot=len(times)
#n_t_plot=len(t_plot)
n_points=len(z)

h = np.zeros((n_t_plot,n_points))
u = np.zeros((n_t_plot,n_points))
v = np.zeros((n_t_plot,n_points))
vel = np.zeros((n_t_plot,n_points))

for i in range(n_t_plot):
    h_aux = stage[i,:]-z[:]
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

fid = open('velocities.csv','w')
fid.write('x[mUTM],y[mUTM],u[m/s],v[m/s],velocity[m/s]\n')

print('Guardando resultados en : velocities.csv')
for t in range(len(times)):
    fid.write('Time: '+repr(times[t])+'\n')
    for i in range(n_points):
        fid.write(repr(x[i])+','+repr(y[i])+','+repr(u[t,i])+','+repr(v[t,i])+','+repr(vel[t,i])+'\n')

fid.close()

infile.close(); print('Dataset is closed!')








