# Script to generate input files for the PDAF online tutorial for 2 fields;
# 'field' is the same field as in the one-fields case; 'fieldB' is added here
#
# This is a translation of the former Matlab script to Python
#
# Note, that the random numbers are different in Matlab and
# Python, so that the varification outputs need to be consistent.
#
# L. Nerger, 9/2024
# update for netCDF files, L. Nerger, 2/2026

import numpy as np
import netCDF4 as nc

dim_x = 36         # Grid dimension in x-direction
dim_y = 18         # Grid dimension in y-direction
dim_ens = 20       # Maximum ensemble size
dim_step = 100     # Number of time steps
stddev_obs = 0.5   # error standard deviation for observations type A
stddev_obsB = 0.25 # error standard deviation for observations type B
dxobs = 5          # x-Grid spacing for observations type A
dyobs = 4          # y-Grid spacing for observations type A
dxobsB = 6         # x-Grid spacing for observations type B
dyobsB = 5         # y-Grid spacing for observations type B
obsB_offsetx = -2  # x-offset in position of observations type B
obsB_offsety = -1  # y-offset in position of observations type B

dowrite = 1        # 1 to write files

# Locations of observations not placed at grid points (x, y)
obs_interp = np.zeros((11,2))
obs_interp[:,:] = [[3.0, 2.1], 
     [3.4, 6.8],  
     [6.1, 6.8], 
     [8.9, 7.6], 
     [8.9, 14.9], 
     [20.0, 6.4], 
     [20.4, 16.1], 
     [14.1, 10.2], 
     [31.0, 5.2], 
     [31.2, 11.9], 
     [28.9, 14.9]];

# Initialize random number generator
np.random.seed(0)

# Generate true field A

field = np.zeros((dim_y, dim_x, dim_step+1))
for j in range(dim_x):
   for i in range(dim_y):
      field[i,j,0] = np.sin(2*np.pi*((i+1)/dim_y + (j+1)/dim_x))

for step in range(1,dim_step+1):
   for i in range(dim_y-1):
      field[i+1,:,step] = field[i,:,step-1]
   field[0,:,step] = field[-1,:,step-1]


# Generate ensemble states A and their mean

ens = np.zeros((dim_y, dim_x, dim_ens))
for k in range(dim_ens):
   for j in range(dim_x):
      for i in range(dim_y):
         ens[i,j,k] = np.sin(2*np.pi*((i+1)/dim_y + (j+1)/dim_x) + 2*0.5*np.pi*(k+1+5)/dim_ens)

state = np.mean(ens,axis=2)


# Generate true field B

fieldB = np.zeros((dim_y, dim_x, dim_step+1))
for j in range(dim_x):
   for i in range(dim_y):
      fieldB[i,j,0] = np.sin(4*np.pi*((i+1)/dim_y - (j+1)/dim_x))

for step in range(1,dim_step+1):
   for i in range(dim_y-1):
      fieldB[i+1,:,step] = fieldB[i,:,step-1]
   fieldB[0,:,step] = fieldB[-1,:,step-1]

# Write truth file for both fields
if dowrite==1:
   ncfile = nc.Dataset('true.nc',mode='w')
   xdim = ncfile.createDimension('dim_x', dim_x)
   ydim = ncfile.createDimension('dim_y', dim_y)
   timedim = ncfile.createDimension('step', dim_step)
   trueA_ini = ncfile.createVariable('trueA_ini',np.float64, ('dim_x', 'dim_y',))
   trueA = ncfile.createVariable('trueA',np.float64, ('step', 'dim_x', 'dim_y',))
   trueB_ini = ncfile.createVariable('trueB_ini',np.float64, ('dim_x', 'dim_y',))
   trueB = ncfile.createVariable('trueB',np.float64, ('step', 'dim_x', 'dim_y',))

   trueA_ini[:,:] = np.transpose(field[:,:,0])
   trueB_ini[:,:] = np.transpose(fieldB[:,:,0])
   for step in range(1,dim_step+1):
      trueA[step-1,:,:] = np.transpose(field[:,:,step])
      trueB[step-1,:,:] = np.transpose(fieldB[:,:,step])
   ncfile.close()


# Generate ensemble states B and their mean

ensB = np.zeros((dim_y, dim_x, dim_ens))
for k in range(dim_ens):
   for j in range(dim_x):
      for i in range(dim_y):
         ensB[i,j,k] = np.sin(4*np.pi*((i+1)/dim_y - (j+1)/dim_x) + 4*0.5*np.pi*(k+1+5)/dim_ens)

# Compute ensemble mean = initial state estimate

stateB = np.mean(ensB,axis=2)

# Write ensemble file for both fields
if dowrite==1:
   ncfile = nc.Dataset('ens.nc',mode='w')
   xdim = ncfile.createDimension('dim_x', dim_x)
   ydim = ncfile.createDimension('dim_y', dim_y)
   timedim = ncfile.createDimension('member', dim_step)
   trueA_ini = ncfile.createVariable('meanstateA',np.float64, ('dim_x', 'dim_y',))
   trueA = ncfile.createVariable('ensA',np.float64, ('member', 'dim_x', 'dim_y',))
   trueB_ini = ncfile.createVariable('meanstateB',np.float64, ('dim_x', 'dim_y',))
   trueB = ncfile.createVariable('ensB',np.float64, ('member', 'dim_x', 'dim_y',))

   trueA_ini[:,:] = np.transpose(state[:,:])
   trueB_ini[:,:] = np.transpose(stateB[:,:])
   for member in range(dim_ens):
      trueA[member,:,:] = np.transpose(ens[:,:,member])
      trueB[member,:,:] = np.transpose(ensB[:,:,member])
   ncfile.close()


# Observations A

obs_error = np.zeros((dim_y, dim_x, dim_step+1))
full_obs = np.zeros((dim_y, dim_x, dim_step+1))
obs_error = stddev_obs * np.random.randn(dim_y, dim_x, dim_step+1) 

full_obs[:,:,:] = field[:,:,:] + obs_error

obs = np.zeros((dim_y, dim_x, dim_step+1)) 
obs[:,:,:] = -999

for step in range(1,dim_step+1):
   for j in range(dxobs-1,dim_x,dxobs):
      for i in range(dyobs-1,dim_y,dyobs):
         obs[i,j,step] = full_obs[i,j,step]

if dowrite==1:
   ncfile = nc.Dataset('obsA.nc',mode='w')
   xdim = ncfile.createDimension('dim_x', dim_x)
   ydim = ncfile.createDimension('dim_y', dim_y)
   timedim = ncfile.createDimension('step', dim_step)
   true = ncfile.createVariable('obsA',np.float64, ('step', 'dim_x', 'dim_y',))

   for step in range(dim_step):
      true[step,:,:] = np.transpose(obs[:,:,step])
   ncfile.close()


# Observations B

obs_errorB = np.zeros((dim_y, dim_x, dim_step+1))
full_obsB = np.zeros((dim_y, dim_x, dim_step+1))
obs_errorB = stddev_obsB * np.random.randn(dim_y, dim_x, dim_step+1) 

full_obsB[:,:,:] = fieldB[:,:,:] + obs_errorB


obsB = np.zeros((dim_y, dim_x, dim_step+1)) 
obsB[:,:,:] = -999

for step in range(1,dim_step+1):
   for j in range(dxobsB-1+obsB_offsetx,dim_x,dxobsB):
      for i in range(dyobsB-1+obsB_offsety,dim_y,dyobsB): 
         obsB[i,j,step] = full_obsB[i,j,step]

if dowrite==1:
   ncfile = nc.Dataset('obsB.nc',mode='w')
   xdim = ncfile.createDimension('dim_x', dim_x)
   ydim = ncfile.createDimension('dim_y', dim_y)
   timedim = ncfile.createDimension('step', dim_step)
   true = ncfile.createVariable('obsB',np.float64, ('step', 'dim_x', 'dim_y',))

   for step in range(dim_step):
      true[step,:,:] = np.transpose(obsB[:,:,step])
   ncfile.close()
