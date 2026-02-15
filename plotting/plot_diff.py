#!/usr/bin/env python3

# A script to plot the 2D field from the online_2D_serialmodel tutorial.
# Requires Python 3, Matplotlib and Numpy.

# Usage: ./plot_field.py <filename1> <filename2>

import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import netCDF4 as nc

def read_and_plot(filename1, filename2):
    ncfile1 = nc.Dataset(filename1)
    ncfile2 = nc.Dataset(filename2)
    field1 = ncfile1['field'][0,:,:]
    field2 = ncfile2['field'][0,:,:]
    field1 = field1.reshape(36,18)
    field2 = field2.reshape(36,18)
    rmse = 0;
    for i in range(36):
       for j in range(18):
	       rmse = rmse + (field1[i,j]-field2[i,j])**2
    rmse = np.sqrt(1/(18*36)*(rmse));
    print('RMSE: ', rmse)
    plt.imshow(np.transpose(field1-field2), origin='lower',interpolation='none',cmap='PiYG')
    plt.colorbar(shrink=0.6)
    plt.show()

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('filename1')
    parser.add_argument('filename2')
    args = parser.parse_args()
    try:
        read_and_plot(args.filename1, args.filename2)
    except OSError as err:
        print(err)
