#!/usr/bin/env python3

# A script to plot the 2D field from the netcdf output
# of the advanced tutorial code
#
# Requires Python 3, netCDF4, Matplotlib and Numpy.

# Usage: ./plot_obs.py <filename> <step>

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import argparse as ap
import netCDF4 as nc

def read_and_plot(filename, step):
    ncfile = nc.Dataset(filename)
    istep = int(step)
    var = ncfile['obs'][istep-1,:,:]
    var = var.reshape(36,18)
    var_m = ma.masked_outside(var,-10.0,10.0)
    plt.imshow(np.transpose(var_m), origin='lower',interpolation='none')
    plt.title(filename)
    plt.colorbar(shrink=0.6)
    plt.show()

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('step')
    args = parser.parse_args()
    try:
        read_and_plot(args.filename, args.step)
    except OSError as err:
        print(err)
