#!/usr/bin/env python3

# A script to plot the 2D field from the netcdf output
# of the advanced tutorial code
#
# Requires Python 3, netCDF4, Matplotlib and Numpy.

# Usage: ./plot_ncfield.py <filename> <step> <field>

import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import netCDF4 as nc

def read_and_plot(filename):
    ncfile = nc.Dataset(filename)
    var = ncfile['field'][0,:,:]
    var = var.reshape(36,18)
    plt.imshow(np.transpose(var), origin='lower',interpolation='none')
    plt.title(filename)
    plt.colorbar(shrink=0.6)
    plt.show()

if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()
    try:
        read_and_plot(args.filename)
    except OSError as err:
        print(err)
