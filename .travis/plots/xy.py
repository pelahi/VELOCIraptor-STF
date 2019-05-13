import sys

import h5py
from matplotlib import pyplot

# Parameters from command-line, type-safe
input_fname, ds_x, ds_y, output_fname, x_log, y_log = sys.argv[1:]
x_log, y_log = int(x_log), int(y_log)

with h5py.File(input_fname) as f:
    x, y = f[ds_x][()], f[ds_y][()]

xlabel, ylabel = ds_x.strip('/'), ds_y.strip('/')
pyplot.plot(x, y, marker='.', linestyle=' ')
pyplot.xlabel(xlabel)
pyplot.ylabel(ylabel)
if x_log:
    pyplot.xscale('log')
if y_log:
    pyplot.yscale('log')
pyplot.savefig(output_fname)
