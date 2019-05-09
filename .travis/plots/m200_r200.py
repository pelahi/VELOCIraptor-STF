import sys

import h5py
from matplotlib import pyplot

input_fname = sys.argv[1]
output_fname = sys.argv[2]
with h5py.File(input_fname) as f:
    m200, r200 = f['/Mass_200mean'].value, f['/R_200mean'].value
pyplot.xscale('log'); pyplot.yscale('log')
pyplot.plot(r200, m200, marker='.', linestyle=' ')
pyplot.savefig(output_fname)
