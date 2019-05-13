import sys

import h5py
from matplotlib import pyplot as plt
import numpy as np

# Parameters from command-line, type-safe
input_fname, dataset, output_fname, do_log, bins = sys.argv[1:]
do_log, bins = int(do_log), int(bins)

with h5py.File(input_fname) as f:
    m200 = f[dataset][()]
label = dataset.strip('/')
if do_log:
    m200 = np.log10(m200)
    label = 'log(%s)' % label
plt.hist(m200, bins=bins)
plt.xlabel(label)
plt.ylabel('Count')
plt.savefig(output_fname)

