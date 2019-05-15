import sys

import h5py
from matplotlib import pyplot as plt
import numpy as np

# Parameters from command-line, type-safe
input_fname, dataset, output_fname, do_log, bins = sys.argv[1:]
do_log, bins = int(do_log), int(bins)

with h5py.File(input_fname) as f:
    x = f[dataset][()]
label = dataset.strip('/')
if do_log:
    positives = np.where(x > 0)
    x[positives] = np.log10(x[positives])
    label = 'log(%s)' % label
plt.hist(x, bins=bins)
plt.xlabel(label)
plt.ylabel('Count')
plt.savefig(output_fname)

