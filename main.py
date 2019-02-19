import numpy as np
import os
import csv
from mpi4py import MPI
import multiprocessing as mp
from functools import partial
import glob
from decimal import Decimal
import re


def split_data(size, data):

    """Split the data on the MPI ranks"""

    n = len(data) // size
    return [data[i:i + n] for i in range(0, len(data), n)]


def run_exe(omega, eta, num_iter, i, j):

    """Run the executable generated by the C++ file"""

    os.system("./cvm {} {} {} {} {}".format(omega, eta, num_iter, i, j))


# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()  # Identification number of node
size = comm.Get_size()  # Number of nodes

# Define data
omega = np.linspace(start=-2., stop=4., num=int(40/0.05), dtype=float)

# Split the data in the zeroth node
if rank == 0:
    data_split = split_data(size, omega)
else:
    data_split = None

# Scatter data from zeroth node onto other nodes and do the calculation in each node
data_in_node = comm.scatter(data_split, root=0)

# Split running of exe in each rank on CPUs
p = mp.Pool(mp.cpu_count())
run_exe_ = partial(run_exe, eta=0.05, num_iter=3, i=1, j=1)
p.map(run_exe_, data_in_node)


correlator = []

# Make sure that all ranks have finished
comm.gather(0, root=0)

# Sort the filenames numerically
filenames = sorted([(element[9:-4]) for element in glob.glob('data/*.txt')], key=float)

if rank == 0:

    for filename in filenames:

        name = 'data/cvm_{}.txt'.format(filename)
        f = open('data/cvm_{}.txt'.format(filename))
        correlator.append(float(f.readline()))
        f.close()
        os.remove(name)

        # Save data in file
    f = open('data/cvm_data.csv', 'w')
    out = csv.writer(f, delimiter=' ')
    out.writerows(zip(omega, correlator))
    f.close()

