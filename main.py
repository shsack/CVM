import numpy as np
import os
import csv

omega = np.linspace(start=0, stop=3, num=2, dtype=float)

for om in omega:
    os.system("./cvm {} 0.1 10".format(om))

correlator = []

for om in omega:

    name = "data/cvm_{}.txt".format(om)
    f = open(name)
    correlator.append(float(f.readline()))
    f.close()
    os.remove(name)

print(correlator)

f = open('data/cvm_data.csv', 'w')
out = csv.writer(f, delimiter=' ')
out.writerows(zip(omega, correlator))
f.close()