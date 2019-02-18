import numpy as np
import os

omega = np.linspace(start=0, stop=3, num=10)

for om in omega:
    os.system("./cvm {} 0.1 10".format(om))