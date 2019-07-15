import matplotlib.pyplot as plt
import pandas as pd

file = 'data/heisenberg_data.csv'
data = pd.read_csv(file, delimiter=' ', header=None, dtype=float)
plt.plot(data[0], data[1], label='cvm')
plt.xlabel('omega')
plt.ylabel('spectral function')
plt.legend()
plt.show()