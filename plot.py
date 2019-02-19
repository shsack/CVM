import matplotlib.pyplot as plt
import pandas as pd

file = 'data/cvm_data.csv'
data = pd.read_csv(file, delimiter=' ', header=None, dtype=float)
plt.plot(data[0], data[1])
plt.show()