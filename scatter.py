import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("flatten_delta.csv")
#df = pd.read_csv("flatten_rho.csv")
#df = pd.read_csv("c1.csv")
mask = df.delta.values < 1.6
#mask = df.rho.values < 1.6
#plt.scatter(df.rho.values[mask], df.ratio.values[mask])
#plt.scatter(df.c1.values, df.c2.values)
plt.scatter(df.delta.values, df.ratio.values)
#xs = [x for x in np.arange(0.4,1.0,0.01)]
#plt.plot(xs,xs)
plt.show()
