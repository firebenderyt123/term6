import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# read the data from the Excel file
df = pd.read_excel('output.xlsx')

# extract the columns
t = df['t']
omega = np.rad2deg(df['omega'])
fi = np.rad2deg(df['fi'])

# plot the data
plt.plot(t, omega, label='Значення кута')
plt.plot(t, fi, label='Значення кутової швидкості')
plt.xlabel('Час, сек')

# show the plot
plt.legend()
plt.show()
