import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file = 'out.xlsx'

# Read the CSV file into a DataFrame
df = pd.read_excel(file)

# Convert fi, omega, and mu to radians
df['fi'] = np.rad2deg(df['fi'])
df['omega'] = np.rad2deg(df['omega'])
df['mu'] = np.rad2deg(df['mu'])

# Create three plots for fi, omega, and mu
fig, axs = plt.subplots(3, 1, figsize=(8, 10))

axs[0].plot(df['t'], df['fi'])
axs[0].set_ylabel('fi')

axs[1].plot(df['t'], df['omega'])
axs[1].set_ylabel('omega')

axs[2].plot(df['t'], df['mu'])
axs[2].set_ylabel('mu')
axs[2].set_xlabel('Time')

plt.show()
