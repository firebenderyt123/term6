import pandas as pd
import matplotlib.pyplot as plt

# read the data from the Excel file
df = pd.read_excel('output.xlsx')

# extract the columns
t = df['t']
h = df['h']
v = df['v']
h__ = df['h^']
v__ = df['v^']

# plot the data
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(t, h)
axs[0, 0].set_title('Висота')
axs[0, 1].plot(t, v)
axs[0, 1].set_title('Швидкість')
axs[1, 0].plot(t, h__)
axs[1, 0].set_title('Очікувана висота')
axs[1, 1].plot(t, v__)
axs[1, 1].set_title('Очікувана швидкість')

# show the plot
plt.show()
