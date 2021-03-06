#!/usr/bin python3

# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

import matplotlib.pyplot as plt
import numpy as np
import plotly.plotly as py  # tools to communicate with Plotly's server

numpy_hist = plt.figure()

plt.hist([1, 2, 1], bins=[0, 1, 2, 3])

plot_url = py.plot_mpl(numpy_hist, filename='numpy-bins')
