#!/usr/bin/env python3

import os, sys, re
import matplotlib
import numpy as np
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

n = int(sys.argv[1])

def run_vcf_parser(n):
	execution = 0
	try:
		x = list(np.random.randn(n))
		print(x)
		y = .4 * x + np.random.randn(100000) + 5

		fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)

		# We can set the number of bins with the `bins` kwarg
		axs[0].hist(x, bins=20)
		axs[1].hist(y, bins=20) 
	except Exception as e:
		print(str(e))
		execution = 1
	return execution 

run_vcf_parser(n)
sys.exit(0)
