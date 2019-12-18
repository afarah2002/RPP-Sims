y_right = [ \

]

y_left = [ \

]

import numpy as np
import matplotlib.pyplot as plt
from pylab import *

x = energy_LIST = list(np.logspace(-6., 3., num=500, endpoint=True, base=10)) # eV
if __name__ == '__main__':

	cutoff = 200
	x = energy_LIST[cutoff:] 
	y_r = y_right[cutoff:]
	y_l = y_left[cutoff:]

	xlog = [log10(i) for i in x]
	ylog_right = [log10(i) for i in y_r]
	ylog_left = [log10(i) for i in y_l]

	m_right, b_right = polyfit(xlog, ylog_right, 1)
	m_left, b_left = polyfit(xlog, ylog_left, 1)

	print m_right, "\n", b_right, "\n", m_left, "\n", b_left

	plt.figure()
	plt.xlabel("Energy (MeV)", fontsize=18)
	plt.xticks(fontsize=14)
	plt.ylabel("Optimal B/E ratio (T/MeV)", fontsize=18)
	plt.yticks(fontsize=14)
	plt.title("Optimal B/E ratio (T/MeV) vs. e+ Energy (MeV)", fontsize=24)
	# plt.plot(xlog, np.multiply(ylog_right, 1000), 'r', xlog, (np.multiply(xlog, m_right)+b_right)*1000, '--k')
	# plt.plot(xlog, np.multiply(ylog_left, 1000), 'y', xlog, (np.multiply(xlog, m_left)+b_left)*1000, '--k')
	# plt.plot(xlog, ylog_right)
	plt.xscale('log')	
	# plt.yscale('log')
	plt.plot(x, np.multiply(y_r, 1), label='right')
	plt.plot(x, np.multiply(y_l, 1), label='left')	
	# plt.ticklabel_format(axis='both', style='sci', scilimits=(-7,0))
	plt.legend()
	plt.show()