import sys

import numpy as np
from sklearn.neighbors.kde import KernelDensity
import matplotlib
matplotlib.use("cairo")
from matplotlib import pyplot as plt

def estimate_distribution(samples, h=0.1, n_points=100):
	kde = KernelDensity(bandwidth=h)
	min_xs = min(samples)
	max_xs = max(samples)
	samples = samples[:, np.newaxis]
	kde.fit(samples)
	xs = np.linspace(min_xs, max_xs, n_points)
	ys = kde.score_samples(xs[:, np.newaxis])
	print xs.shape, ys.shape
	return xs, ys

def plot_log_probs(plot_flname, subset):
	plt.hold(True)
	for loci in xrange(subset.shape[1]):
		samples = subset[:, loci]
		xs, ys = estimate_distribution(samples)
		plt.plot(xs, ys)
	plt.xlabel("Log Prob", fontsize=16)
	plt.ylabel("Frequency", fontsize=16)
	plt.savefig(plot_flname, DPI=200)


if __name__ == "__main__":
	npy_flname = sys.argv[1]
	start = int(sys.argv[2])
	end = int(sys.argv[3])
	plot_flname = sys.argv[4]

	phi_values = np.load(npy_flname)
	subset = phi_values[:, start:end]

	plot_log_probs(plot_flname, subset)





