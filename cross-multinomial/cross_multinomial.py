import sys
import numpy as np
import numpy.random as npr
from calculate_phist import read_counts
from calculate_phist import normalize_haplotypes
from scipy.special import gammaln
import matplotlib.pyplot as plt

def log_factorial(n):
	return gammaln(n+1)

def log_multinomial(xs, ps):
	n = np.sum(xs)
	log_prob = log_factorial(n) - np.sum(log_factorial(xs)) + np.sum(xs * np.log(ps + 0.0000000000001))
	return log_prob

def locus_prob(locus_obs_counts, locus_freq):
		log_prob = 0.0
		n_pop = locus_obs_counts.shape[0]
		for p1 in xrange(n_pop):
			for p2 in xrange(n_pop):
				log_prob += log_multinomial(locus_obs_counts[p1], locus_freq[p2])
		return log_prob

def probability(observed_counts):
	observed_frequencies = normalize_haplotypes(observed_counts)
	n_loci = observed_counts.shape[0]
	locus_probabilities = np.zeros(n_loci)
	for locus in xrange(n_loci):
		prob = locus_prob(observed_counts[locus, :, :], observed_frequencies[locus, :, :])
		locus_probabilities[locus] = prob
	return locus_probabilities

def plot_log_prob(flname, log_probs):
	plt.clf()
	plt.hold(True)
	plt.hist(log_probs, bins=30)
	plt.xlabel("Log Probability", fontsize=16)
	plt.xlim([min(log_probs), 0.0])
	plt.ylabel("Occurrences (Loci)", fontsize=16)
	plt.savefig(flname, DPI=200)

def main(occur_fl, plot_basename):
	observed_counts = read_counts(occur_fl)
	print observed_counts.shape
	locus_log_probs = probability(observed_counts)
	print locus_log_probs.shape
	plot_log_prob(plot_basename + "_log_prob.pdf", locus_log_probs)

if __name__ == "__main__":
	occur_fl = sys.argv[1]
	plot_basename = sys.argv[2]

	main(occur_fl, plot_basename)

	