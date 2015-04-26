import sys
import numpy as np
import numpy.random as npr
from sklearn.neighbors.kde import KernelDensity
from scipy.special import gammaln
import matplotlib.pyplot as plt
from calculate_phist import read_counts
from calculate_phist import normalize_haplotypes

def log_factorial(n):
	return gammaln(n+1)

def log_multinomial(xs, ps):
	n = np.sum(xs)
	log_prob = log_factorial(n) - np.sum(log_factorial(xs)) + np.sum(xs * np.log(ps + 0.0000000000001))
	return log_prob

class KDE_MCMC_Sampler(object):
	def __init__(self, observed_counts):
		"""
		Observed counts is 3D matrix of pop, locus, haplotype
		"""
		self.observed_counts = observed_counts
		self.individual_counts = observed_counts.sum(axis=2)
		self.observed_frequencies = normalize_haplotypes(observed_counts)

		self.n_loci, self.n_pop, self.n_haplotypes = self.observed_counts.shape
		
		# from bamova
		self.DWEIGHT = 1.0
		self.DADD = 0.00001
		self.SMALL_NUM = 0.0000000000001

		print "initializing frequencies"
		self.freq = np.zeros((self.n_loci, self.n_haplotypes))
		for l in xrange(self.n_loci):
			self.freq[l, :] = self.sample_locus_freq(self.observed_frequencies[l, 0, :])

	def sample_locus_freq(self, freq):
		alphas = self.DWEIGHT * freq + self.DADD + self.SMALL_NUM

		return npr.dirichlet(alphas)

	def locus_prob(self, locus_obs_counts, locus_freq):
		log_prob_sum = 0.0
		for p in xrange(self.n_pop):
			log_prob_sum += log_multinomial(locus_obs_counts[p], locus_freq)
		return log_prob_sum


	def step(self):
		total_log_prob = 0.0
		for l in xrange(self.n_loci):
			locus_indiv_counts = self.individual_counts[l, :]
			locus_obs_counts = self.observed_counts[l, :, :]

			log_prob = self.locus_prob(locus_obs_counts, self.freq[l, :])

			proposed_locus_freq = self.sample_locus_freq(self.freq[l, :])

			proposed_log_prob = self.locus_prob(locus_obs_counts, proposed_locus_freq)
				
			log_prob_ratio = proposed_log_prob - log_prob
			log_r = np.log(npr.random())

			if proposed_log_prob >= log_prob or log_r <= log_prob_ratio:
				self.freq[l, :] = proposed_locus_freq
				log_prob = proposed_log_prob

			total_log_prob += log_prob

		locus_prob = []
		for l in xrange(self.n_loci):
			log_prob = self.locus_prob(locus_obs_counts, self.freq[l, :])
			locus_prob.append(log_prob)

		return self.freq, total_log_prob, locus_prob

def plot_log_prob(flname, log_probs):
	plt.clf()
	plt.hold(True)
	plt.hist(log_probs, bins=30)
	plt.xlabel("Log Probability", fontsize=16)
	plt.xlim([min(log_probs), 0.0])
	plt.ylabel("Occurrences (Loci)", fontsize=16)
	plt.savefig(flname, DPI=200)

def simulate(occur_fl, n_steps, plot_basename):
	print "reading occurrences"
	observed_counts = read_counts(occur_fl)
	individual_counts = observed_counts.sum(axis=2)
	observed_frequencies = normalize_haplotypes(observed_counts)

	sampler = KDE_MCMC_Sampler(observed_counts)

	locus_log_prob = []
	for i in xrange(n_steps):
		freq, log_prob, locus_log_prob = sampler.step()
		print "step", i, "log prob", log_prob

	plot_log_prob(plot_basename + "_log_prob.pdf", locus_log_prob)

if __name__ == "__main__":
	occur_fl = sys.argv[1]
	n_steps = int(sys.argv[2])
	plot_basename = sys.argv[3]

	simulate(occur_fl, n_steps, plot_basename)

	