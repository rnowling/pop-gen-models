import sys

from sklearn.lda import LDA
import matplotlib.pyplot as plt
import numpy as np

def read_variants(flname):
	fl = open(flname)
	markers = []
	individuals = []
	population_ids = []
	population = -1
	for ln in fl:
		if "Marker" in ln:
			if len(individuals) == 0:
				continue

			marker = dict()
			marker["individuals"] = np.array(individuals)
			marker["population_labels"] = np.array(population_ids)
			markers.append(marker)
			population = -1
			population_ids = []
			individuals = []
		elif "Population" in ln:
			population += 1
		else:
			individual = map(float, ln.strip().split())
			individuals.append(individual)
			population_ids.append(population)

	if len(individuals) != 0:
		marker = dict()
		marker["individuals"] = np.array(individuals)
		marker["population_labels"] = np.array(population_ids)
		markers.append(marker)
	fl.close()
	return markers

def plot_scores(markers, flname):
	plt.clf()
	scores = []
	for i, marker in enumerate(markers):
		try:
			lda = LDA()
			lda.fit(marker["individuals"], marker["population_labels"])
			scores.append(lda.score(marker["individuals"], marker["population_labels"]))
		except:
			scores.append(0.0)

	plt.hist(scores, bins=np.arange(0.0, 1.0, 0.01))

	plt.xlabel("Score", fontsize=18)
	plt.ylabel("Occurrences", fontsize=18)

	plt.savefig(flname, DPI=200)

def plot_lda_projection(marker, flname):
	lda = LDA()
	lda.fit(marker["individuals"], marker["population_labels"])
	print lda.score(marker["individuals"], marker["population_labels"])
	proj = lda.transform(marker["individuals"])
	n_samples, n_components = proj.shape

	plt.scatter(proj, marker["population_labels"])
	plt.xlabel("Component 0", fontsize=18)
	plt.ylabel("Population Labels", fontsize=18)

	plt.savefig(flname, DPI=200)


if __name__ == "__main__":
	variants_fl = sys.argv[1]
	#variant_id = int(sys.argv[2])
	plot_flname = sys.argv[2]

	variants = read_variants(variants_fl)

	print len(variants)

	#plot_lda_projection(variants[variant_id], plot_flname)
	plot_scores(variants, plot_flname)

