import sys

from sklearn.naive_bayes import MultinomialNB as MNB
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
			mnb = MNB()
			mnb.fit(marker["individuals"], marker["population_labels"])
			scores.append(mnb.score(marker["individuals"], marker["population_labels"]))
		except:
			scores.append(0.0)

	plt.hist(scores, bins=np.arange(0.0, 1.0, 0.01))

	plt.xlabel("Score", fontsize=18)
	plt.ylabel("Occurrences", fontsize=18)

	plt.savefig(flname, DPI=200)


if __name__ == "__main__":
	variants_fl = sys.argv[1]
	plot_flname = sys.argv[2]

	variants = read_variants(variants_fl)

	print len(variants)

	plot_scores(variants, plot_flname)

