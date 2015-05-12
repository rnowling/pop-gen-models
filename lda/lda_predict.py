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

def predict_scores(markers, threshold=0.05):
	scores = []
	for i, marker in enumerate(markers):
		try:
			lda = LDA()
			lda.fit(marker["individuals"], marker["population_labels"])
			scores.append((lda.score(marker["individuals"], marker["population_labels"], i)))
		except:
			scores.append((0.0, i))
	scores.sort()
	scores.reverse()

	cutoff_idx = int(threshold * len(scores))

	return scores[:cutoff_idx]

def write_scores(scores, flname):
	fl = open(flname, "w")
	for loci, score in scores:
		fl.write("%s %s\n" % (loci, score))
	fl.close()

if __name__ == "__main__":
	variants_fl = sys.argv[1]
	scores_flname = sys.argv[2]

	variants = read_variants(variants_fl)

	scores = predict_scores(variants)
	write_scores(scores, scores_flname)

