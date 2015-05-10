import sys

import numpy as np

def predict_loci(phi_values, cutoff_percent):
	average_phi_values = np.mean(axis=1)

	sortable = []
	for loci, phi in enumerate(average_phi_values):
		sortable.append((phi, loci))
	sortable.sort()
	sortable.reverse()

	cutoff_idx = int(len(sortable) * cutoff_percent)

	return sortable[:cutoff_idx]

def write_predicted(flname, predicted_loci):
	fl = open(flname, "w")
	for phi, loci in predicted_loci:
		fl.write("%s %s\n" % (loci, phi))
	fl.close()

if __name__ == "__main__":
	npy_flname = sys.argv[1]
	cutoff_percent = float(sys.argv[2])
	predicted_flname = sys.argv[3]

	phi_values = np.load(npy_flname)
	predicted_loci = predict_loci(phi_values)
	write_predicted(predicted_flname, predicted_loci)
	

	