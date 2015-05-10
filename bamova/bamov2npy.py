import sys
import numpy as np

def read_phi(flname, n_steps, n_loci):
	sampled_phis = np.zeros((n_steps, n_loci))
	fl = open(flname)
	current_iter_idx = 0 # index used for storage
	last_iter_idx = 0 # index used to identify when we finish a step
	for ln in fl:
		cols = ln.strip().split(",")
		iter_idx = int(cols[0])
		locus_idx = int(cols[1])
		phi = float(cols[2])
		if last_iter_idx != iter_idx:
			last_iter_idx = iter_idx
			current_iter_idx += 1
		sampled_phis[current_iter_idx, locus_idx] = phi
	fl.close()

	return sampled_phis


if __name__ == "__main__":
	bamova_phi_output_flname = sys.argv[1]
	n_steps = int(sys.argv[2])
	n_loci = int(sys.argv[3])
	npy_flname = sys.argv[4]

	matrix = read_phi(bamova_phi_output_flname, n_steps, n_loci)

	np.save(matrix, npy_flname)