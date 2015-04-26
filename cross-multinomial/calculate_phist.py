import sys
import numpy as np
import matplotlib.pyplot as plt

def ssdwpfunc(individuals, frequencies):
    """
    Returns the sums squares deviation within populations from the population frequencies.

    individuals[pop] = counts
    frequencies[pop, haplotype] = freq
    """
    ssdwp = 0.0

    n_pop = frequencies.shape[0]
    n_haplotypes = frequencies.shape[1]    

    for pop_idx in xrange(n_pop):
        gene_copies = individuals[pop_idx] * 2 # diploid
        pop_freq = frequencies[pop_idx, :]
        pop_ssd = (np.outer(pop_freq, pop_freq).sum() - np.inner(pop_freq, pop_freq)) / 2.0
        ssdwp += gene_copies * pop_ssd

    return ssdwp

def ssdtotalfunc(individuals, frequencies):
    """
    Calculates the total sum squared deviation for a locus.

    individuals[pop] = counts
    frequencies[pop, haplotype] = freq
    """
    
    # total number of genes across all populations for a given locus
    locus_gene_copies = 2.0 * individuals.sum() # diploid
    total_freq = np.sum(frequencies * individuals[:, np.newaxis], axis=0) / individuals.sum()

    ssd = locus_gene_copies * (np.outer(total_freq, total_freq).sum() - np.inner(total_freq, total_freq)) / 2.0

    return ssd

def onecalcphi(individuals, frequencies):
    """
    individuals[pop] = individuals
    frequencies[pop][haplotype] = individuals
    """
    n_gene_copies = individuals.sum() * 2 # diploid
    n_pop = individuals.shape[0]

    # calculate the sums squared deviation within populations
    ssdwp = ssdwpfunc(individuals, frequencies)

    # sums squared deviations total at the locus
    ssdtotal = ssdtotalfunc(individuals, frequencies)

    # degrees of freedom for between populations
    dfb = n_pop - 1

    # degrees of freedom for total locus
    dfw = n_gene_copies - n_pop

    if dfw == 0:
        return 0.0

    # mean squared deviation within populations
    msdwp = ssdwp / dfw

    # mean squared deviation among populations
    msdap = (ssdtotal - ssdwp)/dfb

    # Calculate the variation among populations
    varAP = (msdap - msdwp)/(float(n_gene_copies)/n_pop)

    if (varAP + msdwp) == 0.0:
        return 0.0

    # PHIst is the proportion of the variation partitioned among populations
    phi_st = varAP/(msdwp + varAP)
    
    assert not(np.isnan(phi_st))

    return phi_st

def calc_locus_phi_st(individuals, frequencies):
    n_loci = individuals.shape[0]
    phi_st = np.zeros(n_loci)
    for locus_i in xrange(n_loci):
        phi_st[locus_i] = onecalcphi(individuals[locus_i, :], frequencies[locus_i, :, :])
    return phi_st

def read_counts(flname):
    fl = open(flname)
    vector = []
    populations = []
    for ln in fl:
        if "Marker" in ln:
            if len(populations) != 0:
                vector.append(populations)
                populations = []
            continue
        cols = ln.split()
        pop_locus_counts = map(float, cols[2:])
        populations.append(pop_locus_counts)
    vector.append(populations)
    fl.close()

    return np.array(vector)

def normalize_haplotypes(counts):
    # sum total haplotype counts for each
    # population-locus combination
    total = np.sum(counts, axis=2)
    frequencies = counts / total[:, :, None]
    return frequencies


def write_phi(flname, phi_values):
    fl = open(flname, "w")
    for i, phi in enumerate(phi_values):
        fl.write("%s,%s\n" % (i, phi))
    fl.close()


if __name__ == "__main__":
    occur_fl = sys.argv[1]
    out_fl = sys.argv[2]

    counts = read_counts(occur_fl)
    frequencies = normalize_haplotypes(counts)
    individuals = counts.sum(axis=2)
    phi_sts = calc_locus_phi_st(individuals, frequencies)
    write_phi(out_fl, phi_sts)

