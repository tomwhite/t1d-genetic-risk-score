import csv
import os


def load_23andme(filename):
    variants = {}
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.split()
                variants[parts[0]] = {'rsid': parts[0], 'chromosome': parts[1], 'position': parts[2], 'genotype': parts[3]}
    return variants


def load_imputed(filename):
    variants = {}
    if os.path.isfile(filename):
        # format definition: http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html
        with open(filename, 'r') as file:
            for line in file:
                parts = line.split()
                rsid, a, b = parts[1], parts[3], parts[4]
                genotypes = (a + a, a + b, b + b)
                probs = (float(parts[5]), float(parts[6]), float(parts[7]))
                ind = probs.index(max(probs))
                variants[rsid] = {'genotype': genotypes[ind], 'info': probs[ind]}
                print(rsid, variants[rsid])
    return variants


def load_analysis(filename):
    grs_snps = {}
    with open(filename) as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # skip the headers
        for row in reader:
            grs_snps[row[0]] = {'snp': row[0], 'weight': row[3], 'effect_allele': row[4]}
    return grs_snps


def allele_counts(genotype, effect_allele, non_effect_allele=None):
    ac = genotype.count(effect_allele)
    nac = genotype.count(non_effect_allele) if non_effect_allele else 2 - ac
    if ac + nac != 2:
        reverse_genotype = genotype.replace('A', 't').replace('C', 'g').replace('G', 'c').replace('T', 'a').upper()
        ac = reverse_genotype.count(effect_allele)
        nac = reverse_genotype.count(non_effect_allele)
    return ac, nac
