import csv
import math
import os
import sys

# Calculate the genetic risk score (GRS) for T1D from the paper
# Frequency and phenotype of type 1 diabetes in the first six decades of life: a cross-sectional, genetically stratified survival analysis from UK Biobank
# https://www.thelancet.com/journals/landia/article/PIIS2213-8587%2817%2930362-5/fulltext?elsca1=tlxpr#sec1

# The definition of GRS is from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5642867/
# 'A GRS is the sum across SNPs of the number of risk increasing alleles (0, 1 or 2) at that SNP multiplied by the
# ln(odds ratio) for each allele divided by the number of alleles. This assumes that each of the risk alleles has a
# log-additive effect on T1D risk.'

# Note that the HLA part of the score is calculated separately

if len(sys.argv) < 2:
    print("Usage: python t1d-grs.py <23andme-file>")
    sys.exit(1)

genome_23andme_file = sys.argv[1]

variants = {}

# add imputed data

if os.path.isfile("imputed-snps-5-types.gen"):
    # format definition: http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html
    with open("imputed-snps-5-types.gen", 'r') as file:
        for line in file:
            parts = line.split()
            rsid, a, b = parts[1], parts[3], parts[4]
            genotypes = (a + a, a + b, b + b)
            probs = (float(parts[5]), float(parts[6]), float(parts[7]))
            ind = probs.index(max(probs))
            variants[rsid] = {'genotype': genotypes[ind], 'info': probs[ind]}
            #print(rsid, variants[rsid])

# load the 23andme data
with open(genome_23andme_file, 'r') as file:
    for line in file:
        if not line.startswith('#'):
            parts = line.split()
            variants[parts[0]] = {'rsid': parts[0], 'chromosome': parts[1], 'position': parts[2], 'genotype': parts[3]}


def allele_counts(genotype, effect_allele, non_effect_allele):
    ac = genotype.count(effect_allele)
    nac = genotype.count(non_effect_allele)
    if ac + nac != 2:
        reverse_genotype = genotype.replace('A', 't').replace('C', 'g').replace('G', 'c').replace('T', 'a').upper()
        ac = reverse_genotype.count(effect_allele)
        nac = reverse_genotype.count(non_effect_allele)
    return ac, nac


def calculate_grs(snps_file):
    # Read in the snps to use for the score and store in a dict
    grs_snps = {}
    with open(snps_file) as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # skip the headers
        for row in reader:
            grs_snps[row[0]] = {'snp': row[0], 'odds_ratio': row[2], 'weight': row[3], 'effect_allele': row[4], 'non_effect_allele': row[5]}

    total_snps = 0
    total_snps_used = 0
    genetic_risk_score = 0.0
    max_possible_grs = 0.0
    total_odds_ratio = 1.0

    for rsid, snp_info in grs_snps.iteritems():
        total_snps += 1
        variant = variants.get(rsid, None)
        odds_ratio = float(snp_info['odds_ratio'])
        weight = math.log(odds_ratio)
        if weight >= 0:
            max_possible_grs += weight * 2
        else:
            max_possible_grs += -weight * 2
        if variant:
            allele_count, non_allele_count = allele_counts(variant['genotype'], snp_info['effect_allele'], snp_info['non_effect_allele'])
            total_snps_used += 1
            if weight >= 0:
                score = weight * allele_count
            else:
                score = -weight * non_allele_count
            info = variant.get('info', 1.0)
            genetic_risk_score += score * info
            for _ in range(allele_count):
                total_odds_ratio = total_odds_ratio * odds_ratio
            #print(rsid, snp_info, variant, allele_count, non_allele_count, score)
        else:
            print("%s not found in variant data" % rsid)

    if total_snps_used < total_snps:
        print("There were %s missing SNPs. Please run imputation (find-missing-snps.sh)." % (total_snps - total_snps_used))
        sys.exit(1)

    return total_snps, total_snps_used, genetic_risk_score, max_possible_grs, total_odds_ratio


for snps_file in ('grs-5types-said.csv', 'grs-5types-sidd.csv', 'grs-5types-sird.csv', 'grs-5types-mod.csv', 'grs-5types-mard.csv'):
    total_snps, total_snps_used, genetic_risk_score, max_possible_grs, total_odds_ratio = calculate_grs(snps_file)
    print(snps_file)
    print("Total SNPs: %s" % total_snps)
    print("Total SNPs used: %s" % total_snps_used)
    print("Genetic risk score: %.3f" % genetic_risk_score)
    print("Max possible genetic risk score: %.3f" % max_possible_grs)
    print("Genetic risk score (percent): %.1f" % (100.0 * genetic_risk_score / max_possible_grs))
    print("Total odds ratio: %.2f" % total_odds_ratio)
    print("")

