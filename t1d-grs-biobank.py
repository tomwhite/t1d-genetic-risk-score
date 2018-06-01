import csv
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
    print("Usage: python t1d-grs-biobank.py <23andme-file>")
    sys.exit(1)

genome_23andme_file = sys.argv[1]

variants = {}

# add imputed data

if os.path.isfile("imputed-snps-biobank.gen"):
    # format definition: http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html
    with open("imputed-snps-biobank.gen", 'r') as file:
        for line in file:
            parts = line.split()
            rsid, a, b = parts[1], parts[3], parts[4]
            genotypes = (a + a, a + b, b + b)
            probs = (float(parts[5]), float(parts[6]), float(parts[7]))
            ind = probs.index(max(probs))
            variants[rsid] = {'genotype': genotypes[ind], 'info': probs[ind]}
            print(rsid, variants[rsid])

# load the 23andme data
with open(genome_23andme_file, 'r') as file:
    for line in file:
        if not line.startswith('#'):
            parts = line.split()
            variants[parts[0]] = {'rsid': parts[0], 'chromosome': parts[1], 'position': parts[2], 'genotype': parts[3]}

# Read in the snps to use for the score and store in a dict
grs_snps = {}
with open('analyses/grs-biobank.csv') as csvfile:
    reader = csv.reader(csvfile)
    next(reader, None)  # skip the headers
    for row in reader:
        grs_snps[row[0]] = {'snp': row[0], 'weight': row[3], 'effect_allele': row[4]}

del grs_snps['rs4948088'] # excluded from Thomas paper since it is out of Hardy-Weinberg equilibrium, see p124

total_snps = 0
total_snps_used = 0
genetic_risk_score = 0.0

# HLA
total_snps += 2

rs2187668 = variants.get('rs2187668', None)
rs2187668_geno = rs2187668['genotype'] if rs2187668 else None
rs7454108 = variants.get('rs7454108', None)
rs7454108_geno = rs7454108['genotype'] if rs7454108 else None

# See Figure 4.A at https://digitalsauna.wordpress.com/2017/12/02/diabetes-and-me/
if rs2187668 and rs7454108:
    total_snps_used += 2
    if rs2187668_geno == 'AG' and rs7454108_geno == 'CT':
        print("DR3/DR4")
        genetic_risk_score += 3.87
    elif rs2187668_geno == 'AA':
        print("DR3/DR3")
        genetic_risk_score += 3.05
    elif rs7454108_geno == 'CC':
        print("DR4/DR4")
        genetic_risk_score += 3.09
    elif rs7454108_geno == 'CT':
        print("DR4/X")
        genetic_risk_score += 1.95
    elif rs2187668_geno == 'AG':
        print("DR3/X")
        genetic_risk_score += 1.51
    else:
        print("DRX/DRX (no risk alleles)")
else:
    if not rs2187668:
        print("rs2187668 not found in variant data")
    if not rs7454108:
        print("rs7454108 not found in variant data")

# Non HLA
for rsid, snp_info in grs_snps.iteritems():
    total_snps += 1
    variant = variants.get(rsid, None)
    if variant:
        allele_count = variant['genotype'].count(snp_info['effect_allele'])
        total_snps_used += 1
        genetic_risk_score += (float(snp_info['weight']) * allele_count)
        print(rsid, snp_info, variant, allele_count, (float(snp_info['weight']) * allele_count))
    else:
        print("%s not found in variant data" % rsid)

if total_snps_used < total_snps:
    print("There were %s missing SNPs. Please run imputation (find-missing-snps.sh)." % (total_snps - total_snps_used))
    sys.exit(1)

total_alleles = total_snps * 2
genetic_risk_score = genetic_risk_score / float(total_alleles)

print("Total SNPs: %s" % total_snps)
print("Total SNPs used: %s" % total_snps_used)
print("Genetic risk score: %.3f" % genetic_risk_score)

if genetic_risk_score > 0.231:
    print("Consistent with T1D")
else:
    print("Unlikely to be T1D")
