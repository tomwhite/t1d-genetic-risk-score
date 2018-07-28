#!/usr/bin/env bash

# Run imputation to find good guesses for the SNPs missing from 23andme

ALL_1000G_DIR=23andme-impute
TMP_DIR=impute-out
GEN_FILE=imputed-snps-5-types.gen

mkdir -p $OUT
rm -r $GEN_FILE

# rs11920090
impute2 -m $ALL_1000G_DIR/genetic_map_chr3_combined_b37.txt \
  -h $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr3_impute.hap.gz \
  -l $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr3_impute.legend.gz \
  -g $ALL_1000G_DIR/imputed.chr3.gen \
  -int 170e6 175e6 \
  -Ne 20000 \
  -o $TMP_DIR/tmp_impute2 \
  -phase \
  -allow_large_regions
grep 'rs11920090 ' $TMP_DIR/tmp_impute2 >> $GEN_FILE

# rs17271305
impute2 -m $ALL_1000G_DIR/genetic_map_chr15_combined_b37.txt \
  -h $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr15_impute.hap.gz \
  -l $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr15_impute.legend.gz \
  -g $ALL_1000G_DIR/imputed.chr15.gen \
  -int 60e6 65e6 \
  -Ne 20000 \
  -o $TMP_DIR/tmp_impute2 \
  -phase \
  -allow_large_regions
grep 'rs17271305 ' $TMP_DIR/tmp_impute2 >> $GEN_FILE

# rs10401969
impute2 -m $ALL_1000G_DIR/genetic_map_chr19_combined_b37.txt \
  -h $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr19_impute.hap.gz \
  -l $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr19_impute.legend.gz \
  -g $ALL_1000G_DIR/imputed.chr19.gen \
  -int 15e6 20e6 \
  -Ne 20000 \
  -o $TMP_DIR/tmp_impute2 \
  -phase \
  -allow_large_regions
grep 'rs10401969 ' $TMP_DIR/tmp_impute2 >> $GEN_FILE


