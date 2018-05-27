#!/usr/bin/env bash

# Run imputation to find good guesses for the SNPs missing from 23andme

ALL_1000G_DIR=/Users/tom/projects-workspace/23andme-impute
TMP_DIR=impute-out
GEN_FILE=imputed-snps.gen

mkdir -p $OUT
rm -r $GEN_FILE

# rs689
impute2 -m $ALL_1000G_DIR/genetic_map_chr11_combined_b37.txt \
  -h $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr11_impute.hap.gz \
  -l $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr11_impute.legend.gz \
  -g $ALL_1000G_DIR/imputed.chr11.gen \
  -int 0e6 5e6 \
  -Ne 20000 \
  -o $TMP_DIR/tmp_impute2 \
  -phase \
  -allow_large_regions
grep 'rs689 ' $TMP_DIR/tmp_impute2 >> $GEN_FILE

# rs12722495
impute2 -m $ALL_1000G_DIR/genetic_map_chr10_combined_b37.txt \
  -h $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr10_impute.hap.gz \
  -l $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr10_impute.legend.gz \
  -g $ALL_1000G_DIR/imputed.chr10.gen \
  -int 5e6 10e6 \
  -Ne 20000 \
  -o $TMP_DIR/tmp_impute2 \
  -phase \
  -allow_large_regions
grep 'rs12722495 ' $TMP_DIR/tmp_impute2 >> $GEN_FILE

# rs17574546
impute2 -m $ALL_1000G_DIR/genetic_map_chr15_combined_b37.txt \
  -h $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr15_impute.hap.gz \
  -l $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr15_impute.legend.gz \
  -g $ALL_1000G_DIR/imputed.chr15.gen \
  -int 35e6 40e6 \
  -Ne 20000 \
  -o $TMP_DIR/tmp_impute2 \
  -phase \
  -allow_large_regions
grep 'rs17574546 ' $TMP_DIR/tmp_impute2 >> $GEN_FILE

# rs4948088
impute2 -m $ALL_1000G_DIR/genetic_map_chr7_combined_b37.txt \
  -h $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr7_impute.hap.gz \
  -l $ALL_1000G_DIR/ALL_1000G_phase1integrated_v3_chr7_impute.legend.gz \
  -g $ALL_1000G_DIR/imputed.chr7.gen \
  -int 50e6 55e6 \
  -Ne 20000 \
  -o $TMP_DIR/tmp_impute2 \
  -phase \
  -allow_large_regions
grep 'rs4948088 ' $TMP_DIR/tmp_impute2 >> $GEN_FILE