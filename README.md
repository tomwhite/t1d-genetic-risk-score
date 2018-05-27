# Type 1 Diabetes Genetic Risk Score

Code for calculating an individual's T1D genetic risk score from their 23andme data.

A genetic risk score (GRS) is a number that quantifies an individual's risk to developing a
particular disease or condition - in this case T1D. It combines the influence from multiple
genes into a single number. To make sense of the number it needs to be compared with a
statistical distribution of a cohort who are known to have T1D or T2D.

The papers below do two things:
1. identify a set of SNPs for calculating the T1D GRS
2. quantify the distribution of GRS for T1D and T2D cohorts and provide criteria to distinguish between them

The code in this repository will calculate an individual's GRS from their genetic data from 23andme.

It is important to note that since we are dealing with probabilities, risk, and a propensity
to developing diabetes, the results of running the code will not give a definitive diagnosis.
If an individual has diabetes, the results _may_ help understand if T1D or T2D is more likely, but
it should never be a replacement for a professional doctor's opinion.

## Papers

[A Type 1 diabetes genetic risk score can aid discrimination between Type 1 and Type 2 diabetes in young adults][Oram], March 2016, Oram et al.

This paper identifies 30 SNPs to build a GRS for discriminating between T1D and T2D. Cohort size n=223,
comprised of white European individuals between 20 and 40.

[Frequency and phenotype of type 1 diabetes in the first six decades of life: a cross-sectional, genetically stratified survival analysis from UK Biobank][Thomas], February 2018, Thomas et al.

This paper uses the same SNPs from [Oram] (although one is excluded), but uses a much larger population size (n=379,511 individuals, white European,
under 60) to find a way of genetically determining T1D vs. T2D.

## Other Papers

[Feature ranking of type 1 diabetes susceptibility genes improves prediction of type 1 diabetes][Winkler], December 2014, Winkler et al

## Running the tool

### Find missing SNPs from 23andme using imputation

Download and install impute2.

```bash
export PATH=$PATH:~/Downloads/impute_v2.3.2_MacOSX_Intel/
```

[Oram]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5642867/
[Thomas]: https://www.thelancet.com/journals/landia/article/PIIS2213-8587%2817%2930362-5/fulltext?elsca1=tlxpr#sec1
[Winkler]: https://link.springer.com/article/10.1007%2Fs00125-014-3362-1