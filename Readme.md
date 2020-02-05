# ACLR
Autocorrelation LD regression: a tool to efficiently estimate the autocorrelation of latent effects in large genetic data sets


## Introduction

ACLR is a software tool to infer the autocorrelation of latent genetic effects on a human trait or disease. Specifically, it uses data from large cohorts that contain genome-wide genetic variant data as well as trait or disease values for each individual in the cohort. It then returns correlation estimates of the trait effects of these genetic variants depending on the distance from each other along the genome, i.e. the autocorrelation as a function of genomic distance. As is typical for such data sets, the method is primarily aimed at applications where the number of genetic effects (typically tens of millions) is much larger than the number of individuals (typically hundreds of thousands). 

Many previous genetic analysis tools have explicitly or implicitly assumed that effects of nearby genetic variants are uncorrelated. This tool was designed to test and potentially challenge this assumption, which has important implications for finding causal genes in human disease as well as improve disease prediction from genetic data. 


## Summary of statistical method and implementation

Statistically, the method corresponds to estimating the autocorrelation of latent effects in a linear or probit model with the number of features being much larger than the number of data points. The tool regresses functions of marginal effect estimates (estimates from single feature regressions) onto so-called distance-dependent LD scores. These distance-dependent LD scores are defined by the matrix RSR, where R = X<sup>T</sup>X with X the design matrix of the linear/probit model, S a matrix indicating which genetic feature pairs are within the distance of interest on the genome, and R the feature-by-feature sample covariance matrix (in genetics referred to as LD matrix). 

Since the number of features are typically in the tens of millions, calculating and using the full feature-by-feature sample covariance matrix R is not possible. However, since covariances of distant features are close to zero, the tool implements a tailored banded matrix-based approach. Since S matrices are usually sparse, computation time is further reduced by using sparse matrix calculations. Also, while the majority of the matrix operations are called using Python, reading in data from compressed genotype files and building initial data matrices is performed using a specifically developed C++ subroutine, that is part of this tool. This, together with the code being fully parallelizable, makes it possible to run the method on large genetic data sets. 

## Requirements and installation
The Python 2.7 based code only uses common libraries such as Numpy, Pandas, scipy, and scikit-learn. C++ library requirements: Gnu Science Library, Intel Math Kernel Library, zlib, C standard library, POSIX threads, all of which are freely available. Before compiling the code, the locations of these libraries have to be updated in the makefile to the appropriate location. Then use ‘make’ do compile the executable ‘ld_mat_calc’. 

## Usage
First, use "dist_ld_score.py" to calculate distance dependent LD scores based on the provided genotype data and list of genetic variants to be used. Then use "regression.py" to combine these scores with marginal effect estimates of the target trait(s) to calculate effect variance and distance dependent covariance effects as well as block-Jackknife error estimates.


## Command line arguments

#### Required arguments for "dist_ld_score.py":

   --bgen-file: path to compressed genotype data file in BGEN v1.2 format (see https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html)

   --snp-file: table of genetic variants (SNPs) to be used; required columns: "bgen_pos", position in the BGEN file; "bp_pos", base pair position in the genome; "chrom", chromosome number; "rsid", SNP ID; "freq", the SNP population frequency; "alleleA" and "alleleB", major and minor allele of the SNP

   --indi-file: binary file containing a vector of positions (32 bit signed integer) in the BGEN file for each individual to be used in the analysis
   
   --sum-stat-file: marginal effects (summary statistics) file in LDSC format (https://github.com/bulik/ldsc/wiki)
   
 
#### Optional arguments for "dist_ld_score.py": 
   
   --out: set name of output csv files; default is "ld_scores.csv"
   
   --batch-num and --batch-num-total: only process input data from a subset of SNPs for parallelization; e.g. "--batch-num 2 --batch-num-total 6" means that only the second of 6 equally sized subsets of all SNPs get processed
   
   --annot-file: zipped functional annotation data file to be used for advanced analyses (please contact author for specific instructions)
   
#### Required arguments for "regression.py":
   
   --ld-score-file: output file from "dist_ld_score.py"; default "ld_scores.csv"
   
   --sum-stat-file: marginal effect (summary statistics) file in LDSC format (https://github.com/bulik/ldsc/wiki)
   
   --prior-env-var: prior variance of environmental noise compared to the variance of a single SNP effect
   
   --trait-num: number of traits used
   
   
#### Optional arguments for "regression.py":

   --ld-score-file: output file from "dist_ld_score.py"; default "ld_scores.csv"
   
   --out-file and --err-file: name of output files for covariance estimates and error estimates respectively; default is "result.txt" and "error.txt"
   
   --jk-block-num: number of Jackknife blocks used for error estimation; default is 200

   --annot-num and --base-var-annot: number functional annotations and base variance annotations to be used for advanced analyses (please contact author for specific instructions)


## Contact
Author: Armin Schoech. Please email arminschoech@g.harvard.edu for comments and questions. 

## License
The ACLR software tool is free under the GNU General Public License v3.0 (GPLv3). 
