#!/usr/bin/env python

'''
Regress products of trait summary statistics on calculated LD scores to obtain
estimates of the covariance between neighboring SNPs
'''

import effect_corr.regression as regression
import numpy as np
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--ld-score-file', type=str, default='ld_scores.csv',
                    help='File with distance dependent LD scores')
parser.add_argument('--sum-stat-file', type=str, required=True,
                    help='File containing summary statistics, i.e. marginal '
                    'SNP effect estimates')
parser.add_argument('--annot-num', type=int, default=1,
                    help='Number of functional annotations used')
parser.add_argument('--base-var-annot', type=int, default=1,
                    help='Number of base variance annotations')
parser.add_argument('--prior-env-var', type=float, required=True,
                    help='Prior environmental variance relative to the SNP '
                    'effect variance')
parser.add_argument('--out-file', type=str, default='result.txt',
                    help='Output file with variance and covariance estimates')
parser.add_argument('--err-file', type=str, default='error.txt',
                    help='File with block jackknife errors of estimates')
parser.add_argument('--jk-block-num', type=int, default=200,
                    help='Number of blocks used in block Jackknife error '
                    'estimation')
parser.add_argument('--trait-num', type=int, required=True,
                    help='Number of traits to be analyzed')


if __name__ == '__main__':
    args = parser.parse_args()
    second_weighting = False
    no_corr_effects = False
    perform_jk = True
    single_chrom = None
    comb_annot = False
    dist_bin_num = 3
    max_reg_snp_dist = 0
    annot_bins = args.annot_num
    base_var_bins = args.base_var_annot
    out_file = args.out_file
    err_file = args.err_file

    reg_tab = pd.read_csv(args.ld_score_file)
    if max_reg_snp_dist is not None:
        reg_snp_dist = np.absolute(reg_tab.bp_pos_1 - reg_tab.bp_pos_2)
        reg_tab = reg_tab[reg_snp_dist <= max_reg_snp_dist]
    if single_chrom is not None:
        reg_tab = reg_tab[reg_tab.bgen_pos_1 // 10**7 == single_chrom]
    if comb_annot:
        reg_tab = regression.combine_funct_annot(reg_tab, annot_bins,
                                                 base_var_bins, annot_bins + 6)
        annot_bins = 1
    sum_stat_tab = pd.read_csv(args.sum_stat_file)
    snp_num = sum_stat_tab.shape[0]
    reg_row_num = reg_tab.shape[0]
    prior_est = np.concatenate((np.ones(args.base_var_annot),
                                np.zeros(args.annot_num - args.base_var_annot),
                                np.zeros(dist_bin_num),
                                np.array([args.prior_env_var])))
    reg_num = reg_tab.filter(regex='regressor', axis=1).shape[1]
    if prior_est.size - 1 != reg_num:
        raise ValueError('Number or regressors (' + str(reg_num) +
                         ') differs from number of prior variances (' +
                         str(prior_est.size - 1) + ').')
    if no_corr_effects:
        prior_est = np.array([prior_est[0], prior_est[-1]])
        out_file = out_file[:-4] + '_no_corr.txt'
        err_file = err_file[:-4] + '_no_corr.txt'

    print('Trait number: ' + str(args.trait_num))

    if no_corr_effects:
        reg_tab = regression.remove_corr_regressors(reg_tab)

    if perform_jk:
        coeffs, coeff_var = regression.jackknife_analysis(reg_tab,
                                                          sum_stat_tab,
                                                          prior_est,
                                                          args.trait_num,
                                                          args.jk_block_num)
    else:
        coeffs = regression.multi_trait_analysis(reg_tab, sum_stat_tab,
                                                 prior_est, args.trait_num,
                                                 second_weighting)
    np.savetxt(out_file, coeffs)
    print('Mean estimate:      ' + str(coeffs.mean(axis=0)[annot_bins:]))
    print('Error in mean:      ' + str(coeffs.std(axis=0)[annot_bins:]
                                       / args.trait_num**0.5))
    print('Std. of estimates:  ' + str(coeffs.std(axis=0)[annot_bins:]))
    if False:
        print('Mean std. estimate: ' + str((coeff_var**0.5).mean(axis=0)))
        print('Std. in std. est.:  ' + str((coeff_var**0.5).std(axis=0)))
        print('Err in mean std est:' +
              str((coeff_var**0.5).std(axis=0) / args.trait_num**0.5))
        np.savetxt(err_file, coeff_var**0.5)
    print('Analysis completed.')
