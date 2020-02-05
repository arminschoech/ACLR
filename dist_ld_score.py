#!/usr/bin/env python

''' Tool to calculate distance dependent LD scores'''

import effect_corr.snps as snps
import effect_corr.genotypes as genotypes
import effect_corr.ld_scores as ld_scores
import numpy as np
import math
import datetime
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--bgen-file', type=str, required=True,
                    help='BGEN file with compressed genotype data')
parser.add_argument('--snp-file', type=str, required=True,
                    help='SNP table file with SNP related information')
parser.add_argument('--indi-file', type=str, required=True,
                    help='File with binary vector of numbers (int32) '
                    'indicating positions in the BGEN file for each individual'
                    ' used in analysis.')
parser.add_argument('--sum-stat-file', type=str, required=True,
                    help='Table of SNPs for which summary statistics are used;'
                    ' default is to use all SNPs in snp-file')
parser.add_argument('--annot-file', default=None, type=str,
                    help='Zipped functional annotation data for all SNPs; '
                    'default is no annotation data')
parser.add_argument('--batch-num-total', type=int, required=True,
                    help='Total number of batches the calculation is divided '
                    'into (enables parallelization)')
parser.add_argument('--batch-num', type=int, required=True,
                    help='Number of batch calculated in current run; has to '
                    'between 0 and (batch-num-total - 1)')
parser.add_argument('--out', default='ld_scores.csv', type=str,
                    help='Name of output csv file; default: "ld_scores.csv"')


if __name__ == '__main__':
    args = parser.parse_args()
    bgen_version = 1.2
    lead_snp_num = 1000  # number of SNPs read in in one batch
    max_reg_snp_dist = 0
    max_imp_snp_dist = 10**6
    max_imp_snp_num = None
    max_imp_snp_num_pre_filter = None
    maf_cutoff = 0.001
    subset_fac = 1
    distance_bins = np.array([0, 1, 100, 1000, 10**16])
    ld_bins = None  # np.array([-1.0, 0.3, 1.0])
    const_snp_dist = 0  # if non-zero, SNPs have constant bp distance
    reg_snps_from_file = True
    simulated_genotypes = False
    use_snp_var_weights = False
    sparse_calc = True
    change_snp_bp_pos_to_blocks = False
    bgen_path = args.bgen_file
    snp_file_path = args.snp_file
    sumstats_path = args.sum_stat_file
    annot_file = args.annot_file
    section_num = args.batch_num
    section_tot = args.batch_num_total
    out_file_path = args.out

    # writing individual file path to file for C++ subroutine
    indi_path_file = open('indi_path.txt', 'w')
    indi_path_file.write(args.indi_file)
    indi_path_file.close()
    individual_num = np.fromfile(args.indi_file, dtype='int32').size

    print('Reading in SNP information...')
    snp_table = snps.SnpTable(snp_file_path, lead_snp_num, max_reg_snp_dist,
                              max_imp_snp_dist, max_imp_snp_num,
                              max_imp_snp_num_pre_filter)
    # snp_table.add_snp_vars_from_file(variance_file)
    snp_table.filter_snps(maf_cutoff, subset_fac)
    # snp_table.calc_snp_var_weights()
    if annot_file == 'annotation.csv':
        snp_table.add_snp_var_annotations(annot_file)
    else:
        snp_table.add_snp_var_annotations_baselineLF(annot_file)
    if const_snp_dist:
        snp_table.set_const_snp_dist(const_snp_dist)
    snp_table.read_in_regression_snps(sumstats_path,
                                      read_from_file=reg_snps_from_file)
    snp_table.remove_mhc_snps()
    if change_snp_bp_pos_to_blocks:
        snp_table.change_snp_bp_pos_to_blocks()
    snp_table.set_distance_bins(distance_bins)
    snp_table.set_ld_bins(ld_bins)

    bgen_link = genotypes.BgenLink(bgen_path, bgen_version)

    print('Starting calculating LD scores...')
    
    # calculating which batches of SNPs to analyze in current section
    # (batches are sets of SNPs processed in one LD matrix)
    batch_num_tot = int(math.ceil(snp_table.reg_snp_num / float(lead_snp_num)))
    print('Total number of batches: ' + str(batch_num_tot))
    batch_num = int(math.ceil(batch_num_tot / float(section_tot)))
    print('Number of batches per section: ' + str(batch_num))
    lead_snps_per_section = batch_num * lead_snp_num
    first_lead_snp = section_num * lead_snps_per_section
    final_lead_snp = (section_num + 1) * lead_snps_per_section
    if section_num == section_tot - 1:  # last section might have less batches
        final_lead_snp = snp_table.reg_snp_num
    print('First lead of current section: ' + str(first_lead_snp))
    print('Final lead of current section: ' + str(final_lead_snp))

    out_file = open(out_file_path, 'a')
    out_file.truncate(0)
    first_iteration = True
    snp_table.next_lead_snp = first_lead_snp
    while snp_table.next_lead_snp < final_lead_snp:
        snp_table.calc_next_analysis_snps()
        print('Calculating LD scores up to SNP ' + str(snp_table.next_lead_snp)
              + ', at ' +
              datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        if simulated_genotypes:
            ld_mat, ld_mat_reg = bgen_link.load_ld_mat('geno.txt')
        else:
            ld_mat, ld_mat_reg = (bgen_link.calc_banded_ld_mats(
                snp_table.imp_snps, snp_table.reg_snps, max_imp_snp_dist))
        if ld_bins is not None:
            ld_mat_approx = bgen_link.calc_approx_ld_mat(snp_table.imp_snps)
            fac = np.diag(ld_mat_approx)**(-0.5)
            ld_mat_approx = (fac.reshape((-1, 1)) * ld_mat_approx
                             * fac.reshape((1, -1)))
        else:
            ld_mat_approx = None

        print(ld_mat_reg.shape)
        fac = np.diag(ld_mat_reg)**(-0.5)
        ld_mat_reg_norm = (fac.reshape((-1, 1)) * ld_mat_reg *
                           fac.reshape((1, -1)))
        print('    LD matrices ready at ' +
              datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        design_mats = snp_table.calc_design_mats(ld_mat_approx)
        if use_snp_var_weights:
            snp_weights = snp_table.get_snp_var_weights()
            if np.isnan(snp_weights).any():
                raise ValueError('Effect variance file does not include values'
                                 ' for all SNPs')
            for i in xrange(len(design_mats)):
                design_mats[i] = design_mats[i].astype(float)
                design_mats[i] *= snp_weights.reshape((-1, 1))
                design_mats[i] *= snp_weights.reshape((1, -1))
        print('    Design matrices read at ' +
              datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        ld_score_mats = ld_scores.calc_ld_score_mats(ld_mat, design_mats,
                                                     sparse_calc=sparse_calc)
        print('    LD score matrices read at ' +
              datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        reg_table = ld_scores.calc_ld_score_table(ld_score_mats, ld_mat_reg,
                                                  ld_mat_reg_norm,
                                                  max_reg_snp_dist,
                                                  snp_table.lead_snps,
                                                  snp_table.reg_snps,
                                                  individual_num)
        print('    Regression table read at ' +
              datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        if first_iteration:
            reg_table.to_csv(out_file, index=False)
        else:
            reg_table.to_csv(out_file, index=False, header=False)
        first_iteration = False
        print('Analysis of current batch finished at '
              + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    print('Analysis finished at '
          + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    out_file.close()
