'''
Calculates LD scores as from LD matrices and distance design matrices.
LD scores can then be used together with summary statistics to calculate
correlation estimates.
'''

import numpy as np
import pandas as pd
import scipy.sparse as sparse


def calc_ld_score_mats(ld_mat_imp, design_mats, sparse_calc=False):
    '''Calculates LD scores given distance design matrices; returns a list of
    MxM matrices of LD scores, one for each design matrix, where M is the SNP
    number'''
    ld_score_mats = []
    diag_mat_num = 0  # counts diagonal design matrices, have to first in list!
    for i in xrange(len(design_mats)):
        if len(design_mats[i].shape) == 2:  # if design mat i isn't a diag. mat
            if sparse_calc & (i != len(design_mats) - 1):
                ld_score_mats.append(ld_mat_imp.T
                                     .dot(sparse.csr_matrix(design_mats[i])
                                          .dot(ld_mat_imp)))
            elif sparse_calc and (i == len(design_mats) - 1):  # last matrix
                # since the last matrix is not sparse we use the fact that
                # the sum of all matrices is the 1 matrix
                col_sum_vec = ld_mat_imp.sum(axis=0)
                last_mat = (col_sum_vec.reshape((-1, 1)) *
                            col_sum_vec.reshape((1, -1)))
                # same as  ld_mat.T.dot(np.ones(.,.)).dot(ld_mat)
                last_mat -= sum(ld_score_mats[diag_mat_num:])
                # remove LD score matrices except from diagonal design mats
                last_mat -= ld_mat_imp.T.dot(ld_mat_imp)
                # remove RT.I.R
                ld_score_mats.append(last_mat)
            else:
                ld_score_mats.append(ld_mat_imp.T
                                     .dot(design_mats[i])
                                     .dot(ld_mat_imp))
        elif len(design_mats[i].shape) == 1:
            diag_mat_num += 1
            ld_score_mats.append((ld_mat_imp.T
                                  * design_mats[i].reshape((1, -1)))
                                 .dot(ld_mat_imp))
        else:
            raise ValueError('Design matrices have to be 1 or 2 dimensional.')
        print('        Finished LD score matrix ' + str(i))
    return ld_score_mats


def _calc_snp_pairs(lead_snps, reg_snps, max_reg_snp_dist):
    '''given set of lead SNPs and reg SNPs, calculates set of SNP pair indices
    that are not more than a maximum distance apart to be used in analysis'''
    lead_index_offset = reg_snps.shape[0] - lead_snps.shape[0]
    # = number of SNPs that are only partner SNPs and never lead SNPs
    snp_indices = []
    for i in xrange(lead_snps.shape[0]):
        max_pos = lead_snps.bp_pos[i]
        min_pos = lead_snps.bp_pos[i] - max_reg_snp_dist
        partners = reg_snps.index[reg_snps.bp_pos
                                  .between(min_pos, max_pos)].values
        leaders = np.full(partners.size, i, dtype=int) + lead_index_offset
        snp_indices.append(np.vstack((leaders, partners)).T)
    snp_indices = np.concatenate(snp_indices)
    lead_indices = snp_indices[:, 0]
    partner_indices = snp_indices[:, 1]
    return lead_indices, partner_indices


def calc_ld_score_table(ld_score_mats, ld_mat_reg, ld_mat_reg_norm,
                        max_reg_snp_dist, lead_snps, reg_snps, individual_num):
    '''Calculates pandas table of LD scores from LD score matrices and LD
    matrices together with SNP specific information to be used as regression
    table given summary statistics'''
    lead, partner = _calc_snp_pairs(lead_snps, reg_snps, max_reg_snp_dist)
    lead_table = reg_snps[['bgen_pos', 'rsid', 'bp_pos']].loc[lead].copy()
    partner_table = (reg_snps[['bgen_pos', 'rsid', 'bp_pos']]
                     .loc[partner].copy())
    lead_table.columns = lead_table.columns + '_1'
    partner_table.columns = partner_table.columns + '_2'
    ld_score_table = pd.concat((lead_table.reset_index(drop=True),
                                partner_table.reset_index(drop=True)), axis=1)

    for i in xrange(len(ld_score_mats)):
        ld_score_table['regressor_' + str(i)] = ld_score_mats[i][lead, partner]

    ld_score_table['environment'] = ld_mat_reg[lead, partner] / individual_num

    # additional terms for regression weighting
    ld_score_table['LD_1'] = ld_mat_reg[lead, lead] / individual_num
    # R_ii term in derviation
    ld_score_table['LD_2'] = ld_mat_reg[partner, partner] / individual_num
    # R_jj term in derivation
    for k in xrange(len(ld_score_mats)):
        ld_score_table['LD_score_1_' + str(k)] = (
            ld_score_mats[k][lead, lead])
    # R_i S_k R_i term in derviation
    for k in xrange(len(ld_score_mats)):
        ld_score_table['LD_score_2_' + str(k)] = (
            ld_score_mats[k][partner, partner])
    # R_j S_k R_j term in derviation

    # off diagonal weighting using normalized LD scores (i.e. overcounting):
    ld_score_table['off_diag_weight'] = (ld_mat_reg_norm.dot(ld_mat_reg_norm)
                                         [lead, partner])

    return ld_score_table
