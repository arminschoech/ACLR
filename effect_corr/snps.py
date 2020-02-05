'''
Stores SNP specific information like rsIDs, base pair, BGEN positions.

Provides sets of SNP indices to be processed next.

Name definitions:
    imp (imputation) SNPs: all imputed SNPs used in analysis
    reg (regression) SNPs: SNPs for which summary statistics are used/available
    lead SNPs: first SNP of sum. stat. pair used for each regression data point
    partner SNPs: second SNP of sum. stat. pair
'''

from __future__ import division
import numpy as np
import pandas as pd
import scipy.stats


def inv_norm_trans(x):
    '''returns rank-based inverse normal transform of input vector'''
    rank = x.argsort().argsort()
    unif_distr = (rank + 1) / float(rank.size + 1)
    norm_distr = scipy.stats.norm.ppf(unif_distr)
    return norm_distr


class SnpTable(object):

    def __init__(self, snp_file_path, lead_snp_num, max_reg_snp_dist,
                 max_imp_snp_dist, max_imp_snp_num=None,
                 max_imp_snp_num_pre_filter=None):
        self.lead_snp_num = lead_snp_num
        self.max_reg_snp_dist = max_reg_snp_dist
        self.max_imp_snp_dist = max_imp_snp_dist
        self.max_imp_snp_num = max_imp_snp_num
        self.next_lead_snp = 0
        self._tab = pd.read_csv(snp_file_path, sep='\t',
                                nrows=max_imp_snp_num_pre_filter)
        self._tab['maf'] = np.minimum(self._tab.freq, 1 - self._tab.freq)
        self._tab.loc[self._tab.rsid == '.', 'rsid'] = (self._tab
                                                        [self._tab == '.'].id)

    def filter_snps(self, maf_cutoff, subset_fac=1, rm_mult_allelic=True):
        '''filter SNPs according to MAF and potentially other metrics'''
        self._tab = self._tab[self._tab.maf > maf_cutoff]
        if rm_mult_allelic:
            # removing all but the highest freq. allele for mult. allelic sites
            self._tab = self._tab.sort_values('maf', ascending=False)
            self._tab = self._tab.drop_duplicates(subset='bp_pos')
            self._tab = self._tab.sort_values('bgen_pos')
            # TODO: find a better solution for multiallelic loci?
        subset = np.arange(0, self._tab.shape[0], subset_fac)
        self._tab = self._tab.iloc[subset, :]
        self._tab = self._tab[:self.max_imp_snp_num]
        self.imp_snp_num = self._tab.shape[0]
        if self.imp_snp_num < self.max_imp_snp_num:
            print('Warning: SNP number is smaller than the maximum number set')

    def read_in_regression_snps(self, sumstats_path, read_from_file=True):
        '''
        Reads regression SNPs from file or takes subset of imp SNPs.
        Adds reg_snp column indicating if it is a regression SNP.
        Add separate pandas table but only include reg SNPs
        '''
        if read_from_file:
            reg_snp_table = pd.read_csv(sumstats_path, sep='\t')
            self._tab = self._tab.merge(reg_snp_table[['SNP']], how='left',
                                        left_on='rsid', right_on='SNP',
                                        indicator='reg_snp')
            self._tab = self._tab.drop('SNP', axis=1)
        else:
            subset_frac = 1
            subset = np.arange(self._tab.shape[0] // subset_frac) * subset_frac
            reg_snp_table = self._tab.iloc[subset]
            self._tab = self._tab.merge(reg_snp_table[['rsid']], how='left',
                                        on='rsid', indicator='reg_snp')

        self._tab.reg_snp = (self._tab.reg_snp == 'both')
        self._reg_tab = self._tab[self._tab.reg_snp].copy()
        self._reg_tab = self._reg_tab.drop('reg_snp', axis=1)
        self._reg_tab = self._reg_tab.reset_index(drop=True)
        self.reg_snp_num = self._reg_tab.shape[0]

    def remove_mhc_snps(self):
        '''Remove SNP from the MHC region on chromosome 6. This region is known
        for frequent and strong non-linear interactions between SNP effects
        and is therefore left out from the analysis'''
        ini_reg_snp_num = self._reg_tab.shape[0]
        in_mhc = ((self._reg_tab.chrom == 6) &
                  (self._reg_tab.bp_pos > 25.5 * 10**6) &
                  (self._reg_tab.bp_pos < 33.5 * 10**6))
        self._reg_tab = self._reg_tab[~in_mhc]
        self._reg_tab = self._reg_tab.reset_index(drop=True)
        self.reg_snp_num = self._reg_tab.shape[0]
        snps_removed = ini_reg_snp_num - self.reg_snp_num
        if snps_removed != 0:
            print(str(snps_removed) +
                  ' regression SNPs were masked in the MHC region')

    def calc_next_analysis_snps(self):
        '''
        calculates lead, reg, and imp SNPs to be used in the analysis of the
        next SNP block and saves them in a pandas table each
        '''
        min_snp_num = self.next_lead_snp
        max_snp_num = self.next_lead_snp + self.lead_snp_num
        self.lead_snps = self._reg_tab[min_snp_num:max_snp_num]
        self.lead_snps = self.lead_snps.copy().reset_index(drop=True)

        max_bp_pos = self._reg_tab.bp_pos[self.next_lead_snp] - 1
        min_bp_pos = max_bp_pos - self.max_reg_snp_dist + 1
        partner_snps = (self._reg_tab
                        [self._reg_tab.bp_pos.between(min_bp_pos, max_bp_pos)])
        self.reg_snps = pd.concat((partner_snps, self.lead_snps), axis=0)
        self.reg_snps = self.reg_snps.copy().reset_index(drop=True)

        min_bp_pos = self.reg_snps.bp_pos.iloc[0] - self.max_imp_snp_dist
        max_bp_pos = self.reg_snps.bp_pos.iloc[-1] + self.max_imp_snp_dist
        self.imp_snps = (self._tab
                         [self._tab.bp_pos.between(min_bp_pos, max_bp_pos)])
        self.imp_snps = self.imp_snps.copy().reset_index(drop=True)
        self.next_lead_snp = max_snp_num  # Note: keeps track of next SNPs
        print('Number of imp SNPs: ' + str(self.imp_snps.shape[0]))
        print('Number of reg SNPs: ' + str(self.reg_snps.shape[0]))

    def set_distance_bins(self, distance_bins):
        '''sets boundaries of internal distance bins given input'''
        self._dist_bins_min = np.zeros(distance_bins.size - 1, dtype='int64')
        self._dist_bins_max = np.zeros(distance_bins.size - 1, dtype='int64')
        for i in xrange(distance_bins.size - 1):
            self._dist_bins_min[i] = distance_bins[i]
            self._dist_bins_max[i] = distance_bins[i + 1] - 1
            # min and max values themselves are included in the bin

    def set_ld_bins(self, ld_bins):
        '''sets boundaries of internal ld bins given input'''
        if ld_bins is None:
            return  # for non-LD stratified analyses
        if (ld_bins[0] != -1.0) or (ld_bins[-1] != 1.0):
            raise ValueError('LD bin bounds have to start at -1 and end at +1')
        self._ld_bins_min = np.zeros(ld_bins.size - 1)
        self._ld_bins_max = np.zeros(ld_bins.size - 1)
        for i in xrange(ld_bins.size - 1):
            self._ld_bins_min[i] = ld_bins[i]
            self._ld_bins_max[i] = ld_bins[i + 1]
        self._ld_bins_max[-1] = 1.1
        # min values are included max values not
        # last max is changed to 1.1 so that 1.0 is included

    def set_const_snp_dist(self, const_snp_dist):
        '''alters SNP base pair positions for them to be all equidistant;
        this is useful for some simulation analysis'''
        snp_num = self._tab.shape[0]
        self._tab.bp_pos = np.arange(0, const_snp_dist * snp_num,
                                     const_snp_dist)

    def calc_design_mats(self, ld_mat_approx=None):
        '''Calculates list of design matrices that indicate if SNP pair is a
        certain distance apart and hence belongs to specific distance bins'''
        pos = self.imp_snps.bp_pos.values
        dist = np.absolute(pos.reshape((-1, 1)) - pos)
        design_mats = []
        if self.annot_num > 1:
            for i in xrange(self.annot_num):
                design_mats.append(self.imp_snps['annotation_' + str(i)]
                                   .values.astype('float32'))
        else:
            design_mats.append(np.ones(self.imp_snps.shape[0],
                                       dtype='float32'))
        for i in xrange(1, self._dist_bins_min.size):
            if ld_mat_approx is None:
                design_mats.append(((dist >= self._dist_bins_min[i]) &
                                    (dist <= self._dist_bins_max[i]))
                                   .astype('float32'))
            else:
                dist_mat = ((dist >= self._dist_bins_min[i]) &
                            (dist <= self._dist_bins_max[i]))
                for j in xrange(self._ld_bins_min.size):
                    ld_design_mat = ((ld_mat_approx >= self._ld_bins_min[j]) &
                                     (ld_mat_approx < self._ld_bins_max[j]))
                    design_mats.append((dist_mat * ld_design_mat)
                                       .astype('float32'))
                    print('    Dist bin ' + str(i) + ', LD bin ' + str(j) +
                          ': ' + str((design_mats[-1] > 0.1).mean()))
        return design_mats

    def get_imp_snp_tab(self):
        '''get SNP table necessary for simulating phenotypes'''
        return self._tab

    def get_reg_snp_tab(self):
        '''get reg SNP table to calculate summary stats from trait values'''
        return self._reg_tab

    def calc_lld(self, ld_mat, quantile_num=10):
        '''Calculate level of LD (LLD) as defined in Gazal et al. 2017 Nat
        Genet; LD scores are binned by minor allele frequency and then rank-
        based inverse normal tansformed within each bin'''
        fac = np.diag(ld_mat)**(-0.5)
        ld_mat_norm = fac.reshape((-1, 1)) * ld_mat * fac.reshape((1, -1))
        self._tab['ld_score'] = np.diag(ld_mat_norm.T.dot(ld_mat_norm))
        maf = self._tab.maf
        snp_num = maf.size
        lld = np.ones(snp_num) * 99.0
        quantile_step = 0.5 / quantile_num
        for i in xrange(quantile_num):
            from_maf = i * quantile_step
            to_maf = (i + 1) * quantile_step
            include = (maf > from_maf) & (maf <= to_maf)
            lld[include] = inv_norm_trans(self._tab.ld_score[include])
        if (lld == 99.0).any():
            raise ValueError('LLD could not be calculated for all SNPs.')
        self._tab['lld'] = lld

    def add_snp_vars_from_file(self, snp_var_file):
        '''Merge external file that contains estimated effect variances for
        each SNP into the main SNP table'''
        var_df = pd.read_csv(snp_var_file)
        var_df = var_df.drop_duplicates(subset='BP')
        # TODO: find a better solution for multiallelic loci?
        var_df.snpvar = np.maximum(var_df.snpvar, 0.0)
        var_df = var_df[['BP', 'snpvar']]
        var_df = var_df.rename(columns={'BP': 'bp_pos', 'snpvar': 'snp_var'})
        self._tab = self._tab.merge(var_df, how='inner', on='bp_pos')
        self.imp_snp_num = self._tab.shape[0]
        self._tab.snp_var /= 2 * self._tab.freq * (1.0 - self._tab.freq)

    def calc_snp_var_weights(self):
        '''Calculate SNP weights from estimated effect variances'''
        self._tab['snp_var_weights'] = (self._tab.snp_var
                                        / self._tab.snp_var.mean())**0.5

    def get_snp_var_weights(self, all_snps=False):
        '''return SNP weights for all or only imputation SNPs'''
        if all_snps:
            weights = self._tab.snp_var_weights.values
        else:
            weights = self.imp_snps.snp_var_weights.values
        return weights

    def add_snp_var_annotations(self, annot_file):
        '''Read in SNP annotations from file and merge into main SNP file;
        this function only is used for simple test files, see full function
        'add_snp_var_annotations_baselineLF' below'''
        if annot_file is None:
            self.annot_num = 1
        else:
            annot_df = pd.read_csv(annot_file)
            annot_df = pd.get_dummies(annot_df, columns=['annotation'])
            self.annot_num = annot_df.shape[1] - 1
            print('Adding ' + str(self.annot_num) + ' functional groups.')
            self._tab = self._tab.merge(annot_df, how='left', on='bgen_pos')
            if self._tab.annotation_0.isnull().values.any():
                raise ValueError('Annotations not available for all SNPs.')

    def get_snp_var_from_annot(self, var_components, convert_to_per_allele):
        '''calculate the expected effect variance of SNPs given annotations'''
        last_col = 'annotation_' + str(self.annot_num - 1)
        annot_df = self._tab.loc[:, 'annotation_0': last_col]
        snp_vars = annot_df.values.dot(var_components)
        if convert_to_per_allele:  # convert per-std. to per-allele effect var.
            freq = self._tab.freq.values
            snp_vars /= 2 * freq * (1.0 - freq)
        return snp_vars

    def add_snp_var_annotations_baselineLF(self, annot_file):
        '''Read in SNP annotations from file and merge into main SNP file'''
        if annot_file is None:
            self.annot_num = 1
        else:
            annot_df = pd.read_csv(annot_file, sep='\t')
            col_subset = 'BP|MAFbin|Pred|LLD|Recomb|Nucleo|Selection|CpG'
            # col_subset = 'BP|MAFbin'
            annot_df = annot_df.filter(regex=col_subset, axis=1)
            self.annot_num = annot_df.shape[1] - 1
            col_names = ['bp_pos']
            for i in xrange(self.annot_num):
                col_names.append('annotation_' + str(i))
            annot_df.columns = col_names
            annot_df = annot_df.drop_duplicates('bp_pos')
            # TODO: find a more proper way to deal with multi-allelic SNPs
            print('Adding ' + str(self.annot_num) + ' functional groups.')
            self._tab = self._tab.merge(annot_df, how='left', on='bp_pos')
            if self._tab.annotation_0.isnull().values.any():
                missing_annot = self._tab.annotation_0.isnull()
                self._tab[missing_annot].to_csv('missing_annot.csv')
                print('Annotations missing for ' + str(missing_annot.sum()) +
                      ' out of ' + str(self._tab.shape[0]) + ' SNPs.')
                self._tab = self._tab[~missing_annot]
            # remove columns that do not vary, i.e. do not differentiate SNPs
            removed_annot = []
            for i in xrange(self.annot_num):
                col = self._tab['annotation_' + str(i)]
                if col.min() == col.max():
                    self._tab = self._tab.drop('annotation_' + str(i), axis=1)
                    removed_annot.append(i)
            print('Annotations ' + str(removed_annot) + ' have been removed' +
                  ' due to lack of variation within the annotations.')
            counter = 0
            for i in xrange(self.annot_num):
                if i not in removed_annot:
                    old = 'annotation_' + str(i)
                    new = 'annotation_' + str(counter)
                    self._tab = self._tab.rename(columns={old: new})
                    counter += 1
            self.annot_num -= len(removed_annot)

    def change_snp_bp_pos_to_blocks(self, block_size=3, bp_dist=100):
        '''changes the base pair positions of all SNPs such that blocks of SNPs
        are all closer than bp_dist but further away from all other SNPs'''
        if block_size > 9:
            raise ValueError('bp position blocks size > 9 are not supported')
        block_num = self.imp_snp_num // block_size + 1
        first_digit = np.mod(np.arange(block_num * block_size, dtype=int),
                             block_size)
        other_digits = (np.arange(block_num * block_size, dtype=int)
                        // block_size)
        # make sure that the overall distance stays roughly the same
        original_dist = self._tab.bp_pos.values[-1] - self._tab.bp_pos.iloc[0]
        factor = original_dist // other_digits[-1]
        other_digits *= factor
        bp_pos = first_digit + other_digits
        self._tab.loc[:, 'bp_pos'] = bp_pos[:self.imp_snp_num]
        self._reg_tab = self._reg_tab.drop('bp_pos', axis=1)
        self._reg_tab = self._reg_tab.merge(self._tab[['bgen_pos', 'bp_pos']],
                                            on='bgen_pos', how='left')

    def get_reg_snp_num(self):
        '''return number of regression SNPs'''
        return self.reg_snp_num
