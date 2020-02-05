'''
Stores link to the BGEN file used.

Selects and stores subset of individuals used in the analysis.

Returns dosage genotypes, frequencies, imput. variances for given set of SNPs.
'''

import numpy as np
import subprocess


class BgenLink(object):

    def __init__(self, bgen_path, bgen_version, max_individual_num=None):
        self.bgen_path = bgen_path
        self.bgen_version = bgen_version
        bgen_path_file = open('bgen_path.txt', 'w')
        bgen_path_file.write(self.bgen_path)
        bgen_path_file.write('\n')
        bgen_path_file.close()

    def calc_banded_ld_mats(self, imp_snp_tab, reg_snp_tab, max_distance):
        '''calculated banded LD matrices, both a imp SNP x reg SNP matrix and a
        reg SNP x reg SNP matrix, using a fast C++ subroutine to read in
        genotypes from BGEN file'''
        imp_pos = imp_snp_tab.bgen_pos.values
        imp_file = open('imp_snps.bin', 'w')
        imp_pos.astype('int32').tofile(imp_file)
        imp_file.close()

        reg_pos = reg_snp_tab.bgen_pos.values
        reg_file = open('reg_snps.bin', 'w')
        reg_pos.astype('int32').tofile(reg_file)
        reg_file.close()

        if self.bgen_version == 1.1:
            cpp_process = subprocess.call('./ld_mat_calc', shell=True)
        elif self.bgen_version == 1.2:
            cpp_process = subprocess.call('./ld_mat_calc_500K', shell=True)
        if cpp_process != 0:
            raise ValueError('C++ subprocess did not finish correctly.')

        ld_mat = np.fromfile('ld_mat.bin', dtype='float32')
        ld_mat = ld_mat.reshape((ld_mat.size // reg_pos.size, reg_pos.size))
        ld_mat = ld_mat[:imp_pos.size, :]
        # subprocess might return too long matrices with last rows invalid

        banding_mat = np.absolute(imp_snp_tab.bp_pos.values.reshape((-1, 1))
                                  - reg_snp_tab.bp_pos.values)
        banding_mat = (banding_mat < max_distance).astype('int')
        ld_mat = ld_mat * banding_mat

        relative_pos = np.argwhere(np.isin(imp_pos, reg_pos))[:, 0]
        # position of reg SNPs within the vector of imp SNPs
        ld_mat_reg = ld_mat[relative_pos, :]
        return ld_mat, ld_mat_reg

    def load_ld_mat(self, file_path):
        '''calculates LD matrices from genotypes loaded from a text file'''
        gen_mat = np.loadtxt(file_path)
        ld_mat = gen_mat.T.dot(gen_mat) / gen_mat.shape[0]
        return ld_mat, ld_mat

    def calc_approx_ld_mat(self, imp_snp_tab):
        '''calculated imp SNP x imp SNP LD matrix, using a fast C++ subroutine
        to read in genotypes from BGEN file; usually only a subset of
        individuals is used, so LD matrix only approximates the sample LD'''
        imp_pos = imp_snp_tab.bgen_pos.values
        imp_file = open('imp_snps.bin', 'w')
        imp_pos.astype('int32').tofile(imp_file)
        imp_file.close()

        if self.bgen_version == 1.1:
            raise ValueError('Approx. LD calc. not supported in BGEN v1.1')
        elif self.bgen_version == 1.2:
            cpp_process = subprocess.call('./approxld_500K', shell=True)
        if cpp_process != 0:
            raise ValueError('C++ subprocess did not finish correctly.')

        ld_mat = np.fromfile('ld_mat.bin', dtype='float32')
        ld_mat = ld_mat.reshape((ld_mat.size // imp_pos.size, imp_pos.size))
        ld_mat = ld_mat[:imp_pos.size, :]
        # subprocess might return too long matrices with last rows invalid

        return ld_mat
