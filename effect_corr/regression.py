'''
Regress products of trait summary statistics on distance-dependent LD scores to
obtain estimates of the covariance between neighboring SNPs
'''

import numpy as np
import pandas as pd
import sklearn.linear_model


def single_trait_analysis(reg_tab, sum_stat_tab, prior_est, i):
    '''Regressing summary statistics of a single traits onto LD scores'''
    col_names = sum_stat_tab.columns
    sum_stat_tab.columns = col_names + '_1'
    reg_tab = reg_tab.merge(sum_stat_tab, on='bgen_pos_1', how='left')
    sum_stat_tab.columns = col_names + '_2'
    reg_tab = reg_tab.merge(sum_stat_tab, on='bgen_pos_2', how='left')
    reg_tab['sum_stat_prod'] = (reg_tab.sum_stats_1 * reg_tab.sum_stats_2)

    if prior_est is not None:
        reg_weights = calc_reg_weights_fixed(reg_tab, prior_est)
        dist_bin_num = prior_est.size - 1
        for i in xrange(dist_bin_num):
            reg_tab['regressor_' + str(i)] *= reg_weights
        reg_tab['environment'] *= reg_weights
        reg_tab['sum_stat_prod'] *= reg_weights

    regression = sklearn.linear_model.LinearRegression(fit_intercept=False)
    regression.fit(reg_tab.loc[:, 'regressor_0':'environment'],
                   reg_tab.sum_stat_prod)
    reg_tab = reg_tab.drop(['sum_stats_1', 'sum_stats_2', 'sum_stat_prod'],
                           axis=1)
    return regression.coef_


def calc_true_mean(reg_tab, true_sigma):
    '''Calculate the expectation of each row in the regression, in other words
    X*beta, to check if this matches simulation results'''
    reg_row_num = reg_tab.shape[0]
    true_mean = np.zeros(reg_row_num)
    for k in xrange(true_sigma.size - 1):
        true_mean += true_sigma[k] * reg_tab['regressor_' + str(k)]
    true_mean += true_sigma[-1] * reg_tab['environment']
    return true_mean.values


def calc_reg_weights_off_diag(reg_tab):
    '''Calculate additional regression weights for overcounting (apart from
    heteroscedastic weights) to improve power; currently not used'''
    return 1.0  # / reg_tab.off_diag_weight.values


def calc_reg_weights_fixed(reg_tab, prior_est):
    ''' calculate regression weights with fixed prior estimates for all
    traits'''
    prior_est = np.maximum(0, prior_est)
    # TODO: think more about what to do with negative estimates;
    # as they can lead to negative variance estimates which doesnt work
    dist_bin_num = prior_est.size - 1  # last entry is noise term

    ld_score_1 = np.zeros(reg_tab.shape[0])  # R_i Sigma R_i
    for k in xrange(dist_bin_num):
        ld_score_1 += reg_tab['LD_score_1_' + str(k)] * prior_est[k]
    ld_score_1 += reg_tab['LD_1'] * prior_est[-1]

    ld_score_2 = np.zeros(reg_tab.shape[0])  # R_j Sigma R_j
    for k in xrange(dist_bin_num):
        ld_score_2 += reg_tab['LD_score_2_' + str(k)] * prior_est[k]
    ld_score_2 += reg_tab['LD_2'] * prior_est[-1]

    ld_score = np.zeros(reg_tab.shape[0])  # R_i Sigma R_j
    for k in xrange(dist_bin_num):
        ld_score += reg_tab['regressor_' + str(k)] * prior_est[k]
    ld_score += reg_tab['environment'] * prior_est[-1]

    heterosced_weights = (ld_score_1 * ld_score_2 + ld_score**2).values**(-0.5)
    off_diag_weights = calc_reg_weights_off_diag(reg_tab)
    return heterosced_weights * off_diag_weights


def remove_corr_regressors(reg_tab):
    '''Remove regressors that estimate SNP effect covariances to estimate the
    effect variances only'''
    reached_regressors = False
    col_names = reg_tab.columns.values
    for col in col_names:
        if col == 'environment':
            break
        if reached_regressors:
            reg_tab = reg_tab.drop(col, axis=1)
        if col == 'regressor_0':
            reached_regressors = True
    return reg_tab


def multi_trait_analysis(reg_tab, sum_stat_tab, prior_est, trait_num,
                         second_weighting):
    '''Regressing summary statistics of multiple traits onto LD scores'''
    print('Regression using prior estimates to calculate weights...')
    # regression with initial weights
    coeff_list = []
    for i in xrange(trait_num):
        if (i % 100 == 0):
            print('Trait number ' + str(i))
        tab_i = sum_stat_tab[['bgen_pos', 'sum_stats_' + str(i)]].copy()
        tab_i.columns = ['bgen_pos', 'sum_stats']
        # coeff_list.append(single_trait_analysis(reg_tab, tab_i))
        coeff_list.append(single_trait_analysis(reg_tab, tab_i, prior_est, i))
    if second_weighting:
        print('Regression using current estimates to calculate weights...')
        # regression with improved weights
        coeff_list_2 = []
        for i in xrange(trait_num):
            tab_i = sum_stat_tab[['bgen_pos', 'sum_stats_' + str(i)]].copy()
            tab_i.columns = ['bgen_pos', 'sum_stats']
            improved_est = coeff_list[i]
            # coeff_list.append(single_trait_analysis(reg_tab, tab_i))
            coeff_list_2.append(single_trait_analysis(reg_tab, tab_i,
                                                      improved_est, i))
    else:
        coeff_list_2 = coeff_list

    coeffs = np.vstack(coeff_list_2)
    return coeffs


def meta_analyze_traits(reg_tab, sum_stat_tab, prior_est, trait_num,
                        second_weighting):
    '''Uses 'meta_analysis_regression' once and optionally a second time with
    improved weights'''
    print('Regression using prior estimates to calculate weights...')
    # regression with initial weights
    coeffs = meta_analysis_regression(reg_tab, sum_stat_tab, prior_est,
                                      trait_num)
    if second_weighting:
        print('Regression using current estimates to calculate weights...')
        # regression with improved weights
        for i in xrange(trait_num):
            tab_i = sum_stat_tab[['bgen_pos', 'sum_stats_' + str(i)]].copy()
            tab_i.columns = ['bgen_pos', 'sum_stats']
            # improved_est = coeffs
            # coeff_list.append(single_trait_analysis(reg_tab, tab_i))
            # coeff_list_2.append(single_trait_analysis(reg_tab, tab_i,
            # improved_est, i))
    return coeffs


def meta_analysis_regression(reg_tab, sum_stat_tab, prior_est, trait_num):
    reg_tab['sum_stat_prod'] = 0.0
    for i in xrange(trait_num):
        col_name = 'sum_stats_' + str(i)
        reg_tab = reg_tab.merge(sum_stat_tab[['bgen_pos', col_name]],
                                left_on='bgen_pos_1', right_on='bgen_pos',
                                how='left')
        reg_tab = reg_tab.drop('bgen_pos', axis=1)
        reg_tab = reg_tab.rename(columns={col_name: 'first'})
        reg_tab = reg_tab.merge(sum_stat_tab[['bgen_pos', col_name]],
                                left_on='bgen_pos_2', right_on='bgen_pos',
                                how='left')
        reg_tab = reg_tab.drop('bgen_pos', axis=1)
        reg_tab = reg_tab.rename(columns={col_name: 'second'})
        reg_tab.sum_stat_prod += reg_tab['first'] * reg_tab['second']
        reg_tab = reg_tab.drop('first', axis=1)
        reg_tab = reg_tab.drop('second', axis=1)
    reg_tab.sum_stat_prod /= trait_num

    if prior_est is not None:
        reg_weights = calc_reg_weights_fixed(reg_tab, prior_est)
        dist_bin_num = prior_est.size - 1
        for i in xrange(dist_bin_num):
            reg_tab['regressor_' + str(i)] *= reg_weights
        reg_tab['environment'] *= reg_weights
        reg_tab['sum_stat_prod'] *= reg_weights

    regression = sklearn.linear_model.LinearRegression(fit_intercept=False)
    regression.fit(reg_tab.loc[:, 'regressor_0':'environment'],
                   reg_tab.sum_stat_prod)
    reg_tab = reg_tab.drop(['sum_stat_prod'], axis=1)
    return regression.coef_


def jackknife_analysis(reg_tab, sum_stat_tab, prior_est, trait_num,
                       jk_block_num):
    '''Regressing summary statistics of multiple traits onto distance-dependent
    LD scores and using the Jackknife estimator on blocks of SNPs to estimate
    the uncertainty in the estimates'''
    print('Regression using prior estimates to calculate weights...')
    bgen_pos = sum_stat_tab.bgen_pos
    first_col = 'sum_stats_0'
    last_col = 'sum_stats_' + str(trait_num - 1)
    sum_stat_tab = sum_stat_tab.loc[:, first_col:last_col]
    sum_stat_tab['bgen_pos'] = bgen_pos
    col_names = sum_stat_tab.columns

    sum_stat_tab.columns = col_names + '_1'
    reg_tab = reg_tab.merge(sum_stat_tab, on='bgen_pos_1', how='left',
                            indicator=True)
    if not (reg_tab._merge == 'both').all():
        raise ValueError('Some SNPs in reg table have unknown sum stats')
    reg_tab = reg_tab.drop('_merge', axis=1)

    sum_stat_tab.columns = col_names + '_2'
    reg_tab = reg_tab.merge(sum_stat_tab, on='bgen_pos_2', how='left',
                            indicator=True)
    if not (reg_tab._merge == 'both').all():
        raise ValueError('Some SNPs in reg table have unknown sum stats')
    reg_tab = reg_tab.drop('_merge', axis=1)

    sum_stat_tab.columns = col_names

    first_col_1 = 'sum_stats_0_1'
    last_col_1 = 'sum_stats_' + str(trait_num - 1) + '_1'
    first_col_2 = 'sum_stats_0_2'
    last_col_2 = 'sum_stats_' + str(trait_num - 1) + '_2'
    prod_names = np.core.defchararray.add('sum_stats_prod_',
                                          np.arange(trait_num).astype(str))
    sum_stats_prod_tab = pd.DataFrame(
        reg_tab.loc[:, first_col_1:last_col_1].values *
        reg_tab.loc[:, first_col_2:last_col_2].values, columns=prod_names)
    reg_tab = reg_tab.iloc[:, :(-2 * trait_num)]
    reg_tab = pd.concat((reg_tab, sum_stats_prod_tab), axis=1)
    del(sum_stats_prod_tab)

    X = reg_tab.loc[:, 'regressor_0':'environment'].values
    init_col = 'sum_stats_prod_0'
    last_col = 'sum_stats_prod_' + str(trait_num - 1)
    Y = reg_tab.loc[:, init_col:last_col].values

    if prior_est is not None:
        reg_weights = calc_reg_weights_fixed(reg_tab, prior_est)
        reg_weights = reg_weights.reshape((-1, 1))
        X *= reg_weights
        Y *= reg_weights
    XX = X.T.dot(X)
    XY = X.T.dot(Y)
    coeffs = np.linalg.solve(XX, XY).T

    coeffs_jk = np.zeros(list(coeffs.shape) + [jk_block_num])
    jk_block_size = X.shape[0] // jk_block_num
    for i in xrange(jk_block_num):
        from_row = i * jk_block_size
        to_row = (i+1) * jk_block_size
        X_jk = X[from_row:to_row, :]
        Y_jk = Y[from_row:to_row, :]
        XX_jk = X_jk.T.dot(X_jk)
        XY_jk = X_jk.T.dot(Y_jk)
        coeffs_jk[:, :, i] = np.linalg.solve(XX - XX_jk, XY - XY_jk).T

    return coeffs_jk.mean(axis=2), coeffs_jk.var(axis=2) * jk_block_num


def combine_funct_annot(reg_tab, annot_bins, base_var_bins, reg_num):
    '''Sum LD scores from multiple functional annotations; this allows
    regression that ignores these annotations and returns results as if they
    were not provided; convenient when functional annotation LD scores already
    exist'''
    for prefix in ['regressor_', 'LD_score_1_', 'LD_score_2_']:
        final_annot = prefix + str(base_var_bins - 1)
        annot_tab = reg_tab.loc[:, (prefix + '0'): final_annot]
        reg_tab[prefix + '0'] = annot_tab.sum(axis=1)
        for i in xrange(1, annot_bins):
            reg_tab = reg_tab.drop(prefix + str(i), axis=1)
        for i in xrange(annot_bins, reg_num):
            reg_tab = reg_tab.rename(columns=
                                     {(prefix + str(i)):
                                      (prefix + str(i - annot_bins + 1))})
    return reg_tab
