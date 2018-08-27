#!/usr/bin/env python

import ConfigParser
import argparse
import ast
import copy
import itertools
import logging
import os
import re
import sys

import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.metrics import mean_absolute_error

import moff

# debug

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


# filtering _outlier
def MahalanobisDist(x, y):
    covariance_xy = np.cov(x, y, rowvar=0)
    inv_covariance_xy = np.linalg.inv(covariance_xy)
    xy_mean = np.mean(x), np.mean(y)
    x_diff = np.array([x_i - xy_mean[0] for x_i in x])
    y_diff = np.array([y_i - xy_mean[1] for y_i in y])
    diff_xy = np.transpose([x_diff, y_diff])
    md = []
    for i in range(len(diff_xy)):
        md.append(np.sqrt(
            np.dot(np.dot(np.transpose(diff_xy[i]), inv_covariance_xy), diff_xy[i])))
    return md


# remove outlier
def MD_removeOutliers(x, y, width):
    MD = MahalanobisDist(x, y)
    threshold = np.mean(MD) * float(width)  # adjust 1.5 accordingly
    nx, ny, outliers = [], [], []
    for i in range(len(MD)):
        if MD[i] <= threshold:
            nx.append(x[i])
            ny.append(y[i])
        else:
            outliers.append(i)  # position of removed pair
    return (np.array(nx), np.array(ny), np.array(outliers))


# combination of rt predicted by each single model
def combine_model(x, model, err, weight_flag):
    x = x.values
    tot_err = np.sum(np.array(err)[np.where(~np.isnan(x))])

    app_sum = 0
    app_sum_2 = 0
    for ii in range(0, len(x)):
        if ~  np.isnan(x[ii]):
            if int(weight_flag) == 0:
                app_sum = app_sum + (model[ii].predict(x[ii])[0][0])
            else:
                app_sum_2 = app_sum_2 + \
                    (model[ii].predict(x[ii])[0][0] *
                     (float(err[ii]) / float(tot_err)))

                # " output weighted mean
    if int(weight_flag) == 1:
        return float(app_sum_2)
    else:
        # output not weight
        return float(app_sum) / float(np.where(~ np.isnan(x))[0].shape[0])


'''
def train_gp(data_A,data_B,c=None):
    bins = np.linspace(data_B.min()-2, data_B.max()+1,20)
    digitized = np.digitize(data_B, bins)
    size_bin = [  digitized[digitized == i].shape[0]  for i in range(1, len(bins))]
    tt_x = np.concatenate([random.sample(np.where(digitized == i)[0], 20) for i in range(1, len(bins))])
    test =  np.setdiff1d(range(data_A.shape[0]),tt_x)
    #ff = ym_test_predicted[:, 0] - np.sqrt(np.diag(y_cov))
    #dd = ym_test_predicted[:, 0] + np.sqrt(np.diag(y_cov))
    # plt.fill_between(data_B[test,0], ff  , dd , alpha=0.5,color='k')
    #size_train= int (data_A.shape[0] * 0.10)
    #rows = random.sample(range(data_A.shape[0]),size_train)
    #data_A= data_A[rows,:]
    #data_B= data_B[rows,:]
    # data_B is x
    # data_A if Y
    #kern = GPy.kern.Linear(1)
    #kernel = 3 * DotProduct() + WhiteKernel(noise_level=0.1,noise_level_bounds=(1e-3, 1e-0))

    kernel = 0.5 * RBF(length_scale=1, length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(noise_level=0.3,noise_level_bounds=(1e-3, 2e-0))
    m = GaussianProcessRegressor(kernel=kernel, alpha=0.1, normalize_y=False, n_restarts_optimizer=1).fit(data_B[tt_x], data_A[tt_x])

    ym_train_predicted, y_cov_train = m.predict(data_B[tt_x], return_std=False, return_cov=True)
    ym_test_predicted, y_cov_test = m.predict(data_B[test], return_std=False, return_cov=True)

    #ff =  np.sqrt(np.diag(y_cov_test))
    #dd = np.sqrt(np.diag(y_cov_test))

    ff = ym_test_predicted[:, 0] - (1.96 * np.sqrt(np.diag(y_cov_test)))
    dd = ym_test_predicted[:, 0] + ( 1.96 * np.sqrt(np.diag(y_cov_test)) )

    ##printing modell
    plt.figure(figsize=(15, 14), dpi=100)
    #plt.scatter(data_B[test],ym_test_predicted,marker='*',c='red',s=15 )

    plt.fill_between(data_B[test, 0], ff, dd, alpha=0.9, color='r')
    plt.scatter(data_B[test],data_A[test],marker='<',c='black',s=25,label='True RT' )
    plt.scatter(data_B[test], ym_test_predicted, marker='*', c='blue', s=25,label='predicted RT')
    plt.legend(loc="best", scatterpoints=1, prop={'size': 18})

    plt.title('GP with Rbf_+ whiteNoise on test set: ' + c )
    plt.savefig( 'D:\\workspace\\ionstar_dataset\\mbr_output\\' + c + '__model.png' )

    log.critical(' Size error training : %i sec ', tt_x.shape[0])
    log.critical(' Mean absolute error training : %4.4f sec',mean_absolute_error(data_A[tt_x], ym_train_predicted))
    log.critical(' Mean absolute error test_set (not sampled point) : %4.4f sec',mean_absolute_error(data_A[test], ym_test_predicted))

    return m, ym_train_predicted,mean_absolute_error(data_A[tt_x], ym_train_predicted)



def combine_model_GP(x, model, err, weight_flag):

    ra_flag= 0


    #tot_err =  1- ( (np.array(err)[np.where(~np.isnan(x))]) / np.max(np.array(err)[np.where(~np.isnan(x))]))
    tot_err = np.sum(np.array(err)[np.where(~pd.isnull(x))])
    #print tot_err
    #print x
    app_sum = 0
    app_sum_2 = 0
    app_var =0
    for ii in range(0, len(x)):


        if ~  np.isnan(x[ii]):
            #sklearn
            pred, y_cov = model[ii].predict(x[ii].reshape(1, 1), return_std=False, return_cov=True)
            var = np.sqrt(np.diag(y_cov))
            #pred, var = model[ii].predict(x[ii].reshape(1, 1), include_likelihood=True)
            #ress = model[ii].predict_quantiles(x[ii].reshape(1, 1))

            #print ' %i Input Rt  %4.4f  Predicted: %4.4f Var %4.4f  Interval at 95 %4.4f <--> %4.4f  ' % (
            #ii, x[ii], float(pred), float(var), pred - (1.965 * var) ,  pred + (1.965 * var)  )
            #print 'intervall width length %4.4f' % abs((pred - (1.965 * var)) - (pred + (1.965 * var)))


            if int(weight_flag) == 0:
                app_sum = app_sum + (pred)
                app_var += var
            else:
                # print ii,model[ii].predict(x[ii])[0][0]
                w = (float(err[ii]) / float(tot_err))
                # w= tot_err[ii]
                # print ii ,'weighted', (model[ii].predict(x[ii])[0][0] * w ),w
                app_sum_2 = app_sum_2 + (pred * w)
                app_var += var
    # " output weighted mean

    if  float(np.where(~ np.isnan(x))[0].shape[0]) == 0.0 :
        return pd.Series({'time_pred': -1 , 'uncertainty_win': -1})

    if int(weight_flag) == 1:
        f_p = app_sum_2
        mean_var = float(app_var) / float(np.where(~ np.isnan(x))[0].shape[0])
    else:
        mean_var = float(app_var) / float(np.where(~ np.isnan(x))[0].shape[0])
        f_p = float(app_sum) / float(np.where(~ pd.isnull(x))[0].shape[0])

    log.critical( '   -->  Final Aggr Predicted: %4.4f Var %4.4f  Interval at 95 %4.4f <--> %4.4f  ' % (  f_p , float(mean_var), f_p - (1.965 * mean_var), f_p + (1.965 * mean_var)))
    log.critical ('   -->  Final intervall width length  %4.4f  |  #_inputs %i' %  (abs(    (f_p - (1.965 * mean_var)) - (f_p + (1.965 * mean_var)) ),int(np.where(~ np.isnan(x))[0].shape[0])) )

    return pd.Series({'time_pred': f_p, 'uncertainty_win': 1.965 * mean_var})

'''

# run the mbr in moFF : input  ms2 identified peptide   output csv file with the matched peptides added


def run_mbr(args):
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    log.addHandler(ch)

    if args.loc_in is None:
        # the user uses --inputtsv option
        if not (args.loc_out is None):
            # if the user use --output_folder the mbr folder will be created there
            output_dir = os.path.join(args.loc_out, 'mbr_output')
        else:
            # if the user does not use  --output_folder the mbr folder will be created on moFF path location
            output_dir = os.path.join('mbr_output')
            print os.path.abspath(output_dir)

    else:
        # the user use the --inputF option
        if os.path.exists(os.path.join(args.loc_in)):
            # if '/' in  str(args.loc_in):
            output_dir = os.path.join(args.loc_in, 'mbr_output')
        else:
            exit(os.path.join(args.loc_in) +
                 ' EXIT input folder path is not well specified --> / missing or wrong path')

            # if not (os.path.isdir(args.loc_in)):
            #   exit(str(args.loc_in) + '-->  input folder does not exist ! ')

            # if str(args.loc_in) == '':
            #    output_dir = 'mbr_output'
            # else:
            #    if os.path.exists(os.path.join(args.loc_in)):
            # if '/' in  str(args.loc_in):
    # output_dir = os.path.join(args.loc_in, 'mbr_output')
    #    else:
    #        exit(os.path.join(args.loc_in) + ' EXIT input folder path not well specified --> / missing ')

    if not (os.path.isdir(output_dir)):

        log.critical("Created MBR output folder in : %s ",
                     os.path.abspath(output_dir))
        os.makedirs(output_dir)
    else:
        log.critical("MBR Output folder in : %s ", os.path.abspath(output_dir))
    # set log to file
    w_mbr = logging.FileHandler(os.path.join(
        output_dir, args.log_label + '_mbr_.log'), mode='w')
    w_mbr.setLevel(logging.INFO)
    log.addHandler(w_mbr)

    moff_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    config = ConfigParser.RawConfigParser()
    config.read(os.path.join(moff_path, 'moff_setting.properties'))

    # it s always placed in same folder of moff_mbr.py
    # read input
    # comment better
    # name of the input file
    exp_set = []
    # list of the input dataframe
    exp_t = []
    # list of the output dataframe
    exp_out = []
    # lsit of input datafra used as help
    exp_subset = []
    # list of the name of the mbr output
    exp_out_name = []

    if args.loc_in is None:
        for id_name in args.tsv_list:
            exp_set.append(id_name)
    else:
        for item in os.listdir(args.loc_in):
            log.critical(item)
            if os.path.isfile(os.path.join(args.loc_in, item)):
                if os.path.join(args.loc_in, item).endswith('.' + args.ext):
                    exp_set.append(os.path.join(args.loc_in, item))

                # sample optiion is valid only if  folder iin option is valid
    if (args.sample is not None) and (args.loc_in is not None):
        exp_set_app = copy.deepcopy(exp_set)
        for a in exp_set:
            if re.search(args.sample, a) is None:
                exp_set_app.remove(a)
        exp_set = exp_set_app

    if (exp_set == []) or (len(exp_set) == 1):
        exit(
            'ERROR input files not found or just one input file selected . check the folder or the extension given in input')

    for a in exp_set:
        log.critical('Reading file: %s ', a)
        exp_subset.append(a)
        data_moff = pd.read_csv(a, sep="\t", header=0)
        data_moff['rt'] = data_moff['rt'] / 60
        list_name = data_moff.columns.values.tolist()
        # get the lists of PS  defaultcolumns from properties file
        list_ps_def = ast.literal_eval(
            config.get('moFF', 'ps_default_export_v1'))
        # here it controls if the input file is a PS export; if yes it maps the input in right moFF name
        if moff.check_ps_input_data(list_name, list_ps_def) == 1:
            log.critical(
                'Detected input file from PeptideShaker export..: %s ', a)
            # map  the columns name according to moFF input requirements
            data_moff, list_name = moff.map_ps2moff(
                data_moff, 'col_must_have_mbr')
            log.critical(
                'Mapping columns names into the  the moFF requested column name..: %s ', a)
            # print data_moff.columns
        if moff.check_columns_name(list_name, ast.literal_eval(config.get('moFF', 'col_must_have_mbr')), log) == 1:
            exit('ERROR minimal field requested are missing or wrong')
        data_moff['matched'] = 0
        data_moff['mass'] = data_moff['mass'].map('{:.4f}'.format)

        data_moff['code_unique'] = data_moff['mod_peptide'].astype(
            str)  # + '_' + data_moff['mass'].astype(str)
        data_moff = data_moff.sort_values(by='rt')
        exp_t.append(data_moff)
        exp_out.append(data_moff)

    log.critical('Read input --> done ')
    # parameter of the number of query
    # set a list of filed mandatory
    # ['matched','peptide','mass','mz','charge','prot','rt']
    n_replicates = len(exp_t)
    exp_set = exp_subset
    aa = range(0, n_replicates)
    out = list(itertools.product(aa, repeat=2))
    # just to save all the model
    # list of the model saved
    model_save = []
    # list of the error in min/or sec
    model_err = []
    # list of the status of the model -1 means model not available for low points in the training set
    model_status = []
    # add matched columns
    list_name.append('matched')
    # final status -1 if one of the output is empty
    out_flag = 1
    # input of the methods

    log.info('Outlier Filtering is %s  ', 'active' if int(
        args.out_flag) == 1 else 'not active')
    log.info('Number of replicates %i,', n_replicates)
    log.info('Pairwise model computation ----')

    if args.rt_feat_file is not None:
        log.critical(
            'Custom list of peptide used  provided by the user in %s', args.rt_feat_file)
        # log.info('Custom list of peptide used  provided by the user in %s', args.rt_feat_file)
        shared_pep_list = pd.read_csv(args.rt_feat_file, sep='\t')
        shared_pep_list['mass'] = shared_pep_list['mass'].map('{:.4f}'.format)
        shared_pep_list['code'] = shared_pep_list['peptide'].astype(
            str) + '_' + shared_pep_list['mass'].astype(str)
        list_shared_pep = shared_pep_list['code']
        log.info('Custom list of peptide contains  %i ',
                 list_shared_pep.shape[0])

    for jj in aa:
        log.info('matching  in %s', exp_set[jj])

        for i in out:
            if i[0] == jj and i[1] != jj:
                if args.rt_feat_file is not None:
                    # use of custom peptide
                    comA = exp_t[i[0]][exp_t[i[0]]['code_unique'].isin(list_shared_pep)][
                        ['code_unique', 'peptide', 'prot', 'rt']]
                    comB = exp_t[i[1]][exp_t[i[1]]['code_unique'].isin(list_shared_pep)][
                        ['code_unique', 'peptide', 'prot', 'rt']]
                    comA = comA.groupby('code_unique', as_index=False).mean()
                    comB = comB.groupby('code_unique', as_index=False).mean()
                    common = pd.merge(
                        comA, comB, on=['code_unique'], how='inner')
                else:
                    # use of shared peptdes.
                    log.info('  Matching  %s peptide in   searching in %s ',
                             exp_set[i[0]], exp_set[i[1]])
                    list_pep_repA = exp_t[i[0]]['code_unique'].unique()
                    list_pep_repB = exp_t[i[1]]['code_unique'].unique()
                    log.info('Peptide unique (mass + sequence) %i , %i ',
                             list_pep_repA.shape[0],
                             list_pep_repB.shape[0])
                    set_dif_s_in_1 = np.setdiff1d(list_pep_repB, list_pep_repA)
                    add_pep_frame = exp_t[i[1]][exp_t[i[1]]['code_unique'].isin(
                        set_dif_s_in_1)].copy()
                    pep_shared = np.intersect1d(list_pep_repA, list_pep_repB)
                    log.info(
                        '  Peptide (mass + sequence)  added size  %i ', add_pep_frame.shape[0])
                    log.info('  Peptide (mass + sequence) )shared  %i ',
                             pep_shared.shape[0])
                    comA = exp_t[i[0]][exp_t[i[0]]['code_unique'].isin(pep_shared)][
                        ['code_unique', 'peptide', 'prot', 'rt']]
                    comB = exp_t[i[1]][exp_t[i[1]]['code_unique'].isin(pep_shared)][
                        ['code_unique', 'peptide', 'prot', 'rt']]
                    # filtering using the variance added 17_08
                    flag_var_filt = False
                    if flag_var_filt:
                        dd = comA.groupby('code_unique', as_index=False)
                        top_res = dd.agg(['std', 'mean', 'count'])
                        # print np.nanpercentile(top_res['rt']['std'].values,[5,10,20,30,50,60,80,90,95,97,99,100])
                        th = np.nanpercentile(top_res['rt']['std'].values, 60)
                        comA = comA[~ comA['code_unique'].isin(
                            top_res[top_res['rt']['std'] > th].index)]
                        # data B '
                        dd = comB.groupby('code_unique', as_index=False)

                        top_res = dd.agg(['std', 'mean', 'count'])
                        # print comB.shape
                        # print np.nanpercentile(top_res['rt']['std'].values,[5,10,20,30,50,60,80,90,95,97,99,100])
                        th = np.nanpercentile(top_res['rt']['std'].values, 60)
                        comB = comB[~ comB['code_unique'].isin(
                            top_res[top_res['rt']['std'] > th].index)]

                    comA = comA.groupby('code_unique', as_index=False).mean()
                    comB = comB.groupby('code_unique', as_index=False).mean()
                    common = pd.merge(
                        comA, comB, on=['code_unique'], how='inner')
                if common.shape[0] <= 10 and args.rt_feat_file is not None:
                    model_status.append(-1)
                    continue
                # filtering outlier option
                else:
                    if int(args.out_flag) == 1:
                        filt_x, filt_y, pos_out = MD_removeOutliers(common['rt_y'].values, common['rt_x'].values,
                                                                    args.w_filt)
                        data_B = filt_x
                        data_A = filt_y
                        data_B = np.reshape(data_B, [filt_x.shape[0], 1])
                        data_A = np.reshape(data_A, [filt_y.shape[0], 1])
                        log.info('Outlier founded %i  w.r.t %i',
                                 pos_out.shape[0], common['rt_y'].shape[0])
                    else:
                        data_B = common['rt_y'].values
                        data_A = common['rt_x'].values
                        data_B = np.reshape(data_B, [common.shape[0], 1])
                        data_A = np.reshape(data_A, [common.shape[0], 1])

                    log.info(' Size trainig shared peptide , %i %i ',
                             data_A.shape[0], data_B.shape[0])
                    clf = linear_model.RidgeCV(alphas=np.power(
                        2, np.linspace(-30, 30)), scoring='neg_mean_absolute_error')
                    clf.fit(data_B, data_A)
                    clf_final = linear_model.Ridge(alpha=clf.alpha_)
                    clf_final.fit(data_B, data_A)
                    # save the model
                    model_save.append(clf_final)
                    model_err.append(mean_absolute_error(
                        data_A, clf_final.predict(data_B)))
                    log.info(' Mean absolute error training : %4.4f sec',
                             mean_absolute_error(data_A, clf_final.predict(data_B)))
                    '''
                    # GP version
                    model_gp, predicted_train, error = train_gp(data_A, data_B,c= str(i[0])+'_'+str(i[1]))
                    #print i[1], comA.shape, error

                    model_err.append(error)
                    model_save.append(model_gp)
                    model_status.append(1)
                    '''
    if np.where(np.array(model_status) == -1)[0].shape[0] >= (len(aa) / 2):
        log.error(
            'MBR aborted :  mbr cannnot be run, not enough shared pepetide among the replicates ')
        exit('ERROR : mbr cannnot be run, not enough shared pepetide among the replicates')

    log.info('Combination of the  model  --------')
    log.info('Weighted combination  %s : ', 'Weighted' if int(
        args.w_comb) == 1 else 'Unweighted')

    diff_field = np.setdiff1d(exp_t[0].columns, [
                              'matched', 'mod_peptide', 'peptide', 'mass', 'mz', 'charge', 'prot', 'rt'])

    for jj in aa:
        pre_pep_save = []
        log.critical('Predict rt for the exp.  in %s ', exp_set[jj])
        c_rt = 0
        for i in out:
            if i[0] == jj and i[1] != jj:
                log.info('Matching peptides found  in  %s ', exp_set[i[1]])
                list_pep_repA = exp_t[i[0]]['peptide'].unique()
                list_pep_repB = exp_t[i[1]]['peptide'].unique()
                set_dif_s_in_1 = np.setdiff1d(list_pep_repB, list_pep_repA)
                add_pep_frame = exp_t[i[1]][exp_t[i[1]]
                                            ['peptide'].isin(set_dif_s_in_1)].copy()
                add_pep_frame = add_pep_frame[[
                    'peptide', 'mod_peptide', 'mass', 'mz', 'charge', 'prot', 'rt']]
                # add_pep_frame['code_unique'] = '_'.join([add_pep_frame['peptide'], add_pep_frame['prot'], add_pep_frame['mass'].astype(str), add_pep_frame['charge'].astype(str)])
                add_pep_frame['code_unique'] = add_pep_frame['mod_peptide'] + '_' + \
                    add_pep_frame['prot'] + '_' + '_' + \
                    add_pep_frame['charge'].astype(str)
                add_pep_frame = add_pep_frame.groupby('code_unique', as_index=False)[
                    'peptide', 'mod_peptide', 'mass', 'charge', 'mz', 'prot', 'rt'].aggregate(max)
                add_pep_frame = add_pep_frame[[
                    'code_unique', 'peptide', 'mod_peptide', 'mass', 'mz', 'charge', 'prot', 'rt']]
                list_name = add_pep_frame.columns.tolist()
                list_name = [w.replace('rt', 'rt_' + str(c_rt))
                             for w in list_name]
                add_pep_frame.columns = list_name
                pre_pep_save.append(add_pep_frame)
                c_rt += 1
                # print 'input columns',pre_pep_save[0].columns
        if n_replicates == 2:
            test = pre_pep_save[0]
        else:
            test = reduce(
                lambda left, right: pd.merge(left, right, on=[
                                             'code_unique', 'peptide', 'mod_peptide', 'mass', 'mz', 'charge', 'prot'], how='outer'),
                pre_pep_save)

        # aggregated by code_unique,  to avoid duplicates
        test = test.groupby('code_unique', as_index=False).aggregate(max)
        test.drop('code_unique', axis=1, inplace=True)

        test['time_pred'] = test.iloc[:, 6: (6 + (n_replicates - 1))].apply(
            lambda x: combine_model(x, model_save[(jj * (n_replicates - 1)):((jj + 1) * (n_replicates - 1))],
                                    model_err[(jj * (n_replicates - 1)):((jj + 1) * (n_replicates - 1))], args.w_comb), axis=1)
        # test[['time_pred','uncertainty_win']] = test.ix[:, 6: (6 + (n_replicates - 1))].apply(lambda x: combine_model_GP(x, model_save[(jj * (n_replicates - 1)):((jj + 1) * (n_replicates - 1))],
        #                            model_err[(jj * (n_replicates - 1)):((jj + 1) * (n_replicates - 1))], args.w_comb), axis=1)

        test['matched'] = 1

        # still to check better
        if test[test['time_pred'] <= 0].shape[0] >= 1:
            log.info(' -- Predicted negative RT: those peptide will be deleted')
            test = test[test['time_pred'] > 0]

        list_name = test.columns.tolist()
        list_name = [w.replace('time_pred', 'rt') for w in list_name]
        test.columns = list_name

        # test = test[['peptide','mod_peptide', 'mass', 'mz', 'charge', 'prot', 'rt', 'matched','uncertainty_win']]

        test = test[['peptide', 'mod_peptide', 'mass',
                     'mz', 'charge', 'prot', 'rt', 'matched']]
        # just put nan with the missing values
        for field in diff_field.tolist():
            test[field] = np.nan  # -1

        log.info('Before adding %s contains %i ',
                 exp_set[jj], exp_t[jj].shape[0])
        exp_out[jj] = pd.concat(
            [exp_t[jj], test], join='outer', axis=0, sort=False)
        log.info('After MBR %s contains:  %i  peptides',
                 exp_set[jj], exp_out[jj].shape[0])
        log.critical('matched features   %i  MS2 features  %i ', exp_out[jj][exp_out[jj]['matched'] == 1].shape[0],
                     exp_out[jj][exp_out[jj]['matched'] == 0].shape[0])
        exp_out[jj].to_csv(
            path_or_buf=os.path.join(output_dir, os.path.split(exp_set[jj])[1].split('.')[0] + '_match.txt'), sep='\t',
            index=False)
        exp_out_name.append(os.path.join(output_dir, os.path.split(
            exp_set[jj])[1].split('.')[0] + '_match.txt'))
        if exp_out[jj].shape[0] > 0:
            out_flag = 1 * out_flag
        else:
            out_flag = -1 * out_flag

    w_mbr.close()
    log.removeHandler(w_mbr)
    return out_flag, exp_out_name


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='moFF match between run input parameter')

    parser.add_argument('--loc_in', dest='loc_in', action='store',
                        help='specify the folder of the input MS2 peptide files  REQUIRED]', required=True)

    parser.add_argument('--sample', dest='sample', action='store',
                        help='specify which replicate files are used fot mbr [regular expr. are valid] ',
                        required=False)

    parser.add_argument('--ext', dest='ext', action='store', default='txt',
                        help='specify the exstension of the input file (txt as default value)', required=False)

    parser.add_argument('--log_label', dest='log_label', default='moFF', action='store',
                        help='a label name for the log file (moFF_mbr.log as default log file name)', required=False)

    parser.add_argument('--w_filt', dest='w_filt', action='store', default=2,
                        help='width value of the filter (k * mean(Dist_Malahobi , k = 2 as default) ', required=False)

    parser.add_argument('--out_flag', dest='out_flag', action='store', default=1,
                        help='filter outlier in each rt time allignment (active as default)', required=False)

    parser.add_argument('--w_comb', dest='w_comb', action='store', default=0,
                        help='weights for model combination combination : 0 for no weight (default) 1 weighted devised by model errors.',
                        required=False)

    parser.add_argument('--rt_feat_file', dest='rt_feat_file', action='store',
                        help='specify the file that contains the features to use in the match-between-run RT prediction ',
                        required=False)

    args = parser.parse_args()

    run_mbr(args)
