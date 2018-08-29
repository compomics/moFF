# /usr/bin/env python

import bisect
import glob
import logging
import multiprocessing
import os as os
import shlex
import subprocess
import sys
# import time
import traceback
from collections import Counter
from itertools import chain
from sys import platform as _platform

import numpy as np
import pandas as pd
import pymzml
import simplejson as json
from brainpy import isotopic_variants
from pyteomics.mass import std_aa_comp
from scipy.stats import spearmanr

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

"""
moFF: this module contains all core utilities such as apex computation, peak apex computation etc..
"""


TXIC_PATH = os.environ.get('TXIC_PATH', './')


def set_logger(name_file):
    if len(log.handlers) == 0:
        ch = logging.StreamHandler()
        ch.setLevel(logging.ERROR)
        log.addHandler(ch)
        fh = logging.FileHandler(name_file, mode='a')
        fh.setLevel(logging.DEBUG)
        # formatter = logging.Formatter('%(message)s')
        # fh.setFormatter(formatter)
        log.addHandler(fh)


def detach_handler():
    handlers = log.handlers[:]
    for handler in handlers:
        handler.close()
        log.removeHandler(handler)


def clean_json_temp_file(loc_output):
    for f in glob.glob(loc_output + "/*.json"):
        os.remove(f)
    return 1


def compute_peptide_matrix(loc_output, log, tag_filename):
    """
    Computation of the export summary intensities peptides
    :param loc_output:
    :param log:
    :param tag_filename:
    :return:
    """
    name_col = []
    name_col.append('prot')
    d = []
    if not glob.glob(loc_output + '/*_moff_result.txt'):
        return False
    for name in glob.glob(loc_output + '/*_moff_result.txt'):

        if 'match_' in os.path.basename(name):
            name_col.append(
                'sumIntensity_' + os.path.basename(name).split('_match_moff_result.txt')[0])
        else:
            name_col.append(
                'sumIntensity_' + os.path.basename(name).split('_moff_result.txt')[0])
        data = pd.read_csv(name, sep="\t")

        '''
        Other possibile quality controll filter
        data = data[ data['lwhm'] != -1]
        data = data[data['rwhm'] != -1 ]
        '''

        data = data[data['intensity'] != -1]
        data.sort_values('rt', ascending=True, inplace=True)
        log.critical('Collecting moFF result file : %s   --> Retrived peptide peaks after filtering:  %i',
                     os.path.basename(name), data.shape[0])
        # cleaning peptide fragmented more than one time. we keep the earliest one
        data.drop_duplicates(subset=[
                             'prot', 'peptide', 'mod_peptide', 'mass', 'charge'], keep='first', inplace=True)
        d.append(data[['prot', 'peptide', 'mod_peptide', 'mass',
                       'charge', 'rt_peak', 'rt', 'intensity']])

    intersect_share = reduce(np.union1d, ([x['peptide'].unique() for x in d]))
    index = intersect_share

    df = pd.DataFrame(index=index, columns=name_col)
    df = df.fillna(0)
    for i in range(0, len(d)):
        grouped = d[i].groupby('peptide', as_index=True)['prot', 'intensity']
        # print grouped.agg({'prot':'max', 'intensity':'sum'}).columns
        df.ix[:, i + 1] = grouped.agg({'prot': 'max',
                                       'intensity': 'sum'})['intensity']
        df.ix[np.intersect1d(df.index, grouped.groups.keys()), 0] = grouped.agg({'prot': 'max', 'intensity': 'sum'})[
            'prot']
    # print df.head(5)
    df.reset_index(level=0, inplace=True)
    df = df.fillna(0)
    df.rename(columns={'index': 'peptide'}, inplace=True)
    log.critical('Writing peptide_summary intensity file')
    df.to_csv(os.path.join(loc_output, "peptide_summary_intensity_" +
                           tag_filename + ".tab"), sep='\t', index=False)
    return True


def save_moff_apex_result(result):
    """
    Collect all CPU results in a data frame

    :param result:
    :return:
    """
    try:
        xx = []
        for df_index in result:
            if result[df_index].get()[1] == -1:
                exit('Raw file not retrieved: wrong path or upper/low case mismatch')
            else:
                xx.append(result[df_index].get()[0])

        final_res = pd.concat(xx)
        if 'index' in final_res.columns:
            final_res.drop('index', axis=1, inplace=True)

    except Exception as e:
        traceback.print_exc()
        print
        raise e
    return (final_res)


def map_ps2moff(data, type_mapping):
    data.drop(data.columns[[0]], axis=1, inplace=True)
    data.columns = data.columns.str.lower()
    if type_mapping == 'col_must_have_mbr':
        data.rename(columns={'sequence': 'peptide', 'modified sequence': 'mod_peptide', 'measured charge': 'charge',
                             'theoretical mass': 'mass', 'protein(s)': 'prot', 'm/z': 'mz'}, inplace=True)
    if type_mapping == 'col_must_have_apex':
        data.rename(columns={'sequence': 'peptide', 'modified sequence': 'mod_peptide', 'measured charge': 'charge', 'theoretical mass': 'mass',
                             'protein(s)': 'prot', 'm/z': 'mz'}, inplace=True)
    return data, data.columns.values.tolist()




def check_ps_input_data(input_column_name, list_col_ps_default):
    """
     Control if the input data is complaint with PS input file

    :param input_column_name:
    :param list_col_ps_default:
    :return:
    """
    input_column_name.sort()
    list_col_ps_default.sort()
    if list_col_ps_default == input_column_name:
        # detected a default PS input file
        return 1
    else:
        # not detected a default PS input file
        return 0


def check_columns_name(col_list, col_must_have, log):
    """
    Controls if the  current input  file informations are complaint with the minimun set of informations need by moFF

    :param col_list:
    :param col_must_have:
    :param log:
    :return:
    """
    for c_name in col_must_have:
        if not (c_name in col_list):
            # fail
            log.critical('This information is missing : %s ', c_name)
            return 1
        # succes
    return 0


def scan_mzml(name):
    # when I am using thermo raw and --raw_repo option used
    if name is None:
        return (-1, -1)
    if 'MZML' in name.upper():

        rt_list = []
        runid_list = []
        run_temp = pymzml.run.Reader(name)
        for spectrum in run_temp:
            if spectrum['ms level'] == 1:
                rt_list.append(spectrum['scan start time'])
                runid_list.append(spectrum['id'])

        return (rt_list, runid_list)
    else:
        # in case of raw file  I put to -1 -1 thm result
        return (-1, -1)


def mzML_get_all(temp, tol, loc, run, rt_list1, runid_list1):
    app_list = []
    for index_ms2, row in temp.iterrows():

        data, status = pyMZML_xic_out(loc, float(
            tol / (10 ** 6)), row['ts'], row['te'], row['mz'], run, runid_list1, rt_list1)
        # status is evaluated only herenot used anymore
        if status != -1:
            app_list.append(data)
        else:
            app_list.append(pd.DataFrame(columns=['rt', 'intensity']))
    return app_list


def pyMZML_xic_out(name, ppmPrecision, minRT, maxRT, MZValue, run, runid_list, rt_list):
    timeDependentIntensities = []
    minpos = bisect.bisect_left(rt_list, minRT)
    maxpos = bisect.bisect_left(rt_list, maxRT)

    for specpos in range(minpos, maxpos):
        specid = runid_list[specpos]
        spectrum = run[specid]
        if spectrum['scan start time'] > maxRT:
            break
        if spectrum['scan start time'] > minRT and spectrum['scan start time'] < maxRT:
            lower_index = bisect.bisect(
                spectrum.peaks, (float(MZValue - ppmPrecision * MZValue), None))
            upper_index = bisect.bisect(
                spectrum.peaks, (float(MZValue + ppmPrecision * MZValue), None))
            maxI = 0.0
            for sp in spectrum.peaks[lower_index: upper_index]:
                if sp[1] > maxI:
                    maxI = sp[1]
            if maxI > 0:
                timeDependentIntensities.append(
                    [spectrum['scan start time'], maxI])

    if len(timeDependentIntensities) > 5:
        return (pd.DataFrame(timeDependentIntensities, columns=['rt', 'intensity']), 1)
    else:
        return (pd.DataFrame(timeDependentIntensities, columns=['rt', 'intensity']), -1)


def check_log_existence(file_to_check):
    """
    Controls the presence of a log file

    :param file_to_check:
    :return:
    """
    if os.path.isfile(file_to_check):
        os.remove(file_to_check)
        return True
    else:
        return False


def check_output_folder_existence(loc_output):
    """
    Controls the presence of a directory. if not, it makes it

    :param loc_output:
    :return:
    """

    if not os.path.exists(loc_output):
        os.mkdir(loc_output)
        return 1
    else:
        return 0


def compute_log_LR(data_xic, index, v_max, disc):
    """
    Computation shape peak metrics log_L_R

    :param data_xic:
    :param index:
    :param v_max:
    :param disc:
    :return:
    """
    log_time = [-1, -1]
    c_left = 0
    find_5 = False
    stop = False
    while c_left <= (index - 1) and not stop:
        if not find_5 and (data_xic.ix[(index - 1) - c_left, 1] <= (disc * v_max)):
            find_5 = True
            log_time[0] = data_xic.ix[(index - 1) - c_left, 0] * 60
            stop = True
        c_left += 1
    find_5 = False
    stop = False
    r_left = 0
    while ((index + 1) + r_left < data_xic.shape[0]) and not stop:
        if not find_5 and data_xic.ix[(index + 1) + r_left, 1] <= (disc * v_max):
            find_5 = True
            log_time[1] = data_xic.ix[(index + 1) + r_left, 0] * 60
            stop = True
        r_left += 1
    return log_time


def compute_peak_simple(x, xic_array, log, mbr_flag, h_rt_w, s_w, s_w_match, offset_index, moff_pride_flag, rt_match_peak, count_match, filt_flag):
    """

    Apex computation method

    :param x:
    :param xic_array:
    :param log:
    :param mbr_flag:
    :param h_rt_w:
    :param s_w:
    :param s_w_match:
    :param offset_index:
    :param moff_pride_flag:
    :param rt_match_peak:
    :param count_match:
    :param filt_flag:
    :return:
    """
    if count_match != -1:
        c = x.prog_xic_index
    else:
        c = x.name
    data_xic = xic_array[c]
    if rt_match_peak > -1:
        time_w = rt_match_peak
    else:
        time_w = x['rt']
        if not moff_pride_flag :
            # NOT moff pride data
            # dealling with rt in minutes
            # standar cases rt must be in second
            time_w = time_w / 60
    # print time_w, x['rt'] , moff_pride_flag, rt_match_peak,time_w, s_w,s_w_match
    if not mbr_flag :
        temp_w = s_w
    else:
        # row['matched'])
        if x['matched'] == 1:
            temp_w = s_w_match
        else:
            temp_w = s_w
    if data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].shape[0] >= 1:
        # data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].to_csv('thermo_testXIC_'+str(c)+'.txt',index=False,sep='\t')
        ind_v = data_xic.index
        pp = data_xic[data_xic["intensity"] == data_xic[(data_xic['rt'] > (
            time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))]['intensity'].max()].index
        pos_p = ind_v[pp]
        if pos_p.values.shape[0] > 1:
            log.info('error, no apex found')
            return pd.Series({'intensity': -1, 'rt_peak': -1, 'lwhm': -1, 'rwhm': -1, '5p_noise': -1, '10p_noise': -1, 'SNR': -1, 'log_L_R': -1, 'log_int': -1})
        val_max = data_xic.ix[pos_p, 1].values
    else:
        if filt_flag == 1:
            if 'matched'in x.axes[0].tolist():
                log.info('peptide %r -->  MZ: %4.4f RT: %4.4f matched (yes=1/no=0): %i Peak not detected  Xic shape %r ',
                         x['mod_peptide'], x['mz'], time_w, x['matched'], data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].shape[0])
            else:
                log.info('peptide %r -->  MZ: %4.4f RT: %4.4f matched (yes=1/no=0): %i Peak not detected  Xic shape %r ',
                         x['mod_peptide'], x['mz'], time_w, 0, data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].shape[0])
        # log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f ', (offset_index +c +2), x['mz'], time_w)
        # log.info("\t LW_BOUND window  %4.4f", time_w - temp_w)
        # log.info("\t UP_BOUND window %4.4f", time_w + temp_w)
        # log.info("\t WARNINGS: Peak not detected  Xic shape %r ", data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].shape[0])
        return pd.Series({'intensity': -1, 'rt_peak': -1,
                          'lwhm': -1,
                          'rwhm': -1,
                          '5p_noise': -1,
                          '10p_noise': -1,
                          'SNR': -1,
                          'log_L_R': -1,
                          'log_int': -1})
    pnoise_5 = np.percentile(data_xic[(data_xic['rt'] > (
        time_w - h_rt_w)) & (data_xic['rt'] < (time_w + h_rt_w))]['intensity'], 5)
    pnoise_10 = np.percentile(data_xic[(data_xic['rt'] > (
        time_w - h_rt_w)) & (data_xic['rt'] < (time_w + h_rt_w))]['intensity'], 10)
    # find the lwhm and rwhm
    time_point = compute_log_LR(data_xic, pos_p[0], val_max, 0.5)
    if (time_point[0] * time_point[1] == 1) or (time_point[0] * time_point[1] < 0):
        # Try a second time FWHM  computation with 0.7 * max intensity
        time_point = compute_log_LR(data_xic, pos_p[0], val_max, 0.7)

    if time_point[0] == -1 or time_point[1] == -1:
        # keep the shape measure to -1 in case on txo point are -1
        log_L_R = -1
    else:
        log_L_R = np.log2(
            abs(time_w - time_point[0]) / abs(time_w - time_point[1]))

    if pnoise_5 == 0 and pnoise_10 > 0:
        SNR = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_10)
    else:
        if pnoise_5 != 0:
            SNR = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_5)
        else:
            log.info('\t 5 percentile is %4.4f (added 0.5)', pnoise_5)
            SNR = 20 * \
                np.log10(data_xic.ix[pos_p, 1].values / (pnoise_5 + 0.5))

    return pd.Series({'intensity': val_max[0], 'rt_peak': data_xic.ix[pos_p, 0].values[0] * 60,
                      'lwhm': time_point[0],
                      'rwhm': time_point[1],
                      '5p_noise': pnoise_5,
                      '10p_noise': pnoise_10,
                      'SNR': SNR[0],
                      'log_L_R': log_L_R,
                      'log_int': np.log2(val_max)[0]})


# def estimate_parameter( df , name_file, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output,  rt_list , id_list, moff_pride_flag ,ptm_map,   log,sample_size, quantile_value, match_filter_flag   ):

def estimate_parameter(df, name_file, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output, rt_list, id_list, moff_pride_flag, ptm_map, sample_size, quantile_value, match_filter_flag, log_file,num_CPU):
    """
    Compute the quality metrics for the filtering during the estimation part

    :param df:
    :param name_file:
    :param raw_name:
    :param tol:
    :param h_rt_w:
    :param s_w:
    :param s_w_match:
    :param loc_raw:
    :param loc_output:
    :param rt_list:
    :param id_list:
    :param moff_pride_flag:
    :param ptm_map:
    :param sample_size:
    :param quantile_value:
    :param match_filter_flag:
    :param log_file:
    :param num_CPU:
    :return:
    """
    set_logger(log_file)
    myPool = multiprocessing.Pool(num_CPU)
    sample = df[df['matched'] == 0].sample(frac=sample_size)
    log.critical(
        'quality measures estimation  using %r  MS2 ident. peptides randomly sampled' % sample.shape[0])
    data_split = np.array_split(sample,num_CPU)
    result = {}
    offset = 0
    # run matchinf filtering for
    #for result in data_split:
    for df_index in range(0, len(data_split)):
        result[df_index] = myPool.apply_async(apex_multithr, args=(data_split[df_index], name_file, raw_name, tol, h_rt_w, s_w, s_w_match,
                                                                   loc_raw, loc_output, offset, rt_list, id_list, moff_pride_flag, ptm_map, 1, -1, -1, match_filter_flag, log_file))
        offset += len(data_split[df_index])
    myPool.close()
    myPool.join()
    ms2_data = save_moff_apex_result(result)
    # log.critical ('Estimated distribution rank correlation exp. int. vs theor. int. %r %r %r ' %( ms2_data[ ms2_data['rankcorr'] != -1]['rankcorr'].quantile(0.25), ms2_data[ms2_data['rankcorr'] != -1]['rankcorr'].quantile(0.50),   ms2_data[ ms2_data['rankcorr'] != -1]['rankcorr'].quantile(0.75) )  )
    log.critical('MAD retention time along all isotope %r',
                 ms2_data[ms2_data['RT_drift'] != -1]['RT_drift'].describe())
    log.critical('Estimated distribition ratio exp. int. left isotope vs. monoisotopic isotope %r ',
                 ms2_data[ms2_data['delta_log_int'] != -1]['delta_log_int'].describe())
    error_relInr = ms2_data[ms2_data['Erro_RelIntensity_TheoExp'] != -1]['Erro_RelIntensity_TheoExp'].quantile(quantile_value)
    rt_drift = ms2_data[ms2_data['RT_drift'] != -1]['RT_drift'].quantile(quantile_value)
    ratio_log_int = ms2_data[ms2_data['delta_log_int'] != -1]['delta_log_int'].quantile(quantile_value)
    return (rt_drift, error_relInr, ratio_log_int)


def compute_match_peak_quality_measure(input_data, moff_pride_flag, log):
    """
    Compute filter quality metrics

    :param input_data:
    :param moff_pride_flag:
    :param log:
    :return:
    """
    sum_intensity =  input_data['intensity'].sum()
    mad_diff_int = np.mean(abs((input_data['intensity'] / sum_intensity) - (
        input_data['ratio_iso'] / input_data['ratio_iso'].sum())))
    rank_spearman = spearmanr(
        (input_data['intensity'] / sum_intensity ), input_data['ratio_iso'])[0]
    mad_rt = np.mean(abs(input_data['rt_peak'] - input_data['rt_peak'].mean()))
    return (mad_diff_int, rank_spearman, mad_rt)


def estimate_on_match_peak(x, input_data, estimate_flag, moff_pride_flag, log, thr_q2, err_ratio_int, xic_data, mbr_flag, h_rt_w, s_w, s_w_match, offset_index):
    """
    Estimation of filter quality measures based on sampling of the MS2 identified peptides.

    :param x:
    :param input_data:
    :param estimate_flag:
    :param moff_pride_flag:
    :param log:
    :param thr_q2:
    :param err_ratio_int:
    :param xic_data:
    :param mbr_flag:
    :param h_rt_w:
    :param s_w:
    :param s_w_match:
    :param offset_index:
    :return:
    """
    test = input_data.loc[input_data['original_ptm'] == x.name, :].copy()
    test.reset_index(inplace=True)
    # print 'local df inside estimate ', input_data.columns
    test.iloc[0:1, 13:22] = test.iloc[0:1, :].apply(lambda x: compute_peak_simple(
        x, xic_data, log, mbr_flag, h_rt_w, s_w, s_w_match, offset_index, moff_pride_flag, -1, 1, 0), axis=1)
    # print 'output -->>  ',input_data.iloc[:,12:22]
    # print  input_data.iloc[0, input_data.columns.get_indexer(['log_L_R'])].all() != -1
    if (test.iloc[0, test.columns.get_indexer(['log_L_R'])]).any() != -1:
        #if not moff_pride_flag :
         #   new_point = test.iloc[0, test.columns.get_indexer(['rt_peak'])]
        #else:
            # to minute - second
            # moffpride data -> convert again in second
        new_point = test.iloc[0,
                                  test.columns.get_indexer(['rt_peak'])] / 60
        test.iloc[1:4, 13:22] = test.iloc[1:4, :].apply(lambda x: compute_peak_simple(
            x, xic_data, log, mbr_flag, h_rt_w, 0.3, 0.3, offset_index, moff_pride_flag, new_point[0], 1, 0), axis=1)
        if (test.iloc[0:3, test.columns.get_indexer(['log_L_R'])] != -1).all()[0]:
            mad_diff_int, rank_spearman, mad_rt = compute_match_peak_quality_measure(
                test.iloc[0:3, :], moff_pride_flag, log)
            if (test.iloc[3, test.columns.get_indexer(['log_L_R'])]).all() == -1:
                return pd.Series({'Erro_RelIntensity_TheoExp': mad_diff_int, 'rankcorr': rank_spearman, 'RT_drift': mad_rt, 'delta_rt': -1, 'delta_log_int': -1})
            else:
                delta_rt_wrong_iso = abs(
                    test.at[3, 'rt_peak'] - test.iloc[0:3, test.columns.get_indexer(['rt_peak'])].mean()[0])
                delta_log_int = test.at[3, 'log_int'] / test.at[0, 'log_int']
                # print pd.Series({'Erro_RelIntensity_TheoExp': mad_diff_int, 'rankcorr': rank_spearman,'RT_drift': mad_rt ,'delta_rt': delta_rt_wrong_iso ,'delta_log_int': delta_log_int})
                return pd.Series({'Erro_RelIntensity_TheoExp': mad_diff_int, 'rankcorr': rank_spearman, 'RT_drift': mad_rt, 'delta_rt': delta_rt_wrong_iso, 'delta_log_int': delta_log_int})
        else:
            return pd.Series({'Erro_RelIntensity_TheoExp': -1, 'rankcorr': -1, 'RT_drift': -1, 'delta_rt': -1, 'delta_log_int': -1})
    else:
        return pd.Series({'Erro_RelIntensity_TheoExp': -1, 'rankcorr': -1, 'RT_drift': -1, 'delta_rt': -1, 'delta_log_int': -1})


def filtering_match_peak(x, input_data, estimate_flag, moff_pride_flag, log, thr_q2, err_ratio_int, xic_data, mbr_flag, h_rt_w, s_w, s_w_match, offset_index):
    """
    Filtering of the matched peptides based on the isotopic envelope and quality measures estimated
    :param x:
    :param input_data:
    :param estimate_flag:
    :param moff_pride_flag:
    :param log:
    :param thr_q2:
    :param err_ratio_int:
    :param xic_data:
    :param mbr_flag:
    :param h_rt_w:
    :param s_w:
    :param s_w_match:
    :param offset_index:
    :return:
    """
    # print 'inside filtering routine ...'
    # log.info('matched peptide  --> %r  mZ: %4.4f RT: %4.4f ', x.mod_peptide , x.mz, x.rt)
    test = input_data.loc[input_data['original_ptm'] == x.name, :].copy()
    test.reset_index(inplace=True)
    test.iloc[0:1, 13:22] = test.iloc[0:1, :].apply(lambda x: compute_peak_simple(
        x, xic_data, log, mbr_flag, h_rt_w, s_w, s_w_match, offset_index, moff_pride_flag, -1, 1, 0), axis=1)
    if (test.iloc[0, test.columns.get_indexer(['log_L_R'])]).all() != -1:
        #if not moff_pride_flag :
        #    new_point = test.iloc[0, test.columns.get_indexer(['rt_peak'])]
        #else:
            # to minute - second
            # moffpride data -> convert again in second
        # from the ssecond isotope always convert to minute case : if new_point is provided
        new_point = test.iloc[0, test.columns.get_indexer(['rt_peak'])] / 60
        test.iloc[1:4, 13:22] = test.iloc[1:4, :].apply(lambda x: compute_peak_simple(
            x, xic_data, log, mbr_flag, h_rt_w, 0.3, 0.3, offset_index, moff_pride_flag, new_point[0], 1, 0), axis=1)
        # check isotope 2-3
        if (test.iloc[1:3, test.columns.get_indexer(['log_L_R'])] != -1).all()[0]:
            mad_diff_int, rank_spearman, mad_rt = compute_match_peak_quality_measure(
                test.iloc[0:3, :], moff_pride_flag, log)
            if (mad_rt < thr_q2 and rank_spearman > 0.8):
                # check isotope -1
                if (test.iloc[3, test.columns.get_indexer(['log_L_R'])]).all() != -1:
                    delta_rt_wrong_iso = abs(
                        test.at[3, 'rt_peak'] - test.iloc[0:3, test.columns.get_indexer(['rt_peak'])].mean()[0])
                    delta_log_int = test.at[3,
                                            'log_int'] / test.at[0, 'log_int']
                    if (delta_rt_wrong_iso < thr_q2) and (delta_log_int > err_ratio_int):
                        # filter  overlapping peptide isotope
                        log.info(' %r mz: %4.4f RT: %4.4f --> Not valid isotope envelope  overlapping detected -->  --  MAD RT  %r  -- rankCorr %r ',
                                 x.mod_peptide, x.mz, x.rt, mad_rt, rank_spearman)
                        return pd.Series({'intensity': -1, 'rt_peak': -1, 'lwhm': -1, 'rwhm': -1, '5p_noise': -1, '10p_noise': -1, 'SNR': -1, 'log_L_R': -1, 'log_int': -1})
                    else:
                        log.info('%r mz: %4.4f RT: %4.4f --> Valid isotope envelope detected after overlapping check -->  --  MAD RT  %r  -- rankCorr %r ',
                                 x.mod_peptide, x.mz, x.rt, mad_rt, rank_spearman)
                        return test.loc[test['ratio_iso'].idxmax(axis=1), ['10p_noise', '5p_noise', 'SNR', 'intensity', 'log_L_R', 'log_int', 'lwhm', 'rt_peak', 'rwhm']]
                else:
                    log.info(' %r mz: %4.4f RT: %4.4f  --> Valid isotope envelope detected and no overlaping detected -->  --  MAD RT  %r  -- rankCorr %r ',
                             x.mod_peptide, x.mz, x.rt, mad_rt, rank_spearman)
                    return test.loc[test['ratio_iso'].idxmax(axis=1), ['10p_noise', '5p_noise', 'SNR', 'intensity', 'log_L_R', 'log_int', 'lwhm', 'rt_peak', 'rwhm']]
            else:
                # not pass the thr. control
                log.info(' %r mz: %4.4f RT: %4.4f  --> Not valid isotope envelope detected  -->  --  MAD RT  %r  -- rankCorr %r ',
                         x.mod_peptide, x.mz, x.rt, mad_rt, rank_spearman)
                return pd.Series({'intensity': -1, 'rt_peak': -1, 'lwhm': -1, 'rwhm': -1, '5p_noise': -1, '10p_noise': -1, 'SNR': -1, 'log_L_R': -1, 'log_int': -1})
        else:
            # I have only the 1st valid isotope peak  but not the second or third
            log.info(' %r mz: %4.4f RT: %4.4f --> not enough isotope peaks detected  ',
                     x.mod_peptide, x.mz, x.rt)
            return pd.Series({'intensity': -1, 'rt_peak': -1, 'lwhm': -1, 'rwhm': -1, '5p_noise': -1, '10p_noise': -1, 'SNR': -1, 'log_L_R': -1, 'log_int': -1})
    else:
        log.info(' %r mz: %4.4f RT: %4.4f  --> first isotope peak not detected ',
                 x.mod_peptide, x.mz, x.rt)
        return pd.Series({'intensity': -1, 'rt_peak': -1, 'lwhm': -1, 'rwhm': -1, '5p_noise': -1, '10p_noise': -1, 'SNR': -1, 'log_L_R': -1, 'log_int': -1})


def apex_multithr(data_ms2, name_file, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output, offset_index, rt_list, id_list, moff_pride_flag, ptm_map, estimate_flag, rt_drift, err_ratio_int, match_filter_flag, log_file):
    """

    General apex method used both for filtering and not filtering usage

    :param data_ms2:
    :param name_file:
    :param raw_name:
    :param tol:
    :param h_rt_w:
    :param s_w:
    :param s_w_match:
    :param loc_raw:
    :param loc_output:
    :param offset_index:
    :param rt_list:
    :param id_list:
    :param moff_pride_flag:
    :param ptm_map:
    :param estimate_flag:
    :param rt_drift:
    :param err_ratio_int:
    :param match_filter_flag:
    :param log_file:
    :return:
    """
    set_logger(log_file)
    # setting flag and ptah
    flag_mzml = False
    flag_windows = False
    mbr_flag = False

    # set platform
    if _platform in ["linux", "linux2", 'darwin']:
        flag_windows = False
    elif _platform == "win32":
        flag_windows = True

    if loc_output != '':
        if not (os.path.isdir(loc_output)):
            os.makedirs(loc_output)
            log.info("created output folder: ", loc_output)

    if '_match' in name_file:
        # in case of mbr , here i dont have evaluate the flag mbr
        start = name_file.find('_match')
        # extract the name of the file
        name_file = name_file[0:start]

    if loc_raw is not None:
        if flag_windows:
            loc = os.path.join(loc_raw, name_file.upper() + '.RAW')

        else:
            # raw file name must have capitals letters :) this shloud be checked
            # this should be done in moe elegant way

            loc = os.path.normcase(os.path.join(loc_raw, name_file + '.RAW'))

            if not (os.path.isfile(loc)):
                loc = os.path.join(loc_raw, name_file + '.raw')

    else:
        # mzML work only with --raw_list option
        loc = raw_name
        if 'MZML' in raw_name.upper():
            flag_mzml = True

    if not (os.path.isfile(loc)):
        log.critical(
            'ERROR: Wrong path or wrong raw file name included: %s' % loc)
        return (None, -1)

    # index_offset = data_ms2.columns.shape[0] - 1
    data_ms2["intensity"] = -1
    data_ms2["rt_peak"] = -1
    data_ms2["lwhm"] = -1
    data_ms2["rwhm"] = -1
    data_ms2["5p_noise"] = -1
    data_ms2["10p_noise"] = -1
    data_ms2["SNR"] = -1
    data_ms2["log_L_R"] = -1
    data_ms2["log_int"] = -1
    data_ms2["rt_peak"] = data_ms2["rt_peak"].astype('float64')
    data_ms2['intensity'] = data_ms2['intensity'].astype('float64')
    data_ms2['lwhm'] = data_ms2['lwhm'].astype('float64')
    data_ms2["rwhm"] = data_ms2['rwhm'].astype('float64')
    data_ms2["5p_noise"] = data_ms2['5p_noise'].astype('float64')
    data_ms2["10p_noise"] = data_ms2['10p_noise'].astype('float64')
    data_ms2["SNR"] = data_ms2['SNR'].astype('float64')
    data_ms2["log_L_R"] = data_ms2['log_L_R'].astype('float64')
    data_ms2["log_int"] = data_ms2['log_int'].astype('float64')
    if estimate_flag == 1:
        # add extra filed if I am in a estimate mode
        data_ms2["Erro_RelIntensity_TheoExp"] = -1
        data_ms2["rankcorr"] = -1
        data_ms2["RT_drift"] = -1
        data_ms2["delta_rt"] = -1
        data_ms2["delta_log_int"] = -1
    # set mbr_flag
    if 'matched' in data_ms2.columns:
        if (data_ms2['matched'] == 1).all():
            # case valid in case of  filtering
            mbr_flag = True
        else:
            if (data_ms2['matched'] == 0).all():
                 # case valiD for estimation
                mbr_flag = False
            else:
                # case valid in case of not filtering
                 mbr_flag = True
    # get txic path: assumes txic is in the same directory as moff.py
    txic_executable_name = "txic_json.exe"
    txic_path = os.path.join(os.path.dirname(
        os.path.realpath(sys.argv[0])), txic_executable_name)
    # for all the input peptide in data_ms2
    try:
        if match_filter_flag == 1:
            all_isotope_df = build_matched_modification(
                data_ms2, ptm_map, tol, moff_pride_flag, h_rt_w)
            xic_data = get_xic_data(flag_mzml, flag_windows, all_isotope_df[[
                                    'mz', 'tol', 'ts', 'te']], loc_output, name_file, txic_path, loc, 1)
            # new filtering
            # not needed
            all_isotope_df['prog_xic_index'] = range(0, len(xic_data))
            all_isotope_df['original_ptm'] = np.repeat(data_ms2.index, 4)
            all_isotope_df["intensity"] = -1
            all_isotope_df["rt_peak"] = -1
            all_isotope_df["lwhm"] = -1
            all_isotope_df["rwhm"] = -1
            all_isotope_df["5p_noise"] = -1
            all_isotope_df["10p_noise"] = -1
            all_isotope_df["SNR"] = -1
            all_isotope_df["log_L_R"] = -1
            all_isotope_df["log_int"] = -1
            if estimate_flag == 0:
                data_ms2[['10p_noise', '5p_noise', 'SNR', 'intensity', 'log_L_R', 'log_int', 'lwhm', 'rt_peak', 'rwhm']] = data_ms2.apply(lambda x: filtering_match_peak(
                    x, all_isotope_df, estimate_flag, moff_pride_flag, log, rt_drift, err_ratio_int, xic_data, mbr_flag, h_rt_w, s_w, s_w_match, offset_index), axis=1)
            else:
                data_ms2[['Erro_RelIntensity_TheoExp', 'RT_drift', 'delta_log_int', 'delta_rt', 'rankcorr']] = data_ms2.apply(lambda x: estimate_on_match_peak(
                    x, all_isotope_df, estimate_flag, moff_pride_flag, log, rt_drift, err_ratio_int, xic_data, mbr_flag, h_rt_w, s_w, s_w_match, offset_index), axis=1)
            if estimate_flag != 1:
                data_ms2 = data_ms2[(data_ms2[['10p_noise', '5p_noise', 'SNR', 'intensity',
                                               'log_L_R', 'log_int', 'lwhm', 'rt_peak', 'rwhm']] != -1).all(1)]
        else:
            # not match  filter
            temp = data_ms2[['mz', 'rt']].copy()  # strange cases
            temp.ix[:, 'tol'] = int(tol)
            if moff_pride_flag == 1:
                temp['ts'] = (data_ms2['rt']) - h_rt_w
                temp['te'] = (data_ms2['rt']) + h_rt_w
            else:
                temp['ts'] = (data_ms2['rt'] / 60) - h_rt_w
                temp['te'] = (data_ms2['rt'] / 60) + h_rt_w
            temp.drop('rt', 1, inplace=True)
            xic_data = get_xic_data(
                flag_mzml, flag_windows, temp, loc_output, name_file, txic_path, loc, 0)
            data_ms2.reset_index(inplace=True)
            data_ms2[['10p_noise', '5p_noise', 'SNR', 'intensity', 'log_L_R', 'log_int', 'lwhm', 'rt_peak', 'rwhm']] = data_ms2.apply(
                lambda x: compute_peak_simple(x, xic_data, log, mbr_flag, h_rt_w, s_w, s_w_match, offset_index, moff_pride_flag, -1, -1, 1), axis=1)

    except Exception as e:
        traceback.print_exc()
        print
        raise e
    return (data_ms2, 1)


def build_matched_modification(data, ptm_map, tol, moff_pride_flag, h_rt_w):
    """
    Computation of th. isotopic envelope tanking into account PSM modification
    :param data:
    :param ptm_map:
    :param tol:
    :param moff_pride_flag:
    :param h_rt_w:
    :return:
    """
    all_isotope_df = pd.DataFrame(
        columns=['peptide', 'mz', 'ratio_iso', 'tol', 'rt', 'matched', 'ts', 'te'])
    for row in data.itertuples():
        # get the sequence
        # for MQ sequence is (mod_tag )
        # for PS sequence is  <mod_tag>
        mq_mod_flag = False
        if mq_mod_flag:
            if not ('(' in row.mod_peptide) and mq_mod_flag:
                #  only fixed mod
                comps = Counter(
                    list(chain(*[list(std_aa_comp[aa].elements()) for aa in row.peptide])))
                comps["H"] += 2
                comps["O"] += 1
                fix_mod_count = row.peptide.count('C')
                if fix_mod_count > 0:
                    comps["H"] += (ptm_map['cC']['deltaChem']
                                   [0] * fix_mod_count)
                    comps["C"] += (ptm_map['cC']['deltaChem']
                                   [1] * fix_mod_count)
                    comps["N"] += (ptm_map['cC']['deltaChem']
                                   [2] * fix_mod_count)
                    comps["O"] += (ptm_map['cC']['deltaChem']
                                   [3] * fix_mod_count)
            else:
                comps = Counter(
                    list(chain(*[list(std_aa_comp[aa].elements()) for aa in row.peptide])))
                for ptm in ptm_map.keys():
                    ptm_c = row.mod_peptide.count(ptm)
                    if ptm_c >= 1:
                        comps["H"] += (ptm_map[ptm]['deltaChem'][0] * ptm_c)
                        comps["C"] += (ptm_map[ptm]['deltaChem'][1] * ptm_c)
                        comps["N"] += (ptm_map[ptm]['deltaChem'][2] * ptm_c)
                        comps["O"] += (ptm_map[ptm]['deltaChem'][3] * ptm_c)
            # add eventually fixed mod/
                fix_mod_count = row.mod_peptide.count('C')
                if fix_mod_count > 0:
                    comps["H"] += (ptm_map['cC']['deltaChem']
                                   [0] * fix_mod_count)
                    comps["C"] += (ptm_map['cC']['deltaChem']
                                   [1] * fix_mod_count)
                    comps["N"] += (ptm_map['cC']['deltaChem']
                                   [2] * fix_mod_count)
                    comps["O"] += (ptm_map['cC']['deltaChem']
                                   [3] * fix_mod_count)
                comps["H"] += 2
                comps["O"] += 1
        else:
            # fixed and variable mod are both in the sequence
            comps = Counter(
                list(chain(*[list(std_aa_comp[aa].elements()) for aa in row.peptide])))
            if '<' in row.mod_peptide :
                # check only if modificatio are present.
                # for the future use dthe tag_mod_sequence_delimiter use in moFF_setting
                for ptm in ptm_map.keys():
                    ptm_c = row.mod_peptide.count(ptm)
                # ptm_c =  sum(ptm in s for s in row.mod_peptide)
                    if ptm_c >= 1:
                        comps["H"] += (ptm_map[ptm]['deltaChem'][0] * ptm_c)
                        comps["C"] += (ptm_map[ptm]['deltaChem'][1] * ptm_c)
                        comps["N"] += (ptm_map[ptm]['deltaChem'][2] * ptm_c)
                        comps["O"] += (ptm_map[ptm]['deltaChem'][3] * ptm_c)
            comps["H"] += 2
            comps["O"] += 1

        theoretical_isotopic_cluster = isotopic_variants(
            comps, charge=int(round(row.mass / float(row.mz))), npeaks=3)
        mz_iso = [peak.mz for peak in theoretical_isotopic_cluster]
        delta = mz_iso[0] - mz_iso[1]
        mz_iso.append(mz_iso[0] + delta)
        ratio_iso = [peak.intensity for peak in theoretical_isotopic_cluster]
        ratio_iso.append(-1)
        isotopic_df = pd.DataFrame({'mz': mz_iso, 'ratio_iso': ratio_iso})

        isotopic_df.loc[:, 'exp_mz'] = row.mz
        isotopic_df.loc[:, 'peptide'] = row.mod_peptide
        isotopic_df.loc[:, 'tol'] = int(tol)
        isotopic_df.loc[:, 'rt'] = row.rt
        isotopic_df.loc[:, 'matched'] = 1
        if moff_pride_flag :
            #moffpridedata  rt is in minutes
            isotopic_df['ts'] = (row.rt) - h_rt_w
            isotopic_df['te'] = (row.rt) + h_rt_w
        else:
            # not moffpridedata rt in second
            isotopic_df['ts'] = (row.rt / 60) - h_rt_w
            isotopic_df['te'] = (row.rt / 60) + h_rt_w

        all_isotope_df = pd.concat(
            [all_isotope_df, isotopic_df], join='outer', axis=0, sort=False)
    all_isotope_df.reset_index(inplace=True)

    return all_isotope_df


def get_xic_data(flag_mzml, flag_windows, data, loc_output, name_file, txic_path, loc, flag_filtering):
    """

    Run the txic_json.xe library to get all the xic requested by each process for Thermo raw file

    :param flag_mzml:
    :param flag_windows:
    :param data:
    :param loc_output:
    :param name_file:
    :param txic_path:
    :param loc:
    :param flag_filtering:
    :return:
    """
    if not flag_mzml:
        # txic-28-9-separate-jsonlines.exe
        if not flag_windows:
            # Linux  to avoid cmd  string too long  and its error. the thresold is mainly base on  from empirical evaluation.
            if len(data.to_json(orient='records')) >= 50000:
                with open(os.path.join(loc_output, multiprocessing.current_process().name + '_' + name_file + '.json'), 'w') as f:
                    f.write(data.to_json(orient='records'))
                args_txic = shlex.split("mono " + txic_path + " -jf " + os.path.join(
                    loc_output, multiprocessing.current_process().name + '_' + name_file + '.json') + " -f " + loc, posix=False)
            else:
                    #  small amount of char. in the request
                args_txic = shlex.split(
                    "mono " + txic_path + " -j " + data.to_json(orient='records') + " -f " + loc, posix=True)
        else:
            # Windows to avoid cmd  string too long  and its error. the thresold is mainly base on  from empirical evaluation.
            if len(data.to_json(orient='records')) >= 10000:
                with open(os.path.join(loc_output, multiprocessing.current_process().name + '_' + name_file + '.json'), 'w') as f:
                    f.write(data.to_json(orient='records'))
                args_txic = shlex.split(txic_path + " -jf " + os.path.join(
                    loc_output, multiprocessing.current_process().name + '_' + name_file + '.json') + " -f " + loc, posix=False)
            else:
                #  small amount of char. in the request
                args_txic = shlex.split(
                    txic_path + " -j " + data.to_json(orient='records') + " -f " + loc, posix=False)
        # start_timelocal = time.time()
        p = subprocess.Popen(args_txic, stdout=subprocess.PIPE)
        output, err = p.communicate()
        xic_data = []
        for l in range(0, data.shape[0]):
            temp = json.loads(output.split('\n')[l].decode("utf-8"))
            xic_data.append(pd.DataFrame(
                {'rt': temp['results']['times'], 'intensity': temp['results']['intensities']}, columns=['rt', 'intensity']))
    else:
        # this does not seem to work yet.
        pass
        # run_temp = pymzml.run.Reader(raw_name)
        # xic_data = mzML_get_all(temp, tol, loc, run_temp, rt_list, id_list)

    return xic_data
