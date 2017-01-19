#!/usr/bin/env python

import ConfigParser
import StringIO
import argparse
import ast
import bisect
import logging
import os as os
import shlex
import subprocess
import sys
import time
from sys import platform as _platform

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

"""
 input
   - MS2 ID file
   - tol
   - half rt time window in minute
 output
   - list of intensities..+
"""

TXIC_PATH = os.environ.get('TXIC_PATH', './')


def map_ps2moff(data):
    data.drop(data.columns[[0]], axis=1, inplace=True)
    data.columns = data.columns.str.lower()
    data.rename(
        columns={'sequence': 'peptide', 'measured charge': 'charge', 'theoretical mass': 'mass', 'protein(s)': 'prot',
                 'm/z': 'mz'}, inplace=True)
    return data, data.columns.values.tolist()


'''
input list of columns
list of column names from PS default template loaded from .properties
'''


def check_ps_input_data(input_column_name, list_col_ps_default):
    res = [i for e in list_col_ps_default for i in input_column_name if e in i]
    if len(res) == len(input_column_name):
        # detected a default PS input file
        return 1
    else:
        # not detected a default PS input file
        return 0


def check_columns_name(col_list, col_must_have):
    for c_name in col_must_have:
        if not (c_name in col_list):
            # fail
            print 'The following filed name is missing or wrong: ', c_name
            return 1
    # succes
    return 0


def pyMZML_xic_out(name, ppmPrecision, minRT, maxRT, MZValue):
    run = pymzml.run.Reader(name, MS1_Precision=ppmPrecision, MSn_Precision=ppmPrecision)
    timeDependentIntensities = []
    for spectrum in run:
        if spectrum['ms level'] == 1 and spectrum['scan start time'] > minRT and spectrum['scan start time'] < maxRT:
            lower_index = bisect.bisect(spectrum.peaks, (float(MZValue - ppmPrecision * MZValue), None))
            upper_index = bisect.bisect(spectrum.peaks, (float(MZValue + ppmPrecision * MZValue), None))
            print lower_index, upper_index
            maxI = 0.0
            for sp in spectrum.peaks[lower_index: upper_index]:
                if sp[1] > maxI:
                    maxI = sp[1]
            if maxI > 0:
                timeDependentIntensities.append([spectrum['scan start time'], maxI])

    if len(timeDependentIntensities) > 5:
        return (pd.DataFrame(timeDependentIntensities, columns=['rt', 'intensity']), 1)
    else:
        return (pd.DataFrame(timeDependentIntensities, columns=['rt', 'intensity']), -1)


def run_apex(file_name, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output):
    # OS detect
    flag_windows = False
    if _platform in ["linux", "linux2", 'darwin']:
        flag_windows = False
    elif _platform == "win32":
        flag_windows = True

    # flag_for matching
    mbr_flag = 0
    config = ConfigParser.RawConfigParser()
    # get the  running path of moff
    moff_path = os.path.dirname(sys.argv[0])

    # it s always placed in same folder of moff.py

    config.read(os.path.join(moff_path, 'moff_setting.properties'))

    # case of moff_all more than one subfolderi
    name = os.path.basename(file_name).split('.')[0]
    if '_match' in name:
        # in case of mbr , here i dont have evaluate the flag mbr
        start = name.find('_match')
        # extract the name of the file
        name = name[0:start]

    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    log.addHandler(ch)

    if loc_output != '':
        if not (os.path.isdir(loc_output)):
            os.makedirs(loc_output)
            log.info("created output folder: ", loc_output)

        # outputname : name of the output
        # it should be ok also in linux
        outputname = os.path.join(loc_output, name + "_moff_result.txt")
        fh = logging.FileHandler(os.path.join(loc_output, name + '__moff.log'), mode='w')
    else:
        outputname = name + "_moff_result.txt"
        fh = logging.FileHandler(os.path.join(name + '__moff.log'), mode='w')

    fh.setLevel(logging.INFO)
    log.addHandler(fh)
    flag_mzml = False
    if loc_raw is not None:
        if flag_windows:
            loc = os.path.join(loc_raw, name + '.RAW')
        else:
            # raw file name must have capitals letters :) this shloud be checked
            loc = os.path.join(loc_raw, name.upper() + '.RAW')
    else:
        ## only with specific file name I usu mzML file :::///
        ## I have already the full location
        loc = raw_name
        if ('MZML' in raw_name.upper()):
            flag_mzml = True

    if os.path.isfile(loc):
        log.info('raw file exist')
    else:
        exit('ERROR: Wrong path or wrong file name included: %s' % loc)

    log.critical('moff Input file: %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ' % (file_name, tol, h_rt_w, s_w))
    log.critical('RAW file  :  %s' % (loc))
    log.info('Output_file in :  %s', outputname)
    log.info('RAW file location :  %s', loc)

    # read data from file
    data_ms2 = pd.read_csv(file_name, sep="\t", header=0)
    if not 'matched' in data_ms2.columns:
        # check if it is a PS file ,
        list_name = data_ms2.columns.values.tolist()
        # get the lists of PS  defaultcolumns from properties file
        list = ast.literal_eval(config.get('moFF', 'ps_default_export'))
        # here it controls if the input file is a PS export; if yes it maps the input in right moFF name
        if check_ps_input_data(list_name, list) == 1:
            # map  the columns name according to moFF input requirements
            data_ms2, list_name = map_ps2moff(data_ms2)
    if check_columns_name(data_ms2.columns.tolist(), ast.literal_eval(config.get('moFF', 'col_must_have_x'))) == 1:
        exit('ERROR minimal field requested are missing or wrong')

    index_offset = data_ms2.columns.shape[0] - 1

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

    # set mbr_flag
    if 'matched' in data_ms2.columns:
        mbr_flag = 1
        log.info('Apex module has detected mbr peptides')
        log.info('moff_rtWin_peak for matched peptide:   %4.4f ', s_w_match)
    c = 0
    log.critical('Starting apex .........')

    start_time = time.time()
    for index_ms2, row in data_ms2.iterrows():
        # log.info('peptide at line: %i',c)
        mz_opt = "-mz=" + str(row['mz'])
        #### PAY ATTENTION HERE , we assume that input RT is in second
        ## RT in Thermo file is in minutes
        #### if it is not the case change the following line
        time_w = row['rt'] / 60
        if mbr_flag == 0:
            log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f ', c, row['mz'], time_w)
            temp_w = s_w
        else:
            log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f matched(y/n): %i', c, row['mz'], time_w,
                     row['matched'])
            if row['matched'] == 1:
                temp_w = s_w_match
            else:
                temp_w = s_w
        if row['rt'] == -1:
            log.warning('rt not found. Wrong matched peptide in the mbr step line: %i', c)
            c += 1
            continue

        # convert rt to sec to min
        try:
            if flag_mzml:
                # mzml raw file
                # transform the tollerance in ppm
                data_xic, status = pyMZML_xic_out(loc, float(tol / (10 ** 6)), time_w - h_rt_w, time_w + h_rt_w,
                                                  row['mz'])

                if status == -1:
                    log.warning("WARNINGS: XIC not retrived line: %i", c)
                    log.warning('MZ: %4.4f RT: %4.4f Mass: %i', row['mz'], row['rt'], index_ms2)
                    c += 1
                    continue
            else:
                #  Thermo RAW file
                if flag_windows:
                    os.path.join('folder_name', 'file_name')
                    args_txic = shlex.split(
                        os.path.join(moff_path, "txic.exe") + " " + mz_opt + " -tol=" + str(tol) + " -t " + str(
                            time_w - h_rt_w) + " -t " + str(time_w + h_rt_w) + " " + loc, posix=False)
                else:
                    args_txic = shlex.split(TXIC_PATH + "txic " + mz_opt + " -tol=" + str(tol) + " -t " + str(
                        time_w - h_rt_w) + " -t " + str(
                        time_w + h_rt_w) + " " + loc)

                p = subprocess.Popen(args_txic, stdout=subprocess.PIPE)
                output, err = p.communicate()
                data_xic = pd.read_csv(StringIO.StringIO(output.strip()), sep=' ', names=['rt', 'intensity'], header=0)
            if data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].shape[0] >= 1:
                ind_v = data_xic.index
                pp = data_xic[data_xic["intensity"] ==
                              data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))][
                                  'intensity'].max()].index
                pos_p = ind_v[pp]
                if pos_p.values.shape[0] > 1:
                    log.warning(" RT gap for the time windows searched. Probably the ppm values is too small %i", c)
                    continue
                val_max = data_xic.ix[pos_p, 1].values
            else:
                log.info("LW_BOUND window  %4.4f", time_w - temp_w)
                log.info("UP_BOUND window %4.4f", time_w + temp_w)
                log.info(data_xic[(data_xic['rt'] > (time_w - +0.60)) & (data_xic['rt'] < (time_w + 0.60))])
                log.info("WARNINGS: moff_rtWin_peak is not enough to detect the max peak line : %i", c)
                log.info('MZ: %4.4f RT: %4.4f Mass: %i', row['mz'], row['rt'], index_ms2)
                c += 1
                continue
            pnoise_5 = np.percentile(
                data_xic[(data_xic['rt'] > (time_w - (h_rt_w / 2))) & (data_xic['rt'] < (time_w + (h_rt_w / 2)))][
                    'intensity'], 5)
            pnoise_10 = np.percentile(
                data_xic[(data_xic['rt'] > (time_w - (h_rt_w / 2))) & (data_xic['rt'] < (time_w + (h_rt_w / 2)))][
                    'intensity'], 10)
        except (IndexError, ValueError, TypeError):
            log.warning(" size is not enough to detect the max peak line : %i", c)
            log.info('MZ: %4.4f RT: %4.4f index: %i', row['mz'], row['rt'], index_ms2)
            continue
            c += 1
        except pd.parser.CParserError:
            log.warning("WARNINGS: XIC not retrived line: %i", c)
            log.warning('MZ: %4.4f RT: %4.4f Mass: %i', row['mz'], row['rt'], index_ms2)

            c += 1
            continue
        else:
            log_time = [-1, -1]
            c_left = 0
            find_5 = False
            stop = False
            while c_left < (pos_p - 1) and not stop:

                if not find_5 and (data_xic.ix[(pos_p - 1) - c_left, 1].values <= (0.5 * val_max)):
                    find_5 = True
                    log_time[0] = data_xic.ix[(pos_p - 1) - c_left, 0].values * 60
                    stop = True
                c_left += 1
            find_5 = False
            stop = False
            r_left = 0
            while ((pos_p + 1) + r_left < len(data_xic)) and not stop:
                if not find_5 and data_xic.ix[(pos_p + 1) + r_left, 1].values <= (0.50 * val_max):
                    find_5 = True
                    log_time[1] = data_xic.ix[(pos_p + 1) + r_left, 0].values * 60
                    stop = True
                r_left += 1

            data_ms2.ix[index_ms2, (index_offset + 1)] = val_max
            data_ms2.ix[index_ms2, (index_offset + 2)] = data_xic.ix[pos_p, 0].values * 60
            data_ms2.ix[index_ms2, (index_offset + 3)] = log_time[0]
            data_ms2.ix[index_ms2, (index_offset + 4)] = log_time[1]
            data_ms2.ix[index_ms2, (index_offset + 5)] = pnoise_5
            data_ms2.ix[index_ms2, (index_offset + 6)] = pnoise_10
            # conpute log_L_R SNR and log intensities
            if (pnoise_5 == 0 and pnoise_10 > 0):
                data_ms2.ix[index_ms2, (index_offset + 7)] = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_10)
            else:
                data_ms2.ix[index_ms2, (index_offset + 7)] = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_5)
            # WARNING  time - log_time 0 / time -log_time 1
            data_ms2.ix[index_ms2, (index_offset + 8)] = np.log2(
                abs(data_ms2.ix[index_ms2, index_offset + 2] - log_time[0]) / abs(
                    data_ms2.ix[index_ms2, index_offset + 2] - log_time[1]))
            data_ms2.ix[index_ms2, (index_offset + 9)] = np.log2(val_max)
            c += 1

    # save  result 
    log.critical('..............apex terminated')
    log.critical('Writing result in %s' % (outputname))
    log.info("--- Running time (measured when start the loop)  %s seconds ---" % (time.time() - start_time))
    # print time.time() - start_time
    data_ms2.to_csv(path_or_buf=outputname, sep="\t", header=True, index=False)
    fh.close()
    log.removeHandler(fh)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='moFF input parameter')

    parser.add_argument('--inputtsv', dest='name', action='store',
                        help='specify the input file with the MS2 peptides/features',
                        required=True)

    parser.add_argument('--inputraw', dest='raw_list', action='store',
                        help='specify directly raw file', required=False)

    parser.add_argument('--tol', dest='toll', action='store', type=float,
                        help='specify the tollerance parameter in ppm', required=True)

    parser.add_argument('--rt_w', dest='rt_window', action='store', type=float, default=3,
                        help='specify rt window for xic (minute). Default value is 3 min', required=False)

    parser.add_argument('--rt_p', dest='rt_p_window', action='store', type=float, default=0.2,
                        help='specify the time windows for the peak ( minute). Default value is 0.1 ', required=False)

    parser.add_argument('--rt_p_match', dest='rt_p_window_match', action='store', type=float, default=0.4,
                        help='specify the time windows for the matched  peak ( minute). Default value is 0.4 ',
                        required=False)

    parser.add_argument('--raw_repo', dest='raw', action='store', help='specify the raw file repository folder',
                        required=False)

    parser.add_argument('--output_folder', dest='loc_out', action='store', default='', help='specify the folder output',
                        required=False)

    args = parser.parse_args()

    if (args.raw_list is None) and (args.raw is None):
        exit('you must specify and raw files  with --inputraw (file name) or --raw_repo (folder)')
    if (args.raw_list is not None) and (args.raw is not None):
        exit('you must specify raw files using only one options --inputraw (file name) or --raw_repo (folder) ')

    file_name = args.name
    tol = args.toll
    h_rt_w = args.rt_window
    s_w = args.rt_p_window
    s_w_match = args.rt_p_window_match
    # if args.raw_list is not None:
    #    raw_list = args.raw_list
    # else:
    #    raw_list = None
    loc_raw = args.raw
    loc_output = args.loc_out
    # " init here the logger
    run_apex(file_name, args.raw_list, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output)
