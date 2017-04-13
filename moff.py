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
import multiprocessing



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


def save_moff_apex_result(list_df, result, folder_output, name):
    xx = []
    for df_index in range(0,len(list_df)):
        if result[df_index].get()[1] == -1:
            exit ('Raw file not retrieved: wrong path or upper/low case mismatch')
        else:
            xx.append( result[df_index].get()[0] )
    
    final_res = pd.concat(xx)
    #print final_res.shape
    # print os.path.join(folder_output,os.path.basename(name).split('.')[0]  + "_moff_result.txt")
    final_res.to_csv(os.path.join(folder_output, os.path.basename(name).split('.')[0] + "_moff_result.txt"), sep="\t",
                     index=False)

    return (1)



def map_ps2moff(data):
    data.drop(data.columns[[0]], axis=1, inplace=True)
    data.columns = data.columns.str.lower()
    data.rename(
        columns={'sequence': 'peptide','modified sequence':'mod_peptide' ,'measured charge': 'charge', 'theoretical mass': 'mass', 'protein(s)': 'prot',
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

def check_log_existence(file_to_check):
    if os.path.isfile(file_to_check):
        os.remove(file_to_check)




def apex_multithr(data_ms2,name_file, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output,offset_index):

    #setting logger for multiprocess
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    log.addHandler(ch)

    #setting flag and ptah
    moff_path = os.path.dirname(sys.argv[0])
    flag_mzml = False
    flag_windows = False
    mbr_flag = 0

    # set platform
    if _platform in ["linux", "linux2", 'darwin']:
        flag_windows = False
    elif _platform == "win32":
        flag_windows = True


    # check output log file in right location
    if loc_output != '':
        if not (os.path.isdir(loc_output)):
            os.makedirs(loc_output)
            log.info("created output folder: ", loc_output)

    # to be checked if it is works ffor both caseses
    fh = logging.FileHandler(os.path.join(loc_output, name_file + '__moff.log'), mode='a')

    fh.setLevel(logging.INFO)
    log.addHandler(fh)

    # check mbr input file
    if '_match' in name_file:
        # in case of mbr , here i dont have evaluate the flag mbr
        start = name_file.find('_match')
        # extract the name of the file
        name_file = name_file[0:start]

    if loc_raw is not None:
        if flag_windows:
           loc  = os.path.join(loc_raw, name_file.upper()+ '.RAW')

        else:
            # raw file name must have capitals letters :) this shloud be checked
            # this should be done in moe elegant way
            
            loc  = os.path.normcase(os.path.join(loc_raw, name_file + '.RAW'))
            
            if not (os.path.isfile(loc)):
                loc  = os.path.join(loc_raw, name_file + '.raw')

    else:
        #mzML work only with --inputraw option
        loc  = raw_name
        if ('MZML' in raw_name.upper()):
            flag_mzml = True

    #print loc 
    if os.path.isfile(loc):
        log.info('raw file exist')
    else:
        #exit('ERROR: Wrong path or wrong raw file name included: %s' % loc  )
        log.info('ERROR: Wrong path or wrong raw file name included: %s' % loc  )
        return (None,-1)


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
        #log.critical('Apex module has detected mbr peptides')
        #log.info('moff_rtWin_peak for matched peptide:   %4.4f ', s_w_match)

    c =0
    for index_ms2, row in data_ms2.iterrows():

        mz_opt = "-mz=" + str(row['mz'])
        #### PAY ATTENTION HERE , we assume that input RT is in second
        ## RT in Thermo file is in minutes
        #### if it is not the case change the following line
        time_w = row['rt'] / 60
        if mbr_flag == 0:
            #print 'xx here mbr_flag == 0'
            log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f',(offset_index +c +2), row['mz'], time_w)
            temp_w = s_w
        else:
            log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f matched (yes=1/no=0): %i',(offset_index +c +2), row['mz'], time_w,row['matched'])
                    # row['matched'])
            if row['matched'] == 1:
                temp_w = s_w_match
            else:
                temp_w = s_w
        if row['rt'] == -1:
            log.warning('rt not found. Wrong matched peptide in the mbr step line: %i', (offset_index +c +2))
            c += 1
            continue
        try:
            if flag_mzml:
                # mzml raw file
                # transform the tollerance in ppm
                data_xic, status = pyMZML_xic_out(loc, float(tol / (10 ** 6)), time_w - h_rt_w, time_w + h_rt_w,row['mz'])

                if status == -1:
                    log.warning("WARNINGS: XIC not retrived line: %i", c)
                    log.warning('MZ: %4.4f RT: %4.4f Mass: %i', row['mz'], row['rt'], index_ms2)
                    c += 1
                    continue
            else:
                #  Thermo RAW file
                if flag_windows:
                    os.path.join('folder_name', 'file_name')
                    args_txic = shlex.split(os.path.join(moff_path, "txic.exe") + " " + mz_opt + " -tol=" + str(tol) + " -t " + str( time_w - h_rt_w) + " -t " + str(time_w + h_rt_w) + " " + loc, posix=False)
                else:
                    args_txic = shlex.split(TXIC_PATH + "txic " + mz_opt + " -tol=" + str(tol) + " -t " + str(time_w - h_rt_w) + " -t " + str(time_w + h_rt_w) + " " + loc )
                
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
                    continue
                val_max = data_xic.ix[pos_p, 1].values
            else:
                log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f ', (offset_index +c +2), row['mz'], time_w)
                log.info("\t LW_BOUND window  %4.4f", time_w - temp_w)
                log.info("\t UP_BOUND window %4.4f", time_w + temp_w)
                #cosche log.info(data_xic[(data_xic['rt'] > (time_w - +0.60)) & (data_xic['rt'] < (time_w + 0.60))])
                log.info("\t WARNINGS: moff_rtWin_peak is not enough to detect the max peak ")
                c += 1
                continue
            pnoise_5 = np.percentile(
                data_xic[(data_xic['rt'] > (time_w - (h_rt_w / 2))) & (data_xic['rt'] < (time_w + (h_rt_w / 2)))][
                    'intensity'], 5)
            pnoise_10 = np.percentile(
                data_xic[(data_xic['rt'] > (time_w - (h_rt_w / 2))) & (data_xic['rt'] < (time_w + (h_rt_w / 2)))][
                    'intensity'], 10)
        except (IndexError, ValueError, TypeError):
            log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f ', (offset_index +c +2), row['mz'], time_w)
            log.warning("\t size is not enough to detect the max peak line : %i", c)
            continue
            c += 1
        except pd.parser.CParserError:
            log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f ', (offset_index +c +2),  row['mz'], time_w)
            log.warning("\t WARNINGS: XIC not retrived line:")
            #log.warning('MZ: %4.4f RT: %4.4f Mass: %i', row['mz'], row['rt'], index_ms2)

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
                if pnoise_5 != 0:
                    data_ms2.ix[index_ms2, (index_offset + 7)] = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_5)
                else:
                    log.info('\t 5 percentile is %4.4f (added 0.5)', pnoise_5)
                    data_ms2.ix[index_ms2, (index_offset + 7)] = 20 * np.log10(data_xic.ix[pos_p, 1].values / (pnoise_5 +0.5))
            # WARNING  time - log_time 0 / time -log_time 1
            data_ms2.ix[index_ms2, (index_offset + 8)] = np.log2(
                abs(data_ms2.ix[index_ms2, index_offset + 2] - log_time[0]) / abs(
                    data_ms2.ix[index_ms2, index_offset + 2] - log_time[1]))
            data_ms2.ix[index_ms2, (index_offset + 9)] = np.log2(val_max)
            c += 1


    return  (data_ms2,1)


def main_apex_alone():
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

    parser.add_argument('--rt_p', dest='rt_p_window', action='store', type=float, default=0.4,
                        help='specify the time windows for the peak ( minute). Default value is 0.4 ', required=False)

    parser.add_argument('--rt_p_match', dest='rt_p_window_match', action='store', type=float, default=0.8,
                        help='specify the time windows for the matched  peak ( minute). Default value is 0.8 ',
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

    loc_raw = args.raw
    loc_output = args.loc_out
    # set stream option for logger
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    log.addHandler(ch)

    config = ConfigParser.RawConfigParser()
    config.read(os.path.join(os.path.dirname(sys.argv[0]), 'moff_setting.properties'))

    df = pd.read_csv(file_name, sep="\t")
    ## check and eventually tranf for PS template
    if not 'matched' in df.columns:
        # check if it is a PS file ,
        list_name = df.columns.values.tolist()
        # get the lists of PS  defaultcolumns from properties file
        list = ast.literal_eval(config.get('moFF', 'ps_default_export'))
        # here it controls if the input file is a PS export; if yes it maps the input in right moFF name
        if check_ps_input_data(list_name, list) == 1:
            # map  the columns name according to moFF input requirements
            data_ms2, list_name = map_ps2moff(df)
    ## check if the field names are good
    if check_columns_name(df.columns.tolist(), ast.literal_eval(config.get('moFF', 'col_must_have_apex'))) == 1:
        exit('ERROR minimal field requested are missing or wrong')

    log.critical('moff Input file: %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ' % (file_name, tol, h_rt_w, s_w))
    if args.raw_list is None:
        log.critical('RAW file from folder :  %s' % loc_raw)
    else:
        log.critical('RAW file  :  %s' % args.raw_list)

    log.critical('Output file in :  %s', loc_output)

    data_split = np.array_split(df, multiprocessing.cpu_count())

    log.critical('Starting Apex for .....')
    #print 'Original input size', df.shape
    name = os.path.basename(file_name).split('.')[0]

    ##check the existencce of the log file before to go to multiprocess
    check_log_existence(os.path.join(loc_output, name + '__moff.log'))

    myPool = multiprocessing.Pool(multiprocessing.cpu_count())

    result = {}
    offset = 0
    start_time = time.time()
    for df_index in range(0, len(data_split)):
        result[df_index] = myPool.apply_async(apex_multithr, args=(
        data_split[df_index], name, args.raw_list, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output, offset))
        offset += len(data_split[df_index])

    myPool.close()
    myPool.join()

    log.critical('...apex terminated')
    log.critical('...apex module execution time %4.4f (sec)' , time.time() - start_time)
    save_moff_apex_result(data_split, result, loc_output, file_name)


if __name__ == '__main__':
    main_apex_alone()
