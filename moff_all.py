#!/usr/bin/env python

import argparse
import ast
import configparser
import json
import logging.config
import multiprocessing
import os
import sys
import time

import numpy as np
import pandas as pd

import moff
import moff_mbr
#import moff_mbr_BF
import moff_anchor_peptide

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

"""moFF: entry point where to run moFF and all its functions """

if __name__ == '__main__':

    multiprocessing.freeze_support()

    parser_1 = argparse.ArgumentParser(
        description='moFF match between run and apex module input parameter', add_help=False)

    parser_1.add_argument('--config_file', dest='config_file', action='store',
                          help='specify a moFF parameter file ', required=False)
    args, remaining_argv = parser_1.parse_known_args()
    if args.config_file:
        config = configparser.SafeConfigParser(allow_no_value=True)
        config.read([args.config_file])
        moFF_parameters = dict(config.items("moFF_parameters"))
        # check if loc_in  is set in the input file
        if not ('loc_in' in moFF_parameters.keys() and 'raw_repo' in moFF_parameters.keys()):
            moFF_parameters['tsv_list'] = moFF_parameters['tsv_list'].split(
                ' ')
        if not ('raw_repo' in moFF_parameters.keys()):
            moFF_parameters['raw_list'] = moFF_parameters['raw_list'].split(
                ' ')
        if not ('tol' in moFF_parameters.keys()):
            exit('you must specify the tollerance in the configuration file ')
        moFF_parameters['tol'] = float(moFF_parameters['tol'])
        moFF_parameters['xic_length'] = float(moFF_parameters['xic_length'])
        moFF_parameters['rt_peak_win'] = float(moFF_parameters['rt_peak_win'])
        moFF_parameters['rt_peak_win_match'] = float(
            moFF_parameters['rt_peak_win_match'])
        moFF_parameters['peptide_summary'] = True if moFF_parameters['peptide_summary'] == 'True' else False
        moFF_parameters['w_comb'] = True if moFF_parameters['w_comb'] == 'True' else False
        moFF_parameters['out_flag'] = True if moFF_parameters['out_flag'] == 'True' else False
        moFF_parameters['w_filt'] = float(moFF_parameters['w_filt'])
        moFF_parameters['quantile_thr_filtering'] = float(moFF_parameters['quantile_thr_filtering'])
        moFF_parameters['cpu_num'] = int(moFF_parameters['cpu'])
        moFF_parameters['sample_size'] = float(moFF_parameters['sample_size'])
        moFF_parameters['match_filter'] = True if moFF_parameters['match_filter'] != '' else False
    args_1, remaining_argv = parser_1.parse_known_args()

    parser = argparse.ArgumentParser(parents=[parser_1],
                                     description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter, )

    parser = argparse.ArgumentParser(
        description='moFF match between run and apex module input parameter')
    parser.add_argument('--loc_in', dest='loc_in', action='store',
                        help='specify the folder of the input MS2 peptide list files', required=False)

    parser.add_argument('--tsv_list', dest='tsv_list', action='store', nargs='*',
                        help='specify the mzid file as a list', required=False)

    parser.add_argument('--raw_list', dest='raw_list', action='store',
                        nargs='*', help='specify the raw file as a list', required=False)

    parser.add_argument('--sample', dest='sample', action='store',
                        help='specify witch replicated to use for mbr reg_exp are valid', required=False)

    parser.add_argument('--ext', dest='ext', action='store', default='txt',
                        help='specify the file extentention of the input like. Default value: txt', required=False)

    parser.add_argument('--log_label', dest='log_label', action='store', default='moFF',
                        help='a label name to use for the log file. Default value: moFF', required=False)

    parser.add_argument('--w_filt', dest='w_filt', action='store', default=2,
                        help='width value of the filter  k * mean(Dist_Malahobis). Default value: 2', required=False)

    parser.add_argument('--out_flag', dest='out_flag', action='store_true', default=True,
                        help='if set, outliers for rt time allignment are filtered. Default value: True',
                        required=False)

    parser.add_argument('--w_comb', dest='w_comb', action='store_true', default=False,
                        help='if set, RT model combination is weighted using traing model errors: Default value: False',
                        required=False)
    parser.add_argument('--tol', dest='tol', action='store', default=10, type=float,
                        help='specify the tollerance  parameter in ppm. Default value: 10', required=False)

    parser.add_argument('--xic_length', dest='xic_length', action='store', type=float, default=3,
                        help='specify rt window for xic (minutes). Default value: 3', required=False)

    parser.add_argument('--rt_peak_win', dest='rt_peak_win', action='store', type=float, default=1,
                        help='specify the time windows for the peak (minutes). Default value: 1', required=False)

    parser.add_argument('--rt_peak_win_match', dest='rt_peak_win_match', action='store', type=float, default=1.2,
                        help='specify the time windows for the matched peptide peak (minutes). Default value: 1.2 ',
                        required=False)

    parser.add_argument('--raw_repo', dest='raw_repo', action='store',
                        help='specify the raw file repository ', required=False)

    parser.add_argument('--loc_out', dest='loc_out', action='store', default='', help='specify the folder output',
                        required=False)

    parser.add_argument('--rt_feat_file', dest='rt_feat_file', action='store',
                        help='specify the file that contains the features to use in the match-between-run RT prediction ',
                        required=False)

    parser.add_argument('--peptide_summary', dest='peptide_summary', action='store_true', default=False,
                        help='if set, export a peptide intesity summary tab-delited file. Default value: False',
                        required=False)

    parser.add_argument('--tag_pepsum', dest='tag_pepsum', action='store', type=str, default='moFF_run',
                        help='a tag text used for peptide summary file name (peptide_summary_intensity_ + tag + .tab ). Default value: moFF_run ',
                        required=False)
    parser.add_argument('--match_filter', dest='match_filter', action='store_true', default=False,
                        help='If set, filtering on the matched peak is activated. Default value: False', required=False)
    parser.add_argument('--ptm_file', dest='ptm_file', action='store', default='ptm_setting.json',
                        help='name of json ptm file. default file ptm_setting.json ', required=False)
    parser.add_argument('--quantile_thr_filtering', dest='quantile_thr_filtering', action='store', type=float,
                        default=0.75,
                        help='quantile value used to compute the filtering threshold for the matched peak .Default value: 0.75',
                        required=False)
    parser.add_argument('--sample_size', dest='sample_size', action='store', type=float, default=0.20,
                        help='percentage of MS2 peptide used to estimated the threshold. Default value: 0.20',
                        required=False)

    parser.add_argument('--mbr', dest='mbr', action='store', type=str, default='on',
                        help='select the moFF workflow: on to run mbr + apex , off to run only apex, only to run obnly mbr. Default value: on   ',
                        required=False)

    parser.add_argument('--cpu', dest='cpu_num', action='store', type=int, default=0,
                        help='number of cpu. as default value it will detect automaticaly the CPU number in your machine.',
                        required=False)

    if args.config_file:
        # load from config file and load the remaining parametes
        parser.set_defaults(**moFF_parameters)
        args = parser.parse_args(remaining_argv)
    else:
        # normal case for the input parsing
        args = parser.parse_args()

    # init global logger
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    log.addHandler(ch)

    if args.tol is None:
        exit('you must specify the tollerance in ppm ')
    if (args.tsv_list is None) and (args.loc_in is None) and (args.raw_list is None) and (args.raw_repo is None):
        exit('you must specify the input and raw files ')
    if (args.tsv_list is not None) and (args.loc_in is not None) and (args.raw_list is not None) and (
            args.raw_repo is not None):
        exit('you must specify the input and raw files or using: --tsv_list and --raw_list or --loc_in and --raw_repo ')
    else:
        if ((args.tsv_list is None) and (args.raw_list is not None)) or (
                (args.tsv_list is not None) and (args.raw_list is None)):
            exit(
                'Missing information: using --tsv_list you must specify the raw file with --raw_list ')
        if ((args.loc_in is None) and (args.raw_repo is not None)) or (
                (args.loc_in is not None) and (args.raw_repo is None)):
            exit(
                'Missing information: using --loc_in you must specify the raw file with --raw_repo ')

    if args.loc_out != '':
        if not (os.path.isdir(args.loc_out)):
            os.makedirs(args.loc_out)
            log.critical("created output folder  %r", args.loc_out)

    # fixed variable number of split and also number of CPU presence in the macine
    # change this variable  with repset to the machine setting of the user

    if args.cpu_num > 0:
        num_CPU = args.cpu_num
    else:
        num_CPU = multiprocessing.cpu_count()

    # only mbr
    if 'only' in args.mbr:
        log.critical('starting matching between run module (mbr)')
        #res_state, output_list_loc = moff_anchor_peptide.run_anchor(args)
        res_state, output_list_loc = moff_mbr.run_mbr(args)
        if res_state == -1:
            exit('An error is occurred during the writing of the mbr file')
        else:
            log.critical('end matching between run module (mbr)')
            exit()

    if 'on' in args.mbr:
        log.critical('Matching between run module (mbr)')
        res_state, output_list_loc = moff_mbr.run_mbr(args)
        # --- debug version-- just to run skip the mbr in for special cases
        # res_state= 1
        # output_list_loc =[]
        # for item in os.listdir(args.loc_in):
        #    #log.critical(item)
        #    if os.path.isfile(os.path.join(args.loc_in, item)):
        #        if os.path.join(args.loc_in, item).endswith('.' + args.ext):
        #            mbr_list_loc.append(os.path.join(args.loc_in, item))
        if res_state == -1:
            exit('An error is occurred during the writing of the mbr file')
        if args.tsv_list is not None:
            # input list of raw and tsv file
            if len(args.tsv_list) != len(args.raw_list):
                exit(
                    'Error:  number of the input files is different from the number of raw files')
            # in case list of file as input , mbr_output is written in local folder
            folder = os.path.join('mbr_output')
        else:
            folder = os.path.join(args.loc_in, 'mbr_output')

        log.critical('Apex module... ')

    if 'off' in args.mbr:
        # put everython in mbr_loc
        output_list_loc = []
        if not (args.loc_in is None):
            for item in os.listdir(args.loc_in):
                # log.critical(item)
                if os.path.isfile(os.path.join(args.loc_in, item)):
                    if os.path.join(args.loc_in, item).endswith('.' + args.ext):
                        output_list_loc.append(os.path.join(args.loc_in, item))
        else:
            output_list_loc = args.tsv_list

    for c, file_name in enumerate(output_list_loc):
        name = os.path.basename(file_name).split('.')[0]
        moff.check_log_existence(os.path.join(
            args.loc_out, name + '__moff.log'))
        fh = logging.FileHandler(os.path.abspath(
            os.path.join(args.loc_out, name + '__moff.log')), mode='a')
        fh.setLevel(logging.DEBUG)

        log.addHandler(fh)

        log_file = os.path.join(args.loc_out, name + '__moff.log')
        tol = args.tol
        h_rt_w = args.xic_length
        s_w = args.rt_peak_win
        s_w_match = args.rt_peak_win_match

        if args.tsv_list is not None:
            # raw_list contains the current raws file provided by args.raw_list option
            raw_list = args.raw_list[c]
        else:
            raw_list = None

        loc_raw = args.raw_repo
        loc_output = args.loc_out

        config = configparser.RawConfigParser()
        config.read(os.path.join(os.path.dirname(
            os.path.realpath(sys.argv[0])), 'moff_setting.properties'))
        df = pd.read_csv(file_name, sep="\t")
        # add same safety checks len > 1
        # Flag for pride pipeline, or to set from second to minute as input rt time scale
        moff_pride_flag = False
        if moff.check_ps_input_data(df.columns.tolist(), ast.literal_eval(config.get('moFF', 'moffpride_format'))) == 1:
            # if it is a moff_pride data I do not check aany other requirement
            log.critical('moffPride input detected')
            moff_pride_flag = True
        else:
            if not 'matched' in df.columns:
                # check if it is a PS file ,
                list_name = df.columns.values.tolist()
                # get the lists of PS  defaultcolumns from properties file
                list = ast.literal_eval(config.get(
                    'moFF', 'ps_default_export_v1'))
                # here it controls if the input file is a PS export; if yes it maps the input in right moFF name
                if moff.check_ps_input_data(list_name, list) == 1:
                    # map  the columns name according to moFF input requirements
                    if not args.peptide_summary:
                        data_ms2, list_name = moff.map_ps2moff(
                            df, 'col_must_have_apex')
                    else:
                        data_ms2, list_name = moff.map_ps2moff(
                            df, 'col_must_have_mbr')
                # check if the field names are     good, in case of pep summary we need same req as in  mbr
            if args.peptide_summary:
                if moff.check_columns_name(df.columns.tolist(),
                                           ast.literal_eval(config.get('moFF', 'col_must_have_mbr')), log) == 1:
                    exit('ERROR minimal field requested are missing or wrong')
            else:
                if moff.check_columns_name(df.columns.tolist(),
                                           ast.literal_eval(config.get('moFF', 'col_must_have_apex')), log) == 1:
                    exit('ERROR minimal field requested are missing or wrong')

        # check if filtering is UP and the input data is not suitable for mbr filtering
        if 'off' in args.mbr and args.match_filter:
            if not 'matched' in df.columns:
                exit(
                    'mbr peptide not detect in the input file, filtering of mbr peptides is not possible. Please set --match_filter to 0 and run again.')
            if not ('mod_peptide' in df.columns):
                exit(
                    'mod_peptide sequence is not present your the input file, filtering of mbr peptides is not possible. Please check your infput file or parameter settings')
        log.critical('Starting Apex for %s ...', file_name)
        log.critical('moff Input file: %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ' % (
            file_name, tol, h_rt_w, s_w))
        if args.raw_list is None:
            log.critical('RAW file from folder :  %s' % loc_raw)
        else:
            log.critical('RAW file  :  %s' % args.raw_list)
        log.critical('Output file in :  %s', loc_output)
        # load the ptm file IF
        # mbr on with filtering  UP
        # mbr off with filtering flag UP (already check if inputdata contains matched field.)
        if 'matched' in df.columns and args.match_filter:
            log.critical('Apex module has detected mbr peptides')
            with open(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), args.ptm_file)) as data_file:
                ptm_map = json.load(data_file)

        name = os.path.basename(file_name).split('.')[0]
        #  IF raw_list contains mzML file -->  I m going to  read the file,
        #  one time just to save all the scan  Id and their RT.
        rt_list, id_list = moff.scan_mzml(raw_list)

        # control id the folder exist
        moff.check_output_folder_existence(loc_output)

        # control if exist the same log file : avoid appending output
        # moff.check_log_existence(os.path.join(loc_output, name + '__moff.log'))

        if args.match_filter:
            with open(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), args.ptm_file)) as data_file:
                ptm_map = json.load(data_file)
            start_time = time.time()

            moff.set_logger(log_file)
            log.critical('starting estimation of quality measures..')
            # run estimation_parameter
            rt_drift, not_used_measure, error_ratio = moff.estimate_parameter(
                df, name, raw_list, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output, rt_list, id_list, moff_pride_flag,
                ptm_map, args.sample_size, args.quantile_thr_filtering, args.match_filter, log_file, num_CPU)
            log.critical(
                'quality threhsold estimated : MAD_retetion_time %r  Ratio Int. FakeIsotope/1estIsotope: %r' % (
                rt_drift, error_ratio))
            log.critical('starting apex quantification of MS2 peptides..')
            log.info('log of MS2 identified peptide not retrived :  ..')
            moff.clean_json_temp_file(loc_output)
            myPool = multiprocessing.Pool(num_CPU)
            data_split = np.array_split(
                df[df['matched'] == 0], num_CPU)
            result = {}
            offset = 0
            for df_index in range(0, len(data_split)):
                result[df_index] = myPool.apply_async(moff.apex_multithr, args=(
                data_split[df_index], name, raw_list, tol, h_rt_w, s_w, s_w_match,
                loc_raw, loc_output, offset, rt_list, id_list, moff_pride_flag, ptm_map, 0, rt_drift, error_ratio, 0,
                log_file))
                offset += len(data_split[df_index])
            # save ms2 resulr
            ms2_data = moff.save_moff_apex_result(result)
            log.critical('end  apex quantification of MS2 peptides..')
            log.critical(
                'starting quantification with matched peaks using the quality filtering...')
            log.critical('initial # matched peaks: %r',
                         df[df['matched'] == 1].shape)
            moff.clean_json_temp_file(loc_output)
            log.info('Log Matched Peptides filtered :')
            data_split = np.array_split(
                df[df['matched'] == 1], num_CPU)
            result = {}
            offset = 0
            for df_index in range(0, len(data_split)):
                result[df_index] = myPool.apply_async(moff.apex_multithr, args=(
                data_split[df_index], name, raw_list, tol, h_rt_w, s_w, s_w_match,
                loc_raw, loc_output, offset, rt_list, id_list, moff_pride_flag, ptm_map, 0, rt_drift, error_ratio,
                args.match_filter, log_file))
                offset += len(data_split[df_index])
            myPool.close()
            myPool.join()
            log.critical('end apex quantification matched peptide ')
            log.critical('Computational time (sec):  %4.4f ' %
                         (time.time() - start_time))
            matched_peak = moff.save_moff_apex_result(result)
            log.critical('after filtering matched peak #%r ',
                         matched_peak.shape[0])
            # concat the ms2 res  + mateched result
            final_res = pd.concat([ms2_data, matched_peak])
            # save result
            final_res.to_csv(os.path.join(loc_output, os.path.basename(
                name).split('.')[0] + "_moff_result.txt"), sep="\t", index=False)
            moff.clean_json_temp_file(loc_output)
        else:
            moff.set_logger(log_file)
            log.critical(
                'starting  peptide quantification (ms2 / matched ) ..')
            myPool = multiprocessing.Pool(num_CPU)
            data_split = np.array_split(df, num_CPU)
            result = {}
            offset = 0
            log.info('log of MS2 identified peptide not retrived ')
            start_time = time.time()
            for df_index in range(0, len(data_split)):
                result[df_index] = myPool.apply_async(moff.apex_multithr,
                                                      args=(data_split[df_index], name, raw_list, tol, h_rt_w,
                                                            s_w, s_w_match, loc_raw, loc_output, offset, rt_list,
                                                            id_list, moff_pride_flag, None, 0, -1, -1, 0, log_file))
                offset += len(data_split[df_index])
            myPool.close()
            myPool.join()
            log.critical('end apex quantification (ms2 / matched ) peptides')
            log.critical('computational time (sec):  %4.4f ' %
                         (time.time() - start_time))
            start_time_2 = time.time()
            result = moff.save_moff_apex_result(result)
            result.to_csv(os.path.join(loc_output, os.path.basename(name).split(
                '.')[0] + "_moff_result.txt"), sep="\t", index=False)
            moff.clean_json_temp_file(loc_output)

        fh.close()
        log.removeHandler(fh)
        moff.detach_handler()

    moff.clean_json_temp_file(loc_output)
    if args.peptide_summary:
        state = moff.compute_peptide_matrix(args.loc_out, log, args.tag_pepsum)
        if not state:
            log.critical(
                'Error during the computation of the peptide intensity summary file: Check the output folder that contains the moFF results file')
