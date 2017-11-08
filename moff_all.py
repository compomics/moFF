#!/usr/bin/env python
import argparse
import logging
import multiprocessing
import os
import time

import numpy as np
import pandas as pd

import moff
import moff_mbr

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)










if __name__ == '__main__':

	multiprocessing.freeze_support()

	parser = argparse.ArgumentParser(description='moFF match between run and apex module input parameter')

	parser.add_argument('--inputF', dest='loc_in', action='store',
						help='specify the folder of the input MS2 peptide list files ', required=False)

	parser.add_argument('--inputtsv', dest='tsv_list', action='store', nargs='*' ,
						help='specify the mzid file as a list ', required=False)


	parser.add_argument('--inputraw', dest='raw_list', action='store',  nargs='*' ,
						help='specify the raw file as a list ', required=False)

	parser.add_argument('--sample', dest='sample', action='store',
						help='specify witch replicated to use for mbr reg_exp are valid ', required=False)

	parser.add_argument('--ext', dest='ext', action='store', default='txt',
						help='specify the file extentention of the input like ', required=False)

	parser.add_argument('--log_file_name', dest='log_label', action='store', default='moFF',
						help='a label name to use for the log file', required=False)

	parser.add_argument('--filt_width', dest='w_filt', action='store', default=2,
						help='width value of the filter  k * mean(Dist_Malahobis)', required=False)

	parser.add_argument('--out_filt', dest='out_flag', action='store', default=1,
						help='filter outlier in each rt time allignment', required=False)

	parser.add_argument('--weight_comb', dest='w_comb', action='store', default=0,
						help='weights for model combination combination : 0 for no weight  1 weighted devised by trein err of the model.',
						required=False)

	# parser.add_argument('--input', dest='name', action='store',help='specify input list of MS2 peptides ', required=True)

	parser.add_argument('--tol', dest='toll', action='store', type=float, help='specify the tollerance  parameter in ppm',
						required=True)

	parser.add_argument('--rt_w', dest='rt_window', action='store', type=float, default=3,
						help='specify rt window for xic (minute). Default value is 3 min', required=False)

	parser.add_argument('--rt_p', dest='rt_p_window', action='store', type=float, default=1,
						help='specify the time windows for the peak ( minute). Default value is 1 minute ', required=False)

	parser.add_argument('--rt_p_match', dest='rt_p_window_match', action='store', type=float, default=1,
						help='specify the time windows for the matched peptide peak ( minute). Default value is 1.2 minute ',
						required=False)

	parser.add_argument('--raw_repo', dest='raw', action='store', help='specify the raw file repository ', required=False)

	parser.add_argument('--output_folder', dest='loc_out', action='store', default='', help='specify the folder output',
						required=False)

	parser.add_argument('--rt_feat_file', dest='rt_feat_file', action='store',
						help='specify the file that contains the features to use in the match-between-run RT prediction ',
						required=False)

	parser.add_argument('--peptide_summary', dest='pep_matrix', action='store',type=int,default= 0,
						help='sumarize all the peptide intesity in one tab-delited file ',
						required=False)

	parser.add_argument('--tag_pep_sum_file', dest='tag_pepsum', action='store',type=str,default= 'moFF_run', help='a tag that is used in the peptide summary file name',required=False)

	args = parser.parse_args()

	## init globa logger
	ch = logging.StreamHandler()
	ch.setLevel(logging.ERROR)
	log.addHandler(ch)

	if (args.tsv_list is None) and  (args.loc_in is None) and  (args.raw_list is None) and (args.raw is None) :
		exit('you must specify the input and raw files ')
	if (args.tsv_list is not None) and  (args.loc_in is not None) and  (args.raw_list is not None) and (args.raw is not None) :
		 exit('you must specify the input and raw files or unsing: --inputtsv and --rawlist or --inputF and --rawrepo ')
	else:
		if ((args.tsv_list is None ) and (args.raw_list is not None) ) or ((args.tsv_list is not  None ) and (args.raw_list is  None) ):
			exit('Missing information: using --inputtsv you must specify the raw file with --inputraw ')
		if ((args.loc_in is None ) and (args.raw is not None) ) or ((args.loc_in is not  None ) and (args.raw is  None) ) :
			exit('Missing information: using --inputF you must specify the raw file with --raw_repo ')


	log.critical('Matching between run module (mbr)')


	# fixed variable number of split and also number of CPU presence in the macine
	# change this variable  with repset to the machine setting of the user
	num_CPU=multiprocessing.cpu_count()


	res_state,mbr_list_loc = moff_mbr.run_mbr(args)
	if res_state == -1:
		exit('An error is occurred during the writing of the mbr file')
	if args.tsv_list is not None:
		# input list of raw and tsv file
		if len(args.tsv_list) != len(args.raw_list) :
			exit('Error:  number of the input files is different from the number of raw files' )
		# in case list of file as input , mbr_output is written in local folder
		folder = os.path.join('mbr_output')
	else:
		folder = os.path.join(args.loc_in, 'mbr_output')

	log.critical('Apex module... ')
	c=0
	start_time_total = time.time()
	for file_name in mbr_list_loc:
		tol = args.toll
		h_rt_w = args.rt_window
		s_w = args.rt_p_window
		s_w_match = args.rt_p_window_match
		if args.tsv_list is not None:
		## list of the raw file and their path
			raw_list = args.raw_list[c]
		else:
			raw_list = None

		loc_raw = args.raw
		loc_output = args.loc_out


		## add multi thredign option
		df = pd.read_csv(file_name,sep="\t")

		data_split= np.array_split(df, num_CPU)
		'''
		# small workaround to prevent max input line in Linux
		if data_split[0].shape[0] > 2500 :
			# increase the number of splitting a bit more than the CPUs in order to get splice smaller than 2500
			data_split = np.array_split(df,  (num_CPU + 8) )
		'''
		log.critical('Starting Apex for %s ...',file_name)
		log.critical('moff Input file: %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ' % (file_name, tol, h_rt_w, s_w))
		if args.raw_list is None:
			log.critical('RAW file from folder :  %s' % loc_raw)
		else:
			log.critical('RAW file  :  %s' % args.raw_list)
		log.critical('Output file in :  %s', loc_output)
		if 'matched' in df.columns:
			log.critical('Apex module has detected mbr peptides')

		#print 'Original input size', df.shape
		name = os.path.basename(file_name).split('.')[0]

		#  IF raw_list contains mzML file -->  I m going to  read the file, one time just to save all the scan  Id and their RT.
		rt_list, id_list = moff.scan_mzml(raw_list)

		#control id the folder exist
		moff.check_output_folder_existence(loc_output)
		#control if exist the same log file : avoid appending output
		moff.check_log_existence(os.path.join(loc_output, name + '__moff.log'))
        # this flag must be set to 0. it is 1 only in case moFF-Pride date and only in tha pex module
		moff_pride_flag= 0
		myPool = multiprocessing.Pool(num_CPU)

		start_time= time.time()
		result = {}
		offset = 0
		for df_index in range(0,len(data_split)):

			result[df_index] = myPool.apply_async(moff.apex_multithr,args = (data_split[df_index],name,raw_list, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output, offset,rt_list , id_list , moff_pride_flag))
			offset += len(data_split[df_index])

		myPool.close()
		myPool.join()
		#print ' TIME multi thre. terminated', time.time() - start_time
		log.critical('...apex terminated in  %4.4f sec', time.time() - start_time )
		moff.save_moff_apex_result (data_split, result, loc_output, file_name  )
		log.critical('TOTAL time for apex %4.4f sec', time.time() - start_time_total)
		c+=1

	moff.clean_json_temp_file(loc_output)

	if args.pep_matrix == 1 :
		state = moff.compute_peptide_matrix(args.loc_out,log,args.tag_pepsum)
		if state == -1 :
			log.critical ('Error during the computation of the peptide intensity summary file: Check the output folder that contains the moFF results file')
