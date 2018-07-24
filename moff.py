#/usr/bin/env python

import ConfigParser
import argparse
import ast
import bisect
import glob
import logging
import multiprocessing
import os as os
import shlex
import subprocess
import sys
import time
import traceback
from sys import platform as _platform

import numpy as np
import pandas as pd
import pymzml
import simplejson as json

from pyteomics.mass import std_aa_comp,Unimod,calculate_mass
from collections import Counter
from  itertools  import chain
from brainpy import isotopic_variants
from scipy.stats import spearmanr

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


def clean_json_temp_file(loc_output):
	for f in glob.glob(loc_output + "/*.json"):
		os.remove(f)
	return 1
def compute_peptide_matrix(loc_output, log, tag_filename):
	name_col = []
	name_col.append('prot')
	d = []
	if not glob.glob(loc_output + '/*_moff_result.txt'):
		return -1
	for name in glob.glob(loc_output + '/*_moff_result.txt'):

		if 'match_' in os.path.basename(name):
			name_col.append('sumIntensity_' + os.path.basename(name).split('_match_moff_result.txt')[0])
		else:
			name_col.append('sumIntensity_' + os.path.basename(name).split('_moff_result.txt')[0])
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
		data.drop_duplicates(subset=['prot', 'peptide', 'mod_peptide', 'mass', 'charge'], keep='first', inplace=True)
		d.append(data[['prot', 'peptide', 'mod_peptide', 'mass', 'charge', 'rt_peak', 'rt', 'intensity']])

	intersect_share = reduce(np.union1d, ([x['peptide'].unique() for x in d]))
	index = intersect_share

	df = pd.DataFrame(index=index, columns=name_col)
	df = df.fillna(0)
	for i in range(0, len(d)):
		grouped = d[i].groupby('peptide', as_index=True)['prot', 'intensity']
		# print grouped.agg({'prot':'max', 'intensity':'sum'}).columns
		df.ix[:, i + 1] = grouped.agg({'prot': 'max', 'intensity': 'sum'})['intensity']
		df.ix[np.intersect1d(df.index, grouped.groups.keys()), 0] = grouped.agg({'prot': 'max', 'intensity': 'sum'})[
			'prot']
	# print df.head(5)
	df.reset_index(level=0, inplace=True)
	df = df.fillna(0)
	df.rename(columns={'index': 'peptide'}, inplace=True)
	log.critical('Writing peptide_summary intensity file')
	df.to_csv(os.path.join(loc_output, "peptide_summary_intensity_" + tag_filename + ".tab"), sep='\t', index=False)
	return 1


def save_moff_apex_result(list_df, result, folder_output):
	#print len(list_df)
	try:
		xx = []
		for df_index in range(0,len(list_df)):
			if result[df_index].get()[1] == -1:
				exit ('Raw file not retrieved: wrong path or upper/low case mismatch')
			else:
				#print result[df_index].get()[0]
					
				xx.append( result[df_index].get()[0] )

		#print len(xx)

		final_res = pd.concat(xx)
		if 'index' in final_res.columns:
			final_res.drop('index',axis=1,inplace=True )
		
		#final_res.to_csv(os.path.join(folder_output, os.path.basename(name).split('.')[0] + "_moff_result.txt"), sep="\t",index=False)
	except Exception as e :
		traceback.print_exc()
		print
		# print os.path.join(folder_output,os.path.basename(name).split('.')[0]  + "_moff_result.txt")
		raise e
	return (final_res)




def map_ps2moff(data,type_mapping):
	data.drop(data.columns[[0]], axis=1, inplace=True)
	data.columns = data.columns.str.lower()
	if type_mapping == 'col_must_have_mbr':
		data.rename(columns={'sequence': 'peptide', 'modified sequence': 'mod_peptide', 'measured charge': 'charge',
		                     'theoretical mass': 'mass', 'protein(s)': 'prot', 'm/z': 'mz'}, inplace=True)
	if type_mapping == 'col_must_have_apex':
		data.rename(columns={'sequence': 'peptide', 'measured charge': 'charge', 'theoretical mass': 'mass',
		                     'protein(s)': 'prot', 'm/z': 'mz'}, inplace=True)
	return data, data.columns.values.tolist()




'''
input list of columns
list of column names from PS default template loaded from .properties
'''

def check_ps_input_data(input_column_name, list_col_ps_default):
	input_column_name.sort()
	list_col_ps_default.sort()
	if list_col_ps_default == input_column_name:
		# detected a default PS input file
		return 1
	else:
		# not detected a default PS input file
		return 0


def check_columns_name(col_list, col_must_have,log):
	for c_name in col_must_have:
		if not (c_name in col_list):
			# fail
			log.critical('This information is missing : %s ', c_name)
			return  1
		# succes
	return 0


def scan_mzml ( name ):
# when I am using thermo raw and --raw_repo option used
	if name is None:
		return (-1,-1)
	if ('MZML' in name.upper()):

		rt_list = []
		runid_list = []
		run_temp = pymzml.run.Reader( name )
		for spectrum in run_temp:
			if spectrum['ms level'] == 1:
				rt_list.append(spectrum['scan start time'])
				runid_list.append(spectrum['id'])

		return (rt_list,runid_list )
	else:
		# in case of raw file  I put to -1 -1 thm result
		return (-1,-1 )


def  mzML_get_all( temp,tol,loc,run,   rt_list1, runid_list1 ):
	app_list=[]
	for index_ms2, row in temp.iterrows():
		
		data, status=pyMZML_xic_out(loc, float(tol / (10 ** 6)), row['ts'], row['te'], row['mz'],run, runid_list1,rt_list1 )
		# status is evaluated only herenot used anymore
		if status != -1 :
			app_list.append(data)
		else:
			app_list.append(pd.DataFrame(columns=['rt','intensity']))
	return app_list



def pyMZML_xic_out(name, ppmPrecision, minRT, maxRT, MZValue,run, runid_list,rt_list ):
	timeDependentIntensities = []
	minpos = bisect.bisect_left(rt_list, minRT)
	maxpos = bisect.bisect_left(rt_list, maxRT)

	for specpos in range(minpos,maxpos):
		specid = runid_list[specpos]
		spectrum = run[specid]
		if spectrum['scan start time'] > maxRT:
			break
		if spectrum['scan start time'] > minRT and spectrum['scan start time'] < maxRT:
			#print 'in ', specid
			lower_index = bisect.bisect(spectrum.peaks, (float(MZValue - ppmPrecision * MZValue), None))
			upper_index = bisect.bisect(spectrum.peaks, (float(MZValue + ppmPrecision * MZValue), None))
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
		return 1
	else:
		return -1


def check_output_folder_existence(loc_output ):
   if not os.path.exists(loc_output):
	os.mkdir(loc_output)
	return 1
   else:
	return 0




def compute_log_LR (data_xic,index,v_max, disc):
	log_time = [-1, -1]
	c_left = 0
	find_5 = False
	stop = False
	while c_left <= (index - 1) and not stop:
		if not find_5 and ( data_xic.ix[(index - 1) - c_left, 1] <= ( disc * v_max)):
			find_5 = True
			log_time[0] = data_xic.ix[(index - 1) - c_left, 0] * 60
			stop = True
		c_left += 1
	find_5 = False
	stop = False
	r_left = 0
	while ((index + 1) + r_left <  data_xic.shape[0] ) and not stop:
		if not find_5 and data_xic.ix[(index + 1) + r_left, 1] <= (disc * v_max):
			find_5 = True
			log_time[1] = data_xic.ix[(index + 1) + r_left, 0] * 60
			stop = True
		r_left += 1
	return log_time





def compute_peak_simple(x,xic_array,log,mbr_flag, h_rt_w,s_w,s_w_match,offset_index,moff_pride_flag ,rt_match_peak,count_match):
	if count_match != -1:
		c = (count_match * 4 ) + x.name
		#print count_match,x.name,'finale index', c
	else:
		c = x.name
	data_xic = xic_array[c]
	if rt_match_peak > -1 :
		time_w= rt_match_peak
	else:
		time_w= x['rt']
        if moff_pride_flag == 0 :
            # daling with rt in minutes , moffpride input data
            ## standar cases rt must be in second
	 	time_w= time_w /60
	#print time_w, x['rt'] , moff_pride_flag, rt_match_peak,time_w, s_w,s_w_match
	if mbr_flag == 0:
		log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f',(offset_index +c +2), x['mz'], time_w)
		temp_w = s_w
	else:
		log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f matched (yes=1/no=0): %i',(offset_index + c +2), x['mz'], time_w,x['matched'])
					# row['matched'])
		if x['matched'] == 1:
			temp_w = s_w_match
			
		else:
			temp_w = s_w
	if data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].shape[0] >= 1:
		#data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].to_csv('thermo_testXIC_'+str(c)+'.txt',index=False,sep='\t')
		ind_v = data_xic.index
		pp = data_xic[data_xic["intensity"] == data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))]['intensity'].max()].index
		pos_p = ind_v[pp]
		if pos_p.values.shape[0] > 1:
			print 'error'
			return pd.Series({'intensity': -1, 'rt_peak': -1,'lwhm':-1,'rwhm':-1,'5p_noise':-1,'10p_noise':-1,'SNR':-1,'log_L_R':-1,'log_int':-1})
		val_max = data_xic.ix[pos_p, 1].values
	else:
		log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f ', (offset_index +c +2), x['mz'], time_w)
		log.info("\t LW_BOUND window  %4.4f", time_w - temp_w)
		log.info("\t UP_BOUND window %4.4f", time_w + temp_w)
		log.info("\t WARNINGS: moff_rtWin_peak is not enough to detect the max peak ")

		return  pd.Series({'intensity': -1, 'rt_peak': -1,
				   'lwhm':-1,
					'rwhm':-1,
					'5p_noise':-1,
					'10p_noise':-1,
					'SNR':-1,
					'log_L_R':-1,
					'log_int':-1})
	pnoise_5 = np.percentile(data_xic[(data_xic['rt'] > (time_w - h_rt_w )) & (data_xic['rt'] < (time_w + h_rt_w ))]['intensity'], 5 )
	pnoise_10 = np.percentile( data_xic[(data_xic['rt'] > (time_w - h_rt_w )) & (data_xic['rt'] < (time_w + h_rt_w)  )]['intensity'], 10)
	# find the lwhm and rwhm
	time_point = compute_log_LR(data_xic, pos_p[0], val_max, 0.5)
	if (time_point[0]*time_point[1] == 1) or (time_point[0]*time_point[1] < 0):
		# Try a second time FWHM  computation with 0.7 * max intensity
		time_point =  compute_log_LR (data_xic,pos_p[0],val_max,0.7)

	if time_point[0]== -1 or  time_point[1] ==-1:
		# keep the shape measure to -1 in case on txo point are -1
		log_L_R=-1
	else:
			log_L_R= np.log2(abs( time_w  - time_point[0]) / abs( time_w - time_point[1]))
	
	if (pnoise_5 == 0 and pnoise_10 > 0):
				SNR  = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_10)
	else:
		if pnoise_5 != 0:
				SNR = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_5)
		else:
				log.info('\t 5 percentile is %4.4f (added 0.5)', pnoise_5)
				SNR = 20 * np.log10(data_xic.ix[pos_p, 1].values / (pnoise_5 +0.5))
	
	return pd.Series({'intensity': val_max[0], 'rt_peak': data_xic.ix[pos_p, 0].values[0] * 60,
			   'lwhm': time_point[0] ,
				'rwhm': time_point[1] ,
				'5p_noise': pnoise_5,
				'10p_noise':pnoise_10,
				'SNR':SNR[0],
				'log_L_R': log_L_R,
				'log_int': np.log2(val_max)[0] })




def estimate_parameter( df , name_file, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output,  rt_list , id_list, moff_pride_flag ,ptm_map,   log,sample_size, quantile_value  ):
	myPool = multiprocessing.Pool(  multiprocessing.cpu_count()   )
	sample= df[df['matched']==0 ].sample(frac=sample_size)
	log.critical('Estimate parameters using %r MS2 peptides randomly sampled' % sample.shape[0] )
	data_split = np.array_split(sample,  multiprocessing.cpu_count())
	result = {}
	offset = 0
	# run matchinf filtering for 
	for df_index in range(0, len(data_split)):
		result[df_index] = myPool.apply_async(apex_multithr_matched_peak, args=(data_split[df_index], name_file, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output,offset,  rt_list , id_list, moff_pride_flag ,ptm_map   ,1,-1,-1  ))
		offset += len(data_split[df_index])
	myPool.close()
	myPool.join()
	ms2_data= save_moff_apex_result(data_split, result, loc_output)
	log.critical ('Estimated distribution rank correlation exp. int. vs theor. int. %r %r %r '%( ms2_data['rankcorr'].quantile(0.25), ms2_data['rankcorr'].quantile(0.50),   ms2_data['rankcorr'].quantile(0.75))  )
	log.critical ('MAD retention time along all isotope %r', ms2_data['RT_drift'].describe())
	log.critical ('Estimated distribition ratio exp. int. left isotope vs. monoisotopic isotope %r ', ms2_data[ms2_data['delta_log_int'] != -1 ]['delta_log_int'].describe() )
	error_relInr =  ms2_data['Erro_RelIntensity_TheoExp'].quantile(quantile_value)
	rt_drift = ms2_data['RT_drift'].quantile(quantile_value)
	ratio_log_int = ms2_data[ms2_data['delta_log_int'] != -1 ]['delta_log_int'].quantile(quantile_value )
	return (rt_drift , error_relInr, ratio_log_int ) 

def compute_match_peak_quality_measure( input_data, moff_pride_flag,log ):
	mad_diff_int = np.mean( abs( (input_data['intensity'] /   input_data['intensity'].sum()) -  (input_data['ratio_iso'] / input_data['ratio_iso'].sum() )    ) )
	rank_spearman= spearmanr( (input_data['intensity'] /   input_data['intensity'].sum()) ,  input_data['ratio_iso'] ) [0]
	mad_rt = np.mean(abs( input_data['rt_peak']  -  input_data['rt_peak'].mean() ))
	return ( mad_diff_int, rank_spearman, mad_rt )

def  estimate_on_match_peak(x,input_data , estimate_flag,moff_pride_flag,log,  thr_q2,err_ratio_int,xic_data , mbr_flag ,h_rt_w,s_w,s_w_match,offset_index   ):
	print x.name
	input_data = input_data.iloc[x.name: x.name+4 ,:]
	print input_data.shape
	input_data.iloc[0:1,11:20] =  input_data.iloc[0:1,:].apply(lambda x : compute_peak_simple( x,xic_data ,log,mbr_flag ,h_rt_w,s_w,s_w_match,offset_index, moff_pride_flag,-1,c_count  ) , axis=1 )
	if input_data.ix[0,'log_L_R'] != -1 :
		if moff_pride_flag == 0 :
			new_point =  input_data.ix[0,'rt_peak']
		else:
		# to minute
			new_point =  input_data.ix[0,'rt_peak'] / 60
		input_data.iloc[1:4,11:20] =  input_data.iloc[1:4,:].apply(lambda x : compute_peak_simple( x,xic_data ,log,mbr_flag ,h_rt_w,0.5,0.5,offset_index, moff_pride_flag,new_point,c_count  ) , axis=1 )
		if  (input_data.ix[0:2,'log_L_R'] != -1).all():
			mad_diff_int, rank_spearman, mad_rt =  compute_match_peak_quality_measure( input_data.iloc[0:3,:], moff_pride_flag,log )
			#print input_data
			#print mad_diff_int, rank_spearman, mad_rt
			if input_data.ix[3,'log_L_R'] == -1:
			#print 'xxx missing wrong iso'
				return pd.Series({'Erro_RelIntensity_TheoExp': mad_diff_int, 'rankcorr': rank_spearman,'RT_drift': mad_rt ,'delta_rt': -1 ,'delta_log_int': -1})
			else:
				delta_rt_wrong_iso =  abs(input_data.ix[3,'rt_peak'] - input_data.ix[0:3,'rt_peak'].mean())
				delta_log_int =  input_data.ix[3,'log_int'] / input_data.ix[0,'log_int']
			#print 'yyy find wrong iso', delta_log_int , delta_rt_wrong_iso
				return pd.Series({'Erro_RelIntensity_TheoExp': mad_diff_int, 'rankcorr': rank_spearman,'RT_drift': mad_rt ,'delta_rt': delta_rt_wrong_iso ,'delta_log_int': delta_log_int})




def  filtering_match_peak(x,input_data , estimate_flag,moff_pride_flag,log,  thr_q2,err_ratio_int,xic_data , mbr_flag ,h_rt_w,s_w,s_w_match,offset_index   ):
 	print x.name
	input_data = input_data.iloc[x.name: x.name+4 ,:]
	print input_data.shape
	input_data.iloc[0:1,11:20] =  input_data.iloc[0:1,:].apply(lambda x : compute_peak_simple( x,xic_data ,log,mbr_flag ,h_rt_w,s_w,s_w_match,offset_index, moff_pride_flag,-1,c_count  ) , axis=1 )
	#print input_data
	if input_data.ix[0,'log_L_R'] != -1 :
		if moff_pride_flag == 0 :
			new_point =  input_data.ix[0,'rt_peak']
		else:
			# to minute
			new_point =  input_data.ix[0,'rt_peak'] / 60

		#input_data.iloc[1:2,11:20] =  input_data.iloc[1:2,:].apply(lambda x : compute_peak_simple( x,xic_data ,log,mbr_flag ,h_rt_w,0.3,0.3,offset_index, moff_pride_flag,new_point,c_count  ) , axis=1 )
		#input_data.iloc[2:3,11:20] =  input_data.iloc[2:3,:].apply(lambda x : compute_peak_simple( x,xic_data ,log,mbr_flag ,h_rt_w,0.3,0.3,offset_index, moff_pride_flag,new_point,c_count  ) , axis=1 )
		#input_data.iloc[3:4,11:20] =  input_data.iloc[3:4,:].apply(lambda x : compute_peak_simple( x,xic_data ,log,mbr_flag ,h_rt_w,0.3,0.3,offset_index, moff_pride_flag,new_point,c_count  ) , axis=1 )
		
		input_data.iloc[1:4,11:20] =  input_data.iloc[1:4,:].apply(lambda x : compute_peak_simple( x,xic_data ,log,mbr_flag ,h_rt_w,0.5,0.5,offset_index, moff_pride_flag,new_point,c_count  ) , axis=1 )
		# check isotope 2-3
		if (input_data.ix[0:2,'log_L_R'] != -1).all() : 
				
			mad_diff_int, rank_spearman, mad_rt =  compute_match_peak_quality_measure( input_data.iloc[0:3,:], moff_pride_flag,log )
			if  (mad_rt < thr_q2 and rank_spearman > 0.8):
			# check isotope -1 
				if  input_data.ix[3,'log_L_R'] != -1: 
					delta_rt_wrong_iso =  abs(input_data.ix[3,'rt_peak'] - input_data.ix[0:3,'rt_peak'].mean())
					delta_log_int =  input_data.ix[3,'log_int'] / input_data.ix[0,'log_int'] 
					if (delta_rt_wrong_iso  < thr_q2 and delta_log_int  >  err_ratio_int):
						# elimina overlapping peptide isotope
						log.info('%s --> Not valid isotope evelope  overlapping detected -->  --  MAD RT  %r  -- rankCorr %r ', input_data['peptide'].unique()[0], mad_rt , rank_spearman)
						return pd.Series({'intensity': -1, 'rt_peak': -1,'lwhm': -1 ,'rwhm': -1 ,'5p_noise': -1,'10p_noise': -1,'SNR': -1,'log_L_R': -1,'log_int': -1 })
					else:
						log.info('%s --> Valid isotope evelope detected after overlaping checkin -->  --  MAD RT  %r  -- rankCorr %r ', input_data['peptide'].unique()[0], mad_rt , rank_spearman)
						return input_data.loc[input_data['ratio_iso'].idxmax(axis=1), ['10p_noise','5p_noise','SNR','intensity','log_L_R','log_int' ,'lwhm','rt_peak','rwhm']]
				else:
					log.info('%s --> Valid isotope evelope detected and no overlaping detected -->  --  MAD RT  %r  -- rankCorr %r ', input_data['peptide'].unique()[0], mad_rt , rank_spearman)
					return input_data.loc[input_data['ratio_iso'].idxmax(axis=1), ['10p_noise','5p_noise','SNR','intensity','log_L_R','log_int' ,'lwhm','rt_peak','rwhm']]
			else:
				# not pass the thr. leveli
				log.info('%s --> Not valid isotope evelope detected  -->  --  MAD RT  %r  -- rankCorr %r ', input_data['peptide'].unique()[0], mad_rt , rank_spearman)
				return pd.Series({'intensity': -1, 'rt_peak': -1,'lwhm': -1 ,'rwhm': -1 ,'5p_noise': -1,'10p_noise': -1,'SNR': -1,'log_L_R': -1,'log_int': -1 })
		else:
			# I have only the 1st valid isotope peak  but not the second
			log.info('%s --> not enough isotope peak detected only  %r over 3(/4) detected ', input_data['peptide'].unique()[0], input_data[input_data['log_L_R']!= -1].shape[0]  )
			return pd.Series({'intensity': -1, 'rt_peak': -1,'lwhm': -1 ,'rwhm': -1 ,'5p_noise': -1,'10p_noise': -1,'SNR': -1,'log_L_R': -1,'log_int': -1 })
	else:
		#print 'not_find the 1st isotope'
		log.info('%s --> first isotope peak not detected %r over 3(/4) not detected   ', input_data['peptide'].unique()[0]  ,  input_data[input_data['log_L_R']!= -1].shape[0]  )
		return pd.Series({'intensity': -1, 'rt_peak': -1,
                           'lwhm': -1 ,
                                'rwhm': -1 ,
                                '5p_noise': -1,
                                '10p_noise': -1,
                                'SNR': -1,
                                'log_L_R': -1,
                                'log_int': -1 })




def apex_multithr_matched_peak(data_ms2,name_file, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output,offset_index,  rt_list , id_list, moff_pride_flag ,ptm_map,estimate_flag, rt_drift, err_ratio_int ):
	# ---
		#WARNING : you must call thi routine with anoter name this is just for quist debug
	# ---
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
	if estimate_flag==1:
		# add extra filed if I am in a estimate mode
		data_ms2["Erro_RelIntensity_TheoExp"] = -1
		data_ms2["rankcorr"] = -1
		data_ms2["RT_drift"] = -1
		data_ms2["delta_rt"]= -1
		data_ms2["delta_log_int"]=-1
	# set mbr_flag
	if 'matched' in data_ms2.columns:
		mbr_flag = 1
		#log.critical('Apex module has detected mbr peptides')
		#log.info('moff_rtWin_peak for matched peptide:   %4.4f ', s_w_match)

	# get txic path: assumes txic is in the same directory as moff.py
	txic_executable_name="txic_json.exe"
	txic_path = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), txic_executable_name)
	## to export a list of XIc
	all_isotope_df = pd.DataFrame(columns=['peptide','mz','ratio_iso','tol','rt','matched','ts','te'])
	# for all the input peptide in data_ms2
	try :
		for row in data_ms2.itertuples():
		# get the sequence
			#print '--  ---'
			## for MQ sequence is (mod_tag )
			# for PS sequenc eis <mod_tag> 
			
			if not ( '(' in row.mod_peptide):
				#  only fixed mod
				comps = Counter(list(chain(*[list(std_aa_comp[aa].elements()) for aa in row.peptide])))
				comps["H"] += 2
				comps["O"] += 1
				#print 'Final COMPS',comp
				fix_mod_count =  row.peptide.count('C')
				if fix_mod_count > 0:
					comps["H"] += (ptm_map['cC']['deltaChem'][0] * fix_mod_count )
					comps["C"] += (ptm_map['cC']['deltaChem'][1] * fix_mod_count )
					comps["N"] += (ptm_map['cC']['deltaChem'][2] * fix_mod_count)
					comps["O"] += (ptm_map['cC']['deltaChem'][3] * fix_mod_count)
			else: 
			
			#print 'still to do it '
				comps = Counter(list(chain(*[list(std_aa_comp[aa].elements()) for aa in row.peptide])))
				for ptm in ptm_map.keys() :
					ptm_c = row.mod_peptide.count(ptm)
				#ptm_c =  sum(ptm in s for s in row.mod_peptide)
					if ptm_c  >=1:
					#print ptm
						comps["H"] += (ptm_map[ptm]['deltaChem'][0] * ptm_c )
						comps["C"] += (ptm_map[ptm]['deltaChem'][1] * ptm_c)
						comps["N"] += (ptm_map[ptm]['deltaChem'][2] * ptm_c)
						comps["O"] += (ptm_map[ptm]['deltaChem'][3] * ptm_c)
						
				comps["H"] += 2
				comps["O"] += 1
			#print 'Modified COMPS',comps 
			theoretical_isotopic_cluster = isotopic_variants(comps, charge=  int(round (  row.mass / float( row.mz)))   ,npeaks=3)
			mz_iso = [ peak.mz  for peak in theoretical_isotopic_cluster ]
			delta = mz_iso[0]-mz_iso[1]
			mz_iso.append( mz_iso[0] + delta)
			ratio_iso = [ peak.intensity  for peak in theoretical_isotopic_cluster ]
			ratio_iso.append(-1)
			isotopic_df =  pd.DataFrame({'mz':mz_iso,'ratio_iso':ratio_iso})
			
			isotopic_df.ix[:,'exp_mz']= row.mz
			isotopic_df.ix[:,'peptide']= row.mod_peptide
			isotopic_df.ix[:,'tol'] = int( tol)
			isotopic_df.ix[:,'rt'] = row.rt
			isotopic_df.ix[:,'matched']=1
			if moff_pride_flag == 1:
				isotopic_df['ts'] = (row.rt  ) - h_rt_w
				isotopic_df['te'] = (row.rt   ) + h_rt_w
			else:
				isotopic_df['ts'] = (row.rt /60 ) - h_rt_w
				isotopic_df['te'] = (row.rt  /60  ) + h_rt_w
			all_isotope_df = pd.concat([all_isotope_df , isotopic_df ],join='outer',axis=0   )
		all_isotope_df.reset_index( inplace=True )
		#print all_isotope_df
		if not flag_mzml :
			# txic-28-9-separate-jsonlines.exe
			if not flag_windows:
			# Linux  to avoid cmd  string too long  and its error. the thresold is mainly base on  from empirical evaluation.
				#print all_isotope_df[['mz','tol','ts','te']].head(1).to_json(orient='records')
				if len(all_isotope_df.to_json(orient='records')) >= 50000:
					with open(os.path.join(loc_output,multiprocessing.current_process().name  +'_' + name_file + '.json'), 'w') as f:
						f.write(all_isotope_df[['mz','tol','ts','te']].to_json(orient='records'))
				
					args_txic = shlex.split( "mono " + txic_path + " -jf " +  os.path.join(loc_output,multiprocessing.current_process().name  +'_' + name_file  + '.json')  + " -f " + loc, posix=False)
				else:
					#  small amount of char. in the request
					args_txic = shlex.split( "mono " + txic_path + " -j " + all_isotope_df[['mz','tol','ts','te']].to_json( orient='records' ) + " -f " + loc,posix=True )
			else:
				# Windows to avoid cmd  string too long  and its error. the thresold is mainly base on  from empirical evaluation.
				if len(isotopic_df.to_json(orient='records')) >= 10000:
					with open(os.path.join(loc_output,multiprocessing.current_process().name +  '_' + name_file + '.json'), 'w') as f:
						f.write(all_isotope_df[['mz','tol','ts','te']].to_json(orient='records'))
					args_txic = shlex.split(txic_path + " -jf " +  os.path.join(loc_output,multiprocessing.current_process().name + '_' + name_file + '.json')  + " -f " + loc, posix=False)
				else:
					#  small amount of char. in the request
					args_txic = shlex.split(txic_path + " -j " + all_isotope_df[['mz','tol','ts','te']].to_json(orient='records') + " -f " + loc, posix=False)
			start_timelocal = time.time()
			p = subprocess.Popen(args_txic, stdout=subprocess.PIPE)
			output, err = p.communicate()
			xic_data=[]
			for l in range ( 0, all_isotope_df.shape[0]  ) :
				temp = json.loads( output.split('\n')[l].decode("utf-8") )
				xic_data.append(pd.DataFrame( { 'rt' : temp['results']['times'], 'intensity':  temp['results']['intensities'] }   , columns=['rt', 'intensity'] ) )
		else:	
				## to test later
			run_temp = pymzml.run.Reader(raw_name)
			xic_data = mzML_get_all( temp,tol,loc, run_temp ,rt_list , id_list  )
		## new filtering
		print all_isotope_df.shape,data_ms2.shape
		print all_isotope_df.head(12)
		all_isotope_df["intensity"] = -1
		all_isotope_df["rt_peak"] = -1
		all_isotope_df["lwhm"] = -1
		all_isotope_df["rwhm"] = -1
		all_isotope_df["5p_noise"] = -1
		all_isotope_df["10p_noise"] = -1
		all_isotope_df["SNR"] = -1
		all_isotope_df["log_L_R"] = -1
		all_isotope_df["log_int"] = -1
		if estimate_flag==0:
			data_ms2.iloc[c_count, data_ms2.columns.get_indexer(['10p_noise','5p_noise','SNR','intensity','log_L_R','log_int' ,'lwhm','rt_peak','rwhm' ]) ].apply( lambda x: filtering_match_peak( x,all_isotope_df,estimate_flag,   moff_pride_flag,log,   rt_drift, err_ratio_int, xic_data , mbr_flag  , h_rt_w,s_w,s_w_match,offset_index  ,  c_count ),axis=1)
		else:
			data_ms2[['Erro_RelIntensity_TheoExp','rankcorr','RT_drift','delta_rt','delta_log_int']].apply( lambda x: estimate_on_match_peak( x,all_isotope_df,estimate_flag,   moff_pride_flag,log,   rt_drift, err_ratio_int, xic_data , mbr_flag  , h_rt_w,s_w,s_w_match,offset_index  ),axis=1)


		'''
		for c_count, row in enumerate(data_ms2.itertuples()):
			working_df =   all_isotope_df[(all_isotope_df['peptide']==row.mod_peptide) & (all_isotope_df['exp_mz']==row.mz) ].copy()
			working_df["intensity"] = -1
			working_df["rt_peak"] = -1
			working_df["lwhm"] = -1
			working_df["rwhm"] = -1
			working_df["5p_noise"] = -1
			working_df["10p_noise"] = -1
			working_df["SNR"] = -1
			working_df["log_L_R"] = -1
			working_df["log_int"] = -1
			#print working_df.columns
			working_df.reset_index(inplace=True)
			#working_df.drop(['peptide'],axis=1,inplace=True)
			## cal the function
			if estimate_flag==0:
				data_ms2.iloc[c_count, data_ms2.columns.get_indexer(['10p_noise','5p_noise','SNR','intensity','log_L_R','log_int' ,'lwhm','rt_peak','rwhm' ]) ] = filtering_match_peak( working_df,estimate_flag,   moff_pride_flag,log,   rt_drift, err_ratio_int, xic_data , mbr_flag  , h_rt_w,s_w,s_w_match,offset_index  ,  c_count )
			else:
				data_ms2.iloc[c_count, data_ms2.columns.get_indexer(['Erro_RelIntensity_TheoExp','rankcorr','RT_drift','delta_rt','delta_log_int' ]) ] =  estimate_on_match_peak( working_df,estimate_flag,   moff_pride_flag,log,   rt_drift, err_ratio_int,xic_data , mbr_flag  , h_rt_w,s_w,s_w_match,offset_index  ,  c_count )
'''
		if estimate_flag != 1:
			data_ms2 =  data_ms2[ (data_ms2[ ['10p_noise','5p_noise','SNR','intensity','log_L_R','log_int' ,'lwhm','rt_peak','rwhm' ]  ] != -1 ).all(1)]
	except Exception as e:
		traceback.print_exc()
		print
		raise e
	return (data_ms2,1)

def apex_multithr(data_ms2,name_file, raw_name, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output,offset_index,  rt_list , id_list, moff_pride_flag ):
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

        # get txic path: assumes txic is in the same directory as moff.py
        txic_executable_name="txic_json.exe"
	txic_path = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), txic_executable_name)

        ## to export a list of XIc
        try:
                temp=data_ms2[['mz','rt']].copy()
        # strange cases

                temp.ix[:,'tol'] = int( tol)
                if moff_pride_flag == 1:
                        temp['ts'] = (data_ms2['rt']  ) - h_rt_w
                        temp['te'] = (data_ms2['rt']   ) + h_rt_w
                else:
                        temp['ts'] = (data_ms2['rt'] /60 ) - h_rt_w
                        temp['te'] = (data_ms2['rt'] /60  ) + h_rt_w
                temp.drop('rt',1,inplace=True )


                if not flag_mzml :
                        # txic-28-9-separate-jsonlines.exe
                        if not flag_windows:
                                # Linux  to avoid cmd  string too long  and its error. the thresold is mainly base on  from empirical evaluation.
                                if len(temp.to_json(orient='records')) >= 50000:
                                        with open(os.path.join(loc_output,multiprocessing.current_process().name  +'_' + name_file + '.json'), 'w') as f:
                                                f.write(temp.to_json(orient='records'))
                                        args_txic = shlex.split( "mono " + txic_path + " -jf " +  os.path.join(loc_output,multiprocessing.current_process().name  +'_' + name_file  + '.json')  + " -f " + loc, posix=False)
                                else:
                                        #  small amount of char. in the request
                                        args_txic = shlex.split( "mono " + txic_path + " -j " + temp.to_json( orient='records' ) + " -f " + loc,posix=True )
                        else:
                                # Windows to avoid cmd  string too long  and its error. the thresold is mainly base on  from empirical evaluation.
                                if len(temp.to_json(orient='records')) >= 10000:
                                        with open(os.path.join(loc_output,multiprocessing.current_process().name +  '_' + name_file + '.json'), 'w') as f:
                                                f.write(temp.to_json(orient='records'))
                                        args_txic = shlex.split(txic_path + " -jf " +  os.path.join(loc_output,multiprocessing.current_process().name + '_' + name_file + '.json')  + " -f " + loc, posix=False)
                                else:
                                        #  small amount of char. in the request
                                        args_txic = shlex.split(txic_path + " -j " + temp.to_json(orient='records') + " -f " + loc, posix=False)
                        start_timelocal = time.time()
                        p = subprocess.Popen(args_txic, stdout=subprocess.PIPE)
                        output, err = p.communicate()
                        xic_data=[]
                        for l in range ( 0,temp.shape[0] ) :
                                        temp = json.loads( output.split('\n')[l].decode("utf-8") )
                                        xic_data.append(pd.DataFrame( { 'rt' : temp['results']['times'], 'intensity':  temp['results']['intensities'] }   , columns=['rt', 'intensity'] ) )
                else:
                        run_temp = pymzml.run.Reader(raw_name)
                        xic_data =  mzML_get_all( temp,tol,loc, run_temp ,rt_list , id_list  )
                        #10p_noise    5p_noise  SNR     intensity  log_L_R    log_int  lwhm rt_peak  rwhm
                data_ms2.reset_index(inplace=True)
                data_ms2[['10p_noise','5p_noise','SNR','intensity','log_L_R','log_int' ,'lwhm','rt_peak','rwhm']] = data_ms2.apply(lambda x : compute_peak_simple( x,xic_data ,log,mbr_flag ,h_rt_w,s_w,s_w_match,offset_index, moff_pride_flag,-1, -1 ) , axis=1   )
        except Exception as e:
                traceback.print_exc()
                print
                raise e

        return (data_ms2,1)


def main_apex_alone():
	parser = argparse.ArgumentParser(description='moFF input parameter')
	parser.add_argument('--tsv_list', dest='name', action='store',
						help='specify the input file with the MS2 peptides/features', required=True)
	parser.add_argument('--raw_list', dest='raw_list', action='store', help='specify directly raw file', required=False)
	parser.add_argument('--toll', dest='toll', action='store', type=float, help='specify the tollerance parameter in ppm',
						required=True)
	parser.add_argument('--xic_length', dest='xic_length', action='store', type=float, default=3,
						help='specify rt window for xic (minute). Default value is 3 min', required=False)
	parser.add_argument('--rt_peak_win', dest='rt_peak_win', action='store', type=float, default=1,
						help='specify the time windows for the peak ( minute). Default value is 1 minute ', required=False)
	parser.add_argument('--rt_peak_win_match', dest='rt_peak_win_match', action='store', type=float, default=1.2,
						help='specify the time windows for the matched  peak ( minute). Default value is 1.2 minute ',
						required=False)
	parser.add_argument('--raw_repo', dest='raw_repo', action='store', help='specify the raw file repository folder',
						required=False)
	parser.add_argument('--loc_out', dest='loc_out', action='store', default='', help='specify the folder output',
						required=False)
	parser.add_argument('--peptide_summary', dest='peptide_summary', action='store', type=int, default=0, help='summarize all the peptide intesity in one tab-delited file ',required=False)
	parser.add_argument('--match_filter', dest='match_filter', action='store', type=int, default=0, help='filtering on the matched peak .default 0',required=False)
	parser.add_argument('--ptm_file', dest='ptm_file', action='store', default='ptm_setting.json', help='name of json ptm file. default file ptm_setting.json ',required=False)
	parser.add_argument('--quantile_thr_filtering', dest='quantile_thr_filtering', action='store', type=float, default=0.75, help='quantile value used to computed the filtering threshold for the matched peak .default 0.75',required=False)
	parser.add_argument('--sample_size', dest='sample_size', action='store', type=float, default=0.05, help='percentage of MS2 peptide used to estimated the threshold',required=False)
	parser.add_argument('--tag_pepsum', dest='tag_pepsum', action='store', type=str, default='moFF_run', help='a tag that is used in the peptide summary file name', required=False)


	args = parser.parse_args()

	if (args.raw_list is None) and (args.raw_repo is None):
		exit('you must specify and raw files  with --raw_list (file name) or --raw_repo (folder)')
	if (args.raw_list is not None) and (args.raw_repo is not None):
		exit('you must specify raw files using only one options --raw_list (file name) or --raw_repo (folder) ')


	file_name = args.name
	tol = args.toll
	h_rt_w = args.xic_length
	s_w = args.rt_peak_win
	s_w_match = args.rt_peak_win_match

	loc_raw = args.raw_repo
	loc_output = args.loc_out

	# set stream option for logger
	ch = logging.StreamHandler()
	ch.setLevel(logging.ERROR)
	log.addHandler(ch)

	config = ConfigParser.RawConfigParser()
	config.read(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 'moff_setting.properties'))
	df = pd.read_csv(file_name, sep="\t")
	## add same safety checks len > 1
	## check and eventually tranf for PS template
	moff_pride_flag = 1

	if check_ps_input_data(df.columns.tolist(), ast.literal_eval(config.get('moFF', 'moffpride_format'))) == 1:
		# if it is a moff_pride data I do not check aany other requirement
		log.critical('moffPride input detected')
		moff_pride_flag = 1
	else:
		if not 'matched' in df.columns:
			# check if it is a PS file ,
			list_name = df.columns.values.tolist()
			# get the lists of PS  defaultcolumns from properties file
			list = ast.literal_eval(config.get('moFF', 'ps_default_export_v1'))
			# here it controls if the input file is a PS export; if yes it maps the input in right moFF name
			if check_ps_input_data(list_name, list) == 1:
				# map  the columns name according to moFF input requirements
				if args.pep_matrix != 1:
					data_ms2, list_name = map_ps2moff(df,'col_must_have_apex')
				else:
					data_ms2, list_name = map_ps2moff(df, 'col_must_have_mbr')
		## check if the field names are good, in case of pep summary we need same req as in  mbr
		if args.peptide_summary == 1:
			if check_columns_name(df.columns.tolist(), ast.literal_eval(config.get('moFF', 'col_must_have_mbr')),log) == 1 :
				exit('ERROR minimal field requested are missing or wrong')
		else:
			if  check_columns_name(df.columns.tolist(), ast.literal_eval(config.get('moFF', 'col_must_have_apex')),log) == 1   :
				exit('ERROR minimal field requested are missing or wrong')

        # flag to check idf the input are from moffPride file.
        # if so, rt is already in minutes


	log.critical('moff Input file: %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ' % (file_name, tol, h_rt_w, s_w))
	if args.raw_list is None:
		log.critical('RAW file from folder :  %s' % loc_raw)
	else:
		log.critical('RAW file  :  %s' % args.raw_list)

	log.critical('Output file in :  %s', loc_output)


	log.critical('starting Apex module .....')
	name = os.path.basename(file_name).split('.')[0]

	
	rt_list , id_list = scan_mzml ( args.raw_list )
	## sampling
	if args.match_filter == 1 :
	#filtering
		# load ptm file 
		# ptm file MUST be located in the moFF folder
		with open(  os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),args.ptm_file)  ) as data_file:
			ptm_map = json.load(data_file)
		start_time = time.time()
		# run estimation _parameter
		rt_drift, not_used_measure,error_ratio =  estimate_parameter( df, name, args.raw_list, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output,  rt_list , id_list,  moff_pride_flag, ptm_map,log,args.sample_size,args.quantile_thr_filtering  )
		#rt_drift = 6
		#thr1= -1
		#error_ratio = 0.90
		log.critical( 'quality threhsold estimated : MAD_retetion_time  %r  Ratio Int. FakeIsotope/1estIsotope: %r '% ( rt_drift ,error_ratio))
		log.critical( 'starting MS2 peaks..')
		myPool = multiprocessing.Pool(  1  )
		data_split = np.array_split(df[df['matched']==0 ].head(1) , 1  )
		result = {}
		offset = 0
		for df_index in range(0, len(data_split)):
			result[df_index] = myPool.apply_async(apex_multithr, args=(data_split[df_index], name, args.raw_list, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output, offset,  rt_list , id_list,  moff_pride_flag  ))
			offset += len(data_split[df_index])
		# save ms2 resulr
		ms2_data = save_moff_apex_result( data_split, result, loc_output )
		log.critical( 'starting matched peaks...')
		log.critical( 'initial # matched peaks: %r', df[ df['matched']==1].shape )
		data_split = np.array_split(df[ df['matched']==1 ].head(10)   , 1  )
		print df[ df['matched']==1 ][['peptide','prot']].head(1)
		result = {}
		offset = 0
		for df_index in range(0, len(data_split)):
			result[df_index] = myPool.apply_async(apex_multithr_matched_peak, args=(data_split[df_index], name, args.raw_list, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output, offset,  rt_list , id_list,  moff_pride_flag, ptm_map,0 , rt_drift,error_ratio ))
			offset += len(data_split[df_index])

		myPool.close()
		myPool.join()
		log.critical('...apex terminated')
		log.critical( 'Computational time (sec):  %4.4f ' % (time.time() -start_time))
		print 'Time no result collect',  time.time() -start_time
		matched_peak= save_moff_apex_result(data_split, result, loc_output)
		log.critical('after filtering matched peak #%r ',matched_peak.shape[0])
		# concat the ms2 res  + mateched result 
		final_res = pd.concat([ms2_data,matched_peak])
		# save result
		final_res.to_csv(os.path.join(loc_output, os.path.basename(name).split('.')[0] + "_moff_result.txt"), sep="\t",index=False)
	else:
	## no filtering matched peak
		myPool = multiprocessing.Pool(  multiprocessing.cpu_count()   )
		data_split = np.array_split(df , multiprocessing.cpu_count()  )
		result = {}
		offset = 0
		start_time = time.time()
		for df_index in range(0, len(data_split)):
			result[df_index] = myPool.apply_async(apex_multithr, args=(data_split[df_index], name, args.raw_list, tol, h_rt_w, s_w, s_w_match, loc_raw, loc_output, offset,  rt_list , id_list,  moff_pride_flag ))
			offset += len(data_split[df_index])
		myPool.close()
		myPool.join()
		log.critical('...apex terminated')
		log.critical( 'Computational time (sec):  %4.4f ' % (time.time() -start_time))
		print 'Time no result collect',  time.time() -start_time
		start_time_2 = time.time()
		save_moff_apex_result(data_split, result, loc_output, file_name)

	if args.peptide_summary == 1 :
		# # TO DO manage the error with retunr -1 like in moff_all.py  master repo
		state = compute_peptide_matrix(loc_output,log,args.tag_pepsum)
		if state ==-1 :
			log.critical ('Error during the computation of the peptide intensity summary file: Check the output folder that contains the moFF results file')

if __name__ == '__main__':
	main_apex_alone()
