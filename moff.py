import numpy as np
import glob as glob
import pandas as pd
import os as os
import sys
import subprocess
import shlex 
import  logging
import argparse
from StringIO import StringIO

### input###
## - MS2 ID file
## - tol
## - half rt time window in minute
###### output
##  list of intensities..+







def run_apex( file_name, tol, h_rt_w , s_w, s_w_match, loc_raw,loc_output  ):
		
	#file_name=  args.name
	#tol= args.toll
	#h_rt_w = args.rt_window
	#s_w= args.rt_p_window
	#s_w_match= args.rt_p_window_match
	#loc_raw = args.raw
	#loc_output = args.loc_out
	# flag_for matching
	#mbr_flag=0

	name =  file_name.split('/')[1].split('.')[0]

	if loc_output != '':
		if  not (os.path.isdir(loc_output)):
			#print "created output Dir ",output_dir
			os.makedirs(loc_output )

		outputname=  loc_output + '/'  +  name  +  "_moff_result.txt"
		logging.basicConfig(filename= loc_output + '/' +  name +    '__moff.log',filemode='w',level=logging.INFO)
	else:

		outputname = name  +  "_moff_result.txt"
		logging.basicConfig(filename= name +  '__moff.log',filemode='w',level=logging.INFO)

	if loc_raw != None:
		loc = loc_raw +  name+ '.RAW'
	else:
		loc =   name + '.RAW'


	loggin.info('moff Input file %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ',file_name,tol,h_rt_w,s_w)
	logging.info('Output_file in %s', outputname)
	logging.info('RAW file location %s',loc)

	##read data from file 
	data_ms2 = pd.read_csv(file_name,sep= "\t" ,header=0)

	index_offset = data_ms2.columns.shape[0]   -1


	data_ms2["intensity"]=-1
	data_ms2["rt_peak"]=-1
	data_ms2["lwhm"]=-1
	data_ms2["rwhm"]=-1
	data_ms2["5p_noise"]=-1
	data_ms2["10p_noise"]=-1
	data_ms2["SNR"]=-1
	data_ms2["log_L_R"]=-1
	data_ms2["log_int"]=-1
	data_ms2["rt_peak"]=data_ms2["rt_peak"].astype('float64')
	data_ms2['intensity'] = data_ms2['intensity'].astype('float64')
	data_ms2['lwhm'] = data_ms2['lwhm'].astype('float64')
	data_ms2["rwhm"]=  data_ms2['rwhm'].astype('float64')
	data_ms2["5p_noise"]=  data_ms2['5p_noise'].astype('float64')
	data_ms2["10p_noise"]=  data_ms2['10p_noise'].astype('float64')
	data_ms2["SNR"]= data_ms2['SNR'].astype('float64')
	data_ms2["log_L_R"]= data_ms2['log_L_R'].astype('float64')
	data_ms2["log_int"]= data_ms2['log_int'].astype('float64')


	##set mbr_flag
	if 'matched' in data_ms2.columns :
		mbr_flag=1
	c=0
	for index_ms2, row in data_ms2.iterrows():
		logging.info('line: %i',c)
		mz_opt= "-mz="+str(row['mz'])
		if mbr_flag==1:
			s_w = s_w_match
		if row['rt']==-1:
			logging.info('rt not found. Wrong matched peptide in the mbr step line: %i',c)
			c+=1
			continue
			
		##convert rt to sec to min
		time_w= row['rt']/60
		## original s_W values is 0.40
		#=0.10 # time of refinement in minutes about 20 sec
		args = shlex.split("./txic " + mz_opt + " -tol="+ str(tol) + " -t " + str(time_w - h_rt_w) + " -t "+ str(time_w +h_rt_w) +" " + loc   )
		p= subprocess.Popen(args,stdout=subprocess.PIPE)
		output, err = p.communicate()
		#print p
		
		try:
			data_xic = pd.read_csv(StringIO(output.strip()), sep=' ',names =['rt','intensity'] ,header=0 )
			ind_v = data_xic.index
			#logging.info ("XIC shape   %i x 2",  data_xic.shape[0] )
			if data_xic[(data_xic['rt']> (time_w - s_w)) & ( data_xic['rt']< (time_w + s_w) )].shape[0] >=1:
				ind_v = data_xic.index
				pp=data_xic[ data_xic["intensity"]== data_xic[(data_xic['rt']> (time_w - s_w)) & ( data_xic['rt']< (time_w + s_w) )]['intensity'].max()].index
			#print 'pp index',pp
			#print 'Looking for ..:',row['mz'],time_w
			#print 'XIC data retrived:',data_xic.shape
			#print data_xic[ data_xic["intensity"]== data_xic[(data_xic['rt']> (time_w - )) & ( data_xic['rt']< (time_w + s_w) )]['intensity'].max()]
			## non serve forzarlo a in
				pos_p = ind_v[pp]
				#logging.info ("Pos of the max int",pp.values)
				#print pos_p.values.shape[0]
				if pos_p.values.shape[0] > 1 :
					logging.info("WARNINGS: Rt gap for the time windows searched. Probably the ppm values is too small %i", c )		
					continue
				val_max = data_xic.ix[pos_p,1].values 
			else:
				logging.info("LW_BOUND window  %4.4f", time_w - s_w )
				logging.info("UP_BOUND windows %4.4f", time_w + s_w )
				logging.info(data_xic[(data_xic['rt']> (time_w - (+0.60))) & ( data_xic['rt']< (time_w + (s_w+0.60)) )]   )
				logging.info("WARNINGS: moff_rtWin_peak is not enough to detect the max peak line : %i", c )
				logging.info( 'MZ: %4.4f RT: %4.4f Mass: %i',row['mz'] ,row['rt'],index_ms2 )
				c+=1
				continue
			pnoise_5 =  np.percentile(data_xic[(data_xic['rt']> (time_w - 1)) & ( data_xic['rt']< (time_w + 1) )]['intensity'],5)
			pnoise_10 = np.percentile(data_xic[(data_xic['rt']> (time_w - 1)) & ( data_xic['rt']< (time_w + 1) )]['intensity'],10) 
		except (IndexError,ValueError,TypeError):
			logging.info("WARNINGS:  size is not enough to detect the max peak line : %i", c )
			logging.info( 'MZ: %4.4f RT: %4.4f index: %i',row['mz'] ,row['rt'],index_ms2 )
			continue
			c +=1
		except pd.parser.CParserError:
			logging.info( "WARNINGS: XIC not retrived line: %i",c)
			logging.info( 'MZ: %4.4f RT: %4.4f Mass: %i',row['mz'] ,row['rt'],index_ms2 )
			
			c +=1
			continue
		else:
			#logging.info("Intensisty at pos_p-1 %4.4f",data_xic.ix[(pos_p-1),1].values )
			log_time = [-1,-1]
			c_left=0
			find_5=False
			stop=False
			while c_left < (pos_p-1) and stop != True :
			#print c_left
				
				if find_5==False and (data_xic.ix[(pos_p-1)-c_left,1].values <= (0.5 * val_max) ) :
					find_5=True
				#print "LWHM",c_left,data_xic.ix[(pos_p-1)-c_left,1]
				#log_point[0] = np.log2(val_max)/np.log2(data_xic.ix[(pos_p-1)-c_left,1]) 
					log_time[0] = data_xic.ix[(pos_p-1)-c_left,0].values*60
					stop=True
				c_left+=1
			find_5=False
			stop=False
			r_left=0
			while ((pos_p+1)+r_left  < len(data_xic) )  and stop != True :
				if find_5==False and data_xic.ix[(pos_p+1)+r_left,1].values <= (0.50 *val_max):
					find_5=True
				#print "RWHM",r_left,data_xic.ix[(pos_p+1)+r_left,1] 
				#log_point[2] = np.log2(val_max) /np.log2(data_xic.ix[(pos_p+1)+r_left,1])
					log_time[1] = data_xic.ix[(pos_p+1)+r_left,0].values*60
					stop=True
				r_left += 1
			
			data_ms2.ix[index_ms2, (index_offset +1) ]= val_max
			data_ms2.ix[index_ms2, (index_offset +2 ) ]= data_xic.ix[pos_p,0].values*60
			data_ms2.ix[index_ms2, (index_offset +3 )]= log_time[0]
			data_ms2.ix[index_ms2, (index_offset +4 )]= log_time[1]
			data_ms2.ix[index_ms2, (index_offset +5 ) ]= pnoise_5
			data_ms2.ix[index_ms2, (index_offset +6 ) ]= pnoise_10
			#conpute log_L_R SNR and log intensities 
			if (pnoise_5 ==0 and pnoise_10 > 0):
				data_ms2.ix[index_ms2, (index_offset +7 ) ] = 20 * np.log10( data_xic.ix[pos_p,1].values  / pnoise_10 )
			else:
				data_ms2.ix[index_ms2, (index_offset +7 ) ] = 20 * np.log10( data_xic.ix[pos_p,1].values  / pnoise_5 )
			#### WARNING  time - log_time 0 / time -log_time 1
			data_ms2.ix[index_ms2, (index_offset +8 ) ] = np.log2(abs( data_ms2.ix[index_ms2,index_offset +2 ]    - log_time[0] )   / abs( data_ms2.ix[index_ms2,index_offset +2 ]  -  log_time[1] ) ) 
			data_ms2.ix[index_ms2, (index_offset +9 ) ]=  np.log2(  val_max  )

			c+=1
	## save  result 
	data_ms2.to_csv(path_or_buf =outputname ,sep="\t",header=True,index=False )






if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='moFF input parameter')


	parser.add_argument('--input', dest='name', action='store',help='specify input list of MS2 peptides ', required=True)

	parser.add_argument('--tol', dest='toll',action='store',type= float,help='specify the tollerance  parameter in ppm', required=True)

	parser.add_argument('--rt_w', dest='rt_window', action='store',type= float, default=3,help='specify rt window for xic (minute). Default value is 3 min', required=True)

	parser.add_argument('--rt_p', dest='rt_p_window', action='store',type= float, default=0.1,help='specify the time windows for the peak ( minute). Default value is 0.1 ', required=False)

	parser.add_argument('--rt_p_match', dest='rt_p_window_match', action='store',type= float, default=0.4,help='specify the time windows for the matched peptide peak ( minute). Default value is 0.4 ', required=False)

	parser.add_argument('--raw_repo', dest='raw', action='store',help='specify the raw file repository ', required=True)

	parser.add_argument('--output_folder', dest='loc_out', action='store', default='',help='specify the folder output', required=False)

	args = parser.parse_args()
	file_name=  args.name
        tol= args.toll
        h_rt_w = args.rt_window
        s_w= args.rt_p_window
        s_w_match= args.rt_p_window_match
        loc_raw = args.raw
        loc_output = args.loc_out

        run_apex(file_name, tol, h_rt_w , s_w, s-w_match, loc_raw,loc_output )

