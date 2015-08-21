import numpy as np
import glob as glob
import pandas as pd
import os as os
import sys
import subprocess
import shlex 
import  logging
from StringIO import StringIO

### input###
## - MS2 ID file
## - tol
## - half rt time window in minute
###### output
##  list of intensities..+

file_name=str(sys.argv[1])
tol=str(sys.argv[2])
h_rt_w = float(sys.argv[3])
s_w= float(sys.argv[4])
loc_raw =str(sys.argv[5])
loc_output =str(sys.argv[6])


##output name
## just for debug one file at time
##print os.getcwd()

## ncomment these two lines, when you run from the split
#os.chdir(file_name.split('/')[0] + "/")
#file_name= file_name.split('/')[1]

first_file=False
#nn =  file_name.split('.')[0]
#count = file_name.split('.')[1]


## 20080 series
nn =  file_name.split('.')[0].split('_')[0 ] + '_' + file_name.split('.')[0].split('_')[1 ] + '_' + file_name.split('.')[0].split('_')[2] + '_' + file_name.split('.')[0].split('_')[3]
count = file_name.split('.')[0].split('_')[4]


## mam series

#if 'r3' in file_name.split('.')[0] :
#	nn =  file_name.split('.')[0].split('_')[0 ] + '_' + file_name.split('.')[0].split('_')[1 ] + '_' + file_name.split('.')[0].split('_')[2] + '_' + file_name.split('.')[0].split('_')[3] +  '_'  + file_name.split('.')[0].split('_')[4 ]    + '_'   + file_name.split('.')[0].split('_')[5 ] 
#	count = file_name.split('.')[0].split('_')[6]
#else:
#	nn =  file_name.split('.')[0].split('_')[0 ] + '_' + file_name.split('.')[0].split('_')[1 ] + '_' + file_name.split('.')[0].split('_')[2] + '_' + file_name.split('.')[0].split('_')[3] +  '_'  + file_name.split('.')[0].split('_')[4]
#	count = file_name.split('.')[0].split('_')[5]


## orbi series

#if ('_01' in file_name.split('.')[0]) or  ('_02' in file_name.split('.')[0]) :
#        nn =  file_name.split('.')[0].split('_')[0 ] + '_' + file_name.split('.')[0].split('_')[1 ] + '_' + file_name.split('.')[0].split('_')[2] + '_' + file_name.split('.')[0].split('_')[3] +  '_'  + file_name.split('.')[0].split('_')[4 ]    + '_'   + file_name.split('.')[0].split('_')[5 ] + '_'   + file_name.split('.')[0].split('_')[6 ] + '_'   + file_name.split('.')[0].split('_')[7 ]  + '_'   + file_name.split('.')[0].split('_')[8]
#        count = file_name.split('.')[0].split('_')[9]
#else:
#        nn =  file_name.split('.')[0].split('_')[0 ] + '_' + file_name.split('.')[0].split('_')[1 ] + '_' + file_name.split('.')[0].split('_')[2] + '_' + file_name.split('.')[0].split('_')[3] +  '_'  + file_name.split('.')[0].split('_')[4] + '_'   + file_name.split('.')[0].split('_')[5 ] + '_'   + file_name.split('.')[0].split('_')[6] + '_'   + file_name.split('.')[0].split('_')[7]
#        count = file_name.split('.')[0].split('_')[8]


outputname=  loc_output  +  nn  + "_" + count  + "_result.txt"
if  int (  count )==0:
        print "First file"
       	first_file=True


#+ "_" + file_name.split('.')[0].split('_')[4]   

logging.basicConfig(filename= nn +  "_" + count +   '__moff.log',filemode='w',level=logging.INFO)

logging.info('moff Input file %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ',file_name,tol,h_rt_w,s_w)
logging.info('Output_file in %s', outputname)


#outputname=  loc_output  +  nn  + "_" + file_name.split('.')[0].split('_')[5] + "_result.txt"

#print "OUTPUT:", outputname
#irst_file=False
#f  int (   file_name.split('.')[0].split('_')[5] )==0:
#print "primo file"
#first_file=True


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


#parametri fixed abs path  
## RAW per acluni file 
## raw per altri
file_raw =  nn +  ".RAW"
#loc_file="/home/compomics/extra_space/2981/"
###
logging.info(file_raw)
#print loc_raw
##location to get file
loc= loc_raw + file_raw
#rint loc

c=0
for index_ms2, row in data_ms2.iterrows():
	logging.info('line: %i',c)
	# versione per il vecchio formato
	mz_opt= "-mz="+str(row['mz'])
	#mz_opt= "-mz="+str(row['mz'])
	if row['rt']==-1:
		logging.info('rt not found. Wrong matched peptide in the mbr step line: %i',c)
		c+=1
		continue
		
	##convert rt to sec to min
	time_w= row['rt']/60
        ## original s_W values is 0.40
	#=0.10 # time of refinement in minutes about 20 sec
	args = shlex.split("./txic " + mz_opt + " -tol="+ tol + " -t " + str(time_w - h_rt_w) + " -t "+ str(time_w +h_rt_w) +" " + loc   )
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
			logging.info("LW_BOUND finestra per il max %4.4f", time_w - s_w )
			logging.info("UP_BOUND finestra per il max %4.4f", time_w + s_w )
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
	##if c==4:
	##  exit()
##print result 
if first_file:
	data_ms2.to_csv(path_or_buf =outputname ,sep="\t",header=True,index=False )
else:
	data_ms2.to_csv(path_or_buf =outputname ,sep="\t",header=False,index=False )
