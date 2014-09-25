import numpy as np
import glob as glob
import pandas as pd
import os as os
import scipy, pylab
import pandas.io.sql as psql
import sys
import subprocess
import shlex 
import argparse
from StringIO import StringIO

### input###
## - MS2 ID file
## - tol
## - half rt time window in minute
###### output
##  list of intensities..+


parser = argparse.ArgumentParser(description='moFF input parameter')


parser.add_argument('--input', dest='name', action='store',
                   help='specify input list of MS2 peptides ', required=True)

parser.add_argument('--tol', dest='toll',action='store',type= float,
                   help='specify the tollerance  parameter in ppm', required=True)

parser.add_argument('--rt_w', dest='rt_window', action='store',type= float,
                   help='specify rt window for xic (minute)', required=True)

parser.add_argument('--rt_p', dest='rt_p_window', action='store',type= float, default=0.4,
                   help='specify the time windows for the peak ( minute). Default value is 0.4 ', required=False)

parser.add_argument('--raw', dest='raw', action='store',
                   help='specify the raw file', required=True)

parser.add_argument('--folder_output', dest='loc_out', action='store', default='',
                   help='specify the folder output', required=False)

args = parser.parse_args()


file_name= args.name 
tol= args.toll 
h_rt_w =args.rt_window
file_raw= args.raw 
#loc_raw =''   
loc_output = args.loc_out    #str(sys.argv[6])
s_w= args.rt_p_window
# defaul value is 0.4 , time of refinement in minutes about 20 sec



#outputname= loc_output  +  file_name.split('.')[0] +"_"+ file_name.split('.')[1] +"_moff_result.txt"

outputname= loc_output  +  file_name +"_moff_result.txt"


##read data from file 
data_ms2 = pd.read_csv(file_name,sep= "\t" ,header=0)
data_ms2["intensity"]=-1
data_ms2["rt_peak"]=-1
data_ms2["lwhm"]=-1
data_ms2["rwhm"]=-1
data_ms2["SNR"]=-1
data_ms2["log_L_R"]=-1
data_ms2["rt_peak"]=data_ms2["rt_peak"].astype('float64')
data_ms2['intensity'] = data_ms2['intensity'].astype('float64')
data_ms2['lwhm'] = data_ms2['lwhm'].astype('float64')
data_ms2["rwhm"]=  data_ms2['rwhm'].astype('float64')
data_ms2["SNR"]=  data_ms2['SNR'].astype('float64')
data_ms2["log_L_R"]=  data_ms2['log_L_R'].astype('float64')

print 'Load ',file_name, '# Int. to compute ', data_ms2.shape[0]  
###
loc= file_raw
#print loc
c=0
for index_ms2, row in data_ms2.iterrows():
	#print c
        
	mz_opt= "-mz="+str(row['precursor'])
	##convert rt to sec to min
	time_w= row['rt']/60
	print 'Computing the Intensity for  ', str(row['precursor']),str(row['rt'])
	args = shlex.split("./txic " + mz_opt + " -tol="+ str(tol) + " -t " + str(time_w - h_rt_w) + " -t "+ str(time_w +h_rt_w) +" " + loc   )
	#print args
	p= subprocess.Popen(args,stdout=subprocess.PIPE)
	output, err = p.communicate()
	data_xic = pd.read_csv(StringIO(output.strip()), sep=' ',names =['rt','intensity'] ,header=0 )
	try: 
		ind_v = data_xic.index
		pp=data_xic[ data_xic["intensity"] == data_xic[(data_xic['rt']> (time_w - s_w)) & ( data_xic['rt']< (time_w + s_w) )]['intensity'].max()].index
		pos_p = ind_v[int(pp)] 
		val_max = data_xic.ix[pos_p,1]
		pnoise_5 =  np.percentile(data_xic[(data_xic['rt']> (time_w - 1)) & ( data_xic['rt']< (time_w + 1) )]['intensity'],5)
		#pnoise_10 = np.percentile(data_xic[(data_xic['rt']> (time_w - 1)) & ( data_xic['rt']< (time_w + 1) )]['intensity'],10) 
	except ValueError:
		print "WARNING: s_w size is not enough to detect the max peak", "line :", c
		c = c + 1
	else:
	## quality controll take 50% FWHI and 0.5 
		log_time = [-1,-1]
		c_left=0
		find_5=False
		stop=False
		while c_left < (pos_p-1) and stop != True :
		#print c_left
			if find_5==False and data_xic.ix[(pos_p-1)-c_left,1] <= (0.5 * val_max) :
				find_5=True
			#print "LWHM",c_left,data_xic.ix[(pos_p-1)-c_left,1]
			#log_point[0] = np.log2(val_max)/np.log2(data_xic.ix[(pos_p-1)-c_left,1]) 
				log_time[0] = data_xic.ix[(pos_p-1)-c_left,0]*60
				stop=True
			c_left+=1
		find_5=False
		stop=False
		r_left=0
		while ((pos_p+1)+r_left  < len(data_xic) )  and stop != True :
			if find_5==False and data_xic.ix[(pos_p+1)+r_left,1] <= (0.50 *val_max):
				find_5=True
            		#print "RWHM",r_left,data_xic.ix[(pos_p+1)+r_left,1] 
			#log_point[2] = np.log2(val_max) /np.log2(data_xic.ix[(pos_p+1)+r_left,1])
				log_time[1] = data_xic.ix[(pos_p+1)+r_left,0]*60
				stop=True
			r_left += 1
	#"assignment output"
		data_ms2.ix[index_ms2,'intensity']=data_xic.ix[pos_p,1]
		data_ms2.ix[index_ms2,'rt_peak']= data_xic.ix[pos_p,0]*60
		data_ms2.ix[index_ms2,'lwhm']= log_time[0]
		data_ms2.ix[index_ms2,'rwhm']= log_time[1]
		data_ms2.ix[index_ms2,'SNR']=  20 * np.log10(  data_xic.ix[pos_p,1]  /  pnoise_5 )
		data_ms2.ix[index_ms2,'log_L_R']=  np.log2(abs( (data_xic.ix[pos_p,0]*60 ) -  log_time[0]  ) / abs( (data_xic.ix[pos_p,0]*60 )-  log_time[1]  ))
		c+=1
##print result 
data_ms2.to_csv(path_or_buf =outputname ,sep="\t",header=True,index=False )
