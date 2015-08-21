import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from sklearn import linear_model
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score
import logging
import itertools
import os
import argparse
import re

## input file folder
parser = argparse.ArgumentParser(description='moFF match between run input parameter')


parser.add_argument('--inputF', dest='loc_in', action='store',
                   help='specify the folder of the input  MS2 peptide list files ', required=True)

parser.add_argument('--sample', dest='sample', action='store',
                   help='specify the folder of the input  MS2 peptide list files ', required=True)


args = parser.parse_args()
## fixed flag for the outlier dev
filt_outlier=1

paths = [os.path.join(args.loc_in,fn) for fn in next(os.walk(args.loc_in))[2] ]

exp_set=paths
output_dir=  'mbr_output'

if  not (os.path.isdir(output_dir)):
 	print "created output Dir ",output_dir
	os.makedirs(output_dir )
else:
	print "Outpur MBR in :",output_dir	

#exp_set=['R24692_4391_Isabo_20141106_Proteome','R24738_4399_Isabo_20141114_Secretome','R24788_4415_Isabo_20141120_secretome']

#exp_set=[' V23257_list.csv', 'V23258_list.csv']



#"Y:\\Andrea\\Cellysat_dataSet\\"

## read input
exp_t=[]
exp_out=[]
exp_subset=[]

#print paths
for a  in exp_set:
    #F1 = "Y:\\Andrea\\Secretome_DataSet\\" + a + "_MS2.txt"
    #F1 = "Y:\\Andrea\\Cellysat_dataSet\\" + a + "_MS2.txt"
    #print a
    if re.search(args.sample ,a) :
	print 'file selected ',a
	exp_subset.append(a)
    	data_moff = pd.read_csv(a, sep="\t", header=0)
    	list_name= data_moff.columns.values.tolist()
    #'spectrumid','filename','sequence','modified_sequence','exp_mass','accession','charge','rt','mz','db_filename','DB,score'
    #data_moff= data_moff[['sequence','modified_sequence','exp_mass','accession','charge','rt','mz']]
    # columns to add
    #matched == 1 if is an matched peak
    	data_moff['matched']=0
    	data_moff = data_moff.sort('rt')
    	exp_t.append(data_moff)
    	exp_out.append(data_moff)
    ## read



print 'Read input --> done '
## parameter of the number of query



exp_set=exp_subset
aa=range(0,len(exp_t))
print aa
out =list (itertools.product(aa,repeat=2))

# add matched columns
list_name.append('matched')


## RT correction unsing a windows 

def Var_corr (x) :
            pep = x[3] 
	    init_rt  = x[7]
            rt_target = exp_t[i[1]][exp_t[i[1]]['mod_peptide']==pep]['rt'].values       
            local_B = exp_t[i[1]][( exp_t[i[1]][ 'rt'] <=   rt_target[0] + 5)  & ( exp_t[i[1]]['rt'] >=   rt_target[0] -5) ][['mod_peptide','rt']] 
            local_B = local_B.reset_index()
        
            local_B= pd.merge(local_B, common ,on=['mod_peptide'],how='left' )
        
        
            local_B['rt_pred'] = local_B[ ~pd.isnull(local_B['rt_y'])]['rt_y'].apply(lambda x: clf_final.predict(x).tolist()[0][0]  ) 

            local_B['Diff_'] =local_B['rt_x'] -local_B['rt_pred']

            #print local_B
            k_point = local_B[ ~pd.isnull(local_B['rt_y'])].shape[0]
            if k_point == 0:
                return    init_rt
            else:
            #print 'point local used', k_point
                return  ( local_B['Diff_'].mean() / k_point) + init_rt


## filtering _outlier 

def MahalanobisDist(x, y):
    covariance_xy = np.cov(x,y, rowvar=0)
    inv_covariance_xy = np.linalg.inv(covariance_xy)
    xy_mean = np.mean(x),np.mean(y)
    x_diff = np.array([x_i - xy_mean[0] for x_i in x])
    y_diff = np.array([y_i - xy_mean[1] for y_i in y])
    diff_xy = np.transpose([x_diff, y_diff])
    
    md = []
    for i in range(len(diff_xy)):
        md.append(np.sqrt(np.dot(np.dot(np.transpose(diff_xy[i]),inv_covariance_xy),diff_xy[i])))
    return md
# x in dependent variable and Y is independar variable
def MD_removeOutliers(x, y):
    MD = MahalanobisDist(x, y)
    threshold = np.mean(MD) * 2 # adjust 1.5 accordingly 
    nx, ny, outliers = [], [], []
    for i in range(len(MD)):
        if MD[i] <= threshold:
            nx.append(x[i])
            ny.append(y[i])
        else:
            outliers.append(i) # position of removed pair
    return (np.array(nx), np.array(ny), np.array(outliers))



##input of the methods

#logging.basicConfig(filename='Y:\\Andrea\\Secretome_DataSet\\moFF_match2runs\\'+ 'production_log_mat2run.log',filemode='w',level=logging.WARNING)
logging.basicConfig(filename= output_dir + '/' + args.sample + '_' +'log_mat2run.log',filemode='w',level=logging.DEBUG)

for jj in aa:
    sum_values=np.array([0,0,0])
    print 'matching  in ', exp_set[jj]
    #match_couple=[]
    final_out = pd.DataFrame( columns = list_name )
    # columns to add
    #'spectrumid','sequence','modified_sequence','exp_mass','protein','charge','rt','mz','score'
    for i in out:
        if  i[0] == jj and i[1] != jj:
            logging.info( 'matching  %s peptide in   searching in %s ', exp_set[i[0]] ,exp_set[i[1]])
	    list_pep_repA = exp_t[i[0]]['mod_peptide'].unique()
	    list_pep_repB =  exp_t[i[1]]['mod_peptide'].unique()
	    set_dif_s_in_1 = np.setdiff1d(list_pep_repB,list_pep_repA)
	    add_pep_frame = exp_t[i[1]][exp_t[i[1]]['mod_peptide'].isin(set_dif_s_in_1)].copy()
	    #print add_pep_frame.columns
	    pep_shared = np.intersect1d(list_pep_repA,list_pep_repB)
	    logging.info( 'Peptide to add size  %i ', add_pep_frame.shape[0])
	    comA =   exp_t[i[0]][ exp_t[i[0]]['mod_peptide'].isin(pep_shared)][['mod_peptide','prot','rt']]
	    comB =   exp_t[i[1]][ exp_t[i[1]]['mod_peptide'].isin(pep_shared)][['mod_peptide','prot','rt']]
	    comA = comA.groupby('mod_peptide',as_index=False).mean()
	    comB = comB.groupby('mod_peptide',as_index=False).mean()
	    common =pd.merge(comA, comB , on=['mod_peptide'], how='inner')
	    # filtering outlier option
	    if filt_outlier == 1:
	    	filt_x,filt_y , pos_out = MD_removeOutliers(common['rt_y'].values,common['rt_x'].values)
            	data_B=  filt_x
            	data_A=  filt_y
	    	data_B=  np.reshape(data_B,[filt_x.shape[0],1])
            	data_A=  np.reshape(data_A,[filt_y.shape[0],1])
		logging.info( 'Outlier founded %i  w.r.t %i', pos_out.shape[0],common['rt_y'].shape[0] )
	    else:
	    	data_B=  common['rt_y'].values
	    	data_A=  common['rt_x'].values
	    	data_B=  np.reshape(data_B,[common.shape[0],1])
            	data_A=  np.reshape(data_A,[common.shape[0],1])
            logging.info('size common peptide , %i %i ',data_A.shape[0], data_B.shape[0]) 
	    clf = linear_model.RidgeCV(alphas=np.power(2,np.linspace(-30,30)),scoring='mean_absolute_error')
	    clf.fit(data_B,data_A)   
	    logging.info( 'alpha of the CV ridge regression model %4.4f',clf.alpha_)
	    clf_final= linear_model.Ridge(alpha=clf.alpha_)
	    clf_final.fit(data_B,data_A)  
	    logging.info('Mean abs Error on training : %4.4f sec', mean_absolute_error(data_A, clf_final.predict(data_B)) )
	    #print add_pep_frame['rt'].apply(lambda x: clf.predict(x).tolist()[0][0])
	    add_pep_frame.loc['rt'] = add_pep_frame['rt'].apply(lambda x: clf_final.predict(x)[0][0]) 
	    # version with retention time correction
	    #add_pep_frame.loc['rt']= add_pep_frame.apply( Var_corr, axis=1 )	
            final_out = pd.concat([ final_out,  add_pep_frame], join='outer', axis=0)
    
    logging.info( 'final set shape : %i',  final_out.shape[0])
    #print final_out.columns
    final_out['matched']=1
    final_out = final_out.drop_duplicates()
    logging.info('final set after drop duplicates pep : %i  ',final_out.shape[0] )
    # drop  strange row.
    final_out.dropna(axis=0, how='any', thresh=12, subset=None, inplace=True)
    print '-----------------'
    #concatenate with the original dataset
    logging.info('Before adding %s contains %i ', exp_set[jj],exp_t[jj].shape[0])
    exp_out[jj]=pd.concat([exp_t[jj], final_out], join='outer', axis=0)
    logging.info('After MBR %s contains:  %i  peptides', exp_set[jj] ,exp_out[jj].shape[0] )
    logging.info('----------------------------------------------')
## Save new file  new file
for a  in aa:
    num_split = 24
    name = os.path.split(exp_set[a])[1].split('.')[0]
    c=0
    for tmp in np.array_split(exp_out[a], num_split):
    #    print "output ", exp ,"part :" +str(c)
        tmp.to_csv(path_or_buf = output_dir +"/"+ name+ '_' + str(c)+".txt",sep="\t",header=True,index=False )
        c+=1
    #exp_out[a].to_csv(  output_dir  +"/"+name+"_ALL.txt" , sep="\t",index=False )

