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
import ConfigParser
import ast
import copy



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
def MD_removeOutliers(x, y, width):
    MD = MahalanobisDist(x, y)
    threshold = np.mean(MD) * width  # adjust 1.5 accordingly
    nx, ny, outliers = [], [], []
    for i in range(len(MD)):
        if MD[i] <= threshold:
            nx.append(x[i])
            ny.append(y[i])
        else:
            outliers.append(i) # position of removed pair
    return (np.array(nx), np.array(ny), np.array(outliers))

def combine_model( x,model,err,weight_flag):
        x =  x.values
        tot_err = np.sum( np.array(err) [ np.where(~np.isnan(x))])

        app_sum=0
        app_sum_2=0
        for ii in range(0,len(x)):
                if ~  np.isnan(x[ii]):
                        if int(weight_flag)==0:
                                app_sum = app_sum   +  ( model[ii].predict(x[ii])[0][0] )
                        else :
                                app_sum_2 = app_sum_2   +  ( model[ii].predict(x[ii])[0][0] *  ( float(err[ii]) / float(tot_err) ) )

        #" output weighted mean
        if int(weight_flag) == 1 :
                return  float(app_sum_2)
        else :
        # output not weight
                return   float(app_sum)/ float(np.where(~ np.isnan(x)  )[0].shape[0] )

def check_columns_name (col_list,col_must_have ):
	for c_name in col_must_have:
		if  not (c_name in col_list):
			# fail
			return 1
	# succes 
	return 0
## fixed flag for the outlier dev
def run_mbr( args):
	if   not (os.path.isdir(args.loc_in)):
		exit(str(args.loc_in) + '-->  inputF folder does not exist ! ')

	filt_outlier= args.out_flag

	if str(args.loc_in) =='':
		output_dir= 'mbr_output'
	else:
		if '/' in  str(args.loc_in):
			output_dir=  str(args.loc_in) + 'mbr_output'
		else:
			 exit(str(args.loc_in) +' EXIT inputF path not well specified --> / missing ')

	if  not (os.path.isdir(output_dir)):
		#print "Created output folder in ",output_dir
		os.makedirs(output_dir )
	else:
		print "Output folder in :",output_dir	

	
	config = ConfigParser.RawConfigParser()
        # it s always placed in same folder of moff_mbr.py
        config.read('moff_setting.properties')

	## read input
	exp_set=[]
	exp_t=[]
	exp_out=[]
	exp_subset=[]
	#if (args.sample) == None:
	for root, dirs, files  in  os.walk(args.loc_in) :
		for f in files :
			if f.endswith( '.' +  args.ext):
				exp_set.append( os.path.join(root, f))
	#old solutions 	
 	#paths = [os.path.join(args.loc_in,fn) for fn in next(os.walk(args.loc_in))[2] ]
        #exp_set=paths
	if not ( (args.sample) == None ) :
		exp_set_app=copy.deepcopy(exp_set)
		for a  in exp_set: 
			if (re.search(args.sample ,a) == None ) :
				exp_set_app.remove(a )
		exp_set = exp_set_app
	if exp_set ==[] :
                exit('ERROR input files not found. check the folder or the extension given in input')
	for a  in exp_set: 
		print 'Reading file.... ',a
		exp_subset.append(a)
		data_moff = pd.read_csv(a, sep="\t", header=0)
		list_name= data_moff.columns.values.tolist()
		if check_columns_name( list_name,  ast.literal_eval( config.get('moFF', 'col_must_have_x'))  ) ==1 :
			exit ('ERROR minimal field requested are missing or wrong')
		data_moff['matched']=0
		data_moff = data_moff.sort('rt')
		exp_t.append(data_moff)
		exp_out.append(data_moff)

		
			
	print 'Read input --> done '
	## parameter of the number of query
	exit('test ')
	## set a list of filed mandatory 
	#['matched','peptide','mass','mz','charge','prot','rt']

	exp_set=exp_subset
	aa=range(0,len(exp_t))
	out =list (itertools.product(aa,repeat=2))
	## just to save all the model
	model_save=[]
	model_err=[]
	model_status=[]
	# add matched columns
	list_name.append('matched')


	##input of the methods

	logging.basicConfig(filename=   args.loc_in +  '/'   +  args.log_label + '_' +'mbr_.log',filemode='w',level=logging.DEBUG)
	logging.info('Filtering is %s :', args.out_flag )
	logging.info( 'Pairwise model computation ----')
	for jj in aa:
	    sum_values=np.array([0,0,0])
	    print 'matching  in ', exp_set[jj]
	    
	    for i in out:
		if  i[0] == jj and i[1] != jj:
		    logging.info( '  matching  %s peptide in   searching in %s ', exp_set[i[0]] ,exp_set[i[1]])
		    list_pep_repA = exp_t[i[0]]['peptide'].unique()
		    list_pep_repB =  exp_t[i[1]]['peptide'].unique()
		    logging.info( '  Peptide unique  %i , %i ',  list_pep_repA.shape[0], list_pep_repB.shape[0])
		    set_dif_s_in_1 = np.setdiff1d(list_pep_repB,list_pep_repA)
		    add_pep_frame = exp_t[i[1]][exp_t[i[1]]['peptide'].isin(set_dif_s_in_1)].copy()
		    #print add_pep_frame.columns
		    pep_shared = np.intersect1d(list_pep_repA,list_pep_repB)
		    logging.info( '  Peptide to add size  %i ',  add_pep_frame.shape[0])
		    logging.info( '  Peptide shared  %i ',  pep_shared.shape[0])
		    comA =   exp_t[i[0]][ exp_t[i[0]]['peptide'].isin(pep_shared)][['peptide','prot','rt']]
		    comB =   exp_t[i[1]][ exp_t[i[1]]['peptide'].isin(pep_shared)][['peptide','prot','rt']]
		    comA = comA.groupby('peptide',as_index=False).mean()
		    comB = comB.groupby('peptide',as_index=False).mean()
		    common =pd.merge(comA, comB , on=['peptide'], how='inner')
		    if common.shape[0]  <= 30 :
			 #print common.shape
			 model_status.append(-1)
			 continue
		    # filtering outlier option
		    else:
				if int(args.out_flag)  == 1:
					filt_x,filt_y , pos_out = MD_removeOutliers(common['rt_y'].values,common['rt_x'].values,args.w_filt)
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
				logging.info(' size trainig shared peptide , %i %i ',data_A.shape[0], data_B.shape[0]) 
				clf = linear_model.RidgeCV(alphas=np.power(2,np.linspace(-30,30)),scoring='mean_absolute_error')
				clf.fit(data_B,data_A)   
				#logging.info( ' alpha of the CV ridge regression model %4.4f',clf.alpha_)
				clf_final= linear_model.Ridge(alpha=clf.alpha_)
				clf_final.fit(data_B,data_A)
				## save the model 
				model_save.append(clf_final)  
				model_err.append( mean_absolute_error(data_A, clf_final.predict(data_B)))
				logging.info(' Mean abs Error on training : %4.4f min', mean_absolute_error(data_A, clf_final.predict(data_B)) )
				model_status.append(1)

	logging.info( 'Combination of the  model  --------')
	logging.info('Weighted combination  %i : ', int( args.w_comb)  )

	diff_field = np.setdiff1d(exp_t[0].columns , ['matched','peptide','mass','mz','charge','prot','rt']) 
	for jj in   aa:
	    pre_pep_save=[]
	    print 'Predict rt for the exp.  in ', exp_set[jj]
	    for i in out:
		if  i[0] == jj and i[1] != jj:
		    logging.info('matching peptides found  in  %s ',exp_set[i[1]])
		    list_pep_repA = exp_t[i[0]]['peptide'].unique()
		    list_pep_repB =  exp_t[i[1]]['peptide'].unique()
		    set_dif_s_in_1 = np.setdiff1d(list_pep_repB,list_pep_repA)
		    add_pep_frame = exp_t[i[1]][exp_t[i[1]]['peptide'].isin(set_dif_s_in_1)].copy()
		    add_pep_frame= add_pep_frame[['peptide','mass','mz','charge','prot','rt']]
		    add_pep_frame['code_unique'] =  add_pep_frame['peptide']+'_'+  add_pep_frame['prot'] + '_' +  add_pep_frame['mass'].astype(str) + '_' +  add_pep_frame['charge'].astype(str)
		    add_pep_frame=add_pep_frame.groupby('code_unique',as_index=False)['peptide','mass','charge','mz','prot', 'rt'].aggregate(max)
		    add_pep_frame= add_pep_frame[['peptide','mass','mz','charge','prot','rt']]
		    pre_pep_save.append( add_pep_frame) 
	    #test = reduce(lambda left,right: pd.merge(left,right,on=['code_unique','peptide','pep_len','prot','mz','mass','charge'],how='outer'), pre_pep_save)
	    test= pd.merge( pre_pep_save[0] ,    pre_pep_save[1],on=['peptide','charge','prot','mass','mz'],how='outer',suffixes=['_x1','_y1'])
	    #test= pd.merge( test ,    pre_pep_save[2],on=['peptide','charge','prot','mass','mz'],how='outer',suffixes=['_x2','_y2'])
	    #test= pd.merge( test,    pre_pep_save[3],on=['peptide','charge','prot','mass','mz'],how='outer',suffixes=['_x3','_y3'])
	    #test= pd.merge( test,    pre_pep_save[4],on=['peptide','charge','prot','mass','mz'],how='outer',suffixes=['_x4','_y4'])
	    #test= pd.merge( test,    pre_pep_save[5],on=['peptide','charge','prot','mass','mz'],how='outer',suffixes=['_x5','_y5'])
	    #test= pd.merge( test,    pre_pep_save[6],on=['peptide','charge','prot','mass','mz'],how='outer',suffixes=['_x6','_y6'])
	    #test= pd.merge( test,    pre_pep_save[7],on=['peptide','charge','prot','mass','mz'],how='outer',suffixes=['_x7','_y7'])     
	    
	    test['rt']= test.ix[:,5:7].apply(  lambda  x: combine_model(x, model_save[(jj*2):((jj+1)*2)],model_err[ (jj*2):((jj+1)*2)],args.w_comb) , axis=1)
	    test['matched']= 1
	    for field in diff_field.tolist():
		test[field]= -1
	    test.drop('rt_x1', axis=1, inplace=True)
	    test.drop('rt_y1', axis=1, inplace=True)
	    ## print the entire file
	    #test.(path_or_buf= output_dir + '/' + str(exp_set[jj].split('.')[0].split('/')[1]) +'_match.txt',sep='\t',index=False)
	    logging.info('Before adding %s contains %i ', exp_set[jj],exp_t[jj].shape[0])
	    exp_out[jj]=pd.concat([exp_t[jj], test ]  , join='outer', axis=0)
	    logging.info('After MBR %s contains:  %i  peptides', exp_set[jj] ,exp_out[jj].shape[0] )
	    logging.info('----------------------------------------------')
	    exp_out[jj].head(10).to_csv(path_or_buf= output_dir + '/' + str(exp_set[jj].split('.')[0].split('/')[1]) +'_match.txt',sep='\t',index=False)
   	    ## forse logging handeles




if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='moFF match between run input parameter')

	parser.add_argument('--inputF', dest='loc_in', action='store', help='specify the folder of the input MS2 peptide list files ', required=False)

	parser.add_argument('--sample', dest='sample', action='store', help='specify which replicate is used fot mbr reg_exp are valid ', required=False)

	parser.add_argument('--ext', dest='ext', action='store', default='txt', help='specify the exstension of the input file (txt as default value) ', required=False)

	parser.add_argument('--log_file_name', dest='log_label', default='moFF',action='store', help='a label name for the log file (moFF_mbr.log as default log file name) ', required=False)

	parser.add_argument ('--filt_width', dest='w_filt', action='store',default=2,help='width value of the filter  k * mean(Dist_Malahobis)', required=False  )

	parser.add_argument('--out_filt', dest='out_flag', action='store',default=1,help='filter outlier in each rt time allignment',  required=False )

	parser.add_argument('--weight_comb', dest='w_comb', action='store',default=0,help='weights for model combination combination : 0 for no weight  1 weighted devised by model errors.', required=False )


	args = parser.parse_args()

        	
	run_mbr(args)
 
