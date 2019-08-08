#!/usr/bin/env python

import ast
import configparser
import copy
import itertools
import logging
import os
import re
import sys
from functools import reduce
from pyds import MassFunction
import bisect
import numpy as np
import pandas as pd
import scipy
from sklearn import linear_model
from sklearn.model_selection import ShuffleSplit
from sklearn.metrics import mean_absolute_error,r2_score
from sklearn.gaussian_process import GaussianProcessRegressor
from  sklearn.gaussian_process.kernels import RBF,WhiteKernel,DotProduct,ConstantKernel
from sklearn.model_selection import KFold
from scipy.interpolate import interp1d
import statsmodels.api as sm
import random as t_random

from numpy import random
import moff

"""moFF: matching between runs module """

# debug

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def find_ge(a, x,total):
	'Find leftmost item greater than or equal to x'
	i = bisect.bisect_left(a, x)
	if i >= total.shape[1]:
		return -1

	if (total[0,i] <  x) and (total[2,i] > x) :
		return a[i]
	else:
		if  (total[0,i-1] <  x) and (total[2,i-1] > x) :
			return a[i-1]
		else:
			return a[i +1]

def index(a, x):
    #'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def get_conflict_matrix_v2(bba_set,err):
    ## pre-filtering
    go_comb = False

    log.critical('initial number of bba %r',len(bba_set))


    #for jj in iter(cons):
        #[jj in bba_set[i].core() for i in range(1, len(bba_set))]
    lista=[]
    for i in range(0, len(bba_set)):

        lista_app=[]
        cons = bba_set[i].core()
        for j in range(0, len(bba_set)):
            if not (cons.intersection(bba_set[j].core())):
                #se non c e interazion aiungi
                lista_app.append(j)
        lista.append(lista_app)
    final=[]
    for  i in range(0, len(lista)):
        if len(lista[i]) <= len(lista)/2:
            final.append(i)
    if  final:
        n_1 = [x for x in range(0, len(bba_set)) if  (x in final)]
        bba_set = [bba_set[i] for i in n_1]
    '''
    finish_flag=False
    while not finish_flag:
        lista=[]
        cons = bba_set[0].core()
        for i in range(1, len(bba_set)):
            if not (cons.intersection(bba_set[i].core())):

            #if  not (jj in bba_set[i].frame()):
                lista.append(i)
            #lista.append(  [jj in bba_set[i].frame() for i in range(1, len(bba_set))]  )

        if not lista:
            finish_flag = True
        if len(lista) >= 1 :
            log.critical('First round check -- bba to filter : %r %b',(lista,finish_flag) )
            #lista.append(0)
            if (len(lista) >= len(bba_set)-1) :
                # caso dove lil primo é un outlier
                log.critical('running another around')
                n_1 = [x for x in range(0, len(bba_set)) if  (x in lista)]
                bba_set = [bba_set[i] for i in n_1]
            else:
                n_1 = [x for x in range(0, len(bba_set)) if not (x in lista)]
                bba_set = [bba_set[i] for i in n_1]
                finish_flag= True
    '''

    '''
      
    if len(lista) >= 1 and   len(lista) <= (len(bba_set) /2) :
        n_1 = [x for x in range(0, len(bba_set)) if not (x in lista)]
        log.critical('case 1: Elimino solo quelle che hano interazioni con la prima ! ')
        bba_set = [bba_set[i] for i in n_1]
        if len(bba_set) == 1:
            log.critical('Special exit 1 : after pre filtering just one bba left %r',lista)
            return lista, go_comb, bba_set
    if len(lista) >= 1 and   len(lista) > (len(bba_set) /2) :
        log.critical('case 2: Elimino la prima + il resto che ha interazioni con la prima ')
        n_1 = [x for x in range(0, len(bba_set)) if x in lista ]
        bba_set = [bba_set[i] for i in n_1]
    
    if len(lista)== len(bba_set):
        return lista,go_comb,bba_set
    '''
    log.critical('After pre-screening number of bba %r', len(bba_set))
    log.critical('Error training, %r',err)
    conflict_val=[]
    norma_comb=[]
    for i in range(0, len(bba_set)):
        # N -1 set
        n_1 = [ x  for x in range(0,len(bba_set)) if x != i ]
        if not (n_1):
            print ('Strange Cases',n_1)
            print(bba_set)
            exit('--  --')
            #print(lista)
        bba_input_n_1 = [bba_set[i] for i in n_1]
        res = conj_bba_set(bba_input_n_1 )
        norma_comb.append(res)
        if not res:
            # no combination
            conflict_val.append(1 )
        else:
            # yes combination
            conflict_part = [bba_set[i].pl({x}) * res.pl({x}) for x in iter(bba_set[i].frame().intersection(res.frame()))]
            if np.isnan(conflict_part).any():
                log.critical('NaN detect in conflict')
            if not conflict_part:
                conflict_val.append(1)
            else:
                conflict_val.append(1 - max(conflict_part))

    #" add another test
    go_comb = True
    ts = pd.Series(conflict_val)
    log.critical('Conflict values : %s', ts.to_string())
    first_step = list(ts[ts >= 0.80].index)
    if len(first_step) ==  ts.shape[0]:
        log.critical('max Conflict for all bba ')
        go_comb = False
        return [], go_comb, []
    if len(first_step) >= 1:
        # if
        left_bba = list(ts[ts < 0.80].index)
        log.critical('inside filtering bba left to use %r',len(left_bba))
        bba_input2 = [bba_set[i] for i in left_bba]
        log.critical(' after filtering bba left to  %r',len(bba_input2))
        left_comb = conj_bba_set_unnorm(bba_input2)
        conflict_left =  1 - max([left_comb.pl({x})  for x in iter(left_comb.frame())])
        if conflict_left < 0.80:
            final_bba_2comb = bba_input2
            log.critical('ready to combine !! ')
        else:
            final_bba_2comb = bba_input2
            log.critical('Second iteration shoudl be inplemented  !! ')
    else:
        # no filtering by conflict
        left_bba = list(ts.index)
        final_bba_2comb = [bba_set[i] for i in list(ts.index)]

    log.info('Subset of expert combined with Conj. Rule %r', left_bba)
    #log.critical('combination on %r', res)


    return left_bba, go_comb,final_bba_2comb


def get_conflict_matrix(bba_set):
    a = np.zeros(( len(bba_set) ,  len(bba_set )) )
    #list_inf= []
    #list_source_2cut= []
    go_comb=True
    if len(bba_set) ==1:
        return [], a,go_comb
    for i in range(0,len(bba_set)  ):
        for j in range(i+1, len(bba_set)):
            conflict_part = [bba_set[i].pl({x}) * bba_set[j].pl({x}) for x in iter(bba_set[i].frame().intersection(bba_set[j].frame()))]
            if not conflict_part:
                a[i, j]= 1
            else:
                a[i, j]= 1 - max(conflict_part)
    ## to copy
    aa_df = pd.DataFrame(a)
    log.info('%s',aa_df.to_string())
    index_to_drop= aa_df[aa_df == 1].sum(axis=0)[aa_df[aa_df == 1].sum(axis=0) >= int(len(bba_set) * 0.6)].index
    log.critical('Drop index %r', index_to_drop)
    aa_df.drop(index_to_drop, axis=1, inplace=True)
    aa_df.drop(index_to_drop, axis=0, inplace=True)
    aa_df.replace(0,np.nan, inplace= True )
    #summed_conflic = [  ( aa_df.iloc[:,x].sum() + aa_df.iloc[x,:].sum()) /aa_df.shape[0]   for x in aa_df.index ]
    #summed_conflic= np.array([(aa_df.iloc[:, x].sum() + aa_df.iloc[x, :].sum()) / (aa_df.shape[0]- 1) for x in range(0, aa_df.shape[0])])
    min_conflict=  np.array([np.nanmin( [aa_df.iloc[:, x].min(skipna=True) ,aa_df.iloc[x, :].min(skipna=True)] )  for x in range(0, aa_df.shape[0])])
    max_conflict=  np.array([np.nanmax( [aa_df.iloc[:, x].max(skipna=True) ,aa_df.iloc[x, :].max(skipna=True)] )  for x in range(0, aa_df.shape[0])])
    log.info('%r', max_conflict)
    log.info('%r', min_conflict)
    # use of series for a need of non-continuous  index
    ts = pd.Series(min_conflict, index=aa_df.index)
    log.critical('%s', ts.to_string())
    res= list(ts[ts <= 0.8].index)

    if  res  :
        return res , aa_df, go_comb
    else:
        app=[]
        return app, aa_df, go_comb

def conj_bba_set(bba_set):
    app=  bba_set[0]
    for ii in range(1,len(bba_set)):
            app = app.combine_conjunctive(bba_set[ii])
    return app


def conj_bba_set_unnorm(bba_set):
    printa_bba(bba_set)
    log.critical('inside conj unnormalized %r',len(bba_set))
    app=  bba_set[0]
    for ii in range(1,len(bba_set)):
            app = app.combine_conjunctive(bba_set[ii],normalization=False)
    return app



def printa_bba(input):
    if isinstance(input, list):
        for c,j in enumerate(input):
            log.info( '- Exp %i', c)
            for a in j.keys():
                log.info(' -  %r ->  %r ', a,j[a])
    else:
        log.info(' Combined bba  ')
        for a in input.keys():
            log.info(' -  %r ->  %r ', a, input[a])

def conj_combination(bba, log, intervall, radius, alpha):
    # print bba
    m_comb = conj_bba_set(bba)
    log.info('COMBINATION: Conjunctive comb. rule normalized')
    printa_bba(m_comb)
    #log.info('conflict %r', m_comb.local_conflict())
    log.info('COMBINATION: Decision on the belief space ')

    set_res = []
    # if not m_comb.pignistic()
    if not m_comb:
        print('STOP')

    if len(bba) ==1:
        log.critical('COMBINATION just one bba as output , no combination ')
        alpha= .70
        #res = np.array( list(m_comb.keys()))[np.where(list(m_comb.values()) == np.max(list((m_comb.values()))) )] [0]

    #else:
    min_set = 200
    final_set= None
    for x in iter(m_comb.focal()):
        #print(1 - alpha, x, m_comb.pl(x))
        if m_comb.pl(x) >= (1- alpha) :
            #print(1-alpha, x, m_comb.pl(x) )
            if (len(x)  < min_set) and (len(x) > 1) :
                #print ('finding some minimal pl',x, m_comb.pl(x) )
                min_set= len(x)
                final_set = x
    log.critical(final_set)


    if  final_set is None :
        res = np.array( list(m_comb.keys()))[np.where(list(m_comb.values()) == np.max(list((m_comb.values()))) )] [0]
    else:
        res = final_set
    #res = np.array( list(m_comb.keys()))[np.where(list(m_comb.values()) == np.max(list((m_comb.values()))) )] [0]
    # res = np.array(m_comb.pignistic().keys())[np.where(m_comb.pignistic().values() == np.max(m_comb.pignistic().values()))[0]]
    # log.info('pig_trans %r', m_comb.pignistic())
    # print 'res ::', res

    if len(res) == 1:
        log.info('COMBINATION: Final decison on %r', res)
        index = int(list(res)[0])
        output = intervall[1, index]
        lw_bound = intervall[0, index]
        up_bound = intervall[2, index]
    else:
        # " multiple best
        log.info('COMBINATION: Final decison on %r', res)
        for s in res:
            set_res.append(int(s))
        # set_res.append(s)
        lw_bound = intervall[0, sorted(set_res)[0]]
        up_bound = intervall[2, sorted(set_res)[-1:]][0]
        output = ((up_bound - lw_bound) / 2 ) +  lw_bound

    '''
    for s in  m_comb.pignistic().keys():
        #log.info('Pig.prob %.4f',  m_comb.pignistic()[s])
        if m_comb.pignistic()[s] >= max_set:
            max_set= m_comb.pignistic()[s]
            set_res.append(s)
            #print set_res


    log.info( 'choosen interval %s', list(set_res))
    #log.info( 'max prob %.4f', max_set)
    # integer  as hypothesis in the frame
    output = intervall[1, int(list(set_res)[0])]
    lw_bound = intervall[0, int(list(set_res)[0])]
    up_bound = intervall[2, int(list(set_res)[0])]

    '''

    # log.info( 'All intervall  %r',intervall[:,pos_inex_union])
    log.info('COMBINATION: final value %.4f | %.4f  %.4f ', output, lw_bound, up_bound)
    return lw_bound, up_bound, output

def combine_model_GP_consonant_bf( x,  model, log,intervall,r,err):
    debug__mode= True

    #lower_bound_ra, upper_bound_ra, ra_flag = refine_bound_rankAggr(x, ref_rep, shared_pep, data)
    #print 'int RA lower_bound', index(intervall[1, :], find_ge(intervall[1, :],  lower_bound_ra   , intervall))
    #print 'int RA upper_bound ', index(intervall[1, :], find_ge(intervall[1, :],  upper_bound_ra   , intervall))

    app_name = x['code_unique']
    log.info('feature %s',app_name)
    if  'K.AAPDVQLLMNDSQNDQSK.Q + Deamidated (NQ)_16130088__2' in app_name:
        #ESGIIQGDLIAK_1242.6819_2
        print ('stop and check')
    # RT_benchmark
    x = np.array(x[0:10])
    # RT_benchmak2
    #x = np.array(x[0:7])

    #interval_rankAggr = get_combinedInt_fromAggrRank(x, app_name.split('_')[0], ref_rep, shared_pep, data)
        # if interval_rankAggr != -1: than
        #bisect.bisect(intervall[1, :].tolist(), interval_rankAggr[0]) - 1
        #bisect.bisect(intervall[1, :].tolist(), interval_rankAggr[1]) - 1
    log.info ('radius in min %.4f # interval %r  ', r, intervall.shape)
    #pos_inex_union = define_frame_GP(x, model, intervall, 20)
    #out_map, out_map_set = create_dict_frame(pos_inex_union)
    #print pos_inex_union
    # rt_benchameark
    #x = np.array(x[0:10])
    #rt_benchmark2
    #x = np.array(x[0:6])

    bba_input=[]
    #log.info('frame final Theta %r %r', pos_inex_union, pos_inex_union.shape)
    for ii in range(0, len(x) ):
        if ~  np.isnan(x[ii] ):
            pred, y_cov = model[ii].predict(np.array(x[ii]).reshape(-1, 1), return_std=False, return_cov=True)
            var = np.sqrt(np.diag(y_cov))[0]
                #z_value=[0.06,0.13,0.26,0.39,0.53,0.68,0.85,1.04,1.29,1.65,1.965]
                #z_value=[ 0.13,  0.39,0.68,1.04,1.65]
            alpha = np.array([0.0,0.0,0.0,0.0,0.0, 1.0])
            delta_rt=[0.3,0.5,0.7,1,2,2.3] # [1.7,1,0.75,0.5,0.2]
            # delta_rt=[0.2,0.4,0.6,0.8,1]
            # alpha= scipy.stats.norm.sf(abs(z_scores))*2
            log.info('#  %i  model  predicted %r  var %r ', ii, pred, var)
            upper_list = []
            lower_list=[]
            m1 = MassFunction()
            log.info('---   ----')
            for cc, delta_val in enumerate(delta_rt):
                #pred + delta_val
                upper =index(intervall[1, :] ,find_ge(intervall[1, :],  pred + delta_val  , intervall)  )
                lower  = index(intervall[1, :], find_ge(intervall[1, :], pred - delta_val, intervall))
                # print alpha[cc], 1-alpha[cc],  pred - ( z_value[cc]* np.sqrt(var)), pred, pred + (   z_value[cc]  * np.sqrt(var))
                upper_list.append(upper)
                lower_list.append(lower)
                z_score =  ( (pred + delta_val) - pred  )/  np.sqrt(var)
                alpha[cc] = scipy.stats.norm.sf(abs(z_score)) * 2

                #log.info('# %i alpha %r  1-alpha %r lower_limit %r  %r  upper_limit %r ', cc,alpha[cc], 1-alpha[cc],  pred - ( z_score* np.sqrt(var)), pred, pred + (   z_score  * np.sqrt(var))  )
                log.info ('# %i model rt delta (min): %r  interval_lower %r --- interval_upper %r alpha %r  1-alpha %r', cc,delta_val, lower, upper,alpha[cc], 1-alpha[cc])
                _alpha = 1 - alpha
            for cc,alpha_val in enumerate(alpha):
                if alpha_val == 1:
                    lower_cur -= 1
                    upper_cur += 1
                else:
                    lower_cur = lower_list[cc]
                    upper_cur = upper_list[cc]
                    #print alpha[cc], 1-alpha[cc],  pred - ( z_value[cc]* np.sqrt(var)), pred, pred + (   z_value[cc]  * np.sqrt(var))
                    #print 'interval lw :', index(intervall[1, :], find_ge(intervall[1, :], pred - (  z_value[cc] * np.sqrt(var)), intervall)) , 'interval uw :',index(
                    #intervall[1, :], find_ge(intervall[1, :], pred + ( z_value[cc] * np.sqrt(var)), intervall))
                pos_index = np.arange(lower_cur, upper_cur + 1)
                if (cc == 0) or (cc == len(alpha)- 1 ):
                    if ((cc == 0)):
                        m1[[str(a) for a in pos_index.tolist()]] = _alpha[cc]
                    else:
                        if _alpha[cc] < 1 :
                            m1[[str(a) for a in pos_index.tolist()]] = 1 - _alpha[cc]
                else:
                    if frozenset([str(a) for a in  pos_index.tolist() ]) in m1.keys()  :
                        m1[[str(a) for a in pos_index.tolist()]] = m1[[str(a) for a in pos_index.tolist()]]  +  ( _alpha[cc] - _alpha[cc-1])
                    else:
                        m1[[str(a) for a in pos_index.tolist()]] = _alpha[cc] - _alpha[cc-1]
            bba_input.append(m1)

    printa_bba( bba_input)
    '''
    [print(bba_input[0].pl({x}),x) for x in iter(bba_input[0].frame()) ]
    [print(bba_input[0].pl({x}), bba_input[1].pl({x}),x) for x in iter(bba_input[0].frame().intersection(bba_input[1].frame() ))]
    [(bba_input[0].pl({x}), bba_input[4].pl({x}), bba_input[0].pl({x}) * bba_input[4].pl({x}), x) for x in iter(bba_input[0].frame().intersection(bba_input[4].frame()))]
        [(bba_input[0].pl({x}), bba_input[4].pl({x}), bba_input[0].pl({x}) * bba_input[4].pl({x}), x) for x in iter(bba_input[0].frame().intersection(bba_input[4].frame()))]

    '''

    ii_index_2=[]
    if (len(bba_input) == 1 ) :
        lw_bound, up_bound, output = conj_combination(bba_input, log, intervall, r,0.05)
        combination = True
    else:
        ii_index_2, comb_flag, bba_2comb = get_conflict_matrix_v2(bba_input,err)
        #print(ii_index_2)
        #print (bba_2comb)
        if len(ii_index_2) >= 1 and comb_flag :
            lw_bound, up_bound, output = conj_combination(bba_2comb, log, intervall, r,0.05)
            combination = True
        else:
            combination= False
        #print conf_x
        #print ii_index_2
    '''
    if ( (len(ii_index_2) >= 1 )or (len(bba_2comb) == 1 ) ) and (comb_flag):
        # uso la ConJ rule if two or more have  same common intervall
        if  (len(ii_index_2) >= 1 ) :
            # more than one sources
            #log.info('Subset of expert combined with Conj. Rule %r', ii_index_2)
            #bba_input = [bba_input[i] for i in ii_index_2]
            lw_bound, up_bound, output = conj_combination(bba_2comb, log, intervall, r)
        if ( not ii_index_2 ) and (len(bba_2comb) == 1 ):
            # case for just one szources
            lw_bound, up_bound, output = conj_combination(bba_2comb, log, intervall,  r)
        combination = True
        # else:
        #   output = disj_combination (bba_input, log, intervall, out_map_set,pos_inex_union)
    else:
        # empty set for intersection and  nont sincgel sources case
        #lw_bound, up_bound, output = disj_combination(bba_input, log, intervall, None, pos_inex_union)
        combination = False
        # print 'stop here please', app_name, output
    '''
    log.info('----- ----- -----')
    # print output
    # " output basic check control
    if combination:
        return pd.Series({'time_ev': output,
                          'ev_upper': up_bound, 'ev_lower':lw_bound , 'conflict_flag':0})
    else:
        # output not weight
        return pd.Series({'time_ev':-1,'ev_upper':-1, 'ev_lower':-1,'conflict_flag':1})


def	 create_belief_RT_interval (max_rt, min_rt,n_interval):
	# print max_rt, min_rt, float(max_rt-min_rt) / float(20)
	off_set = float(max_rt-min_rt) / float(n_interval)
	#print 'length_interval',off_set
	interval_mat  = np.zeros(shape=(3,n_interval))
	for i in range(0,n_interval):
		#print i , min_rt + (i * off_set )
		interval_mat[0,i] = min_rt + (i * off_set )
		interval_mat[2,i] = (min_rt + (i * off_set )) +  (( float(max_rt-min_rt) / float(n_interval)  ) -0.001)
		interval_mat[1,i] = (interval_mat[0,i] + interval_mat[2,i] ) /2

	#print interval_mat[:,0], interval_mat[:,19]
	return interval_mat,off_set

# filtering _outlier
def MahalanobisDist(x, y):
    """
    Computee the Mahalanobis distance to filter outlier in the RT allignment
    :param x:
    :param y:
    :return:
    """
    covariance_xy = np.cov(x, y, rowvar=0)
    inv_covariance_xy = np.linalg.inv(covariance_xy)
    xy_mean = np.mean(x), np.mean(y)
    x_diff = np.array([x_i - xy_mean[0] for x_i in x])
    y_diff = np.array([y_i - xy_mean[1] for y_i in y])
    diff_xy = np.transpose([x_diff, y_diff])
    md = []
    for i in range(len(diff_xy)):
        md.append(np.sqrt(
            np.dot(np.dot(np.transpose(diff_xy[i]), inv_covariance_xy), diff_xy[i])))
    return md


# remove outlier
def MD_removeOutliers(x, y, width):
    """
    Remove outliers point using MahalanobisDist function
    :param x:
    :param y:
    :param width:
    :return:
    """
    MD = MahalanobisDist(x, y)
    threshold = np.mean(MD) * float(width)  # adjust 1.5 accordingly
    nx, ny, outliers = [], [], []
    for i in range(len(MD)):
        if MD[i] <= threshold:
            nx.append(x[i])
            ny.append(y[i])
        else:
            outliers.append(i)  # position of removed pair
    return np.array(nx), np.array(ny), np.array(outliers)


# combination of rt predicted by each single model
def combine_model(x, model, err, weight_flag):

    #x = x.values
    x = np.array(x[0:10])

    tot_err = np.sum(np.array(err)[np.where(~np.isnan(x))])

    app_sum = 0
    app_sum_2 = 0
    for ii in range(0, len(x)):
        if ~  np.isnan(x[ii]):
            if not weight_flag :
                app_sum = app_sum + (model[ii].predict(x[ii])[0][0])
            else:
                app_sum_2 = app_sum_2 + \
                    (model[ii].predict(x[ii])[0][0] *
                     (float(err[ii]) / float(tot_err)))

                # " output weighted mean
    if not  weight_flag :
        # not weight outpuy
        #return float(app_sum) / float(np.where(~ np.isnan(x))[0].shape[0])
        return pd.Series({'time_pred': float(app_sum) / float(np.where(~ np.isnan(x))[0].shape[0]), 'uncertainty_win': 0})
    else:
        # output weight
        #return float(app_sum_2)
        return pd.Series({'time_pred': float(app_sum_2), 'uncertainty_win': 0})




def train_gp_mono(data_A,data_B,c=None):
    """
    Using GP for retention time alligment
    """
    log.critical('Target data shape %r  Input Data shape %r ', data_A.shape[0], data_B.shape[0] )
    bins = np.linspace(data_B.min(), data_B.max(),10)
    digitized = np.digitize(data_B, bins)
    size_bin = np.array([  digitized[digitized == i].shape[0]  for i in range(1, len(bins)+1)])
    k = 0.9
    w = 1 - ( k * ( (size_bin - size_bin.min()) / (  size_bin.max() - size_bin.min()) ) )
    log.info('Using k: %4.4f  --> Sampling weights from bins size : %r ', k, w  )
    clean_lista = [e for e in  [t_random.sample(set(np.where(digitized == i)[0]), int(size_bin[i - 1] * w[i - 1])) for i in
                    range(1, len(bins) + 1)] if e]
    tt_x = np.concatenate(clean_lista)
    test =  np.setdiff1d(range(data_A.shape[0]),tt_x)
    #ff = ym_test_predicted[:, 0] - np.sqrt(np.diag(y_cov))
    #dd = ym_test_predicted[:, 0] + np.sqrt(np.diag(y_cov))
    # plt.fill_between(data_B[test,0], ff  , dd , alpha=0.5,color='k')
    #size_train= int (data_A.shape[0] * 0.10)
    #rows = random.sample(range(data_A.shape[0]),size_train)
    #data_A= data_A[rows,:]
    #data_B= data_B[rows,:]
    # data_B is x
    # data_A if Y
    #kernel =  ConstantKernel(1.0, (1e-2, 1e2))  * DotProduct() +  WhiteKernel(noise_level=5,noise_level_bounds=(1e-3, 1e5))
    kernel = ConstantKernel(1.0, (1e-2, 1e2)) * RBF(length_scale=1, length_scale_bounds=(1e-5, 1e5)) + WhiteKernel(noise_level=0.2,noise_level_bounds=(1e-5, 1e2))
    m = GaussianProcessRegressor(kernel=kernel, alpha=0.1, normalize_y=False, n_restarts_optimizer=1).fit(data_B[tt_x], data_A[tt_x])

    ym_train_predicted, y_cov_train = m.predict(data_B[tt_x], return_std=False, return_cov=True)
    ym_test_predicted, y_cov_test = m.predict(data_B[test], return_std=False, return_cov=True)
    '''
    #ff =  np.sqrt(np.diag(y_cov_test))
    #dd = np.sqrt(np.diag(y_cov_test))

    ff = ym_test_predicted[:, 0] - (1.96 * np.sqrt(np.diag(y_cov_test)))
    dd = ym_test_predicted[:, 0] + ( 1.96 * np.sqrt(np.diag(y_cov_test)) )

    ##printing modell
    plt.figure(figsize=(15, 14), dpi=100)
    #plt.scatter(data_B[test],ym_test_predicted,marker='*',c='red',s=15 )

    plt.fill_between(data_B[test, 0], ff, dd, alpha=0.5, color='r')
    plt.scatter(data_B[test],data_A[test],marker='h',c='black',s=25,label='True RT' )
    plt.scatter(data_B[test], ym_test_predicted, marker='*', c='blue', s=25,label='predicted RT')
    plt.legend(loc="best", scatterpoints=1,  prop={'size': 18})

    plt.title('GP with Rbf_+ whiteNoise on test set: ' + c )
    plt.savefig( 'D:\\workspace\\ionstar_dataset\\mbr_output\\' + c + '__model.png' )
    '''
    log.info(' MONO Size sampled training set : %i perc/total %4.4f ', tt_x.shape[0], tt_x.shape[0]/data_A.shape[0]  )
    log.info(' MONO Size remain validation set : %i  perc/total %4.4f  ', test.shape[0], test.shape[0]/data_A.shape[0]  )
    log.info(' MONO Mean absolute error training : %4.4f sec',mean_absolute_error(data_A[tt_x], ym_train_predicted))
    log.info(' MONO R^2 training : %4.4f sec',r2_score(data_A[tt_x], ym_train_predicted))
    log.info(' MONO  -- Mean absolute error validation set (not sampled point) : %4.4f sec',mean_absolute_error(data_A[test], ym_test_predicted))
    log.info(' MONO  -- R^2 error validation set (not sampled point) : %4.4f sec',r2_score(data_A[test], ym_test_predicted))

    return m, ym_train_predicted,mean_absolute_error(data_A[tt_x], ym_train_predicted)


def train_loess_gp( data_A,data_B,c=None):
    log.critical('Target data shape %r  Input Data shape %r ', data_A.shape[0], data_B.shape[0])
    rs = ShuffleSplit(n_splits=1, test_size=.80, random_state=0)
    for train_index, test_index in rs.split(data_B):
        lowess = sm.nonparametric.lowess(data_A[train_index].flatten(), data_B[train_index].flatten(), frac=.4)

        lowess_x = list(zip(*lowess))[0]
        lowess_y = list(zip(*lowess))[1]

        kernel = ConstantKernel(1.0, (1e-2, 1e2)) * RBF(length_scale=1, length_scale_bounds=(1e-5, 1e5)) + WhiteKernel(
            noise_level=0.2, noise_level_bounds=(1e-5, 1e2))
        m = GaussianProcessRegressor(kernel=kernel, alpha=0.1, normalize_y=False, n_restarts_optimizer=1).fit(
                                                                                                np.array(lowess_x).reshape(-1,1),
                                                                                                np.array(lowess_y).reshape(-1, 1))
        ym_train_predicted, y_cov_train = m.predict(data_B[train_index], return_std=False, return_cov=True)
        ym_test_predicted, y_cov_test = m.predict(data_B[test_index], return_std=False, return_cov=True)

    log.critical(' Size sampled training set : %i perc/total %4.4f ', train_index.shape[0], train_index.shape[0] / data_A.shape[0])
    log.critical(' Size remain validation set : %i  perc/total %4.4f  ', test_index.shape[0],
             test_index.shape[0] / data_A.shape[0])
    log.critical(' Mean absolute error training no_smoothed: %4.4f sec', mean_absolute_error(data_A[train_index], ym_train_predicted))
    log.critical(' R^2 training : %4.4f sec', r2_score(data_A[train_index], ym_train_predicted))
    log.critical('  -- Mean absolute error validation set (not sampled point) : %4.4f sec',
             mean_absolute_error(data_A[test_index], ym_test_predicted))
    log.critical('  -- R^2 error validation set (not sampled point) : %4.4f sec', r2_score(data_A[test_index], ym_test_predicted))

    return m, ym_train_predicted, mean_absolute_error(data_A[train_index], ym_train_predicted)


def train_gp(data_A,data_B,c=None):
    """
    Using GP for retention time alligment 
    """
    log.critical('Target data shape %r  Input Data shape %r ', data_A.shape[0], data_B.shape[0] )
    bins = np.linspace(data_B.min() , data_B.max()  ,10)
    digitized = np.digitize(data_B, bins)
    size_bin = np.array([  digitized[digitized == i].shape[0]  for i in range(1, len(bins)+1)])
    k=0.9
    w = 1 - ( k * ( (size_bin - size_bin.min()) / (  size_bin.max() - size_bin.min()) ) )
    log.info('Using k: %4.4f  --> Sampling weights from bins size : %r ', k, w  )
    clean_lista= [e for e in [t_random.sample(set(np.where(digitized == i)[0]), int(size_bin[i-1] *  w[i-1]))  for i in range(1, len(bins)+1)] if e]
    tt_x = np.concatenate(clean_lista)
    test =  np.setdiff1d(range(data_A.shape[0]),tt_x)
    #ff = ym_test_predicted[:, 0] - np.sqrt(np.diag(y_cov))
    #dd = ym_test_predicted[:, 0] + np.sqrt(np.diag(y_cov))
    # plt.fill_between(data_B[test,0], ff  , dd , alpha=0.5,color='k')
    #size_train= int (data_A.shape[0] * 0.10)
    #rows = random.sample(range(data_A.shape[0]),size_train)
    #data_A= data_A[rows,:]
    #data_B= data_B[rows,:]
    # data_B is x
    # data_A if Y
    #kernel =  ConstantKernel(1.0, (1e-2, 1e2)) * DotProduct() +  WhiteKernel(noise_level=5,noise_level_bounds=(1e-3, 1e5))
    kernel = ConstantKernel(1.0, (1e-2, 1e2)) * RBF(length_scale=1, length_scale_bounds=(1e-5, 1e5)) + WhiteKernel(noise_level=0.2,noise_level_bounds=(1e-5, 1e2))
    m = GaussianProcessRegressor(kernel=kernel, alpha=0.1, normalize_y=False, n_restarts_optimizer=1).fit(data_B[tt_x], data_A[tt_x])

    ym_train_predicted, y_cov_train = m.predict(data_B[tt_x], return_std=False, return_cov=True)
    ym_test_predicted, y_cov_test = m.predict(data_B[test], return_std=False, return_cov=True)
    '''
    #ff =  np.sqrt(np.diag(y_cov_test))
    #dd = np.sqrt(np.diag(y_cov_test))

    ff = ym_test_predicted[:, 0] - (1.96 * np.sqrt(np.diag(y_cov_test)))
    dd = ym_test_predicted[:, 0] + ( 1.96 * np.sqrt(np.diag(y_cov_test)) )

    ##printing modell
    plt.figure(figsize=(15, 14), dpi=100)
    #plt.scatter(data_B[test],ym_test_predicted,marker='*',c='red',s=15 )

    plt.fill_between(data_B[test, 0], ff, dd, alpha=0.5, color='r')
    plt.scatter(data_B[test],data_A[test],marker='h',c='black',s=25,label='True RT' )
    plt.scatter(data_B[test], ym_test_predicted, marker='*', c='blue', s=25,label='predicted RT')
    plt.legend(loc="best", scatterpoints=1, prop={'size': 18})
    plt.plot(data_B[tt_x],ym_train_predicted )

    plt.title('GP with Rbf_+ whiteNoise on test set: ' + c )
    plt.savefig( 'D:\\workspace\\ionstar_dataset\\mbr_output\\' + c + '__model.png' )
    '''
    log.info(' Size sampled training set : %i perc/total %4.4f ', tt_x.shape[0], tt_x.shape[0] / data_A.shape[0])
    log.info(' Size remain validation set : %i  perc/total %4.4f  ', test.shape[0],
             test.shape[0] / data_A.shape[0])
    log.info(' Mean absolute error training : %4.4f sec',mean_absolute_error(data_A[tt_x], ym_train_predicted))
    log.info(' R^2 training : %4.4f sec',r2_score(data_A[tt_x], ym_train_predicted))
    log.info('  -- Mean absolute error validation set (not sampled point) : %4.4f sec',mean_absolute_error(data_A[test], ym_test_predicted))
    log.info('  -- R^2 error validation set (not sampled point) : %4.4f sec',r2_score(data_A[test], ym_test_predicted))

    return m, ym_train_predicted,mean_absolute_error(data_A[tt_x], ym_train_predicted)

def combine_model_GP_mono(x, model):
    log.critical(' mono model predicting : %r',x.shape)
    pred, y_cov = model.predict( x, return_std=False, return_cov=True)
    var = np.sqrt(np.diag(y_cov))
    return (pred,  1.965 * var )

def combine_model_GP(x, model, err, weight_flag):
    """
    Combination of GP model
    """
    ra_flag = 0
    look_pep = x.mod_peptide
    '''
    lista_shared = exp_t[ii+1][(exp_t[ii+1]['rt'] > pred[0][0]  -1 ) & (exp_t[ii+1]['rt'] < pred[0][0] +1)]['mod_peptide']
    list_2 = np.intersect1d(rt_df.index,lista_shared)

   ply.scatter(rt_df.loc[list_2].columns,rt_df.loc[list_2].mean(axis=0)  )

   ub= rt_df.loc[list_2]['s7'].sort_values().index[0]
   lb= rt_df.loc[list_2]['s7'].sort_values().index[-1]
    '''
    # case RT_dataset  x = np.array(x[7:17])
    # case isofix  x = np.array(x[0:19])
    # ionstar
    #x = np.array(x[0:19])
    #  Rt datasetbenchmark 2
    #x = np.array(x[0:6])
    #  Rt datasetbenchmark 1
    x = np.array(x[0:10])
    log.info('%s -- missing value %i over %i ' , look_pep,  np.where(pd.isnull(x))[0].shape[0] , x.shape[0] )
    log.info('value %r std: %2.2f  mean: %2.2f', x.tolist(), np.nanstd(x,dtype='float'), np.nanmean(x,dtype='float'))
    # tot_err =  1- ( (np.array(err)[np.where(~np.isnan(x))]) / np.max(np.array(err)[np.where(~np.isnan(x))]))
    tot_err = np.sum(np.array(err)[np.where(~pd.isnull(x))])
    # print tot_err
    # print x
    app_sum = 0
    app_sum_2 = 0
    app_var = 0
    lista_shared = []
    for ii in range(0, len(x)):

        if ~  np.isnan(x[ii]):

            # sklearn
            pred, y_cov = model[ii].predict(np.array(x[ii]).reshape(-1, 1), return_std=False, return_cov=True)
            var = np.sqrt(np.diag(y_cov))[0]
            # pred, var = model[ii].predict(x[ii].reshape(1, 1), include_likelihood=True)
            # ress = model[ii].predict_quantiles(x[ii].reshape(1, 1))

            # print ' %i Input Rt  %4.4f  Predicted: %4.4f Var %4.4f  Interval at 95 %4.4f <--> %4.4f  ' % (
            # ii, x[ii], float(pred), float(var), pred - (1.965 * var) ,  pred + (1.965 * var)  )
            # print 'intervall width length %4.4f' % abs((pred - (1.965 * var)) - (pred + (1.965 * var)))


            if not weight_flag:
                app_sum = app_sum + (pred[0][0])
                app_var += var
            else:
                # print ii,model[ii].predict(x[ii])[0][0]
                w = (float(err[ii]) / float(tot_err))
                # w= tot_err[ii]
                # print ii ,'weighted', (model[ii].predict(x[ii])[0][0] * w ),w
                app_sum_2 = app_sum_2 + (pred[0][0] * w)
                app_var += var
    # " output weighted mean

    if weight_flag:
        f_p = app_sum_2
        mean_var = float(app_var) / float(np.where(~ np.isnan(x))[0].shape[0])
    else:
        mean_var = float(app_var) / float(x[~pd.isnull(x)].shape[0])
        f_p = float(app_sum) / float(x[~pd.isnull(x)].shape[0])

    # log.info( '   -->  Final Aggr Predicted: %4.4f Var %4.4f  Interval at 95 %4.4f <--> %4.4f  ' % (  f_p , float(mean_var), f_p - (1.965 * mean_var), f_p + (1.965 * mean_var)))
    # log.info ('   -->  Final intervall width length  %4.4f  |  #_inputs %i' %  (abs(    (f_p - (1.965 * mean_var)) - (f_p + (1.965 * mean_var)) ),int(np.where(~ np.isnan(x))[0].shape[0])) )

    return pd.Series({'time_pred': f_p, 'uncertainty_win': 1.965 * mean_var })



def combine_model_GP_rankAgg(x, model, err, weight_flag,rank_df,rt_df,exp_t,ref_rep):
    """
    Combination of GP model 
    """
    ra_flag= 0
    look_pep = x.mod_peptide
    '''
    lista_shared = exp_t[ii+1][(exp_t[ii+1]['rt'] > pred[0][0]  -1 ) & (exp_t[ii+1]['rt'] < pred[0][0] +1)]['mod_peptide']
    list_2 = np.intersect1d(rt_df.index,lista_shared)
    
   ply.scatter(rt_df.loc[list_2].columns,rt_df.loc[list_2].mean(axis=0)  )
   
   ub= rt_df.loc[list_2]['s7'].sort_values().index[0]
   lb= rt_df.loc[list_2]['s7'].sort_values().index[-1]
    '''
    x = np.array(x[7:17])
    #tot_err =  1- ( (np.array(err)[np.where(~np.isnan(x))]) / np.max(np.array(err)[np.where(~np.isnan(x))]))
    tot_err = np.sum(np.array(err)[np.where(~pd.isnull(x))])
    #print tot_err
    #print x
    app_sum = 0
    app_sum_2 = 0
    app_var =0
    lista_shared=[]
    for ii in range(0, len(x)):


        if ~  np.isnan(x[ii]):



            #sklearn
            pred, y_cov = model[ii].predict(np.array(x[ii]).reshape(-1,1), return_std=False, return_cov=True)
            var = np.sqrt(np.diag(y_cov))[0]
            #pred, var = model[ii].predict(x[ii].reshape(1, 1), include_likelihood=True)
            #ress = model[ii].predict_quantiles(x[ii].reshape(1, 1))

            #print ' %i Input Rt  %4.4f  Predicted: %4.4f Var %4.4f  Interval at 95 %4.4f <--> %4.4f  ' % (
            #ii, x[ii], float(pred), float(var), pred - (1.965 * var) ,  pred + (1.965 * var)  )
            #print 'intervall width length %4.4f' % abs((pred - (1.965 * var)) - (pred + (1.965 * var)))

            if ii >= ref_rep:
                c = ii + 1
            else:
                c = ii


            #lista_shared.append( exp_t[c][(exp_t[c]['rt'] > pred[0][0] - 3) & (exp_t[c]['rt'] < pred[0][0] + 3)][
            #    'mod_peptide'].values)
            ff= np.intersect1d(rt_df.index,exp_t[c][(exp_t[c]['rt'] > x[ii] - 1) & (exp_t[c]['rt'] < x[ii] + 1 )][
                'mod_peptide'].values)
            if  ff.shape[0] >= 4:
                lista_shared.append(ff)
            else:
                lista_shared.append(np.nan)

            if not weight_flag :
                app_sum = app_sum + (pred[0][0])
                app_var += var
            else:
                # print ii,model[ii].predict(x[ii])[0][0]
                w = (float(err[ii]) / float(tot_err))
                # w= tot_err[ii]
                # print ii ,'weighted', (model[ii].predict(x[ii])[0][0] * w ),w
                app_sum_2 = app_sum_2 + (pred[0][0] * w)
                app_var += var
    # " output weighted mean

    if pd.isnull(lista_shared).any() :
        ub = -1
        lb = -1
    else:
        top_x = reduce(np.union1d, lista_shared)
        res= rt_df.loc[top_x]['s1'].sort_values()
        ub= res.head(1).values[0]
        lb= res.tail(1).values[0]
    if  pd.isnull(x).all()  :
        return pd.Series({'time_pred': -1 , 'uncertainty_win': -1})

    if  weight_flag :
        f_p = app_sum_2
        mean_var = float(app_var) / float(np.where(~ np.isnan(x))[0].shape[0])
    else:
        mean_var = float(app_var) / float(x[~pd.isnull(x)].shape[0])
        f_p = float(app_sum) / float(x[~pd.isnull(x)].shape[0])

    #log.info( '   -->  Final Aggr Predicted: %4.4f Var %4.4f  Interval at 95 %4.4f <--> %4.4f  ' % (  f_p , float(mean_var), f_p - (1.965 * mean_var), f_p + (1.965 * mean_var)))
    #log.info ('   -->  Final intervall width length  %4.4f  |  #_inputs %i' %  (abs(    (f_p - (1.965 * mean_var)) - (f_p + (1.965 * mean_var)) ),int(np.where(~ np.isnan(x))[0].shape[0])) )

    return pd.Series({'time_pred': f_p, 'uncertainty_win': 1.965 * mean_var,'rank_up_bound':ub,'rank_low_bound':lb})

    #return pd.Series({'time_pred': f_p, 'uncertainty_win': 1.965 * mean_var})


def rank_matrix(intersect_share,data):
    intersect_share_local = reduce(np.intersect1d, ([x['mod_peptide'].unique() for x in data]))
    rt_val = np.zeros((len(intersect_share_local), 11))

    c= 0
    for a in data:
        rt_val[:,c] =a[a.mod_peptide.isin(intersect_share_local)][['mod_peptide', 'rt']].groupby('mod_peptide', as_index=True)['rt'].mean().values
        #rank_val[:, c]=a[a.mod_peptide.isin(intersect_share)][['mod_peptide', 'rt']].groupby('mod_peptide', as_index=True)[
        #    'rt'].mean().rank().values
        c+=1
    th = np.percentile(rt_val.std(axis=1), 70)
    index = np.where(rt_val.std(axis=1) <= th)
    drt = pd.DataFrame(rt_val[index[0], :], index=intersect_share_local[index[0]],
                       columns=['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10', 's11'])

    rank_val = np.zeros((len(intersect_share_local[index[0]]), 11), dtype=int)
    c = 0
    for a in data:
        rank_val[:, c] = \
        a[a.mod_peptide.isin(intersect_share_local[index[0]])][['mod_peptide', 'rt']].groupby('mod_peptide', as_index=True)[
            'rt'].mean().rank().values
        c += 1
    dd = pd.DataFrame(rank_val, index=intersect_share_local[index[0]],
                      columns=['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10', 's11'])

    return dd,drt


# run the mbr in moFF : input  ms2 identified peptide   output csv file with the matched peptides added


def run_mbr(args):
    """
    Macthing Between Run module.
    :param args:
    :return:
    """
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    log.addHandler(ch)

    if args.loc_in is None:
        # the user uses --inputtsv option
        if not (args.loc_out is None):
            # if the user use --output_folder the mbr folder will be created there
            output_dir = os.path.join(args.loc_out, 'mbr_output')
        else:
            # if the user does not use  --output_folder the mbr folder will be created on moFF path location
            output_dir = os.path.join('mbr_output')
            print(os.path.abspath(output_dir))

    else:
        # the user use the --inputF option
        if os.path.exists(os.path.join(args.loc_in)):
            # if '/' in  str(args.loc_in):
            output_dir = os.path.join(args.loc_in, 'mbr_output')
        else:
            exit(os.path.join(args.loc_in) +
                 ' EXIT input folder path is not well specified --> / missing or wrong path')

            # if not (os.path.isdir(args.loc_in)):
            #   exit(str(args.loc_in) + '-->  input folder does not exist ! ')

            # if str(args.loc_in) == '':
            #    output_dir = 'mbr_output'
            # else:
            #    if os.path.exists(os.path.join(args.loc_in)):
            # if '/' in  str(args.loc_in):
    # output_dir = os.path.join(args.loc_in, 'mbr_output')
    #    else:
    #        exit(os.path.join(args.loc_in) + ' EXIT input folder path not well specified --> / missing ')

    if not (os.path.isdir(output_dir)):

        log.critical("Created MBR output folder in : %s ",
                     os.path.abspath(output_dir))
        os.makedirs(output_dir)
    else:
        log.critical("MBR Output folder in : %s ", os.path.abspath(output_dir))
    # set log to file
    w_mbr = logging.FileHandler(os.path.join(
        output_dir, args.log_label + '_mbr_.log'), mode='w')
    w_mbr.setLevel(logging.INFO)
    log.addHandler(w_mbr)

    moff_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    config = configparser.RawConfigParser()
    config.read(os.path.join(moff_path, 'moff_setting.properties'))

    # it s always placed in same folder of moff_mbr.py
    # read input
    # comment better
    # name of the input file
    exp_set = []
    # list of the input dataframe
    exp_t = []
    # list of the output dataframe
    exp_out = []
    # lsit of input datafra used as help
    exp_subset = []
    # list of the name of the mbr output
    exp_out_name = []

    if args.loc_in is None:
        for id_name in args.tsv_list:
            exp_set.append(id_name)
    else:
        for item in os.listdir(args.loc_in):

            if os.path.isfile(os.path.join(args.loc_in, item)):
                if os.path.join(args.loc_in, item).endswith('.' + args.ext):
                    log.critical(item)
                    exp_set.append(os.path.join(args.loc_in, item))

                # sample optiion is valid only if  folder iin option is valid
    if (args.sample is not None) and (args.loc_in is not None):
        exp_set_app = copy.deepcopy(exp_set)
        for a in exp_set:
            if re.search(args.sample, a) is None:
                exp_set_app.remove(a)
        exp_set = exp_set_app

    if (exp_set == []) or (len(exp_set) == 1):
        exit(
            'ERROR input files not found or just one input file selected . check the folder or the extension given in input')

    min_RT = 100
    max_RT = 0
    for a in exp_set:
        log.critical('Reading file: %s ', a)
        exp_subset.append(a)
        data_moff = pd.read_csv(a, sep="\t", header=0)
        list_name = data_moff.columns.values.tolist()
        # get the lists of PS  defaultcolumns from properties file
        list_ps_def = ast.literal_eval(
            config.get('moFF', 'ps_default_export_v1'))
        # here it controls if the input file is a PS export; if yes it maps the input in right moFF name
        if moff.check_ps_input_data(list_name, list_ps_def) == 1:
            log.critical(
                'Detected input file from PeptideShaker export..: %s ', a)
            # map  the columns name according to moFF input requirements
            data_moff, list_name = moff.map_ps2moff(
                data_moff, 'col_must_have_mbr')
            log.critical(
                'Mapping columns names into the  the moFF requested column name..: %s ', a)
            # print data_moff.columns
        if moff.check_columns_name(list_name, ast.literal_eval(config.get('moFF', 'col_must_have_mbr')), log) == 1:
            exit('ERROR minimal field requested are missing or wrong')
        data_moff['matched'] = 0
        data_moff['mass'] = data_moff['mass'].map('{:.4f}'.format)

        data_moff['code_unique'] = data_moff['mod_peptide'].astype(str)    + '_' \
             + data_moff['prot'].astype(str) +'_' +'_'  \
              +  data_moff['charge'].astype(str)

        data_moff = data_moff.sort_values(by='rt')
        # data_moff['rt'] = data_moff['rt'] # / 60
        exp_t.append(data_moff)
        exp_out.append(data_moff)
        if data_moff['rt'].min() <= min_RT:
            min_RT = data_moff['rt'].min()
        if data_moff['rt'].max() >= max_RT:
            max_RT = data_moff['rt'].max()

    c_data = copy.deepcopy(exp_t)

    #intersect_share = reduce(np.intersect1d, ([x['code_unique'].unique() for x in exp_t]))

    intersect_share = reduce(np.union1d, ([x['code_unique'].unique() for x in exp_t]))

    #'K.VAVLGAAGGIGQALALLLK.T_16131126__3',
    #'K.AMLQDIATLTGGTVISEEIGMELEK.A_16131968__3', 'R.ELLSQYDFPGDDTPIVR.G_16131218__2', 'R.MGHIELASPTAHIWFLK.S_16131818__3',
    #'R.YTLAGTEVSALLGR.M_16131600__2', 'K.DKPEDAVLDVQGIATVTPAIVQACTQDK.Q_16131382__3'

    #intersect_share= np.array(['R.SALQYAASVAGLMITTECMVTDLPK.N_16131968__3','-.MRHPLVMGNWK.L_16131757__2'])

    log.critical('Read input --> done ')
    # parameter of the number of query
    # set a list of filed mandatory
    # ['matched','peptide','mass','mz','charge','prot','rt']
    n_replicates = len(exp_t)
    exp_set = exp_subset
    aa = range(0, n_replicates)
    #aa = range(6, n_replicates)
    out = list(itertools.product(aa, repeat=2))
    # just to save all the model
    # add matched columns
    list_name.append('matched')
    # final status -1 if one of the output is empty
    out_flag = 0
    # input of the methods
    diff_field = np.setdiff1d(exp_t[0].columns, [
        'matched', 'mod_peptide', 'peptide', 'mass', 'mz', 'charge', 'prot', 'rt'])

    log.info('Outlier Filtering is %s  ', 'active' if
        args.out_flag else 'not active')
    log.info('Number of replicates %i,', n_replicates)
    log.info('Pairwise model computation ----')

    if args.rt_feat_file is not None:
        log.critical(
            'Custom list of peptide used  provided by the user in %s', args.rt_feat_file)
        # log.info('Custom list of peptide used  provided by the user in %s', args.rt_feat_file)
        shared_pep_list = pd.read_csv(args.rt_feat_file, sep='\t')
        shared_pep_list['mass'] = shared_pep_list['mass'].map('{:.4f}'.format)
        shared_pep_list['code'] = shared_pep_list['peptide'].astype(
            str) + '_' + shared_pep_list['mass'].astype(str)
        list_shared_pep = shared_pep_list['code']
        log.info('Custom list of peptide contains  %i ',
                 list_shared_pep.shape[0])
    list_inter = 220
    interval, l_int = create_belief_RT_interval(max_RT + 1, min_RT - 1, list_inter)
    kf = KFold(n_splits=4)
    for jj in [3]:
        fold=0
        #200:1000
        for train_index, test_index in kf.split(intersect_share):
            if jj == 10:
                print ('stop')
            X_train, X_test = intersect_share[train_index], intersect_share[test_index]
            log.critical('# testing :%r',X_test.shape)
            for c in range(0,len(c_data)):
                exp_t[c] = c_data[c]
            #exp_t[jj] = c_data[jj]

            # exp_t= c_data
            exp_t[jj] = exp_t[jj][~exp_t[jj].code_unique.isin(X_test)]

            # list of the model saved
            model_save = []
            # list of the error in min/or sec
            model_err = []
            # list of the status of the model -1 means model not available for low points in the training set
            model_status = []
            c_rt = 0
            pre_pep_save = []

            #rank_df,rt_df = rank_matrix( intersect_share,exp_t )

            log.info('matching  in %s', exp_set[jj])
            result = itertools.filterfalse(lambda x: x[0] != jj or x[1] == jj, out)
            train_data = []
            for i in result:
                #if i[0] == jj and i[1] != jj:
                if args.rt_feat_file is not None:
                    # use of custom peptide
                    comA = exp_t[i[0]][exp_t[i[0]]['code_unique'].isin(list_shared_pep)][
                        ['code_unique', 'peptide', 'prot', 'rt']]
                    comB = exp_t[i[1]][exp_t[i[1]]['code_unique'].isin(list_shared_pep)][
                        ['code_unique', 'peptide', 'prot', 'rt']]
                    comA = comA.groupby('code_unique', as_index=False).mean()
                    comB = comB.groupby('code_unique', as_index=False).mean()
                    common = pd.merge(
                        comA, comB, on=['code_unique'], how='inner')
                else:
                    # use of shared peptdes.
                    log.info('  Matching  %s peptide in   searching in %s ',
                             exp_set[i[0]], exp_set[i[1]])
                    list_pep_repA = exp_t[i[0]]['code_unique'].unique()
                    list_pep_repB = exp_t[i[1]]['code_unique'].unique()
                    #log.info('Peptide unique (mass + sequence) %i , %i ',
                     #        list_pep_repA.shape[0],
                      #       list_pep_repB.shape[0])
                    #set_dif_s_in_1 = np.setdiff1d(list_pep_repB, list_pep_repA)
                    add_pep_frame = exp_t[i[1]][exp_t[i[1]]['code_unique'].isin(
                        X_test)].copy()
                    #-- prepare the testing set
                    add_pep_frame = add_pep_frame[[
                        'peptide', 'mod_peptide','code_unique', 'mass', 'mz', 'charge', 'prot', 'rt']]
                    # add_pep_frame['code_unique'] = '_'.join([add_pep_frame['peptide'], add_pep_frame['prot'], add_pep_frame['mass'].astype(str), add_pep_frame['charge'].astype(str)])
                    add_pep_frame['code_unique'] = add_pep_frame['mod_peptide'] + '_' + \
                                                   add_pep_frame['prot'].astype(str) + '_' + '_' + \
                                                   add_pep_frame['charge'].astype(str)
                    add_pep_frame = add_pep_frame.groupby('code_unique', as_index=False)[
                        'peptide', 'mod_peptide', 'mass', 'charge', 'mz', 'prot', 'rt'].aggregate(max)
                    add_pep_frame = add_pep_frame[[
                        'code_unique', 'peptide', 'mod_peptide', 'mass', 'mz', 'charge', 'prot', 'rt']]
                    list_name = add_pep_frame.columns.tolist()
                    list_name = [w.replace('rt',  str(c_rt) )
                                 for w in list_name]
                    add_pep_frame.columns = list_name
                    #if add_pep_frame.shape[0]==0:
                    log.critical('investigate fold %r target run %r input run% r',fold,i[0],i[1])
                    pre_pep_save.append(add_pep_frame)

                    #--------
                    pep_shared = np.intersect1d(list_pep_repA, list_pep_repB)
                    log.info(
                        '  Peptide (mass + sequence)  added size  %i ', add_pep_frame.shape[0])
                    log.info('  Peptide (mass + sequence) )shared  %i ',
                             pep_shared.shape[0])
                    comA = exp_t[i[0]][exp_t[i[0]]['code_unique'].isin(pep_shared)][
                        ['code_unique', 'peptide', 'prot', 'rt']]
                    comB = exp_t[i[1]][exp_t[i[1]]['code_unique'].isin(pep_shared)][
                        ['code_unique', 'peptide', 'prot', 'rt']]
                    # filtering using the variance added 17_08
                    flag_var_filt = False
                    if flag_var_filt:
                        dd = comA.groupby('code_unique', as_index=False)
                        top_res = dd.agg(['std', 'mean', 'count'])
                        # print np.nanpercentile(top_res['rt']['std'].values,[5,10,20,30,50,60,80,90,95,97,99,100])
                        th = np.nanpercentile(top_res['rt']['std'].values, 80)
                        log.critical('Distribution of rt std %r', np.nanpercentile(top_res['rt']['std'].values,[5,20,50,75,80,95]))
                        log.critical('Target data th std used %4.4f', th)
                        comA = comA[~ comA['code_unique'].isin(
                            top_res[top_res['rt']['std'] > th].index)]
                        # data B '
                        dd = comB.groupby('code_unique', as_index=False)

                        top_res = dd.agg(['std', 'mean', 'count'])
                        # print comB.shape
                        # print np.nanpercentile(top_res['rt']['std'].values,[5,10,20,30,50,60,80,90,95,97,99,100])
                        log.critical('Distribution of rt std %r', np.nanpercentile(top_res['rt']['std'].values,[5,20,50,75,80,95]))
                        th = np.nanpercentile(top_res['rt']['std'].values, 80)
                        log.critical('Input data th std used %4.4f', th)
                        comB = comB[~ comB['code_unique'].isin(
                            top_res[top_res['rt']['std'] > th].index)]

                    comA = comA.groupby('code_unique', as_index=False).mean()
                    comB = comB.groupby('code_unique', as_index=False).mean()
                    common = pd.merge(
                        comA, comB, on=['code_unique'], how='inner')

                #common[['rt_x', 'rt_y']].to_csv('D:\\workspace\\ionstar_dataset\\RT_benchmark2\\testdata_Ax_By_'+ str(i[0]) + '_' + str(i[1]) + '.txt', sep="\t")
                if common.shape[0] <= 10 and args.rt_feat_file is not None:
                    model_status.append(-1)
                    continue
                # filtering outlier option
                else:
                    if args.out_flag :
                        filt_x, filt_y, pos_out = MD_removeOutliers(common['rt_y'].values, common['rt_x'].values,
                                                                    args.w_filt)
                        data_B = filt_x
                        data_A = filt_y
                        data_B = np.reshape(data_B, [filt_x.shape[0], 1])
                        data_A = np.reshape(data_A, [filt_y.shape[0], 1])
                        log.info('Outlier founded %i  w.r.t %i',
                                 pos_out.shape[0], common['rt_y'].shape[0])
                    else:
                        data_B = common['rt_y'].values
                        data_A = common['rt_x'].values
                        data_B = np.reshape(data_B, [common.shape[0], 1])
                        data_A = np.reshape(data_A, [common.shape[0], 1])
                    '''
                    log.info(' Size trainig shared peptide , %i %i ',
                             data_A.shape[0], data_B.shape[0])
                    clf = linear_model.RidgeCV(alphas=np.power(
                        2, np.linspace(-30, 30)), scoring='neg_mean_absolute_error')
                    clf.fit(data_B, data_A)
                    clf_final = linear_model.Ridge(alpha=clf.alpha_)
                    clf_final.fit(data_B, data_A)
                    # save the model
                    model_save.append(clf_final)
                    model_err.append(mean_absolute_error(
                        data_A, clf_final.predict(data_B)))
                    log.info(' Mean absolute error training : %4.4f sec',
                             mean_absolute_error(data_A, clf_final.predict(data_B)))
                    model_status.append(1)
                    '''
                    #common['run_id'] = int(i[1])
                    list_name_1=  common.columns.tolist()
                    list_name_1 = [w.replace('rt_y', 'input_' + str(c_rt)) for w in  list_name_1]
                    list_name_12 = [w.replace('rt_x','target_' + str(c_rt)) for w in  list_name_1]
                    common.columns= list_name_12
                    train_data.append(common[['code_unique','target_' + str(c_rt),'input_' + str(c_rt)]])
                    log.critical('model %i ',c_rt)
                    c_rt += 1
                    # GP version
                    model_gp, predicted_train, error = train_gp(data_A, data_B,c= str(i[0])+'_'+str(i[1]))
                    #print i[1], comA.shape, error

                    model_err.append(error)
                    model_save.append(model_gp)
                    model_status.append(1)

            if np.where(np.array(model_status) == -1)[0].shape[0] >= (len(aa) / 2):
                log.error(
                    'MBR aborted :  mbr cannnot be run, not enough shared pepetide among the replicates ')
                exit('ERROR : mbr cannnot be run, not enough shared pepetide among the replicates')

            ##--------- mono model prep
            test_all = reduce(lambda left, right: pd.merge(left, right, on=['code_unique'], how='outer'), train_data)
            #test_all = test_all.fillna(0)
            #test_all.filter(regex='target_*')
            #test_all.filter(regex='input_*')
            #test_all.to_csv('D:\\workspace\\ionstar_dataset\\RT_benchmark2\\testdata_Mono_' + str(i[0]) + '_' + str(i[1]) + '.txt', sep="\t")

            Y_train= test_all.filter(regex='target_*').mean(axis=1).values
            X_train = test_all.filter(regex='input_*').mean(axis=1).values

            fullmodel_gp, fullpredicted_train, error = train_gp_mono(Y_train, X_train.reshape(-1, 1), c=str(i[0]) + '_' + str(i[1]))
            ##----- end mono model part

            log.info('Combination of the  model  --------')
            log.info('Weighted combination  %s : ', 'Weighted' if
            args.w_comb else 'Unweighted')
            if n_replicates == 2:
                test = pre_pep_save[0]
            else:
                test = reduce(
                    lambda left, right: pd.merge(left, right, on=[
                                                 'code_unique', 'peptide', 'mod_peptide', 'mass', 'mz', 'charge', 'prot'], how='outer'),
                    pre_pep_save)
            test = test.reindex(sorted(test.columns), axis=1)
            test = test.groupby('code_unique', as_index=False).aggregate(max)
            ## RT_benchmark_2
            #test = test[
            #    ['0', '1', '2', '3', '4', '5', '6',
            #     'code_unique', 'charge', 'mass', 'mod_peptide', 'mz', 'peptide', 'prot']]

            test = test[
                ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                 'code_unique', 'charge', 'mass', 'mod_peptide', 'mz', 'peptide', 'prot']]

            #test = test[['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17',
            #      '18','code_unique','charge','mass','mod_peptide','mz','peptide','prot']]

            if test.shape[0]==0:
                log.critical('no testing data ...')

            log.critical('size testing data : %r',test.shape)
            #test.drop('code_unique', axis=1, inplace=True)
            #test[['time_pred','uncertainty_win']] = test.iloc[:, 7: (7 + (n_replicates - 1))].apply(
            #    lambda x: combine_model_GP(x, model_save, model_err, args.w_comb,rank_df,rt_df),axis=1)


            #test[['time_pred', 'uncertainty_win']] = test.apply( lambda x:
            #                                    combine_model_GP_rankAggr(x, model_save, model_err, args.w_comb, rank_df, rt_df,exp_t,jj), axis=1)
            # rank bound output
            #test[['time_pred', 'uncertainty_win','rank_up_bound','rank_low_bound']] = test.apply(lambda x:
            #                                                    combine_model_GP(x, model_save, model_err, args.w_comb,
            #                                                                     rank_df, rt_df, exp_t, jj), axis=1)
            #Mono model -1
            #test_data = test.iloc[:,7:17].fillna(0).values
            ## aggregated input MonoModel_v2

            # rt_benchmark2
            #test_data =test.iloc[:, 0:6].mean(axis=1).values

            test_data =test.iloc[:, 0:10].mean(axis=1).values
            #test_data = test.iloc[:, 0:10].mean(axis=1).values
            mono_mean_pred, mono_uncert_win =  combine_model_GP_mono(test_data.reshape(-1,1),fullmodel_gp)

            test[['time_pred', 'uncertainty_win']] = test.apply(lambda x: combine_model_GP(x,
                                                                    model_save,model_err,args.w_comb) ,  axis=1)


            test[['time_ev','ev_upper','ev_lower','conflict_flag']]  =  test.apply(lambda x: combine_model_GP_consonant_bf
                        (x, model_save,log, interval,l_int /2,model_err) ,axis=1)

            #test['time_ev'] = np.nan
            #test['ev_upper'] = np.nan
            #test['ev_lower']  = np.nan
            #test['conflict_flag'] = np.nan

            test['mono_time_pred']= mono_mean_pred
            test['mono_uncertainty_win']=  mono_uncert_win



            #test[['time_pred','uncertainty_win']] = test.iloc[:, 7: (7 + (n_replicates - 1))].apply(
            #       lambda x: combine_model(x, model_save, model_err, args.w_comb),axis=1)
            test['matched'] = 1

            # still to check better
            if test[test['time_pred'] <= 0].shape[0] >= 1:
                log.info(' -- Predicted negative RT: those peptide will be deleted')
                test = test[test['time_pred'] > 0]

            list_name = test.columns.tolist()
            list_name = [w.replace('time_pred', 'rt_pred') for w in list_name]
            test.columns = list_name

            # test = test[['peptide','mod_peptide', 'mass', 'mz', 'charge',
            # 'prot', 'rt', 'matched','uncertainty_win']]
            # 'rank_up_bound','rank_low_bound'
            test = test[['code_unique','peptide', 'mod_peptide', 'mass',
                         'mz', 'charge', 'prot', 'rt_pred', 'uncertainty_win','matched','mono_rt_pred','mono_uncertainty_win','conflict_flag','ev_lower','ev_upper','time_ev']]
            test = test.merge(c_data[jj][c_data[jj].code_unique.isin(X_test)][['prot', 'rt', 'code_unique']],
                       on=['code_unique', 'prot'], how='inner')
            # just put nan with the missing values
            #for field in diff_field.tolist():
            #    test[field] = np.nan
            print (fold,jj)
            if fold == 0:
                output = test
            else:
                output = pd.concat([output,test],axis=0,join='outer')
            fold +=1
            #log.info('Before adding %s contains %i ',
            #         exp_set[jj], exp_t[jj].shape[0])
            #exp_out[jj] = pd.concat(
            #    [exp_t[jj], test], join='outer', axis=0, sort=False)
            #log.info('After MBR %s contains:  %i  peptides',
             #        exp_set[jj], exp_out[jj].shape[0])
            #log.critical('matched features   %i  MS2 features  %i ', exp_out[jj][exp_out[jj]['matched'] == 1].shape[0],
             #            exp_out[jj][exp_out[jj]['matched'] == 0].shape[0])
        # save result
        output.to_csv(
                path_or_buf=os.path.join(output_dir, os.path.split(exp_set[jj])[1].split('.')[0] + '_match.txt'), sep='\t',
                index=False)
        #exp_out[jj].to_csv(
        #        path_or_buf=os.path.join(output_dir, os.path.split(exp_set[jj])[1].split('.')[0] + '_match.txt'), sep='\t',
        #        index=False)
        exp_out_name.append(os.path.join(output_dir, os.path.split(
                exp_set[jj])[1].split('.')[0] + '_match.txt'))
        if exp_out[jj].shape[0] > 0:
            out_flag = 1 * out_flag
        else:
            out_flag = -1 * out_flag


    w_mbr.close()
    log.removeHandler(w_mbr)
    return out_flag, exp_out_name

