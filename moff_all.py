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

import moff_mbr
#import moff
parser = argparse.ArgumentParser(description='moFF match between run input parameter')

parser.add_argument('--inputF', dest='loc_in', action='store', help='specify the folder of the input MS2 peptide list files ', required=False)

parser.add_argument('--sample', dest='sample', action='store', help='specify witch replicated to use for mbr reg_exp are valid ', required=False)

parser.add_argument('--ext', dest='ext', action='store', default='.txt', help='specify the file extentention of the input like ', required=False)

parser.add_argument('--log_file_name', dest='log_label', action='store', help='a label name to use for the log file', required=False)

parser.add_argument ('--filt_width', dest='w_filt', action='store',default=2,help='width value of the filter  k * mean(Dist_Malahobis)', required=False  )

parser.add_argument('--out_filt', dest='out_flag', action='store',default=1,help='filter outlier in each rt time allignment',  required=False )

parser.add_argument('--weight_comb', dest='w_comb', action='store',default=0,help='weights for model combination combination : 0 for no weight  1 weighted devised by trein err of the model.', required=False )


args = parser.parse_args()

moff_mbr.run_mbr(args)

 
