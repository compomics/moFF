moFF
====

 A modest Feature Finder (but still robust) to extract apx MS1 itensity directly from Thermo raw

moFF is written in python and it is based on a Go library that is able to read raw file from Thermo machine

Required library :

Python 2.7
pandas  > 0.14.1
numpy > 1.9.0
argparse > 1.2.1 i
glob
sci-learn >
logging
itertools
os
argparse
re

moFF is composed by two  stand alone modules : 
moff_mbr.py to run the matching between run 
moff.py for the apex intensity

To run both modules moff_all.py is able to launch the mbr and the apex modules for the total  workflow



moFF uses txic to read the .raw file, the execute txic must be located in the same folder where you have moFF.py file.



use  python moff.py -h

to  see all possible input parameters.

  --input NAME           specify an input tab-separated file that contains a list of MS2 peptides
  --tol TOLL             specify the tollerance parameter in ppm
  --rt_w RT_WINDOW       specify  the half rt window for xic (minute) 
  --rt_p RT_P_WINDOW     specify the time windows for the peak ( minute). Default value is 0.4 (optional)
  --raw RAW              specify the raw file
  --foder_output LOC_OUT specify the folder output. As default the output file is written in the same folder of the script (optional)
  
example 

For each MS2 identified peptide (x) , characterized by  RT and precursor mass values (RT_x, mass_x) , it extracts a XIC for the th massses mass_x +/- toll  considering a time window from  RT - rt_w  to RT + rt_w 

The input file must a be a tab-delimited file with an header that contains at least two fields named  as folowing
"precursor"   that contains the precursor mass of the MS2 identified peptides
"rt" that contains the rt time  in seconds of the MS2 identified peptides.

Please notice that the name of the fields is case sensitive

If the input file contains also other information fields such as peptide sequence and protein it does not mater for the script.
As result, moFF just add the following fields to the your origin input file:

"intensity" intensity, taking the highest peak in the XIC
"rt_peak" rt of the highest peak
"lwhm" left width half maximun of the signal in seconds
"rwhm" right width half maximun of the signal in seconds
"SNR" signal-to-noise
"log_L_R" log ratio of lwhm over rwhm (peak shape )


