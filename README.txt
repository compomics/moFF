moFF
====

 A modest Feature Finder (but still robust)  to extract feature in MS1 Data

Required library :

Python 2.7
pandas  (latest version)
numpy (latest version)
argparse

moFF uses txic to read the .raw file, the execute txic must be located in the same folder where you have moFF.py file.



use 
python moff.py -h

to  see all possible parameter.

optional arguments
  --input NAME           specify input list of MS2 peptides
  --tol TOLL             specify the tollerance parameter in ppm
  --rt_w RT_WINDOW       specify rt window for xic (minute)
  --rt_p RT_P_WINDOW     specify the time windows for the peak ( minute). Default value is 0.4
  --raw RAW              specify the raw file
  --foder_output LOC_OUT specify the folder output
  
The imput file must a be a tab delimited file with an header that contains at least two fields  called  as folowing
"precursor"   that contains the precursor mass of the MS2 identified peptides
"rt" that contains the rt time  in seconds of the MS2 identified peptides.

Please notice that the name of fields is case sensitive

If the file contains also other information for each MS2 peptides (peptide sequence, protein, etc..)it does not mater  for moFF,
in the output file moFF just add the following fields tothe origin input file:

"intensity" intensity
"rt_peak" rt of the feature 
"lwhm" left width half maximun of the signal in seconds
"rwhm" right width half maximun of the signal in seconds
"SNR" signal-to-noise
"log_L_R" log ratio of lwhm over rwhm (peak shape )


