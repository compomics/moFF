# moFF #
## A modest Feature Finder (but still robust) to extract apex MS1 intensity directly from Thermo raw file ##

[Introduction](#introduction)
[Requirement](#requirement)
[Sample data](#sample-data)
[Matching Between Runs](#matching-between-runs)
[Apex intensity](#apex-intensity)
[Entire workflow](#entire-workflow)


---

## Introduction ##

moFF is a python tool that quantifies MS1 intensity peak starting from a list of MS2 idenfied peptide 
moFF works directly and only on Thermo Raw file thanks to a Go library that is able to read the data from Thersmos raw file without any kind of conversion in other formats.

moFF is composed by two stand alone modules :
- *moff_mbr.py* :  matching between run
- *moff.py* :  apex intensity

To run  the entire workflow (mbr and apex ) you should  use  *moff_all.py*

[Top of page](#moff)

----

## Requirement ##

Required python libraries :
- Python 2.7
- pandas  > 0.17.
- numpy > 1.10.0
- argparse > 1.2.1 
- scikit-learn > 0.17

moFF uses *txic* to extract the XiC data from the raw files, so  *txic*  must be located in the same folder where you have all moFF scripts.

The txic program is compatibale with  the raw file of all the Orbitrap and triple quadrupole Thermo machines. 
For the moment it does not work with the Thermo Fusion machine.

The input files that contain the list of the MS2 identified peptides (you can use any search engines) must contains the information showed in *moFF_setting.property* for each peptide. The minimun specificic requirements of the input files are:
- tab delimited file
- the header of the input file should contain the following the fields  and columnns names :  
  - 'peptide' : sequence of the peptide
  - 'prot': protein ID 
  - 'rt': retention time of peptide   ( The retention time must be specified in second )
  - 'mz' : mass over charge
  - 'mass' : mass of the peptide
  - 'charge' : charge of the ionized peptide

see the sample input files in the folder *f1_folder* for more information .

[Top of page](#moff)

---

## Sample data  ##

In the folder *f1_folder* you have three input files, that contain the MS2 identified  peptides founded by MASCOT of three runs (three tecnical replicates) from  the CPTAC study 6. 
You can download the relative [raw files]( https://goo.gl/ukbpCI), in order to run the next examples.

---

## Matching Between Runs ##

use :  `python moff_mbr.py -h`
```
  --inputF LOC_IN             specify the folder of the input MS2 peptide files [REQUIRED]
  --sample SAMPLE            specify which replicate files are used fot mbr [regular expr. are valid]
  --ext EXT                  specify the exstension of the input file (txt as default value)
  --log_file_name LOG_LABEL  a label name for the log file (moFF_mbr.log as default log file name)
  --filt_width W_FILT        iwidth value of the filter (k * mean(Dist_Malahobis) , k = 2 as default)
  --out_filt OUT_FLAG        filter outlier in each rt time allignment (active as default)
  --weight_comb W_COMB       weights for model combination combination : 0 for no weight (default) 1 weighted devised by model errors.
```

`python moff_mbr.py --inputF f1_folder/` 

It runs the mbr modules and save the output files in a subfolder  called 'mbr_output' inside the folder given in input.
The mbr module will take all the .txt files in your input folder as replicates. (to select specific files or different extension see below))
In *f1_folder/mbr_output* you will find the same number of the input files, but they will have a new field called 'matched' that specifies which peptides are matched  (1) or the not (0)
The rt field of the matched peptide contains the predicted rt retentioins time.

if your input files inside your working fodler  have another exstension like (.list, etc) you can use :

use : `python --inputF f1_folder/ --ext list ` ( Do not specify '.list' but only 'list')

if you need to select specific input files from your working folder  ( choose  ) , you can use an regular expression as:

use : `python --inputF f1_folder/  --sample *_6A ` (you can also use --ext option if you need)

the mbr will output a log file (moFF_mbr.log as default log file name) with all the details and it is saved inside the  folder given in inout

[Top of page](#moff)

---

## Apex Intensity ##

use  `python moff.py -h`
````
  --input NAME                        specify the input file with the of MS2 peptides
  --tol TOLL                          specify the tollerance parameter in ppm
  --rt_w RT_WINDOW                    specify rt window for xic (minute). Default value is 3 min
  --rt_p RT_P_WINDOW                  specify the time windows for the peak ( minute). Default value is 0.1
  --rt_p_match RT_P_WINDOW_MATCH      specify the time windows for the matched peptide peak ( minute). Default value is 0.4
  --raw_repo RAW                      specify the raw file repository
  --output_folder LOC_OUT             specify the folder output
```
`python moff.mbr --input f1_folder/20080311_CPTAC6_07_6A005.txt  --raw_rep f1_folder/ --tol 1O ` 
 
It run the apex module on the input file , extraxing the apex intesity from the respective raw files in folder specified.
In the output files moFF just adds the following fields to your origin input file:
- "intensity" intensity, taking the highest peak in the XIC
- "rt_peak" rt of the highest peak
- "lwhm" left width half maximun of the signal in seconds
- "rwhm" right width half maximun of the signal in seconds
- "SNR" signal-to-noise
- "log_L_R" log ratio of lwhm over rwhm (peak shape )
- "log_int" log 2 of the intesity 

It generates a .log file (with same name of input file) that contains  detailesd information for each peak retrieved.
This module determines automaticaly if the input file contains matched peptides or not.

WARNING : the raw file names  MUST be the same of the input file otherwise the script give you an error !

use `python moff.mbr --input f1_folder/20080311_CPTAC6_07_6A005.txt  --raw_rep f1_folder/ --tol 1O --output_folder output_moff`
It will save the results in the folder output_moff

[Top of page](#moff)

---

## Entire workflow ##

use `python moff_all.py -h`
```
	--inputF LOC_IN       specify the folder of the input MS2 peptide list files
  	--sample SAMPLE       specify witch replicated use for mbr reg_exp are valid
  	--ext EXT             specify the file extentention of the input like
  	--log_file_name LOG_LABEL a label name to use for the log file
  	--filt_width W_FILT   width value of the filter k * mean(Dist_Malahobis)
  	--out_filt OUT_FLAG   filter outlier in each rt time allignment
  	--weight_comb W_COMB  weights for model combination combination : 0 for no weight 1 weighted devised by trein err of the model.
  	--tol TOLL            specify the tollerance parameter in ppm
  	--rt_w RT_WINDOW      specify rt window for xic (minute). Default value is  3  min
  	--rt_p RT_P_WINDOW    specify the time windows for the peak ( minute). Default value is 0.1
  	--rt_p_match RT_P_WINDOW_MATCH	specify the time windows for the matched peptide peak ( minute). Default value is 0.4
  	--raw_repo RAW        	specify the raw file repository
  	--output_folder LOC_OUT		specify the folder output
```
`python moff_all.py --inputF  f1_folder/   --raw_repo f1_folder/ --output_folder output_moff`

The options are the same of the two modules, the the output mbr files are stores in the folder f1_folder/mbr_output  and the result of the apex module are stored in output_moff. Also the log files are stored in the respective folders

[Top of page](#moff)

---
