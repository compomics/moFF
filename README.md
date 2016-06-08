# moFF #

 * [Introduction](#introduction)
 * [Minimum Requirements](#minimum-requirements)
 * [Input Data](#input-data)
 * [Sample Data](#sample-data)
 * [Match between runs](#match-between-runs)
 * [Apex Intensity](#apex-intensity)
 * [Entire workflow](#entire-workflow)
 * [Output Data](#output-data)

---

## Introduction ##

moFF is a cross-platform  tool that extracts apex MS1 intensity  starting from a list of MS2 idenfied peptide.
moFF works directly  on Thermo Raw file using a Go library that is able to read the data from Thermos raw file without any kind of conversion in other formats.

moFF is composed by two stand alone modules :
- *moff_mbr.py* :  matching between run
- *moff.py* :  apex intensity

To run  the entire workflow (mbr and apex ) you should  use  *moff_all.py*

[Top of page](#moff)

----

## Minimum Requirements ##

Required python libraries :
- Python 2.7
- pandas  > 0.17.
- numpy > 1.10.0
- argparse > 1.2.1 
- scikit-learn > 0.17

moFF uses *txic* to extract the XiC data from the raw files, so  *txic*  must be located in the same folder where you have all moFF scripts.

The txic program is compatibale with  the raw file of all the Orbitrap and triple quadrupole Thermo machines. 
For the moment it does not work with the Thermo Fusion machine.

[Top of page](#moff)

---


##Input data 


The tab-delimited file that contains the MS2 identified peptides must contain the following information for each peptides:
  - 'peptide' : sequence of the peptide
  - 'prot': protein ID 
  - 'rt': retention time of peptide   (The retention time must be specified in second)
  - 'mz' : mass over charge
  - 'mass' : mass of the peptide
  - 'charge' : charge of the ionized peptide

see the sample input files in the folder *f1_folder* for more information .

[Top of page](#moff)

---

## Sample data  ##

The  *f1_folder* contains as  the input files the MS2 identified  peptides founded by MASCOT  for 3 runs of  the CPTAC study 6 (Paulovich, MCP Proteomics, 2010). 
You can download the relative [raw files]( https://goo.gl/ukbpCI) in order to run  moFF on the sample data.

---

## Match between runs ##

use :  `python moff_mbr.py -h`
```
	--inputF              specify the folder of the input MS2 peptide list files
  	--sample	      specify witch replicated to use for mbr reg_exp are valid
  	--ext                 specify the file extention of the input file
  	--log_file_name       a label name to use for the log file
  	--filt_width          width value for  the outlier  filtering 
  	--out_filt            filtering (on/off) of the outlier in the training set
  	--weight_comb         combination weighting : 0 for no weight 1 for a weighted schema
```

`python moff_mbr.py --inputF f1_folder/` 

It runs the mbr modules and save the output files in a subfolder  called 'mbr_output' inside the folder given in input.
The match-between-runs module will take all the .txt files in your input folder as replicates. (to select specific files or different extension see below))
In *f1_folder/mbr_output* you will find the same number of the input files, but they will have a new field called 'matched' that specifies which peptides are matched  (1) or the not (0)

if your input files inside your working fodler have another exstension like (.list, etc) you can use :

use : `python --inputF f1_folder/ --ext list ` ( Do not specify '.list' but only 'list')

if you need to select specific input files from your working folder  ( choose  ) , you can use an regular expression as:

use : `python --inputF f1_folder/  --sample *_6A ` (you can also use --ext option if you need)

The match between runs will produce a log file  with all the details and it will be  saved inside the folder given in input

[Top of page](#moff)

---

## Apex intensity ##

use  `python moff.py -h`
````
  --input NAME        specify the input file with the of MS2 peptides
  --tol               specify the tollerance parameter in ppm
  --rt_w              specify rt windows for xic (minute). Default value is  3  min
  --rt_p     	      specify the time windows used to get  the apex  for the ms2 peptide/feature  ( minute). Default value is 0.2
  --rt_p_match 	      specify the time windows used to get  the apex  for machted features ( minute). Default value is 0.4
  --raw_repo          specify the raw file folder
  --output_folder     specify the folder output
```
`python moff.mbr --input f1_folder/20080311_CPTAC6_07_6A005.txt  --raw_rep f1_folder/ --tol 1O ` 
It will save the results in the folder inside the f1_folder

use `python moff.mbr --input f1_folder/20080311_CPTAC6_07_6A005.txt  --raw_rep f1_folder/ --tol 1O --output_folder output_moff`
It will save the results in the folder output_moff

WARNING : the raw file names  MUST be the same of the input file otherwise the script give you an error !

[Top of page](#moff)

---

## Entire workflow ##

use `python moff_all.py -h`
```
	--inputF              specify the folder of the input MS2 peptide list files
  	--sample	      specify witch replicated to use for mbr reg_exp are valid
  	--ext                 specify the file extention of the input file
  	--log_file_name       a label name to use for the log file
  	--filt_width          width value for  the outlier  filtering 
  	--out_filt            filtering (on/off) of the outlier in the training set
  	--weight_comb         combination weighting : 0 for no weight 1 for a weighted schema
  	--tol         	      specify the tollerance parameter in ppm
  	--rt_w                specify rt windows for xic (minute). Default value is  3  min
  	--rt_p     	      specify the time windows used to get  the apex  for the ms2 peptide/feature  ( minute). Default value is 0.2
  	--rt_p_match 	      specify the time windows used to get  the apex  for machted features ( minute). Default value is 0.4
  	--raw_repo            specify the raw file folder
  	--output_folder       specify the folder output
```
`python moff_all.py --inputF  f1_folder/   --raw_repo f1_folder/ --output_folder output_moff`

The options are the same of the two modules, the the output mbr files are stores in the folder f1_folder/mbr_output  and the result of the apex module are stored in output_moff. Also the log files are stored in the respective folders

[Top of page](#moff)

---
## Output data

The output consists of : 

- Tab delimited file (with the same name of the input raw file) that contains the apex intensity values and some other information (a)
- Log file specific to the apex module (b) or the MBR module (c)

(a) Description of the fields added by moFF in the output file:

Parameter | Meaning
--- | -------------- | 
*rt_peak* | retention time (in seconds) of the discovered apex peak
*SNR*     | signal-to-noise  ratio of the peak intensity.
*log_L_R*'| peak shape. 0 indicates that the peak is centered. Positive or negative values are an indicator for respectively right or left skewness 
*intensity* |  MS1 intensity
*log_int* | log 2 transformed MS1 intensity 
*lwhm* | first rt value where the intensity is at least the 50% of the apex peak intensity on the left side
*rwhm* | first rt value where the intensity is at least the 50% of the apex peak intensity on the right side
*5p_noise* | 5th percentile of the intensity values contained in the XiC. This value is used for the *SNR* computation
*10p_noise* |  10th percentile of the intensity values contained in the XiC.
*code_unique* | this field is concatenation of the peptide sequence and mass values. It is used by moFF during the match-between-runs.
*matched* | this value indicated if the featured has been added by the match-between-run (1) or is a ms2 identified features (0) 

(b) A log file is also provided containing the process output. 

(c) A log file where all the information about all the trained linear model are displayed.

NOTE : The log files and the output files are in the output folder specified by the user. 

[Go to top of page](#moff-gui)
