# moFF #

 * [Introduction](#introduction)
 * [Minimum Requirements](#minimum-requirements)
 * [Input Data](#input-data)
 * [Sample Data](#sample-data)
 * [Match between runs](#match-between-runs)
 * [Apex Intensity](#apex-intensity)
 * [Entire workflow](#entire-workflow)
 * [Post Translation Modification file](#post-translation-modification-file)
 * [Output Data](#output-data)

---

## Introduction ##

moFF is an OS independent tool designed to extract apex MS1 intensity using a set of identified MS2 peptides. It currently uses a Go library to directly extract data from Thermo Raw spectrum files, eliminating the need for conversions from other formats. Moreover, moFF also allows to work directly with mzML files.

moFF is built up from two standalone modules :
- *moff_mbr.py* :  match between run (mbr)
- *moff.py*: apex intensity

NOTE : Please use *moff_all.py* script to run the entire pipeline with both MBR and apex strategies.

The version presented here is a commandline tool that can easily be adapted to a cluster environment. A graphical user interface can be found [here](https://github.com/compomics/moff-gui). The latter is designed to be able to use [PeptideShaker](https://github.com/compomics/peptide-shaker) results as an input format. Please refer to the [moff-GUI](https://github.com/compomics/moff-gui) manual for more information on how to do this.


[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/moff/README.html) 

moFF is also available on bioconda. To install with conda, use the following command: 
```
conda install -c bioconda moff
```
This automatically installs all dependencies. Note that bioconda only supports 64-bit macOS and Linux. 

[Top of page](#moff)

----

## moFF Publication:
  * [Argentini et al. Nature Methods. 2016 12(13):964–966](http://www.nature.com/nmeth/journal/v13/n12/full/nmeth.4075.html).
  * If you use moFF as part of a publication, please include this reference.

---

## Minimum Requirements ##

Required java version :
- Java Runtime Environment (JRE) 8

Required python libraries :
- Python 2.7
- pandas  > 0.20.
- numpy > 1.10.0
- argparse > 1.2.1 
- scikit-learn > 0.18
- pymzML > 0.7.7


Required linux library:
- Mono version 4.2.1

Required windows library:
- .NET Framework 4.6.2


Optional requirements :
-when using PeptideShaker results as a source, a PeptideShaker installation (<http://compomics.github.io/projects/peptide-shaker.html>) needs to be availabe.
 

During processing, moFF makes use of a third party algorithm (txic_json.exe) which allows for the parsing of the Thermo RAW data.


[Top of page](#moff)

---


## Input Data ##

moFF requires two types of input for the quantification procedure :
 - Thermo RAW file or mzML file
 - MS2 identified peptide information

The MS2 identified peptides can be presented as a tab-delimited file containing mimimal (mandatory) annotation for each peptide (a)

(a) The tab-delimited file must contain the following information for all the peptides:
  - 'peptide' : peptide-spectrum-match  sequence
  - 'prot' : protein ID 
  - 'mod_peptide' :  peptide-spectrum-match  sequence that contains also possible modification (i.e `NH2-M<Mox>LTKFESK-COOH` )
  - 'rt': peptide-spectrum-match retention time  (i.e the retention time contained in the mgf file; The retention time must be specified in second)
  - 'mz' : mass over charge
  - 'mass' : mass of the peptide
  - 'charge' : charge of the ionized peptide
 
NOTE 1 : In case the tab-delimited file provided by the user contains fields that are not mentioned here (i.e petides length, search engines score) the algorithm will retain these in the final output. The peptide-spectrum-match sequence with its modications  and the protein id  and  informations are used only in the match-between-run module.

NOTE 2 : Users can also provide the default PSM export provided by PeptideShaker as source material for moFF.


[Top of page](#moff)

---

## Sample data  ##

The  *sample_folder* contains a resultset for 3 runs of the CPTAC study 6 (Paulovich, MCP Proteomics, 2010). These MS2  peptides are identified by  X!Tandem and MSGF+ using SearchGUI and then processed by PeptidesShaker . The [raw files]( https://goo.gl/ukbpCI) for this study are required to apply moFF to the sample data.

---

## Match between runs ##

use :  `python moff_all.py -mbr only `
```
  --loc_in                      the folder where the input files are located
  --sample                      reg exp to filter the input file names (only with --loc_in input option-
  --ext                         file extention of the input file. Default .txt)
  --log_label                   filename for the mbr log file. Default moFF_mbr
  --w_filt                      width value for outlier filtering. Default 3
  --out_flag                    filtering (on/off) of the outlier in the training set. Default 1
  --w_comb                      combination weighting : 0 for no weight 1 for a weighted schema. Default 1
```

`python moff_mbr.py --loc_in sample_folder/ --mbr only `

This command runs the MBR modules. The output will be stored in a subfolder ('mbr_output') inside the specified input folder.
The MBR module will consider all the .txt files present in the specified input folder as replicates (to select specific files or different extension, please refer to the example below).
The files in *sample_folder/mbr_output* will be identical to the input files, but they will have an additional field ('matched') that specifies which peptides have match (1) or not (0). The MBR algorithm also produces a log file in the provided input directory.


### Customizing Match between runs ###

In case of a different extension (.list, etc), please use :

`python moff_mbr.py  --loc_in  sample_folder/ --ext list ` (Provide the extension without the period ('.'))

In case of using only specific input files within the provided directory, please use a regular expression:

`python moff_mbr.py --loc_in sample_folder/  --sample *_6A` (This can be combined with the aforementioned syntax)

You can set all the parameters values in a file and load them using  `--config_file`. For an example see `example_parameter_file.ini`



[Top of page](#moff)

---

## Apex intensity ##

use `python moff_all.py -mbr off `
```
  --loc_in                      the folder containing all input files
  --raw_repo                    the folder containing all the raw files
  --tsv_list                    the input file with for MS2 peptides
  --raw_list                    pecify directly the  raw file
  --toll                        mass tollerance (ppm)
  --xic_length                  rt windows for xic (minutes). Default value is 3  min
  --rt_peak_win                 time windows used to get the apex for the ms2 peptide/feature  (minutes). Default value is 1
  --rt_peak_win_match           time windows used to get the apex for machted features (minutes). Default value is 1.2
  --peptide_summary             flag that allows have as output the peptided summary intensity file. Default is disable(0)
  --tag_pepsum                  tag string that will be part of the  peptided summary intensity file name. Default is moFF_run
  --loc_out                     output folder
  --tag_pepsum                  a tag that is used in the peptide summary file name

  --match_filter                filtering on the matched peak . default 0
  --ptm_file                    modification json ptm file. Default file ptm_setting.json
  --quantile_thr_filtering      quantile value used to computed the filtering threshold for the matched peak . Default is 0.75
  --sample_size                 percentage of MS2 identified peptides used to estimated the threshold
```

You can run the apex module in two ways:

`python moff_all.py  --mbr off --tsv_list sample_folder/20080311_CPTAC6_07_6A005.txt  --raw_list sample_folder/20080311_CPTAC6_07_6A005.RAW --toll 1O --loc_out output_moff --peptide_summary 1 `
in this case you specify more than a file separated by a blanck space

In case you want to run the apex module  on all the files in a folder (all so the raw files shold located in a foder)

`python moff_all.py  --mbr on  --loc_in sample_folder/sample_data/  --raw_repo sample_folder/sample_data/your_raw_folder   --toll 1O --loc_out output_moff --peptide_summary 1 `

You can activate the filtering of the matching peptides setting `--match_filter 1`. In order to do the filtering:
- `--ptm_file` MUST be specified and input files MUST contain a matched field.

This option is usefull in the case you have run the mbr module alone and later you want to run the apex module separately.

WARNING :  in case of  --loc_in  and  --raw_repo  raw file names MUST be the same of the input file otherwise the script gives you an error !

WARNING 1  :  you can not mixed the two input ways ( --loc_in / --raw_repo and --tsv_list / --raw_list  ) otherwise the script gives you an error !

WARNING 2: mzML raw file MUST be only specified using `--tsv_list | --raw_list`. The `--raw_repo` option is not available for mzML files.

NOTE: all the parameters related to the the time windows (xic_lentgh,rt_peak_win, rt_peak_win_match) are basicaly the half of the entire time windows where the apex peak is searched or the XiC is retrieved. For a correct rt windows, we suggest to set the **rt_peak_win** value equal or slighly greater to the __dynamic exclusion duration set in your machine.__
We suggest also to set the rt_peak_win_match  always slightly bigger than tha values used for rt_peak_win


[Top of page](#moff)

---


## Entire workflow ##

use `python moff_all.py -mbr on`
```
  --config_file                 specify a moFF parameter file
  --loc_in                      the folder containing all input files
  --raw_repo                    the folder containing all the raw files
  --tsv_list                    the input file with for MS2 peptides
  --raw_list                    pecify directly the  raw file

  --sample                      reg exp to filter the input file names (only with --loc_in input option-
  --ext                         file extention of the input file. Default .txt)
  --log_label                   filename for the mbr log file. Default moFF_mbr
  --w_filt                      width value for outlier filtering. Default 3
  --out_flag                    filtering (on/off) of the outlier in the training set. Default 1
  --w_comb                      combination weighting : 0 for no weight 1 for a weighted schema. Default 1

  --toll                        mass tollerance (ppm)
  --xic_length                  rt windows for xic (minutes). Default value is 3  min
  --rt_peak_win                 time windows used to get the apex for the ms2 peptide/feature  (minutes). Default value is 1
  --rt_peak_win_match           time windows used to get the apex for machted features (minutes). Default value is 1.2
  --peptide_summary             flag that allows have as output the peptided summary intensity file. Default is disable (0)
  --tag_pepsum                  tag string that will be part of the  peptided summary intensity file name. Default value is moFF_run
  --loc_out                     output folder  default is the input folder, raw_repo)
   --tag_pepsum                  a tag that is used in the peptide summary file name

  --match_filter                filtering on the matched peak . default 0
  --ptm_file                    modification json ptm file. Default file ptm_setting.json
  --quantile_thr_filtering      quantile value used to computed the filtering threshold for the matched peak . Default is 0.75
  --sample_size                 percentage of MS2 identified peptides used to estimated the threshold
```

Like for the apex module, you input  you input data specifing the folder :

`python moff_all.py --mbr all  --loc_in  sample_folder/   --raw_repo sample_folder/ --toll 10  --loc_out output_moff --peptide_summary 1`

OR, specifing a list of input and raw files using:

`python moff_all.py  --mbr all --tsv_list  sample_folder/input_file1.txt sample_folder/input_file2.txt  --raw_list sample_folder/input_file1.raw sample_folder/input_file2.raw --toll 10 --loc_out output_moff --peptide_summary 1 `

The options are identical for both apex and MBR modules.The output for the latter (MBR) is stored in the folder sample_folder/mbr_output, while the former (apex) generates files in the specified output_moff folder.Log files for both algorithms are generated in the respective folders.

In case you activate the filtering of the mached peptides  you have to specify with `--ptm_file` a valid json file that describes the modificatiuon used in your experiment. See section

You can set all the parameters values in a file and load them using `--config_file`. For an example see `example_parameter_file.ini`

WARNING: Using `--tsv_list | --raw_list`  you can not filterted the input file using `--sample --ext` like in the case with `--loc_in | --raw_repo`

WARNING: **mzML raw file  MUST be specified  using `--tsv_list | --raw_list`. The `--raw_repo` option is not available for mzML files.

NOTE: The consideration of retention time window parameters (xic_length,rt_peak_win,rt_peak_win_match) mentioned for apex module are stil valid also for the entire workflow


[Top of page](#moff)


---
## Post Translation Modification file ##

The Post Translation Modificatio must be indicated in json file with the following structure :
```
{
"tagModification": {"deltaChem":[H atom, C atom, N atom ,O atom],"desc":"name unimod : unimod_id"},
}
```

- `"tagModification"` : the tag used in modified sequence for the modification
- `"deltaChem":[H atom, C atom, N atom ,O atom]` : the delta of chemical composition if the modification. The order of the elements is fixed, so pay attention when you add your modification
- `desc` : name of the modification and its unimod id.

For example a ptm file (ptm_setting_ps.json) with Carboxyamidomethylation of Cysteine and Oxidation for PeptideShaker output looks like:
```
{
"<cmm>": {"deltaChem":[3,2,1,1],"desc":"Carboxyamidomethylation C unimod:4"},
"<ox>": {"deltaChem":[0,0,0,1],"desc":"oxidation oxidation unimod:35" }
}
```


WARNING: The handling of the modification is still a working part; so maybe in the future could be changed.



[Top of page](#moff)

---
## Output data ##

The output consists of : 

- a tab delimited file (with the same name of the input raw file) containing the apex intensity values and additional information (a)
- a log file specific to the apex module (b) or the MBR module (c)
- peptide summary intensity file (when peptide summary option is enabled) (d) 

(a) Description of the fields added by moFF in the output file:

Parameter | Meaning
--- | -------------- | 
*rt_peak* | retention time (in seconds) for the discovered apex peak
*SNR*     | signal-to-noise ratio of the peak intensity.
*log_L_R*'| peak shape. 0 indicates that the peak is centered. Positive or negative values are an indicator for respectively right or left skewness 
*intensity* |  MS1 intensity
*log_int* | log2 transformed MS1 intensity 
*lwhm* | first rt value where the intensity is at least the 50% of the apex peak intensity on the left side
*rwhm* | first rt value where the intensity is at least the 50% of the apex peak intensity on the right side
*5p_noise* | 5th percentile of the intensity values contained in the XiC. This value is used for the *SNR* computation
*10p_noise* |  10th percentile of the intensity values contained in the XiC.
*code_unique* | this field is concatenation of the peptide sequence and mass values. It is used by moFF during the match-between-runs.
*matched* | this value indicated if the featured has been added by the match-between-run (1) or is a ms2 identified features (0) 

(b) A log file is also provided containing the process output. 

(c) A log file where all the information about all the trained linear model are displayed.

(d) The peptide summary intensity is a tab delimited file where for each peptide sequence MS1 intensities are summed for all the occurences in each run (aggregated by charge states and modification).

In case you run the entire workflow on an a settings that contains N runs, the size of the file (rows and columns) will be **M x (N+2)**, where M is number of peptides (across all the runs) and N are summed intensity columns plus the peptide sequence and the protein ids. In case of running only the apex module, the size of the file  will be on M x 3 (only one replicate is considered).

If a peptide is shared across several proteins, the protein column will also contains all the shared protein ids usually separed by _;_ or _,_.
In case a peptide is not quantified it has 0 as intensities. The peptide summary intensity could be used for downstream statistical analysis such as in MsQRob


NOTE : The log files and the output files are in the output folder specified by the user.

[Go to top of page](#moff)
