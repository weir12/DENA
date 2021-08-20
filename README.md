
![image]( ./DENA_log.jpg)
![image](https://img.shields.io/badge/Python-3.7-green) ![image](https://img.shields.io/badge/Pytorch-1.8-blue) ![image](https://img.shields.io/badge/Centos-8-yellow)
- [DENA (Deeplearning Explore Nanopore m6A)](#dena-deeplearning-explore-nanopore-m6a)
  - [Author: liangou](#author-liangou)
  - [E-mail:liangou@ips.ac.cn](#e-mailliangouipsaccn)
  - [Getting Started](#getting-started)
  - [0.Prerequisites](#0prerequisites)
  - [1.Input data required](#1input-data-required)
    - [1.Obtain coordinates matching motif in reference](#1obtain-coordinates-matching-motif-in-reference)
      - [New version function](#new-version-function)
  - [2.Signal re-sqguiggle and sequence alignment](#2signal-re-sqguiggle-and-sequence-alignment)
    - [2.1 tombo re-sqguiggle](#21-tombo-re-sqguiggle)
    - [2.2 sequence alianment based on minimap2](#22-sequence-alianment-based-on-minimap2)
  - [3.extract features](#3extract-features)
      - [New version function](#new-version-function-1)
      - [Parameters panel](#parameters-panel)
    - [4.Predict(v3.0)](#4predictv30)
      - [New version function](#new-version-function-2)
      - [Parameters](#parameters)
      - [Requirement](#requirement)
      - [Example](#example)
    - [](#)
    - [TroubleShoot](#troubleshoot)
  - [Utils](#utils)
    - [1.Dimension reduction & Cluster of a dataset](#1dimension-reduction--cluster-of-a-dataset)
    - [2.Absolute difference of mean](#2absolute-difference-of-mean)
#  DENA (Deeplearning Explore Nanopore m6A)
 
Deep learning model used to detect RNA m6a with read level based on the Nanopore direct RNA data.

## Author: liangou
## E-mail:liangou@ips.ac.cn
## Getting Started
 
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.
 
## 0.Prerequisites
Utilizing Conda or virtualenv to create a relatively independent & clean work environment may be a wise choice for using DENA 
Here are What things you need to install(Please confirm one by one):
1. Unix like system(centos,ubuntu,etc)
2. Cuda-supported graphics cards(optional)
3. Python>=3.7.x and Pytorch
4. tombo,minimap2,samtools
## 1.Input data required
1. a batch of fast5 files containing the raw current signals
2. a fastq file which is contain basecalled sequence corresponding fast5 above
3. Appropriate reference sequence(Transcriptome is recommended for RNA data)


Tips:

${variable} : You need to assign it with the your actual value
### 1.Obtain coordinates matching motif in reference
```bash
python3 LSTM_extract.py get_pos --fasta ${fasta_fn}  --motif 'RRACH' --output ./candidate_predict_pos.txt
```
You will get result(candidate_predict_pos.txt) like this
```
AT1G01010.1     17      22      +       AAACC
```
#### New version function

- We no longer need external C++ tools
- Fixed compatibility bugs in FASTA file of some species

## 2.Signal re-sqguiggle and sequence alignment
### 2.1 tombo re-sqguiggle
This step is to obtain a unique mapping between the signal fragment of each base of each reads and the reference sequence
For detailed help, please see https://github.com/nanoporetech/tombo
```
tombo resquiggle --processes {thread} --ignore-read-locks --max-scaling-iterations 5 --rna --basecall-group "Basecall_1D_001" --num-most-common-errors 5 --include-event-stdev --overwrite --signal-length-range 0 500000 {params.fast5} {params.ref}
```

### 2.2 sequence alianment based on minimap2
For detailed help, please see [minimap2](https://github.com/lh3/minimap2) [samtools](https://github.com/samtools/samtools) 

```
minimap2 -a -uf -k10 --sam-hit-only --secondary=no ${transcriptome} basecalls.reverse.fa > basecalls.sam
samtools flagstat basecalls.sam
samtools view -bS -F 2048 -F 16 -F 4 basecalls.sam >basecalls.bam
samtools sort -@ 12 basecalls.bam>basecalls.sort.bam
samtools index basecalls.sort.bam
```


## 3.extract features

#### New version function

- Support for reading BAM files in BRI mode to reduce memory consumption

Install the C ++ libraries and Python wrappers to enable this functionality
[https://github.com/nanoporetech/bripy](https://github.com/nanoporetech/bripy) [https://github.com/jts/bri](https://github.com/jts/bri)

- Flexible window Settings are now supported
- In this step,you need provide two input params for program:fast5_folder(re-squiggle by tombo) and bam file(sorted & index)
#### Parameters panel
```python
	parser.add_argument('--processes',default=24,type=int,
						help=("Number of processes allocated"))
	parser_a = subparsers.add_parser('get_pos',formatter_class=argparse.RawDescriptionHelpFormatter,help='get candidate position')
	parser_a.add_argument('--fasta',required=True, default=None,
						help=("reference fasta"))
	parser_a.add_argument('--motif',  default='RRACH',
						help=("specifies a motif pattern"))
	parser_a.add_argument('--output', default='./candidate_predict_pos.txt',
						help=("output file"))
	parser_a.set_defaults(func=get_pos)
	parser_b = subparsers.add_parser('predict',formatter_class=argparse.RawDescriptionHelpFormatter,help='predict')
	parser_b.add_argument('--fast5',required=True, default=None,
						help=("a directory(has been re-squiggled by tombo) that contains the FAST5 files"))	
	parser_b.add_argument('--corr_grp',default="RawGenomeCorrected_000",
						help=("Analysis slot containing the re-squiggle information"))		
	parser_b.add_argument('--bam',required=True, default=None,	
						help=("BAM file used to extract base-quality(feature)"))
	parser_b.add_argument('--sites',default='./candidate_predict_pos.txt',	
						help=("candidate position are used to extract features of mapped reads"))
	parser_b.add_argument('--label',required=True,	
						help=("The string used to distinguish the sample"))	
	parser_b.add_argument('--windows',required=True,nargs=2,metavar='3',type=int,
						help=("Window drift away from the center of m6A"))
	parser_b.add_argument('--debug',action='store_true',default=False,
						help=("Enable debug mode (output more detailed run log)"))
	parser_b.add_argument('--bri',action='store_true',default=False,
						help=("Enable BRI mode (Reduce RAM consumption of BAM files)"))	
```

```bash
python3 LSTM_extract.py --fast5 ${fast5_fn}  --corr_grp ${RawGenomeCorrected_000} --bam ${bam_fn}  --sites ${candidate_predict_pos.txt} --label ${any meaningful string} --windows 3 3
```

- You will get result(*.tmp) like this
```
>AT4G35300.4_2258_GGACT
37b79f1c-c3c2-4c6f-a25c-65e618b7bb6f    28.0,27.0,24.0,31.0,24.0,74.0,6.0,33.0,27.0,63.0,2.6680189601886424,2.4252261046166588,0.16661375589914051,-0.4574264926055352,0.8548283129364287,2.6957830373199303,2.568959475115271,0.08321765590395497,-0.5128530864579424,0.8695237415728408,0.45954324955726,0.7483187244003591,0.23116851931302654,0.21901705697503512,0.23758996696952467
e944b3ff-156c-409f-95f3-996dfa3d3fd3    26.0,30.0,28.0,25.0,31.0,158.0,6.0,47.0,37.0,27.0,2.3582877438291905,2.354059678502405,-0.22480153918158693,-0.5296594244218555,1.084794467401169,2.5690832257777787,2.4814037210635487,-0.24292374684288695,-0.5435391915773902,1.122371397992982,0.7615589741556634,0.6712493289989668,0.19113118923616254,0.20057491958760532,0.21721818493112435
```
### 4.Predict(v3.0)
Tips :If the input features  **NOT changed**  here is **NO** need to repeat run step 2
#### New version function
add "-d" in cmd for output m6a probability for each read at each site
 Added support for deep learning
**Caution** :Using deep learning model will occupy a lot of computing resources and time costing without GPU
#### Parameters
In this step,you need Provide the following parameters:

1. ${path_features} :Path contain [0-9]*.tmp (generated by step 2)
1. ${path_models} :Path contain *.dat(ensemble learning) or *.pkl (deep learning) 
1. ${path_output} The output path
1. ${prefix_outfile} The prefix of the output file
#### Requirement

- pip install lightgbm
- pip install xgboost
- conda install pytorch torchvision torchaudio cpuonly -c pytorch
#### Example
```bash
python LSTM_predict.py -i ${path_features} -m ${path_models} -o ${path_output} -p ${prefix_outfile} -d 
```

- You will get result **${prefix_outfile}**.tsv in **${path_output}** like this
```
AT1G01010.1     30      AAACA   0       3       0.0
AT1G01010.1     212     AAACA   0       5       0.0
AT1G01010.1     341     AAACA   1       5       0.2
AT1G01020.2     679     AAACA   2       2       1.0
AT1G01030.1     306     AAACA   0       10      0.0
AT1G01030.1     422     AAACA   0       10      0.0
AT1G01030.1     726     AAACA   0       11      0.0
AT1G01030.1     838     AAACA   1       11      0.09090909090909091
AT1G01030.1     876     AAACA   1       11      0.09090909090909091
AT1G01030.1     1233    AAACA   0       11      0.0
```

- if "-d" was added,You will get result **${prefix_outfile}_details**.tsv in **${path_output}** like this
```
AT1G01010.1	30
b129005a-01e7-49f8-bb50-aadf0d57f079	0.235
85a565d6-08c1-4819-a12a-2beacbf63319	0.647
a5fafedc-1539-4cb8-ad6d-72c8699220cd	0.498
```
### 
### TroubleShoot

- Make sure the rules for gene names are consistent among bam file,fast5 files and fasta file
## Utils
### 1.Dimension reduction & Cluster of a dataset

- **python3 ./pca_cluster.py  ${some params}**
```python
	parser.add_argument('--processes',default=24,type=int,
						help=("Number of processes allocated"))
	parser.add_argument('--input',required=True, default=None,
						help=("A directory containing both 'positive_dataset.txt' and 'negative_dataset.txt'"))
	parser.add_argument('--output', default='./dimRe_clust_fig',
						help=("output directory"))
	parser.add_argument('--pos_fast5',required=True, default=None,
						help=("a directory(has been re-squiggled by tombo) that contains the FAST5 files of positive_dataset"))	
	parser.add_argument('--neg_fast5',required=True, default=None,
						help=("a directory(has been re-squiggled by tombo) that contains the FAST5 files of negative_dataset"))		
	parser.add_argument('--corr_grp',default="RawGenomeCorrected_000",
						help=("Analysis slot containing the re-squiggle information"))		
	parser.add_argument('--windows',required=True,nargs=2,metavar='3',type=int,
						help=("Window drift away from the center of m6A"))
	parser.add_argument('--features',default='norm_mean',nargs=+,choices=['length','norm_mean','norm_med','base_q','norm_stdev']
						help=("Input features for dimensionality reduction clustering"))													
	parser.add_argument('--algorithm_DimRe',default="PCA",choices=['PCA']
						help=("Algorithms for dimension reduction"))
	parser.add_argument('--algorithm_cluster',default="kmeans",choices=['kmeans']
						help=("Algorithms for cluster"))
```

- For each candidate site, the corresponding dimensionality reduction result is printed, like this

BLUE:WT sample RED : KO/KD sample 
![](https://cdn.nlark.com/yuque/0/2021/png/644304/1618468161343-928ec141-b00c-4343-ba92-9bac7aff15e6.png#clientId=ub8b101f3-54a8-4&from=paste&height=263&id=u909b4ce7&margin=%5Bobject%20Object%5D&originHeight=1051&originWidth=1384&originalType=binary&size=130003&status=done&style=none&taskId=u3de8527d-cda4-446a-acf8-6133bf728b0&width=346)
### 2.Absolute difference of mean

- Calculate and visualize the difference between the current mean values of the two samples covered on the candidate sites
- **python3 ./box_plot.py  **
- **Cautions!!!:**Please change the custom variables (below) within the script. This script does not accept CMD parameters
```python
#Defining variables
sites_file="/home/weir/m6a_model/DENA/plotter/final_overlep_sites"
pos_fast5="/home/weir/tair_rawdata/elife_rawdata/VIRc"
neg_fast5="/home/weir/tair_rawdata/elife_rawdata/vir-1"
corrected_group='RawGenomeCorrected_000'

In [4]: !head "/home/weir/m6a_model/DENA/plotter/final_overlep_sites"
AT5G67590.1     787
AT5G67590.1     733
AT5G67560.1     1156
AT5G67560.1     1126
AT5G67510.1     677
AT5G67330.1     2164
AT5G67250.1     2329
AT5G67130.1     1714
AT5G67030.1     2469
AT5G67030.1     2308

```

-  result is printed, like this（Data used for drawing is also saved for downstream analysis and visual adjustment）

![](https://cdn.nlark.com/yuque/0/2021/png/644304/1618468669638-aef7bc8b-726c-40d1-b79c-15f5afa5fcba.png#clientId=ub8b101f3-54a8-4&from=paste&height=331&id=u1122a54a&margin=%5Bobject%20Object%5D&originHeight=1323&originWidth=1822&originalType=binary&size=146707&status=done&style=none&taskId=ub6f0d7ec-52b6-4c98-b06a-badb0dd96f6&width=455.5)
