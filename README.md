## Bulk RNA-Sequencing pipeline
Contains scripts for bulk RNA-Seq preprocessing

<!-- badges: start -->
![GitHub last commit](https://img.shields.io/github/last-commit/fouratilab/RNAseq)
<!-- badges: end -->

## Dependencies
- RSeQC (version 3.0.1)
- STAR aligner (version 2.7.0e)
- HTSeq (version 0.11.2)

## Usage
#### Rename raw sequencing files
The pipeline expect raw sequencing files to be named using the convention
*SampleID*_*mate*.fq.gz where:  
*SampleID* can be any string  
*mate* is the digit [1,2]  
ex: S01_1.fq.gz

#### Run RNA-Seq pipeline
The RNA-Seq pipeline will map sequencing reads to a reference
genome
```bash
bash rna.preprocessing_master.sh -d {raw_directory}

arguments:  
d=[d]irectory with raw data (directory; required)  
g=reference [g]enome  
    accepted values: GRCh38, GRCm39  
p=[p]air-end sequencing files  
    if arugment -p is not given, single-end sequencing files are  
    expected
e=[e]mail address  
i=[i]soform transcript/exon counts  
h=print [h]elp
```
