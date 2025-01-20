# Fuscan
Fuscan is an ultra-sensitive and noise-free fusion detector that can be used to analyze paired-end DNA sequencing data. It allows for the personalized selection of interested genomic regions, and predicts their potential breakpoints accurately.
## Getting started
```bash
git clone https://github.com/YJmedLab/Fuscan.git
```
## Requirement
- bwa >= 0.7.18
- BLAT >= v.39
- samtools >= 1.9
- bedtools >= v2.29
- seqkit >= v2.9.0
- java >= 11.0.25
- singularity >= 3.11.3
## Usage
The most basic usage of Fuscan is as follows:
```bash
bin/Fuscan -f <REFERENCE_FASTA> -b <BED> -R1 <R1> -R2 <R2> -o <OUTDIR>
```
If the same bed file is used for long-term, you can utilize `Fuscan_pre` to create a pkl file, thereby avoiding the need to repeatedly analyze the homologous regions.
```bash
bin/Fuscan_pre -f <FASTA> -b <BED> -o <OUTDIR>
bin/Fuscan -f <REFERENCE_FASTA> -p <PKL> -R1 <R1> -R2 <R2> -o <OUTDIR>
```
Fuscan can use `Fuscan_bg` to construct a background with negative samples without gene fusion (such as leukocyte samples from healthy individuals), which improves the specificity of detection.
```bash
bin/Fuscan_bg -i <INPUTDIR> -fp <FUSION_PAIR> -o <OUTDIR>
bin/Fuscan -f <REFERENCE_FASTA> -p <PKL> -bg <BACKGORUND_PKL> -R1 <R1> -R2 <R2> -o <OUTDIR>
```
Based on the detection experience with different sample types, you can provide key gene fusion pairs and set different thresholds of split reads and discordant reads counts for the fusion pairs of interest, as well as for intron, exon, and intergenic regions.
The format for the gene fusion pairs file uses a tab to separate the two genes, with the partner genes separated by commas, as shown below:
```txt
ALK	CLTC,EML4,HIP1,KIF5B,KLC1,MSN,STRN
RET	CCDC6,CUX1,KIF5B,NCOA4,PRKAR1A
ROS1	CCDC6,CD74,CLTC,EZR,GOPC,LRIG3
```
```bash
bin/Fuscan -f <REFERENCE_FASTA> -p <PKL> -bg <BACKGORUND_PKL> -fq fusion_pair.txt -ts 1,1,10,10,15,15,20,20 -R1 <R1> -R2 <R2> -o <OUTDIR>

## Parameters
```txt

An ultra-sensitive and noise-free fusion detector for genomic breakpoints discovery in DNA sequencing data.

options:
  -h, --help            show this help message and exit
  -f , --FASTA          FASTA file of reference sequence
  -b , --BED            BED file of targeted regions of interested genes
  -p , --PKL            PKL file generated from Fuscan_pre
  -R1                   FASTQ or FASTQ.gz file of R1 reads
  -R2                   FASTQ or FASTQ.gz file of R2 reads
  -o , --OUTDIR         Output directory for results [default: current directory]
  -bg , --BG            Background PKL file generated from Fuscan_bg [default: None]
  -fp , --FUSION_PAIR   Tab-separated file of interested fusion pairs [default: None]
  -ts , --THRESHOLD     Threshhold of split reads and discordant reads count for interested fusion pairs, intron, exon and intergenic region [default: 1,1,10,10,15,15,20,20]
  -t , --THREADS        Number of threads to use [default: 1]
  -v, --version         show program's version number and exit
```
## Results
**Improper_ratio**: This refers to the proportion of detected fusions among all improper mapped reads. Generally, the proportion of true fusions is relatively high. The default threshold of Fuscan is 0.05.

**Background**：`BIB`: Both breakpoints are in the background; `2R`: Backgrounds that occur more than twice during the construction of the background; `IB`: Only the breakpoint of the reference genome is located in the background library.
