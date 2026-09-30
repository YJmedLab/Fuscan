# Fuscan
Fuscan is a robust DNA fusion caller for targeted sequencing data in cancer diagnostics. It allows for the personalized selection of interested genomic regions, and predicts their potential breakpoints accurately.
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
- SnpEff >= 4.3t
## Usage
The most basic usage of Fuscan is as follows:
```bash
Fuscan -f <REFERENCE_FASTA> -b <BED> -R1 <R1> -R2 <R2> -o <OUTDIR>
```
If the same bed file is used for long-term, you can utilize `Fuscan_pre` to create a pkl file, thereby avoiding the need to repeatedly analyze the homologous regions.
```bash
Fuscan_pre -f <FASTA> -b <BED> -o <OUTDIR>
Fuscan -f <REFERENCE_FASTA> -p <PKL> -R1 <R1> -R2 <R2> -o <OUTDIR>
```
When reads are already aligned to the same reference as `-f`, you can pass a name-sorted or coordinate-sorted BAM with `-bam`. The pipeline will sort it by read name for improper-read selection; `-R1` and `-R2` must still be the matching FASTQ (or FASTQ.gz) files for the same library, because they are used to extract read sequences for targeted re-alignment.
```bash
Fuscan -f <REFERENCE_FASTA> -p <PKL> -bam <BAM> -R1 <R1> -R2 <R2> -o <OUTDIR>
```
Fuscan can use `Fuscan_bg` to construct a background with negative samples without gene fusion (such as leukocyte samples from healthy individuals), which improves the specificity of detection.
```bash
Fuscan_bg -i <INPUTDIR> -fp <FUSION_PAIR> -o <OUTDIR>
Fuscan -f <REFERENCE_FASTA> -p <PKL> -bg <BACKGROUND_PKL> -R1 <R1> -R2 <R2> -o <OUTDIR>
```
Based on the detection experience with different sample types, you can provide key gene fusion pairs and set different thresholds of split reads and discordant reads counts for the fusion pairs of interest, as well as for intron, exon, and intergenic regions.
The format for the gene fusion pairs file uses a tab to separate the two genes, with the partner genes separated by commas, as shown below:
```txt
ALK	CLTC,EML4,HIP1,KIF5B,KLC1,MSN,STRN
RET	CCDC6,CUX1,KIF5B,NCOA4,PRKAR1A
ROS1	CCDC6,CD74,CLTC,EZR,GOPC,LRIG3
```
```bash
Fuscan -f <REFERENCE_FASTA> -p <PKL> -bg <BACKGROUND_PKL> -fp fusion_pair.txt -ts 1,1,10,10,15,15,20,20 -R1 <R1> -R2 <R2> -o <OUTDIR>
```
Optional `-bp` supplies a blacklist of reference positions (one per line: `chrom<TAB>pos`) that are labeled `BLACK` in the summary when either breakpoint falls on a listed site. `-td` should be a whitespace-separated table of gene symbols and transcription direction (`5>3` or `3>5`); the summary step uses it when writing `results_summary.filtered.txt` so fusion partner order follows strand context (provide for all genes you expect in the output, e.g. those in your panel).

> This tool uses SnpEff for gene annotation. The data and jar package are included in the repository, and the hg19 genome is used by default.
## Parameters
```txt
usage: Fuscan -f <FASTA> -b <BED> -R1 <R1> -R2 <R2> -o <OUTDIR>

A robust DNA fusion caller for targeted sequencing data in cancer diagnostics.

options:
  -h, --help            show this help message and exit
  -f , --FASTA          FASTA file of reference sequence
  -b , --BED            BED file of targeted regions of interested genes
  -p , --PKL            PKL file generated from Fuscan_pre
  -R1                   FASTQ or FASTQ.gz file of R1 reads
  -R2                   FASTQ or FASTQ.gz file of R2 reads
  -bam , --BAM          BAM file of R1 and R2 reads
  -o , --OUTDIR         Output directory for results [default: current
                        directory]
  -bg , --BG            Background PKL file generated from Fuscan_bg [default:
                        None]
  -bp , --BLACK_POS     Blacklisted positions to exclude [default: None]
  -fp , --FUSION_PAIR   Tab-separated file of interested fusion pairs
                        [default: None]
  -td , --TRANSCRIPT_DIR
                        Transcription direction of interested genes [default:
                        None]
  -ts , --THRESHOLD     Threshhold of split reads and discordant reads count
                        for interested fusion pairs, intron, exon and
                        intergenic region [default: 1,1,10,10,15,15,20,20]
  -t , --THREADS        Number of threads to use [default: 1]
  -v, --version         show program's version number and exit
```
## Results
**Improper_ratio**: This refers to the proportion of detected fusions among all improper mapped reads. Generally, the proportion of true fusions is relatively high. The default threshold of Fuscan is 0.05.

**Background**：`BIB`: Both breakpoints are in the background; `2R`: Backgrounds that occur more than twice during the construction of the background; `IB`: Only the breakpoint from the partner gene is located in the background library.

**Other flags in the summary**: `BLACK` indicates a breakpoint overlapping a position from `-bp`; `LOW_DEPTH` and `NO_SA` reflect optional depth and supplementary-alignment filters applied when a background library is provided.

## Citation

If you use Fuscan in your research, please cite:

> Liu Z, Wang S, Chen S, Feng H, Hu X, Zhou P, Shi D. Fuscan: a robust DNA fusion caller for targeted sequencing data in cancer diagnostics. *Bioinformatics Advances*. 2026;6(1):vbag152. doi: [10.1093/bioadv/vbag152](https://doi.org/10.1093/bioadv/vbag152)
