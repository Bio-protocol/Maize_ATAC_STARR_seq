[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# Computational Analysis of Maize Enhancer REgulatory Elements Using ATAC-STARR-seq

The blueprints to development, response to the environment, and cellular function are largely the manifestation of distinct gene expression programs controlled by the spatiotemporal activity of cis-regulatory elements. Although biochemical methods for identifying accessible chromatin – a hallmark of cis-regulatory elements – have been developed, approaches capable of measuring and quantifying cis-regulatory activity are only beginning to be realized. Massively Parallel Reporter Assays coupled to chromatin accessibility profiling presents a high-throughput solution for testing the transcription activating capacity of millions of putatively regulatory DNA sequences. However, clear computational pipelines for analyzing these high-throughput sequencing-based reporter assays are lacking. In this protocol, I layout and rationalize a computational framework for the processing and analysis of Assay for Transposase Accessible Chromatin profiling followed by Self-Transcribed Active Regulatory Region sequencing (ATAC-STARR-seq) data from a recent study in Zea mays. The approach described herein can be adapted to other sequencing-based reporter assays and it largely agnostic to the model organism.   


## Software dependencies
[BWA MEM](http://bio-bwa.sourceforge.net/bwa.shtml)
[SAMtools](http://www.htslib.org)
[BEDtools](https://bedtools.readthedocs.io/en/latest/)
[SRA-toolkit](https://github.com/ncbi/sra-tools)
[fastp](https://github.com/OpenGene/fastp)
[MACS2](https://pypi.org/project/MACS2/)
[UCSC binaries](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
[tabix](http://www.htslib.org/doc/tabix.html)
[IGV](https://software.broadinstitute.org/software/igv)
[MEME](https://meme-suite.org/meme/index.html)
[CrossMap](http://crossmap.sourceforge.net)
[DeepTools](https://deeptools.readthedocs.io/en/develop/index.html)


## Input data
The computational pipeline uses paired-end sequencing data from an ATAC-STARR-seq experiment performed on maize protoplasts (Ricci et al., 2019). The ATAC-STARR-seq experiment consisted of a DNA input (ATAC-seq library) and a mRNA readout (self-transcribed regulatory regions) to identify genomic regions exhibiting transcription-activating regulatory activity. 

1.	Transfected ATAC-seq DNA-input [FASTQ](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10964904)
2.	Transcribed ATAC-seq mRNA [FASTQ](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10964905)


## Procedure
Additional details can be found in [paper](https://bio-protocol.org/default.aspx).

1.	Download data

```
# set variables and download FASTQ files
mkdir FASTQ_files
cd FASTQ_files
fasterq-dump -o B73_maize_DNA_input.fastq SRR10964904
fasterq-dump -o B73_maize_mRNA_output.fastq SRR10964905

# compress fastq files
pigz *.fastq

# NOT RUN
# Tip: gzip can be used as an alternative to pigz (parallel gzip)
# gzip *.fastq

# download reference data
cd ../
mkdir Genome_Reference
cd Genome_Reference
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz

# create indices for reference genome FASTA
gunzip Zm-B73-REFERENCE-NAM-5.0.fa.gz
samtools faidx Zm-B73-REFERENCE-NAM-5.0.fa
bwa index Zm-B73-REFERENCE-NAM-5.0.fa
```

2.	Trim adapters and remove low quality reads

```
# run step 1
sbatch step01_trim_raw_reads.sh
```

3.	Align and process sequenced reads
```
# run step 2
sbatch step02_align_STARR_data.sh
```

4.	Extract fragments
```
# run step 3
sbatch step03_extract_fragments.sh
```

5.	Call peaks
```
# run step 4
sbatch step04_call_peaks.sh
```

6.	Estimate enhancer activity
```
# run step 5
sbatch step05_estimate_enhancer_activity.sh

# estimate enhancer activity
cd BED_files/
Rscript Estimate_Enhancer_Activity.R
```

7.	Filter noisy STARR peaks using empirical FDR
```
# run step 6
sbatch step06_create_control_regions.sh

# create directory to contain analysis
cd ../
mkdir 01_Peak_Analysis
cd 01_Peak_Analysis

# map maximum enhancer activity to putative regulatory regions (wdups)
bedtools map -a ../Peak_data/STARR_merged_peaks.bed -b ../BED_files/B73_maize.enhancer_activity.bdg -o max -c 4 > STARR_merged_peaks.enhancer_activity.bed

# map maximum enhancer activity to control 
bedtools map -a ../Peak_data/STARR_CONTROL.bed -b ../BED_files/B73_maize.enhancer_activity.bdg -o max -c 4 > STARR_CONTROL.enhancer_activity.bed

# run eFDR filter
Rscript eFDR_Filter_STARR_Peaks.R

```

8.	Plot heatmaps
```
# run step 7
sbatch step07_plot_enhancer_activity.sh
```

## Downstream analysis
See the [paper](https://bio-protocol.org/default.aspx) for a downstream analysis use case.

