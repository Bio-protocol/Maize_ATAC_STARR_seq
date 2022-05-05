#!/bin/bash

## submission properties

#SBATCH --partition=batch
#SBATCH --job-name=STEP03_extract_fragments
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00:00
#SBATCH --mem=100gb
#SBATCH --output=LOGS_step03_extract_fragments.%j.log
#SBATCH --error=LOGS_step03_extract_fragments.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules
ml BWA/0.7.17-GCC-8.3.0

# variables
outdir=$PWD/BAM_files
beddir=$PWD/BED_files
dna=B73_maize_DNA_input
rna=B73_maize_mRNA_output

# extract DNA input fragments
echo " extracting fragments from STARR DNA input ..."
samtools sort -n $outdir/$dna.mq10.pp.unique.bam \
	| bedtools bamtobed -bedpe -i - \
	| sort -k1,1 -k2,2n - \
	| cut -f1,2,6 - > $beddir/$dna.fragments.bed

# extract mRNA output fragments
echo " extracting fragments from STARR mRNA output ..."
samtools sort -n $outdir/$rna.mq10.pp.unique.bam \
	| bedtools bamtobed -bedpe -i - \
	| sort -k1,1 -k2,2n - \
	| cut -f1,2,6 - > $beddir/$rna.fragments.bed


