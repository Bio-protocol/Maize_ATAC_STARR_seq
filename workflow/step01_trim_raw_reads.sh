#!/bin/bash

## submission properties

#SBATCH --partition=batch
#SBATCH --job-name=STEP01_trim_raw_reads
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-12:00:00
#SBATCH --mem=50gb
#SBATCH --output=LOGS_step01_trim_raw_reads.%j.log
#SBATCH --error=LOGS_step01_trim_raw_reads.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules
ml fastp/0.23.2

# variables
fastqdir=$PWD/FASTQ_files
dna=B73_maize_DNA_input
rna=B73_maize_mRNA_output
threads=16

# trim and filter DNA input reads
fastp -j $dna.json -h $dna.html -w $threads \
	-i $fastqdir/${dna}_1.fastq.gz -I $fastqdir/${dna}_2.fastq.gz \
	-o $fastqdir/${dna}_1.trim.fastq.gz -O $fastqdir/${dna}_2.trim.fastq.gz

# trim and filter mRNA output reads
fastp -j $rna.json -h $rna.html -w $threads \
	-i $fastqdir/${rna}_1.fastq.gz -I $fastqdir/${rna}_2.fastq.gz \
	-o $fastqdir/${rna}_1.trim.fastq.gz -O $fastqdir/${rna}_2.trim.fastq.gz
