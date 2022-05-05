#!/bin/bash

## submission properties

#SBATCH --partition=schmitz_p
#SBATCH --job-name=STEP06_create_control_regions
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=0-12:00:00
#SBATCH --mem=115gb
#SBATCH --output=LOGS_step06_create_control_regions.%j.log
#SBATCH --error=LOGS_step06_create_control_regions.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules
ml BWA/0.7.17-GCC-8.3.0
ml SAMtools/1.14-GCC-8.3.0


########################################################################################
## create synthetic reads --------------------------------------------------------------
########################################################################################

# estimate number of reads
mRNA_counts=$(samtools view -c ./BAM_files/B73_maize_mRNA_output.raw.bam)
DNA_counts=$(samtools view -c ./BAM_files/B73_maize_DNA_input.raw.bam)

# move into the Genome_References directory
cd ./Genome_Reference

# generate simulated reads matching the mRNA output
wgsim -1 36 -2 36 -d 300 -N $mRNA_counts Zm-B73-REFERENCE-NAM-5.0.fa simulated_STARR_mRNA_r1.fq simualted_STARR_mRNA_r2.fq 

# generate simulated reads matching the DNA input
wgsim -1 36 -2 36 -d 300 -N $DNA_counts Zm-B73-REFERENCE-NAM-5.0.fa simulated_STARR_DNA_r1.fq simualted_STARR_DNA_r2.fq 

# compress synthetic fastq files
pigz *.fq

# move back to project root
cd ../


#########################################################################################
## map synthetic reads ------------------------------------------------------------------
#########################################################################################

# variables
outdir=$PWD/BAM_files
refdir=$PWD/Genome_Reference
ref=$refdir/Zm-B73-REFERENCE-NAM-5.0.fa.fai.sorted
fastqdir=$refdir
mq=10

# align synthetic mRNA output and pipe to samtools for SAM to BAM conversion
bwa mem -M -t 24 $refdir/Zm-B73-REFERENCE-NAM-5.0.fa \
	$fastqdir/simulated_STARR_mRNA_r1.fq.gz \
	$fastqdir/simulated_STARR_mRNA_r2.fq.gz \
	| samtools view -bS - \
	| samtools sort - > $outdir/simulated_STARR_mRNA.raw.bam

# align synthetic DNA input and pipe to samtools for SAM to BAM conversion
bwa mem -M -t 24 $refdir/Zm-B73-REFERENCE-NAM-5.0.fa \
	$fastqdir/simulated_STARR_DNA_r1.fq.gz \
	$fastqdir/simulated_STARR_DNA_r2.fq.gz \
	| samtools view -bS - \
	| samtools sort - > $outdir/simulated_STARR_DNA.raw.bam

# filter synthetic mRNA alignments
echo " filtering synthetic STARR mRNA alignments ..."
samtools view -h -q $mq -f 3 $outdir/simulated_STARR_mRNA.raw.bam \
       | grep -v -E -e '\bXA:Z:' \
       | samtools view -bSh - > $outdir/simulated_STARR_mRNA.mq$mq.pp.unique.bam

# filter synthetic DNA alignments
echo " filtering STARR DNA alignments ..."
samtools view -h -q $mq -f 3 $outdir/simulated_STARR_DNA.raw.bam \
       | grep -v -E -e '\bXA:Z:' \
       | samtools view -bSh - > $outdir/simulated_STARR_DNA.mq$mq.pp.unique.bam

# extract mRNA output fragments
echo " extracting fragments from simulated STARR mRNA output ..."
samtools sort -n $outdir/simulated_STARR_mRNA.mq$mq.pp.unique.bam \
        | bedtools bamtobed -bedpe -i - \
        | sort -k1,1 -k2,2n - \
        | cut -f1,2,6 - > $refdir/simulated_STARR_mRNA.fragments.bed

# extract DNA input fragments
echo " extracting fragments from simulated STARR DNA input ..."
samtools sort -n $outdir/simulated_STARR_DNA.mq$mq.pp.unique.bam \
        | bedtools bamtobed -bedpe -i - \
        | sort -k1,1 -k2,2n - \
        | cut -f1,2,6 - > $refdir/simulated_STARR_DNA.fragments.bed

# merge all fragments (sorting by coordinate at this step may take a while)
cat $refdir/simulated_STARR_mRNA.fragments.bed $refdir/simulated_STARR_DNA.fragments.bed \
	| sort -k1,1 -k2,2n \
	| bedtools merge -i - > $refdir/mappable_genomic_regions.bed

# create controls
peaks=$PWD/Peak_data/STARR_merged_peaks.bed
bedtools shuffle -i $peaks \
	-g $ref \
	-incl $refdir/mappable_genomic_regions.bed \
	-excl $peaks \
	| sort -k1,1 -k2,2n - > $PWD/Peak_data/STARR_CONTROL.bed
