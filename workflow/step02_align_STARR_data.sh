#!/bin/bash

## submission properties

#SBATCH --partition=batch
#SBATCH --job-name=STEP02_align_STARR_data
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=0-12:00:00
#SBATCH --mem=100gb
#SBATCH --output=LOGS_step02_align_STARR_data.%j.log
#SBATCH --error=LOGS_step02_align_STARR_data.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules
ml BWA/0.7.17-GCC-8.3.0
ml SAMtools/1.14-GCC-8.3.0

# variables
outdir=$PWD/BAM_files
refdir=$PWD/Genome_Reference
fastqdir=$PWD/FASTQ_files
dna=B73_maize_DNA_input
rna=B73_maize_mRNA_output
mq=10

# align DNA input and pipe to samtools for SAM to BAM conversion
echo " aligning STARR DNA input ..."
bwa mem -M -t 24 $refdir/Zm-B73-REFERENCE-NAM-5.0.fa \
	$fastqdir/${dna}_1.trim.fastq.gz $fastqdir/${dna}_2.trim.fastq.gz \
	| samtools view -bSh - \
	| samtools sort - > $outdir/$dna.raw.bam

# align RNA output and pipe to samtools for SAM to BAM conversion
echo " aligning STARR mRNA output ..."
bwa mem -M -t 24 $refdir/Zm-B73-REFERENCE-NAM-5.0.fa \
	$fastqdir/${rna}_1.trim.fastq.gz $fastqdir/${rna}_2.trim.fastq.gz \
	| samtools view -bSh - \
	| samtools sort - > $outdir/$rna.raw.bam

# filter DNA input alignments
echo " filtering STARR DNA alignments ..."
samtools view -h -q $mq -f 3 $outdir/$dna.raw.bam \
	| grep -v -E -e '\bXA:Z:' \
	| samtools view -bSh - > $outdir/$dna.$mq.pp.unique.bam

# filter mRNA output alignments
echo " filtering STARR mRNA alignments ..."
samtools view -h -q $mq -f 3 $outdir/$rna.raw.bam \
	| grep -v -E -e '\bXA:Z:' \
	| samtools view -bSh - > $outdir/$rna.$mq.pp.unique.bam


