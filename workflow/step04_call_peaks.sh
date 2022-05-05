#!/bin/bash

## submission properties

#SBATCH --partition=batch
#SBATCH --job-name=STEP04_find_regulatory_regions
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00:00
#SBATCH --mem=50gb
#SBATCH --output=LOGS_step04_find_regulatory_regions.%j.log
#SBATCH --error=LOGS_step04_find_regulatory_regions.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules
ml MACS2/2.2.7.1-foss-2019b-Python-3.7.4

# variables
beddir=$PWD/BED_files
dna=B73_maize_DNA_input
rna=B73_maize_mRNA_output

# generate input files
uniq $beddir/$rna.fragments.bed > $beddir/$rna.fragments.uniq.bed
uniq $beddir/$dna.fragments.bed > $beddir/$dna.fragments.uniq.bed
cat $beddir/$rna.fragments.uniq.bed $beddir/$dna.fragments.uniq.bed | sort -k1,1 -k2,2n - > $beddir/ALL.fragments.uniq.bed

# find regulatory regions
echo " calling regulatory regions without duplicate removal ..."
macs2 callpeak -t $beddir/$rna.fragments.bed \
	-c $beddir/ALL.fragments.uniq.bed \
	--keep-dup all \
	--max-gap 50 \
	--min-length 300 \
	--nolambda \
	-f BEDPE \
	-g 1.6e9 \
	--bdg \
	-n STARR_wdups

# find regulatory regions
echo " calling regulatory regions with duplicate removal ..."
macs2 callpeak -t $beddir/$rna.fragments.uniq.bed \
        -c $beddir/$dna.fragments.uniq.bed \
        --keep-dup all \
	--max-gap 50 \
	--min-length 300 \
	--nolambda \
        -f BEDPE \
        -g 1.6e9 \
        --bdg \
        -n STARR_nodups

# find regulatory regions
echo " calling regulatory regions by aggregating all unique fragments ..."
macs2 callpeak -t $beddir/ALL.fragments.uniq.bed \
	--keep-dup all \
	--max-gap 50 \
	--min-length 300 \
	--nolambda \
	-f BEDPE \
	-g 1.69e9 \
	--bdg \
	-n STARR_ALL

# clean-up output
mv STARR_* Peak_data

# merge peaks
cd Peak_data
cat STARR_wdups_peaks.narrowPeak STARR_nodups_peaks.narrowPeak STARR_ALL_peaks.narrowPeak | sort -k1,1 -k2,2n - | bedtools merge -i - > STARR_merged_peaks.bed
