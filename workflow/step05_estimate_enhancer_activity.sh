#!/bin/bash

## submission properties

#SBATCH --partition=batch
#SBATCH --job-name=STEP05_quantify_regulatory_regions
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00:00
#SBATCH --mem=50gb
#SBATCH --output=LOGS_step05_quantify_regulatory_regions.%j.log
#SBATCH --error=LOGS_step05_quantify_regulatory_regions.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules
ml MACS2/2.2.7.1-foss-2019b-Python-3.7.4

# variables
beddir=$PWD/BED_files
dna=B73_maize_DNA_input
rna=B73_maize_mRNA_output
ref=./Genome_Reference/Zm-B73-REFERENCE-NAM-5.0.fa.fai

# sort the reference
sort -k1,1 -k2,2n $ref > $ref.sorted

# merge RNA/DNA
cat $beddir/$rna.fragments.uniq.bed $beddir/$dna.fragments.uniq.bed \
	| sort -k1,1 -k2,2n - \
	| bedtools genomecov -i - -bga -g $ref.sorted \
	| sort -k1,1 -k2,2n - \
	| cut -f1-3 - > $beddir/Unique_genomic_intervals.bed

# count fragments and add pseudocount
bedtools intersect -a $beddir/Unique_genomic_intervals.bed \
	-b $beddir/$rna.fragments.bed \
	-c -sorted -g $ref.sorted > $beddir/$rna.activity.raw.bed

bedtools intersect -a $beddir/$rna.activity.raw.bed \
        -b $beddir/$dna.fragments.bed \
        -c -sorted -g $ref.sorted > $beddir/B73_maize_mRNA_DNA.activity.raw.bed

# clean up temporary files
rm $beddir/Unique_genomic_intervals.bed
rm $beddir/$rna.activity.raw.bed
