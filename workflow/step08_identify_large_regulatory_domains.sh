#!/bin/bash

## submission properties

#SBATCH --partition=batch
#SBATCH --job-name=STEP07_identify_TFBS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-24:00:00
#SBATCH --mem=100gb
#SBATCH --output=LOGS_step07_identify_TFBS.%j.log
#SBATCH --error=LOGS_step07_identify_TFBS.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules
ml MEME/5.4.1-foss-2019b-Python-3.7.4

# variables
threads=16
ref=../Genome_Reference/Zm-B73-REFERENCE-NAM-5.0.fa
peaks=../01_Peak_Analysis/STARR_merged_peaks.enhancer_activity.eFDR05.bed
controls=../01_Peak_Analysis/STARR_CONTROL.enhancer_activity.bed
motifs=./motif_databases/ARABD/ArabidopsisDAPv1.meme

# extract fasta sequences
#bedtools getfasta -bed $peaks -fi $ref -fo $peaks.fasta
bedtools getfasta -bed $controls -fi $ref -fo $controls.fasta

# identify putative TFBS
#fimo --oc TFBS_peaks $motifs $peaks.fasta
fimo --oc TFBS_controls $motifs $controls.fasta

# reformat fimo output (filtering p-value > 1e-5) using the perl script provided in the github repository (/path/to/script)
perl convertMotifCoord.pl TFBS_peaks/fimo.gff | sed -e 's/_tnt//g' - | sort -k1,1 -k2,2n - > TFBS_peaks.motifs.bed
perl convertMotifCoord.pl TFBS_controls/fimo.gff | sed -e 's/_tnt//g' - | sort -k1,1 -k2,2n - > TFBS_controls.motifs.bed


# annotate motif coverage/counts for STARR and control peaks
bedtools annotate -i ../01_Peak_Analysis/STARR_merged_peaks.enhancer_activity.eFDR05.bed -files TFBS_peaks.motifs.bed -both | sort -k1,1 -k2,2n - > STARR_merged_peaks.enhancer_activity.eFDR05.ann.bed
bedtools annotate -i ../01_Peak_Analysis/STARR_CONTROL.enhancer_activity.bed -files TFBS_controls.motifs.bed -both | sort -k1,1 -k2,2n - > STARR_CONTROL.enhancer_activity.ann.bed

# extract genes
perl -ne 'if($_ =~ /^#/){next;}chomp;my@col=split("\t",$_);if($col[2] eq 'gene'){print"$_\n";}' ../Genome_Reference/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 | sort -k1,1 -k4,4n - > ../Genome_Reference/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.genes.gff3

# classify genomic context of STARR and control peaks
bedtools closest -a STARR_merged_peaks.enhancer_activity.eFDR05.ann.bed -b ../Genome_Reference/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.genes.gff3 -D b > STARR_merged_peaks.enhancer_activity.eFDR05.ann2.bed
bedtools closest -a STARR_CONTROL.enhancer_activity.ann.bed -b ../Genome_Reference/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.genes.gff3 -D b > STARR_CONTROL.enhancer_activity.ann2.bed

# clean up
mv STARR_merged_peaks.enhancer_activity.eFDR05.ann2.bed STARR_merged_peaks.enhancer_activity.eFDR05.ann.bed
mv STARR_CONTROL.enhancer_activity.ann2.bed STARR_CONTROL.enhancer_activity.ann.bed

