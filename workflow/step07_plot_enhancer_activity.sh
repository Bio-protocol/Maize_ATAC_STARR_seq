#!/bin/bash

## submission properties

#SBATCH --partition=schmitz_p
#SBATCH --job-name=STEP07_plot_starr_data
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=0-12:00:00
#SBATCH --mem=100gb
#SBATCH --output=LOGS_step07_plot_STARR_data.%j.log
#SBATCH --error=LOGS_step07_plot_STARR_data.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# load module
module load deepTools/3.5.1-intel-2020b-Python-3.8.6

# load function
getmaps(){

        # input
        ina=$1
        id=$2
        dat=../BED_files/*.bw

        # output
        outa=$id.mat.gz
        outm=$id.mat.txt
        fig=$id.pdf

        # parameters
        threads=48
        window=2000
        bin=20
        cols=YlGnBu

        # create matrix
        computeMatrix reference-point --referencePoint center \
		-S $dat \
                -b $window -a $window \
                -R $ina \
		--missingDataAsZero \
                -o $outa \
		--outFileNameMatrix $outm \
                -p $threads --binSize $bin

	# plot heatmap
	plotHeatmap --matrixFile $outa \
		--colorMap YlGnBu \
		-out $id.heatmap.pdf

}
export -f getmaps

# run for each file
getmaps STARR_merged_peaks.enhancer_activity.bed STARR_peaks
getmaps STARR_merged_peaks.enhancer_activity.eFDR05.bed STARR_peaks_filtered
getmaps STARR_CONTROL.enhancer_activity.bed control_regions
