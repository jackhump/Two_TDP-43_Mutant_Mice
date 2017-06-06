#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

iFolder=/Users/Jack/SAN/IoN_RNAseq/RNA_Maps/data/
clusters=${iFolder}/iCLIP/All_TDP_iCLIP_merged.bed 

# intron beds
F210I_skipped=${iFolder}/F210I_embryonic_brain_se_intron_skipped.bed
F210I_included=${iFolder}/F210I_embryonic_brain_se_intron_included.bed
M323K_skipped=${iFolder}/M323K_adult_brain_se_intron_skipped.bed
M323K_included=${iFolder}/M323K_adult_brain_se_intron_included.bed
controls=${iFolder}/F210I_embryonic_brain_notsig_se_intron_all.bed

# outFolder
outFolder="/Users/Jack/Google Drive/TDP_paper/RNA_maps/results/"

# exon beds
F210I_skipped_exons=${iFolder}/noheader/F210I_embryonic_brain_flank0_se_skipped.bed
F210I_included_exons=${iFolder}/noheader/F210I_embryonic_brain_flank0_se_included.bed
M323K_skipped_exons=${iFolder}/noheader/M323K_adult_brain_flank0_se_skipped.bed
M323K_included_exons=${iFolder}/noheader/M323K_adult_brain_flank0_se_included.bed
control_exons=${iFolder}/noheader/F210I_embryonic_brain_flank0_se_notsig_se_all.bed

# F210I skipped
python whole_intron_cluster_coverage.py $F210I_skipped $clusters ${outFolder}/F210I_skipped

# F210I included 
python whole_intron_cluster_coverage.py $F210I_included $clusters ${outFolder}/F210I_included

# M323K skipped
python whole_intron_cluster_coverage.py $M323K_skipped $clusters ${outFolder}/M323K_skipped

# M323K included
python whole_intron_cluster_coverage.py $M323K_included $clusters ${outFolder}/M323K_included

# controls
python whole_intron_cluster_coverage.py $controls $clusters ${outFolder}/controls



# create F210I RNA Map
Rscript plot_whole_intron_coverage.R --skipped ${outFolder}/F210I_skipped_100_coverage.csv \
				     --included ${outFolder}/F210I_included_100_coverage.csv \
				     --control ${outFolder}/controls_100_coverage.csv \
				     --skipped_exons $F210I_skipped_exons \
				     --included_exons $F210I_included_exons \
				     --control_exons $control_exons \
				     --code F210I \
				     --outFolder $outFolder

# M323K 
Rscript plot_whole_intron_coverage.R --skipped ${outFolder}/M323K_skipped_100_coverage.csv \
                                     --included ${outFolder}/M323K_included_100_coverage.csv \
                                     --control ${outFolder}/controls_100_coverage.csv \
                                     --skipped_exons $M323K_skipped_exons \
                                     --included_exons $M323K_included_exons \
                                     --control_exons $control_exons \
                                     --code M323K \
                                     --outFolder $outFolder



