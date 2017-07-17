#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

iFolder=/Users/Jack/SAN/IoN_RNAseq/F210I_M323K_paper/RNA_Maps/data/

# don't use this! just the F210I!
#clusters=${iFolder}/iCLIP/All_TDP_iCLIP_merged.bed 

#clusters=${iFolder}/iCLIP/F210I_WT_rep23_lowFDR_clusters_strand.bed
clusters=/Users/Jack/google_drive/TDP_paper/RNA_maps/F210I_WT_rep23_lowFDR_clusters_strand.bed

# intron beds
F210I_skipped=${iFolder}/F210I_embryonic_brain_se_intron_skipped.bed
F210I_included=${iFolder}/F210I_embryonic_brain_se_intron_included.bed
M323K_skipped=${iFolder}/M323K_adult_brain_se_intron_skipped.bed
M323K_included=${iFolder}/M323K_adult_brain_se_intron_included.bed
controls=${iFolder}/F210I_embryonic_brain_notsig_se_intron_all.bed

# outFolder
outFolder="/Users/Jack/google_drive/TDP_paper/RNA_maps/results/"

# exon beds
F210I_skipped_exons=${iFolder}/noheader/F210I_embryonic_brain_flank0_se_skipped.bed
F210I_included_exons=${iFolder}/noheader/F210I_embryonic_brain_flank0_se_included.bed
M323K_skipped_exons=${iFolder}/noheader/M323K_adult_brain_flank0_se_skipped.bed
M323K_included_exons=${iFolder}/noheader/M323K_adult_brain_flank0_se_included.bed
control_exons=${iFolder}/noheader/F210I_embryonic_brain_flank0_se_notsig_se_all.bed




# # F210I skipped
# python whole_intron_cluster_coverage.py $F210I_skipped $clusters ${outFolder}/F210I_skipped 100

# # F210I included 
# python whole_intron_cluster_coverage.py $F210I_included $clusters ${outFolder}/F210I_included 100

# # M323K skipped
# python whole_intron_cluster_coverage.py $M323K_skipped $clusters ${outFolder}/M323K_skipped 100

# # M323K included
# python whole_intron_cluster_coverage.py $M323K_included $clusters ${outFolder}/M323K_included 100

# # controls
# python whole_intron_cluster_coverage.py $controls $clusters ${outFolder}/controls 100

# cryptic exons

# skiptic exons

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



