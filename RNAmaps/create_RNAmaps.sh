#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#----- set locations
iFolder=../data/RNAmaps
biomart=../data/biomart_annotations_mouse.tab.gz

# unzip biomart
if [ -e $biomart ]; then
	gunzip $biomart
fi

biomart=${biomart%.gz}


# ------ iCLIP clusters
# all TDP iCLIP:
# F210I WT and Mutant, M323K WT and Mutant, same published iCLIP as in cryptex TDP_paper

clusters=$iFolder/All_TDP_iCLIP_merged.bed.gz


# unzip clusters
if [ -e $clusters ];then
	gunzip $clusters
fi

clusters=${clusters%.gz}

clusterCode=All_TDP_iCLIP




# intron beds
F210I_skipped=${iFolder}/F210I_embryonic_brain_se_intron_skipped.bed
F210I_included=${iFolder}/F210I_embryonic_brain_se_intron_included.bed
M323K_skipped=${iFolder}/M323K_adult_brain_se_intron_skipped.bed
M323K_included=${iFolder}/M323K_adult_brain_se_intron_included.bed
controls=${iFolder}/F210I_embryonic_brain_notsig_se_intron_all.bed

# outFolder
outFolder=results/

# create outFolder

if [ ! -e $outFolder ];then
	mkdir $outFolder
fi

# exon beds
F210I_skipped_exons=${iFolder}/F210I_embryonic_brain_flank0_se_skipped.bed
F210I_included_exons=${iFolder}/F210I_embryonic_brain_flank0_se_included.bed
M323K_skipped_exons=${iFolder}/M323K_adult_brain_flank0_se_skipped.bed
M323K_included_exons=${iFolder}/M323K_adult_brain_flank0_se_included.bed
control_exons=${iFolder}/F210I_embryonic_brain_flank0_se_notsig_se_all.bed

# scripts
covScript=whole_intron_cluster_coverage.py
plotScript=plot_whole_intron_coverage.R
heatmapScript=make_heatmaps.R 


#---------- create coverage of iCLIP peaks for each bed file

# # # F210I skipped
python $covScript $F210I_skipped $clusters ${outFolder}/F210I_skipped_${clusterCode} 100 --stranded

# F210I included 
python $covScript $F210I_included $clusters ${outFolder}/F210I_included_${clusterCode} 100 --stranded

# M323K skipped
python $covScript $M323K_skipped $clusters ${outFolder}/M323K_skipped_${clusterCode} 100 --stranded

# M323K included
python $covScript $M323K_included $clusters ${outFolder}/M323K_included_${clusterCode} 100  --stranded

# controls
python $covScript $controls $clusters ${outFolder}/controls_${clusterCode} 100  --stranded

#------------ create RNA Map

# # create F210I RNA Map
Rscript $plotScript \
					 --skipped ${outFolder}/F210I_skipped_${clusterCode}_100_coverage.csv \
				     --included ${outFolder}/F210I_included_${clusterCode}_100_coverage.csv \
				     --control ${outFolder}/controls_${clusterCode}_100_coverage.csv \
				     --skipped_exons $F210I_skipped_exons \
				     --included_exons $F210I_included_exons \
				     --control_exons $control_exons \
				     --no_scaled_intron \
				     --no_scaled_exon \
				     --code F210I \
				     --outFolder $outFolder \
				     --mode iCLIP \
				     --input $clusterCode \
                                     --annotations $biomart


# M323K 
Rscript $plotScript --skipped ${outFolder}/M323K_skipped_${clusterCode}_100_coverage.csv \
                                     --included ${outFolder}/M323K_included_${clusterCode}_100_coverage.csv \
                                     --control ${outFolder}/controls_${clusterCode}_100_coverage.csv \
                                     --skipped_exons $M323K_skipped_exons \
                                     --included_exons $M323K_included_exons \
                                     --control_exons $control_exons \
                                     --code M323K \
                                     --outFolder $outFolder \
                                     --no_scaled_intron \
				     --no_scaled_exon \
                                     --mode iCLIP \
                                     --input $clusterCode \
                                     --annotations $biomart


#---------------- create heatmaps

Rscript $heatmapScript --data ${outFolder}/M323K_iCLIP_${clusterCode}_coverage_data.Rdata \
 					   --outFolder ${outFolder} \
 					   --annotations $biomart \
 					   --code F210I \
 					   --mode iCLIP \
 					   --input $clusterCode




