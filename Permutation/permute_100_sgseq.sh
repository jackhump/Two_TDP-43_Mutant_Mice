# 100 permutations of SGSeq step2b
support="support_files/m323k_2018_support.tab"

permuteR="permute_support.R"

R=/share/apps/R-3.3.2/bin/R

export LD_LIBRARY_PATH=/share/apps/zlib-1.2.8/lib:/share/apps/zlib-1.2.8/lib:/opt/gridengine/lib/linux-x64:/opt/gridengine/lib/linux-x64:/opt/openmpi/lib:/opt/python/lib


for i in `seq 1 100`; do

	echo "iteration: $i" >> permutation_progress.txt

	iteration=$i
	outFolder=results/permute_${i}
	permuted_support=${support}.permute.${i}.tab

	# permute sample order
	Rscript $permuteR $support $i

	if [ ! -e $outFolder ];then
		mkdir -p $outFolder
	fi

	# get around the need for Rdata 
	realDataFile="/SAN/vyplab/TDP43_RNA/TDP_F210I_M323K/M323K/New_adult_brain/processed/sgseq/2018_4_HET_4_HOM/M323K_New_sgv_novel.RData"
	fakeDataFile=${outFolder}/M323K_permute_${i}_sgv_novel.RData
	# create symlink
	ln -s $realDataFile $fakeDataFile

	# do step2 

	${R}script --vanilla \
		/SAN/vyplab/HuRNASeq/RNASeq_pipeline/SGSeq//sgseq_step2.R \
		--step step2b \
		--support.tab ${permuted_support} \
		 --code M323K_permute_${i} \
		 --output.dir $outFolder \
		 --annotation /SAN/vyplab/HuRNASeq/reference_datasets/RNASeq//Mouse/biomart_annotations_mouse.tab >> R_output.txt



done
