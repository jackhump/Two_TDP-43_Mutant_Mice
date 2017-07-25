# Conservation scoring
# requires bigwigSummary from UCSC

library(data.table)
library(ggplot2)
library(dplyr)
library(optparse)
library(stringr)

options(echo = TRUE)

bigWigSummary <- "/SAN/vyplab/HuRNASeq/UCSC_tools/bigWigSummary"

# compare the per-exon mean conservation score for a group of exons to the scores of all protein-coding exons in the genome.

# skipped and included exons discovered by SGSeq

# to do this read in GENCODE mouse GTF 
# write out all exons of all protein coding transcript of all protein coding genes in BED format.
# do in Python!



# bed.files.list <- c(
#   "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/all_mouse_exons_merged.bed",
# "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/F210I_se_included.bed",
# "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/F210I_se_skipped.bed",
# "/SAN/vyplab/IoN_RNAseq/Kitty/M323K/sgseq/M323K_adult_brain_se_included.bed",
# "/SAN/vyplab/IoN_RNAseq/Kitty/M323K/sgseq/M323K_adult_brain_se_skipped.bed" 
#   )



# for(i in bed.files.list){
#   if( !file.exists(i) ){
#     stop( paste0(i, " doesn't exist"))
#   }

# }



# bed.names <- c("All annotated exons",
#             "F210I included",
#             "F210I skipped",
#             "M323K included",
#             "M323K skipped" )

# species <- "mouse"
# code <- "F210I_M323K"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K"
# graph_title <- "F210I and M323K exon conservation (PhyloP 60 Way)"



#############################
# just the cryptic (F210I) and skiptics (M323K)

bed.files.list <- c(
 "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/all_mouse_exons_merged_1000_random.bed",
 "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/cryptic_skiptic/cryptics_nohead.bed",
 "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/cryptic_skiptic/skiptics_nohead.bed"
  )

bed.names <- c(
  "All annotated exons",
  "F210I extreme",
  "M323K extreme"
  )

species <- "mouse"
code <- "extreme_splicing"
outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/cryptic_skiptic/"
graph_title <- "F210I and M323K extreme splicing \nexon conservation (PhyloP 60 Way)"




conservation.outFolder <- paste0(outFolder,"/conservation/")
if (! file.exists(conservation.outFolder)) dir.create(conservation.outFolder)



if(species == "mouse"){
  phyloP.bw <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/phyloP/mm10.60way.phyloP60way.bw"
}
if(species == "human"){
  phyloP.bw <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/phyloP/hg38.phyloP100way.bw"
}




# PhyloP conservation!
# use multiple PhyloP builds

exon.conservation <- list()
for(i in 1:length(bed.names)){
	bed.out <- bed.files.list[i]
	conservation.tmp.out <- paste0(conservation.outFolder,"/temp.txt")
  cmd <- paste0("cat ",bed.out," | while read chr start end name exonID strand; do ",bigWigSummary, " ", phyloP.bw, " $chr $start $end 1; done  2>&1 | awk \' $0 ~ \"data\"{print \"NA\"}$0 !~ \"data\" {print $0}\' | sed \'/^$/d\'", " > ", conservation.tmp.out)
	system(cmd)
  conservation <- fread( conservation.tmp.out )
  	names(conservation)[1] <- "phyloP.score"

	conservation$exon.type <- bed.names[i]

	conservation.out <- paste0(bed.names[i],".conservation")

  # read in bed file to annotate the conservation scores

  d <- as.data.frame(fread(bed.files.list[i]))

  conservation <- cbind( conservation,d)


	assign(conservation.out,conservation)
	exon.conservation[[i]] <- conservation
}


# prepare tables for graphing
exon.conservation.merge <- as.data.frame(do.call(args = c(exon.conservation, fill = TRUE) , what = rbind, ))

res <- group_by(exon.conservation.merge, exon.type) %>% 
        summarise( 
          n = n() , 
          mean = mean(phyloP.score, na.rm= TRUE), 
          sd = sd(phyloP.score, na.rm = TRUE), 
          sem = sd / sqrt(n) ) %>% 
        mutate(pretty = paste0(signif(mean,3),"±",signif(sem,3)) )

res$pvalue <- ""
for(i in 2:3){
  control <- subset( exon.conservation.merge,  exon.type == res$exon.type[1] )$phyloP.score
  case <-  subset( exon.conservation.merge,  exon.type == res$exon.type[i] )$phyloP.score
  res$p.value[i] <- t.test( control, case)$p.value
}



# save data - not the all exon null though!
conservation.table <- paste0(conservation.outFolder,code,"_exon_conservation_1000_null.tab")

exon.conservation.no.null <- exon.conservation.merge[ exon.conservation.merge$exon.type != bed.names[1],]

write.table(exon.conservation.merge, conservation.table, sep = "\t", row.names = FALSE, quote = FALSE)

#exon.conservation.merge <- fread(conservation.table)

conservation.graph <- paste0(conservation.outFolder,code,"_exon_conservation_plot.pdf")


###### get order of exons the right way round!
exon.conservation.merge$exon.type <- factor(exon.conservation.merge$exon.type, levels = bed.names)
exon.n <- as.numeric(table(exon.conservation.merge$exon.type))
bed.labels <- paste0(gsub(" ","\n",bed.names), "\n(",exon.n,")")

####### make the plot

pdf(conservation.graph)
p <- ggplot(exon.conservation.merge, aes(x = gsub(".","\n",exon.type,fixed =T), y = phyloP.score, fill = exon.type, colour = exon.type)) + 
  geom_violin(colour = "black", fill = "white") +
	geom_boxplot(colour = "black", fill = "skyblue", notch = T, width = 0.1 ) +
	xlab("") + 
	ylab("per exon mean phyloP score") + 
	ylim(c(-1,6)) +
	#scale_fill_manual(values = c("skyblue","firebrick3","skyblue","skyblue")) +
	scale_x_discrete(name=" ",
					limits=bed.names,
                    labels=bed.labels) +
	ggtitle(paste0(graph_title,"\naverage phyloP conservation score per exon")) + 
	theme_bw() + 
	theme(panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
          axis.line = element_line(colour = "black"),
          legend.position="none"
          )
print(p)
dev.off()

#exon.conservation.out <- paste0( conservation.outFolder ,"conservation_table.tab" )
#write.table( exon.conservation.merge, file = exon.conservation.out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# two-sample t test of each sample against the null
exon.types <- names(table(exon.conservation.merge$exon.type))

test <- lapply( exon.types[2:length(exon.types)],
  FUN = function(x) { 
    t.test( exon.conservation.merge[ exon.conservation.merge$exon.type == "All annotated exons",]$phyloP.score,
      exon.conservation.merge[ exon.conservation.merge$exon.type == x,]$phyloP.score )
  }
  )

#test <- t.test(exon.conservation.merge$phyloP.score, cryptic_nulls.conservation$phyloP.score)
test.results <- paste0(conservation.outFolder,code,"_t_test.tab")

capture.output(test, file = test.results )

