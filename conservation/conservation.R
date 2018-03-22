# Conservation scoring
# requires bigwigSummary from UCSC

library(data.table)
library(ggplot2)
library(dplyr)
library(optparse)
library(stringr)

options(echo = TRUE)

# check bigWigSummary is in the $PATH
if( Sys.which("bigWigSummary") == "" ){
	stop("bigWigSummary is not on your PATH")
}


#############################
# just the cryptic (F210I) and skiptics (M323K)

bed.files.list <- c(
	"../data/conservation/all_mouse_exons_merged_1000_random.bed",
	"../data/cryptics.bed",
	"../data/skiptics.bed"
)
bed.names <- c(
  "All annotated exons",
  "F210I extreme",
  "M323K extreme"
  )

species <- "mouse"
code <- "extreme_splicing"
outFolder <- "results/"
graph_title <- "F210I and M323K extreme splicing \nexon conservation (PhyloP 60 Way)"


conservation.outFolder <- outFolder

#conservation.outFolder <- paste0(outFolder,"/conservation/")


# create outFolder

if (! file.exists(conservation.outFolder)) dir.create(conservation.outFolder)



# make sure phyloP bigwig exists - if not then download it!

phyloP.bw <- "../mm10.60way.phyloP60way.bw"

if( ! file.exists(phyloP.bw) ){
	message("downloading phyloP bigwig not found! download it to root folder of repository")
	message( "wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phyloP60way/mm10.60way.phyloP60way.bw " )
}



message("PhyloP conservation!")
# use multiple PhyloP builds

exon.conservation <- list()
for(i in 1:length(bed.names)){
	bed.out <- bed.files.list[i]
	conservation.tmp.out <- paste0(conservation.outFolder,"/temp.txt")
  	cmd <- paste0(
		"cat ",
		bed.out,
		" | while read chr start end name exonID strand; do ",
		"bigWigSummary ",
		phyloP.bw, 
		" $chr $start $end 1; done  2>&1 | awk \' $0 ~ \"data\"{print \"NA\"}$0 !~ \"data\" {print $0}\' | sed \'/^$/d\'",
		" > ", 
		conservation.tmp.out
	)
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
  res$pvalue[i] <- t.test( control, case)$p.value
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

