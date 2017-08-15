# long genes in F210I
library(data.table)
library(ggplot2) #ggplot package was installed
library(dplyr) #dpplyr package was installed
library(gridExtra) #gridExtra was installed
library(stringr)
library(tidyr)

options(echo=TRUE)
# mydata.vector <- c( "../../Downloads/deseq_Fratta_F210I_adult_2016_differential_expression.tab",
#                	"../../Downloads/deseq_Fratta_F210I_embryo_differential_expression.tab",
#                	"../../Downloads/deseq_Cleveland_TDP_differential_expression.tab",
#                	"../../Downloads/deseq_Fratta_F210I_embryo_June_2016_long_differential_expression.tab") #The pathways to four different datasets were created
mydata.vector <- c(

 # F210I embryo brain
"/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/deseq2/F210I_CTL_long_F210I_HOM_long/deseq_F210I_embryo_June_2016_differential_expression.tab",
"/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/deseq2/F210I_CTL_norm_F210I_HOM_norm/deseq_F210I_embryo_June_2016_differential_expression.tab",

  # F210I adult brain
"/SAN/vyplab/IoN_RNAseq/F210I/New_adult_brain/F210I_adult/processed/deseq2/CONTROL_HET/deseq_F210I_adult_differential_expression.tab",

 # F210I embryo head
"/SAN/vyplab/IoN_RNAseq/F210I/Old_embryo_brain/processed/deseq2/Ctl_Het/deseq_F210I_emb_differential_expression.tab",
"/SAN/vyplab/IoN_RNAseq/F210I/Old_embryo_brain/processed/deseq2/Ctl_Hom/deseq_F210I_emb_differential_expression.tab",

  # MEFs
"/SAN/vyplab/IoN_RNAseq/F210I/Embryonic_fibroblasts/processed/deseq2/F210I_CTL_F210I_HET/deseq_F210I_MEF_differential_expression.tab",
"/SAN/vyplab/IoN_RNAseq/F210I/Embryonic_fibroblasts/processed/deseq2/F210I_CTL_F210I_HOM/deseq_F210I_MEF_differential_expression.tab",


	# M323K adult
"/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/deseq2/CONTROL_HET/deseq_M323K_adult_differential_expression.tab",
"/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/deseq2/CONTROL_HOM/deseq_M323K_adult_differential_expression.tab",

 # M323K embryo head

"/SAN/vyplab/IoN_RNAseq/M323K/Old_embryonic_head/processed/deseq2/Ctl_Het/deseq_M323K_emb_differential_expression_keep_dups.tab",
"/SAN/vyplab/IoN_RNAseq/M323K/Old_embryonic_head/processed/deseq2/Ctl_Hom/deseq_M323K_emb_differential_expression_keep_dups.tab",

    # M323K MEFs
"/SAN/vyplab/IoN_RNAseq/F210I/Embryonic_fibroblasts/processed/deseq2/M323K_CTL_M323K_HET/deseq_F210I_MEF_differential_expression.tab",

"/SAN/vyplab/IoN_RNAseq/F210I/Embryonic_fibroblasts/processed/deseq2/M323K_CTL_M323K_HOM/deseq_F210I_MEF_differential_expression.tab",

  # Cleveland TDP
  "/cluster/scratch3/vyp-scratch2/Humphrey_datasets/Cleveland/TDP43/processed/deseq2/CTL_TDP_KD/deseq_Cleveland_TDP_differential_expression.tab",
  # Cleveland FUS
  "/cluster/scratch3/vyp-scratch2/Humphrey_datasets/Cleveland/FUS/processed/deseq2/CTL_FUS/deseq_Cleveland_FUS_differential_expression.tab",
 
 # ENCODE 
  "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/expression/deseq2/control_TDP/deseq_dataset_1_differential_expression.tab",
  
  "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/expression/deseq2/control_TDP/deseq_dataset_2_differential_expression.tab"

 )

 


file.exists(mydata.vector) #Checking the file exists
titles <- c(
    # F210I embryonic brain
    "Mouse embryonic brain TDP-43 F210I HOM\n(normal inserts)",
  	"Mouse embryonic brain TDP-43 F210I HOM\n(long inserts) ",
  	# adult
    "Mouse adult brain\nTDP-43 F210I HET",
  	# head
    "Mouse embryo whole head\nTDP-43 F210I HET",
    "Mouse embryo whole head\nTDP-43 F210I HOM",
  	 # MEFs
    "Mouse embryonic fibroblasts\nTDP-43 F210I HET",
    "Mouse embryonic fibroblasts\nTDP-43 F210I HOM",

    # M323K adult brain
    "Mouse adult brain\nTDP-43 M323K HET",
    "Mouse adult brain\nTDP-43 M323K HOM",
    # head
    "Mouse embryonic head\nTDP-43 M323K HET",
    "Mouse embryonic head\nTDP-43 M323K HOM",
    # MEFs
    "Mouse embryonic fibroblasts\nTDP-43 M323K HET",
    "Mouse embryonic fibroblasts\nTDP-43 M323K HOM",
    
    # Cleveland
    "Mouse adult striatum\nTDP-43 knockdown",
  	"Mouse adult striatum\nFUS knockdown",

    # ENCODE
  	"Human K562 cell mRNA\nTDP-43 knockdown",
  	"Human K562 cell total RNA\nTDP-43 knockdown"
        	) #A vector of titles was created

plots <- list() #An empty list was created so that each plot can be kept after being run
res.list <- list()
i_lengths <- list(
	Mouse = fread("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/long_genes/gencode_mouse_intron_lengths.txt"),
	Human = fread("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/long_genes/gencode_human_intron_lengths.txt")
)


print(length(titles) )
print(length(mydata.vector) )

print(titles)

 # The amount of times the iteration was to be performed was stated, and the length function was used so more data can be added in future without re-writing this section of the script, curly brackets were used to state what needs to be re-iterated with each dataset
for(i in 1:length(mydata.vector)){
	print( mydata.vector[i])
	species <- str_split_fixed(titles[i], " ",2)[,1]
	
	intron_lengths <- i_lengths[[species]]
	
	head(intron_lengths) 
	d <- as.data.frame( fread( mydata.vector[i], header=T) ) # The FUS dataset was reloaded
	 
	signed.p <- ifelse(test = d$log2FoldChange > 0, 
	  					 yes = -log10(d$pvalue), 
	  					 no = -1*(-log10(d$pvalue))) #A logical condition was tested
	 
	signed.z <- ifelse(test = d$log2FoldChange > 0,
	  					yes = qnorm(1 - (d$pvalue / 2)),
	  					no = qnorm(d$pvalue / 2) )	 
	d <- mutate(d, signed.p = signed.p, signed.z = signed.z) # A new column of signed (-log10)p values was created
	  
	d$intron.length <- intron_lengths$V4[ match( d$external_gene_id, intron_lengths$V1 ) ]

	#  d <- mutate(d, gene.length = end_position - start_position) #A new column, gene.length, was created using pre-existing data
	  
	d <- filter(d, baseMean > 0.1 & intron.length > 1000) #All genes with a baseMean >1 remained
	d <- filter(d, !is.na(d$intron.length)) # Any NA values were removed, and all values that weren't NA were kept

	bin.size <- 200

	bin <- c()
	for ( b in 1:ceiling( nrow(d) / bin.size) ){
	  bin <- c( bin, rep(b, bin.size) )
	}
	bin <- bin[1:nrow(d)]
	d <- arrange(d, signed.z)
	ranked.bins.z <- bin
	d <- cbind(d, ranked.bins.z)
	#d <- mutate(d, ranked.bins.p = ntile(signed.p, 50), ranked.bins.z = ntile(signed.z, 50)) #A new column 'ranked.bin' was created with 50 equally sized bins

	res.list[[i]] <- d

	d.mean <- group_by(d, ranked.bins.z) #The d.mean column was grouped according to ranked bins

	d.mean <- summarise(d.mean, mean.intron.length = mean(intron.length), sem.intron.length = sd(intron.length) / sqrt(bin.size)  )

	print(as.data.frame(d.mean))
	
	print(titles[i])
	
	p <- ggplot(d.mean, aes(x=ranked.bins.z , y= mean.intron.length/1000)) +
	  geom_line(colour = "red", size = 0.75) +
	  geom_ribbon( aes( ymin =  (mean.intron.length - sem.intron.length) / 1000,
			    ymax =  (mean.intron.length + sem.intron.length  ) /1000  
				), 
			fill = "red", alpha = 0.3 ) +
	  ggtitle( titles[i] ) +
	  xlab("") +
	  xlab("Downregulated                   Not regulated              Upregulated") +
	  ylab("Mean intron length (kb)") +
	  scale_x_continuous(breaks=NULL) +
	  ylim(c(0,250)) +
	  theme_bw(base_size = 8)


	 
	  plots[[i]] <- p #Setting the ith value of list of plots to p
 
}


pdf("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/long_genes/all_mice_TDP-43_long_gene_plot_Z.pdf") #A png file was created
grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 2)
grid.arrange(plots[[4]], plots[[5]], plots[[6]], plots[[7]]) 
grid.arrange(plots[[8]], plots[[9]], plots[[10]], plots[[11]])
grid.arrange(plots[[12]],plots[[13]],plots[[14]], plots[[15]])
grid.arrange(plots[[16]],plots[[17]], ncol = 2, nrow = 2)
dev.off() #The device was turned off

warnings()

quit()

# create a list of mouse long genes from the first DESeq results
lengths <- read.table(mydata.vector[1], header=T) %>%
			mutate(gene.length = end_position - start_position) %>%
			select(EnsemblID, external_gene_id, gene.length) %>%
			filter(gene.length >= 1E5)

for(i in 1:length(plots)){
	dataset <- res.list[[i]]
	lengths$dataset.log2FC <- dataset$log2FoldChange[match(lengths$EnsemblID,dataset$EnsemblID)]
	names(lengths)[i + 3] <- titles[i]
}

x <- gather(lengths, key = dataset, value = log2FoldChange, -EnsemblID, -external_gene_id, -gene.length) %>%
	 arrange(desc(gene.length))

# just take the genes < 100kb but > 1kb
mid.lengths <- read.table(mydata.vector[1], header=T) %>%
			mutate(gene.length = end_position - start_position) %>%
			select(EnsemblID, external_gene_id, gene.length) %>%
			filter(gene.length < 1E5 & gene.length > 1E3)

for(i in 1:6){
	dataset <- res.list[[i]]
	mid.lengths$dataset.Z <- dataset$signed.z[match(mid.lengths$EnsemblID,dataset$EnsemblID)]
	names(mid.lengths)[i + 3] <- titles[i]
}

long.gene.plot <- ggplot(x, aes(x = dataset, y = external_gene_id) ) + 
				  geom_tile(aes(fill = log2FoldChange), colour = "white") + 
				  scale_fill_distiller(palette = "Spectral")
				  #scale_fill_gradient(low = "red", high = "steelblue")

ggsave("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/long_genes/F210I_compound_long_gene_heatmap.pdf", height = 100, width = 7, units = "in", limitsize = F)
