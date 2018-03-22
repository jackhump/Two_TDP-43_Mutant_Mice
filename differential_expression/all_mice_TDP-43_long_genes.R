# long genes in F210I
#Jack Humphrey and Shannon Edwards

library(data.table)
library(ggplot2) #ggplot package was installed
library(dplyr) #dpplyr package was installed
library(gridExtra) #gridExtra was installed
library(stringr)
library(tidyr)


iFolder <- "../data/"

F210Ires <- paste0(iFolder, "DESeq2/deseq_F210I_embryo_June_2016_differential_expression.tab")
M323Kres <- paste0(iFolder, "DESeq2/deseq_M323K_adult_differential_expression.tab" )

# generated from Polymenidou et al, 2011 - PRJNA141971 
Cleveland <- paste0(iFolder, "DESeq2/deseq_Cleveland_TDP_differential_expression.tab" )

# mouse intron lengths
mouse_intron_lengths <- "../data/DESeq2/gencode_mouse_intron_lengths.txt"

mydata.vector <- c( F210Ires, M323Kres, Cleveland )

file.exists(mydata.vector) #Checking the file exists

titles <- c(
    "Mouse embryonic brain\nTDP-43 F210I HOM",
    "Mouse adult brain\nTDP-43 M323K HOM",
    "Mouse adult striatum\nTDP-43 knockdown") 


plots <- list() #An empty list was created so that each plot can be kept after being run
res.list <- list()
i_lengths <- list(
	Mouse = fread( paste0( "zless ", mouse_intron_lengths ) )
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
	d <- as.data.frame( fread( mydata.vector[i], header=T) ) # The dataset was reloaded
	 
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


pdf("all_mice_TDP-43_long_gene_plot_Z.pdf") #A png file was created
grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 2)
dev.off() #The device was turned off

