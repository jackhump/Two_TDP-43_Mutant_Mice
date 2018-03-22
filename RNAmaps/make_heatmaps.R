library(ggplot2)
library(stringr)
library(optparse)
library(reshape2)

options(echo=TRUE)

opt <- parse_args(
  OptionParser(option_list=list(
    make_option( "--data", type="character", default=NULL, "Rdata file produced by coverage script"),
    make_option( "--code", type="character"),
    make_option( "--mode", type="character"),
    make_option( "--outFolder",type="character", default=NULL, help = "where you want the plot to go"),
    make_option( "--input", type="character", default=NULL, help = "what clusters or motifs"),
    make_option( "--annotations", type="character", "biomart annotations"))
))

data <- opt$data
outFolder <- opt$outFolder
code <- opt$code
mode <- opt$mode
input <- opt$input
biomart <- opt$annotations

n <- 20
# load in Rdata from iCLIP and UG RNA maps
annotations <- read.table(biomart,header=TRUE, stringsAsFactors = FALSE)
# remember - the two coverage data are not in matching row order - need to use exon list to match correctly.
# BUT! there are NA values and duplicates in the exon list - how stupid!

# make heatmaps of the most decorated genes

# for testing!
# mat <- matrix_list[[1]]
# exons <- exon_list[[1]]
# matching <- clip$genes
# n <- 20
# counts <- count_exon_list[[1]]

prepareHeatMap <- function(mat, exons, counts, n, matching=NULL, centre_point ){
    # exon list is in same order as coverage matrix so a straightforward merge is fine
    row.names(exons) <- 1:nrow(exons)
    # remove duplicates
    duplicates <- duplicated(exons)
    exons <- exons[!duplicates,]
    mat <- mat[!duplicates,]
    

    counts <- counts[!duplicates]
    counts <- data.frame(x = counts)
    
    # replace EnsemblIDs with gene symbols
    if( any(grepl("ENS.*G.*0",exons$V4)) ){
      gene_symbol <- annotations$external_gene_name[ match(exons$V4, annotations$EnsemblID)]
      # fix compound names
      if( nrow( exons[is.na(gene_symbol),] ) > 0 ){
        fixed_names <- apply( exons[is.na(gene_symbol),], MAR = 1, FUN = function(x){
          genes <- str_split( x[4], "\\+", 2)[[1]]
          symbols <- annotations$external_gene_name[ match(genes, annotations$EnsemblID)]
          symbols <- paste(symbols, collapse="+")
          return(symbols)
        })
        gene_symbol[ is.na(gene_symbol)] <- as.character(fixed_names)
      }
      exons$V4 <- gene_symbol
    }
    
  # either find the top genes by CLIP cluster number or match a set of genes
  if( is.null(matching)){  
    
    #exons$V4[ match( matching, exons$V4 ) ]
    
    # find the top genes and keep their exons
    mat <- mat[ order(rowSums(mat), decreasing = TRUE), ]
    mat <- mat[1:n,] # pick top
    
    counts <- counts$x[ match( row.names(mat), row.names(counts))]
    
    genes <- exons$V4[ match(row.names(mat), row.names(exons))]
    row.names(mat) <- genes
    
    counts <- data.frame( gene = genes, N = counts)
  
    
  }else{
    # get the rows of the matrix matching the row numbers of the list of genes to match
    mat <- mat[ match( matching, exons$V4 ) , ]
    x <- counts$x[ match( row.names(mat), row.names(counts))]
    
    row.names(mat) <- matching
    counts <- data.frame(gene = matching, N = x)
  }
  #pvalues <- signif(exons$V5[ match(row.names(mat), row.names(exons))],digits=2)
  names(mat) <- 1:ncol(mat)

  mat.melt <- melt(as.matrix(mat))
  small <- subset(mat.melt, value > 0)
  small$Var1 <- factor(small$Var1, levels = rev(row.names(mat)) )
  #print(exons)
  print(counts)
  counts$x <- centre_point
  counts <- subset(counts, N > 0)
  
  return(list(melted=small, genes = genes, counts = counts) )
}


makePlot <- function(data, motifs = NULL, colour, direction){
  title <- paste(code, mode, input, direction)
  p <- ggplot(data$melted, aes(x = Var2, y = Var1 ) ) +
    geom_tile(fill = colour, colour = colour, alpha = 0.2 ) +
    scale_y_discrete("")
  if( !is.null(motifs) ){
    p <- p + 
      geom_tile(data=motifs$melted, aes(x = Var2, y = Var1 ), fill = "gray", colour = "gray", alpha = 1  ) +
      geom_point(data=motifs$counts, aes(x = x, y = gene ), colour = "gray", alpha = 1)
  }
  p <- p +
    
    scale_x_continuous(
      "",
      breaks = x_breaks,
      label = x_labels,
      limits = c(0,total_length+1),
      expand = c(0, 0)
    ) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = NA),
          panel.grid.major.y = element_line(colour = "gray")
          ) +
    geom_vline( data = data.frame(xintercept = x_breaks[c(2:6,8:12)], slope=0),
                aes(xintercept=xintercept),colour = "black", linetype="dashed",alpha=0.5 ) +
    geom_polygon(data=blank_df, aes(x,y, group = group), fill = "white", colour = "NA") +
    geom_point(data=data$counts, aes(x = x, y = gene ), colour = colour, alpha = 0.8) +
    ggtitle(title)
return(p)
}


#load("../iCLIP/results/F210I_coverage_data.Rdata")
#data <- "/Users/Jack/google_drive/TDP_paper/RNA_maps/iCLIP/results//M323K_iCLIP_F210I_WT_coverage_data.Rdata"
load(data)
print(centre_point)
included_iCLIP <- prepareHeatMap(matrix_list[[1]], exons = exon_list[[1]], n, counts = count_exon_list[[1]], centre_point = centre_point)
skipped_iCLIP <- prepareHeatMap(matrix_list[[2]], exons = exon_list[[2]], n, counts = count_exon_list[[2]], centre_point = centre_point)
#load("TGT/F210I_TGT_coverage_data.Rdata")
#included_motifs <- prepareHeatMap(matrix_list[[1]], exons = exon_list[[1]], matching = clip$genes, n = n, counts = count_exon_list[[1]], centre_point = centre_point)
#skipped_motifs <- prepareHeatMap(matrix_list[[2]], exons = exon_list[[2]], matching = clip$genes, n = n, counts = count_exon_list[[2]], centre_point = centre_point)
# fudge - x_breaks are loaded in with the data
blank_df <- data.frame(
  x = c(
    x_breaks[3],
    x_breaks[3],
    x_breaks[4],
    x_breaks[4],
    
    # x_breaks[6],
    # x_breaks[6],
    # x_breaks[8],
    # x_breaks[8],
    
    x_breaks[10],
    x_breaks[10],
    x_breaks[11],
    x_breaks[11]
  ),
  y = rep(c(n+1,0,0,n+1),2),
  group = c(
    1,1,1,1,
    # 2,2,2,2,
    3,3,3,3
  )
)


p1 <- makePlot(included_iCLIP, "red", motifs = NULL, direction = "included")
p2 <- makePlot(skipped_iCLIP, "blue", motifs = NULL, direction = "skipped")
#p3 <- makePlot(included_iCLIP, "red", motifs = included_motifs)
#p4 <- makePlot(skipped_iCLIP, "blue", motifs = skipped_motifs)

fileName <- paste0(outFolder, "/", gsub(" ", "_", code),"_", mode, "_", input,"_heatmap.Rdata")

pdf(fileName,width = 12,height=4)
#grid.arrange(p1,p2, nrow=2)
print(p1)
print(p2)
#print(p3)
#print(p4)
dev.off()


