# Aggregating all the functional analyses together
library(ggplot2)
library(stringr)
library(data.table)
library(dplyr)
library(gridExtra)

iFolder <- "/Users/Jack/google_drive/TDP_paper/"
  
biomartAnnotations <- paste0(iFolder, "/differential_splicing/biomart_annotations_mouse.tab")
biomart <- as.data.frame(fread(biomartAnnotations))

# do separately for F210I and M323K
# F210I
skipped_exon_list <- paste0(iFolder,"Sgseq/F210I_embryonic_brain/F210I_embryonic_brain_flank0_se_skipped.bed")
included_exon_list <- paste0(iFolder,"Sgseq/F210I_embryonic_brain/F210I_embryonic_brain_flank0_se_included.bed")
skipped_exon_annotation <- paste0(iFolder, "intron_annotation/F210I_skipped/cluster_annotations.tab")
included_exon_annotation <- paste0(iFolder, "intron_annotation/F210I_included/cluster_annotations.tab")
deseq_results <- paste0(iFolder,"/differential_expression/deseq_F210I_embryo_June_2016_differential_expression.tab")
cryptic_exon_list <- paste0(iFolder,"cryptics.bed")
skipped_exon_NMD <- "/Users/Jack/project/F210I_M323K/NMD_prediction/F210I/F210I_skipped_exons_individual_exon_results.tab"
included_exon_NMD <- "/Users/Jack/project/F210I_M323K/NMD_prediction/F210I/F210I_included_exons_individual_exon_results.tab"
code <- "F210I"

# M323K
skipped_exon_list <- paste0(iFolder,"Sgseq/M323K_adult_brain/M323K_adult_brain_flank0_se_skipped.bed")
included_exon_list <- paste0(iFolder,"Sgseq/M323K_adult_brain/M323K_adult_brain_flank0_se_included.bed")
skipped_exon_annotation <- paste0(iFolder, "intron_annotation/M323K_skipped/cluster_annotations.tab")
included_exon_annotation <- paste0(iFolder, "intron_annotation/M323K_included/cluster_annotations.tab")
deseq_results <- paste0(iFolder,"/differential_expression/deseq_M323K_adult_differential_expression.tab")
cryptic_exon_list <- paste0(iFolder,"skiptics.bed")
skipped_exon_NMD <- "/Users/Jack/project/F210I_M323K/NMD_prediction/M323K/M323K_skipped_exons_individual_exon_results.tab"
included_exon_NMD <- "/Users/Jack/project/F210I_M323K/NMD_prediction/M323K/M323K_included_exons_individual_exon_results.tab"
code <- "M323K"
# LOAD IN DATA

bringTogether <- function(skipped_exon_list,
                          included_exon_list,
                          skipped_exon_annotation,
                          included_exon_annotation,
                          deseq_results,
                          cryptic_exon_list,
                          skipped_exon_NMD,
                          included_exon_NMD,
                          code){

skipped_exons <- read.table(skipped_exon_list,header=TRUE, stringsAsFactors = FALSE)
included_exons <- read.table(included_exon_list, header=TRUE, stringsAsFactors = FALSE)
skipped_annotation <- read.table(skipped_exon_annotation, header= TRUE, sep = "\t", stringsAsFactors = FALSE)
included_annotation <- read.table(included_exon_annotation, header=TRUE, sep = "\t", stringsAsFactors = FALSE)
deseq <- read.table(deseq_results, header=TRUE, stringsAsFactors = FALSE, fill = TRUE)
cryptics <- read.table(cryptic_exon_list, header=FALSE, stringsAsFactors = FALSE)
names(cryptics) <- c("chr", "start", "end", "geneName", stringsAsFactors = FALSE)
skipped_NMD <- read.table(skipped_exon_NMD,header=TRUE, stringsAsFactors = FALSE)
included_NMD <- read.table(included_exon_NMD, header=TRUE, stringsAsFactors = FALSE)

# create simplified NMD prediction classes
group.names <- c(
  "Transcript switch is neutral",
  "Inclusion destabilises transcript",
  "Skipping destabilises transcript",
  "Inclusion destablises transcript (possibly non-NMD)",
  "exon in UTR",
  "exon evades classification"
)
simple.names <- c(
  "Neutral",
  "Inclusion toxic",
  "Skipping toxic",
  "Inclusion toxic",
  "UTR",
  "NA"
)
simplify_df <- data.frame(
  group.names, simple.names)


# create big table of all the data for skipped and included
mergeAllData <- function( exon_list, annotation, deseqRes, cryptic_list, NMD_list, code){
  d <- exon_list
  # merge extreme
  d$is.extreme <- match( 
    paste0( d$start,d$end, d$geneName), 
    paste0(cryptic_list$start, cryptic_list$end, cryptic_list$geneName) 
  )
  d$is.extreme <- ifelse( !is.na(d$is.extreme), TRUE, FALSE)
  # if( !all(is.na(d$is.extreme))){
  # d[ !is.na( d$is.extreme ), ]$is.extreme <- TRUE
  # }
  # merge annotation
  d$annotation <- annotation$class[ match( d$geneName, annotation$geneID)]
  # merge deseq
  if(code == "M323K"){
    d$geneLog2FC <- deseq$log2FoldChange[ match( d$geneName, deseq$external_gene_id)]
    # take adjusted P
    d$genePvalue <- deseq$padj[ match( d$geneName, deseq$external_gene_id)]
  }else{
    d$geneLog2FC <- deseq$log2FoldChange[ match( d$geneName, deseq$EnsemblID)]
    # adjusted P
    d$genePvalue <- deseq$padj[ match( d$geneName, deseq$EnsemblID)]
  }
  # merge NMD
  d$NMD.prediction <- NMD_list$verdict[ match(
    paste0( d$start,d$end, d$geneName), 
    paste0(NMD_list$central.start, NMD_list$central.end, NMD_list$ensemblID) 
  ) ]
  d$group <- NMD_list$group[ match(
    paste0( d$start,d$end, d$geneName), 
    paste0(NMD_list$central.start, NMD_list$central.end, NMD_list$ensemblID) 
  ) ]
  # UTR prediction
  utr_prediction <- NMD_list$UTR.verdict[ match(
    paste0( d$start,d$end, d$geneName), 
    paste0(NMD_list$central.start, NMD_list$central.end, NMD_list$ensemblID) 
  ) ]
  not_in_CDS <- which( d$NMD.prediction == "Exon not in CDS" ) 
  
  d$NMD.prediction[ not_in_CDS ] <- utr_prediction[ not_in_CDS ] 
  
  #d$NMD.prediction <- ifelse( is.na(d$NMD_prediction))
  d$group[ not_in_CDS ] <- utr_prediction[ not_in_CDS ] 
  #d[ d$NMD.prediction == "Exon not in CDS", ]$group <- "Exon not in CDS"
  d$NMD_prediction_simplified <- simplify_df$simple.names[match(d$group, simplify_df$group.names)]
  d$group <- NULL
  
  #change gene names from EnsemblIDs to gene symbols if F210I
  if( code == "F210I"){
    d$geneName <- biomart$external_gene_name[ match(d$geneName, biomart$EnsemblID)]
  }
  
  return(d)
}

included <- mergeAllData( included_exons, included_annotation, deseq, cryptics, included_NMD, code = code )
skipped <- mergeAllData( skipped_exons, skipped_annotation, deseq, cryptics, skipped_NMD, code = code )
write.table(included, paste0(iFolder, "/NMD_prediction/",code,"_all_included_exons.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE , sep = "\t")
write.table(skipped, paste0(iFolder, "/NMD_prediction/",code,"_all_skipped_exons.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE , sep = "\t")

return(list(included,skipped) )
}

# F210I
F210I_exons <- bringTogether( skipped_exon_list = paste0(iFolder,"Sgseq/F210I_embryonic_brain/F210I_embryonic_brain_flank0_se_skipped.bed"),
                              included_exon_list = paste0(iFolder,"Sgseq/F210I_embryonic_brain/F210I_embryonic_brain_flank0_se_included.bed"),
                              skipped_exon_annotation = paste0(iFolder, "intron_annotation/F210I_skipped/cluster_annotations.tab"),
                              included_exon_annotation = paste0(iFolder, "intron_annotation/F210I_included/cluster_annotations.tab"),
                              deseq_results = paste0(iFolder,"/differential_expression/deseq_F210I_embryo_June_2016_differential_expression.tab"),
                              cryptic_exon_list = paste0(iFolder,"cryptics.bed"),
                              skipped_exon_NMD = "/Users/Jack/project/F210I_M323K/NMD_prediction/F210I/F210I_skipped_exons_individual_exon_results.tab",
                              included_exon_NMD = "/Users/Jack/project/F210I_M323K/NMD_prediction/F210I/F210I_included_exons_individual_exon_results.tab",
                              code = "F210I"
)

# M323K
M323K_exons <- bringTogether(skipped_exon_list = paste0(iFolder,"Sgseq/M323K_adult_brain/M323K_adult_brain_flank0_se_skipped.bed"),
                             included_exon_list = paste0(iFolder,"Sgseq/M323K_adult_brain/M323K_adult_brain_flank0_se_included.bed"),
                             skipped_exon_annotation = paste0(iFolder, "intron_annotation/M323K_skipped/cluster_annotations.tab"),
                             included_exon_annotation = paste0(iFolder, "intron_annotation/M323K_included/cluster_annotations.tab"),
                             deseq_results = paste0(iFolder,"/differential_expression/deseq_M323K_adult_differential_expression.tab"),
                             cryptic_exon_list = paste0(iFolder,"skiptics.bed"),
                             skipped_exon_NMD = "/Users/Jack/project/F210I_M323K/NMD_prediction/M323K/M323K_skipped_exons_individual_exon_results.tab",
                             included_exon_NMD = "/Users/Jack/project/F210I_M323K/NMD_prediction/M323K/M323K_included_exons_individual_exon_results.tab",
                             code = "M323K"
)

# already made tables, read back in 
#F210I
code <- "F210I"
included <- as.data.frame(fread(paste0(iFolder, "/NMD_prediction/",code,"_all_included_exons.txt") ))
code <- "M323K"
skipped <- as.data.frame(fread( paste0(iFolder, "/NMD_prediction/",code,"_all_skipped_exons.txt") ))



#PLOTS

makeScatterPlots <- function(mergedData, colourBy, shapeBy = NULL){
  p <- ggplot(mergedData, aes_string(x = "geneLog2FC", y = "-log10(genePvalue)", colour = colourBy, shape = shapeBy )) + 
    geom_point(size = 2.5) + 
    geom_hline(yintercept = -log10(0.05), linetype = 3)
  return(p)
}

gridPlots <- function(mergedData, title, outFile){
  p1 <- makeScatterPlots(mergedData, "group", "is.extreme")
  p2 <- makeScatterPlots(mergedData, "NMD.prediction")
  p3 <- makeScatterPlots(mergedData, "annotation")
  #grid.arrange(p1,p3,p2, top = title, layout_matrix = rbind( c(1,2), 4) )
  g <- arrangeGrob(p1,p3,p2, top = title, layout_matrix = rbind( c(1,2), 4) )
  ggsave(g, filename =  outFile, width = 20, height = 10 , units = "in" )

}

makeBarPlots <- function( mergedData, drop = NULL,colourBy, fillLabel, title){
  #print(colourBy)
  control <- filter(mergedData, genePvalue > 0.1 | is.na(genePvalue))
  up <- filter( mergedData, geneLog2FC > 0 & genePvalue < 0.1 & colourBy != "exon evades classification" )
  down <- filter( mergedData, geneLog2FC < 0 & genePvalue < 0.1 & colourBy != "exon evades classification" )
  control_group <- group_by_(control, colourBy)
  print(nrow(mergedData))
  print(nrow(up))
  print(nrow(down))
  print(nrow(control_group))
  up_group <- group_by_(up, colourBy)
  down_group <- group_by_(down, colourBy)
  
  control_sum <- summarise(control_group, count = n() / nrow(control) ) %>% mutate(direction = paste0("None\n(", nrow(control), ")") )
  up_sum <- summarise(up_group, count = n() / nrow(up)) %>% mutate( direction = paste0("Upregulated\n(", nrow(up), ")" ) )
  down_sum <- summarise(down_group, count = n() / nrow(down) ) %>% mutate( direction = paste0("Downregulated\n(",nrow(down), ")" ) )

  toPlot <- rbind(control_sum,up_sum,down_sum)
  print(toPlot)
  if(!is.null(drop)){
    toPlot <- rbind(control_sum, down_sum)
  }
  x_limits <- c( paste0("None\n(", nrow(control), ")"),
  paste0("Upregulated\n(", nrow(up), ")" ),
  paste0("Downregulated\n(",nrow(down), ")" ) 
  )
  if (!is.null(drop)){
    x_limits <- c( paste0("None\n(", nrow(control), ")"),
                   paste0("Downregulated\n(",nrow(down), ")" ) 
    )
  }
  

  toPlot$dataset <- title
  p <- ggplot( toPlot, aes_string(x = "direction", y = "count", fill = colourBy ) ) +
    geom_col() +
    scale_fill_discrete(fillLabel) + 
    xlab("Gene direction") + 
    ggtitle(title) +
    theme_bw() +
    scale_x_discrete( limits = x_limits
      ) +
    scale_y_continuous("Proportion",labels=scales::percent) +
    theme_bw()
  return(p)
}

p1 <- makeBarPlots(included[ included$is.extreme, ],
             colourBy = "NMD_prediction_simplified", 
             fillLabel = "Functional\nprediction\n(exon)",
             drop = NULL,
             title = paste0("RRM2mut cryptic exons"))
p2 <- makeBarPlots(skipped[ skipped$is.extreme,],
             colourBy = "NMD_prediction_simplified",
             fillLabel = "Functional\nprediction\n(exon)",
             drop = "Upregulated",
             title = "LCDmut skiptic exons")
p3 <- makeBarPlots(included,
                   colourBy = "NMD_prediction_simplified", 
                   fillLabel = "Functional\nprediction\n(exon)",
                   title = paste0("RRM2mut all included exons"))
p4 <- makeBarPlots(skipped,
                   colourBy = "NMD_prediction_simplified",
                   fillLabel = "Prediction",
                   title = "LCDmut all skipped exons")

pdf("/Users/Jack/google_drive/TDP_paper/NMD_prediction/NMD_prediction_bar_plots.pdf")
grid.arrange(p1,p2,p3,p4)
dev.off()
#makeBarPlots(included, "NMD.prediction")
#makeBarPlots(included, "annotation")

# gridPlots( included, "M323K included exons", 
#            outFile = "/Users/Jack/google_drive/TDP_paper/NMD_prediction/M323K_included_exons_expression_scatters.pdf") 
# 
# gridPlots( skipped, "M323K skipped exons",
#            outFile = "/Users/Jack/google_drive/TDP_paper/NMD_prediction/M323K_skipped_exons_expression_scatters.pdf")
# 
# gridPlots( included, "F210I included exons", 
#            outFile = "/Users/Jack/google_drive/TDP_paper/NMD_prediction/F210I_included_exons_expression_scatters.pdf") 
#   
# gridPlots( skipped, "F210I skipped exons",
#            outFile = "/Users/Jack/google_drive/TDP_paper/NMD_prediction/F210I_skipped_exons_expression_scatters.pdf")
# 

