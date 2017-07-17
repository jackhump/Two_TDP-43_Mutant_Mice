library(ggplot2)
library(stringr)
library(dplyr)
library(data.table)

iFolder <- "/Users/Jack/google_drive/TDP_paper/"
F210Ires <- paste0(iFolder, "differential_expression/deseq_F210I_embryo_June_2016_differential_expression.tab")
M323Kres <- paste0(iFolder, "differential_expression/deseq_M323K_adult_differential_expression.tab" )
crypticExons <- paste0(iFolder,"cryptics.bed")
skipticExons <- paste0(iFolder,"skiptics.bed")
F210I_included <- paste0(iFolder, "intron_annotation/", "F210I_embryonic_brain_flank0_se_included.bed")
F210I_skipped <- paste0(iFolder, "intron_annotation/", "F210I_embryonic_brain_flank0_se_skipped.bed" )
M323K_included <- paste0(iFolder, "intron_annotation/", "M323K_adult_brain_flank0_se_included.bed" )
M323K_skipped <-  paste0(iFolder, "intron_annotation/", "M323K_adult_brain_flank0_se_skipped.bed")

prepareAllAnnotatedComparison <- function(deseqRes, exonList, exonCode, genotypeCode){
  res <- as.data.frame(fread(deseqRes)) 
  cryptics <- as.data.frame(fread(exonList))
  if( genotypeCode == "F210I" ){
    res$is.cryptic <- cryptics$V4[ match(res$EnsemblID, cryptics$V4)]
  }
  else{
    res$is.cryptic <- cryptics$V4[ match(res$external_gene_id, cryptics$V4)]
  }
  crypticRes <- filter(res, !is.na(is.cryptic))
  meanExpr <- min(crypticRes$baseMean)
  clean <- filter(res, baseMean >= meanExpr )
  annotated <- clean[ is.na(clean$is.cryptic),]
  cryptic <- clean[ !is.na(clean$is.cryptic),]
  return( list(annotated, cryptic))
}

prepareOtherSplicingComparison <- function(deseqRes, exonList, nullIncluded, nullSkipped, genotypeCode){
  res <- as.data.frame(fread(deseqRes))
  # read in list of cryptic or skiptic exons
  cryptics <- as.data.frame(fread(exonList))
  if( genotypeCode == "F210I" ){
    res$is.cryptic <- cryptics$V4[ match(res$EnsemblID, cryptics$V4)]
  }else{
    res$is.cryptic <- cryptics$V4[ match(res$external_gene_id, cryptics$V4)]
  }
  includedExons <- as.data.frame(fread(nullIncluded))
  skippedExons <- as.data.frame(fread(nullSkipped))
  allChangedExons <- rbind(includedExons, skippedExons)
  
  if( genotypeCode == "F210I" ){
    res$is.alt.spliced <- allChangedExons$V4[ match(res$EnsemblID, allChangedExons$V4) ]
  }else{
    res$is.alt.spliced <- allChangedExons$V4[ match(res$external_gene_id, allChangedExons$V4) ]
  }
  annotated <- filter(res, !is.na(is.alt.spliced) & is.na(is.cryptic) )
  cryptic <- filter(res, !is.na(is.cryptic) )
  return( list(annotated, cryptic))
}

makeProportionPlot <- function(annotated,cryptic,genotypeCode, exonCode, controlCode){
  results <- data.frame( 
    group = c(
    "annotated UP", 
    "cryptic UP", 
    "annotated DOWN", 
    "cryptic DOWN"
    ),
    proportion = c(
      sum(annotated$log2FoldChange > 0 & annotated$padj < 0.1, na.rm = TRUE ), 
      sum(cryptic$log2FoldChange > 0 & cryptic$padj < 0.1, na.rm = TRUE),
      sum(annotated$log2FoldChange < 0 & annotated$padj < 0.1, na.rm = TRUE) ,
      sum(cryptic$log2FoldChange < 0 & cryptic$padj < 0.1, na.rm = TRUE )
    ),
    n = c(
      nrow(annotated),
      nrow(cryptic),
      nrow(annotated),
      nrow(cryptic)
    )
  )
  binomResults <- apply(results, MAR = 1, FUN = function(y){ binom.test(x = as.numeric( y[2] ), n = as.numeric( y[3] ) ) }$conf.int)
  results$confA <- binomResults[1,]
  results$confB <- binomResults[2,]
  
  # perform two tests comparing the up and the down
  tests <- c("upregulated", NA, "downregulated")
  pvalues <- c()
  for( i in c(1,3) ){
    hitInSample <- results[i+1, 2]
    hitInPop <- results[i, 2]
    failInPop <- results[i, 3] - results[i,2]
    sampleSize <- results[i+1,3]
    pval <- phyper(q = hitInSample-1, m = hitInPop, n = failInPop, k = sampleSize, lower.tail= FALSE)
    print(tests[i])
    print(pval)
    pvalues[i] <- pval
  }
  pvalues <- ifelse(
    pvalues > 0.05, "ns",
    ifelse(
      pvalues <= 0.05 & pvalues >= 0.001, "*", 
      ifelse(
        pvalues < 0.001 & pvalues > 0.0001, "**", "***"
      )
    )
  )
  names(pvalues) <- tests
  
  p <- ggplot(
    results, 
    aes( 
      x = str_split_fixed(group, " ", 2)[,2],
      y = proportion/n, fill = str_split_fixed(group, " ", 2)[,1] ) 
  ) +
    geom_col( position=position_dodge() ) +
    geom_errorbar( aes( ymin = confA, ymax = confB ), width = 0.2,position=position_dodge(.9) ) + 
    scale_y_continuous(labels=scales::percent, limits = c(0,0.5)) +
    ylab("Proportion significant at FDR < 10%") +
    scale_x_discrete(name = "", labels = c("downregulated", "upregulated") ) +
    scale_fill_discrete(name = "",
                        labels = c(
                          paste0(controlCode,"\n(", nrow(annotated), ")" ),
                          paste0( exonCode, "\n(", nrow(cryptic), ")" ) ) 
    ) +
    ggtitle(genotypeCode) +
    # p values
    annotate("text",x = 1, y = 0.5, label = pvalues[3]) +
    annotate("text", x = 2, y = 0.5, label = pvalues[1])
  
return(p)
}

combineAllPlot <- function(annotated, spliced, cryptic, genotypeCode, exon_codes){
  results <- data.frame( 
    group = c(
      "annotated genes\tUP", 
      "genes with cassette exons\tUP", 
      "genes with extreme events\tUP", 
      "annotated genes\tDOWN", 
      "genes with cassette exons\tDOWN", 
      "genes with extreme events\tDOWN"
    ),
    proportion = c(
      sum(annotated$log2FoldChange > 0 & annotated$padj < 0.1, na.rm = TRUE ), 
      sum(spliced$log2FoldChange > 0 & spliced$padj < 0.1, na.rm = TRUE ), 
      sum(cryptic$log2FoldChange > 0 & cryptic$padj < 0.1, na.rm = TRUE),
      sum(annotated$log2FoldChange < 0 & annotated$padj < 0.1, na.rm = TRUE) ,
      sum(spliced$log2FoldChange < 0 & spliced$padj < 0.1, na.rm = TRUE ), 
      sum(cryptic$log2FoldChange < 0 & cryptic$padj < 0.1, na.rm = TRUE )
    ),
    n = c(
      nrow(annotated),
      nrow(spliced),
      nrow(cryptic),
      nrow(annotated),
      nrow(spliced),
      nrow(cryptic)
    )
  )
  binomResults <- apply(results, MAR = 1, FUN = function(y){ binom.test(x = as.numeric( y[2] ), n = as.numeric( y[3] ) ) }$conf.int)
  results$confA <- binomResults[1,]
  results$confB <- binomResults[2,]
  
  tests <- NULL
  pvalues <- rep(NA,5)
  for( i in c(1,4) ){
    # sample doesn't change
    hitInSample <- results[i+2, 2]
    sampleSize <- results[i+2,3]
    # population does
    for( j in c(0,1) ){
      hitInPop <- results[i+j, 2]
      failInPop <- results[i+j, 3] - results[i+j,2]  
      pval <- phyper(q = hitInSample-1, m = hitInPop, n = failInPop, k = sampleSize, lower.tail= FALSE)
      pvalues[i+j] <- pval
      print(pval)
    }
  }
  pvalues <- ifelse(
    pvalues > 0.05, "ns",
    ifelse(
      pvalues <= 0.05 & pvalues >= 0.001, "*", 
      ifelse(
        pvalues < 0.001 & pvalues > 0.0001, "**", "***"
      )
    )
  )
  names(pvalues) <- tests
  print(pvalues)
  print(results)
  p <- ggplot(
    results, 
    aes( 
      x = str_split_fixed(group, "\t", 2)[,2],
      y = proportion/n, fill = str_split_fixed(group, "\t", 2)[,1] ) 
  ) +
    geom_col( position=position_dodge() ) +
    geom_errorbar( aes( ymin = confA, ymax = confB ), width = 0.2,position=position_dodge(.9) ) + 
    scale_y_continuous(labels=scales::percent, limits = c(0,0.5)) +
    ylab("Proportion significant at FDR < 10%") +
    scale_x_discrete(name = "", labels = c("downregulated", "upregulated") ) +
    scale_fill_discrete(name = "Splicing status",
                        #values = exon_codes,
                        labels = c(
                          paste0(exon_codes[1],"\n(", nrow(annotated), ")" ),
                          paste0( exon_codes[2], "\n(", nrow(spliced), ")" ),
                          paste0( exon_codes[3], "\n(", nrow(cryptic), ")" ))
    ) +
    ggtitle(genotypeCode) +
    # p values
    annotate("text",x = 2, y = 0.5, label = pvalues[1]) +
    annotate("text",x = 2.15, y = 0.46, label = pvalues[2]) +
    annotate("text", x = 1, y = 0.5, label = pvalues[4]) +
    annotate("text", x = 1.15, y = 0.46, label = pvalues[5]) +
    annotate("segment", x = 0.7, xend = 1.3, y = 0.48, yend = 0.48 ) + 
    annotate("segment", x = 1.7, xend = 2.3, y = 0.48, yend = 0.48 ) + 
    annotate("segment", x = 1, xend = 1.3, y = 0.44, yend = 0.44 ) +
    annotate("segment", x = 2, xend = 2.3, y = 0.44, yend = 0.44 ) +
    theme_bw()
  
  return(p)
  #return(results)
}


p_list <- list()
# F210I
annotatedComparison <- prepareAllAnnotatedComparison( deseqRes = F210Ires, exonList = crypticExons, genotypeCode = "F210I", exonCode = "cryptics" )
p_list[[1]] <- makeProportionPlot(annotated = annotatedComparison[[1]], cryptic = annotatedComparison[[2]], genotypeCode = "F210I",exonCode = "cryptics",controlCode = "annotated")
splicedComparison <- prepareOtherSplicingComparison(deseqRes = F210Ires, exonList = crypticExons,nullIncluded = F210I_included, nullSkipped = F210I_skipped,genotypeCode = "F210I") 
p_list[[2]] <- makeProportionPlot(annotated = splicedComparison[[1]], cryptic = splicedComparison[[2]], genotypeCode = "F210I",exonCode = "cryptics",controlCode = "other spliced")
p_list[[5]] <- combineAllPlot(annotated = annotatedComparison[[1]], spliced = splicedComparison[[1]], cryptic = splicedComparison[[2]], genotypeCode = "RRM2mut gene expression", exon_codes = c("no cassettes", "cassettes", "extreme") )


# M323K
annotatedComparison <- prepareAllAnnotatedComparison( deseqRes = M323Kres, exonList = skipticExons, genotypeCode = "M323K" )
p_list[[3]] <- makeProportionPlot(annotated = annotatedComparison[[1]], cryptic = annotatedComparison[[2]], genotypeCode = "M323K",exonCode = "skiptics",controlCode = "annotated")
splicedComparison <- prepareOtherSplicingComparison(deseqRes = M323Kres, exonList = skipticExons,nullIncluded = M323K_included, nullSkipped = M323K_skipped,genotypeCode = "M323K") 
p_list[[4]] <- makeProportionPlot(annotated = splicedComparison[[1]], cryptic = splicedComparison[[2]], genotypeCode = "M323K",exonCode = "skiptics",controlCode = "other spliced")
p_list[[6]] <- combineAllPlot(annotated = annotatedComparison[[1]], spliced = splicedComparison[[1]], cryptic = splicedComparison[[2]], genotypeCode = "LCDmut gene expression", exon_codes = c("no cassettes", "cassettes", "extreme") )


pdf(paste0(iFolder, "/differential_expression/", "F210I_differential_expression_plot.pdf"))
grid.arrange(p_list[[1]], p_list[[2]], ncol = 1)
dev.off()

pdf( paste0(iFolder, "/differential_expression/", "M323K_differential_expression_plot.pdf"))
grid.arrange(p_list[[3]], p_list[[4]], ncol = 1)
dev.off()

# for combining with NMD prediction plots:
pdf( paste0(iFolder, "/differential_expression/", "extreme_differential_expression_NMD_prediction.pdf"))
grid.arrange(p_list[[6]], p2)
grid.arrange(p_list[[5]], p1)
dev.off()
