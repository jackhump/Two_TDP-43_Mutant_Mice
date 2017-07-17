library(data.table)
library(dplyr)
library(stringr)
library(gridExtra)
# intronList <- "/SAN/vyplab/IoN_RNAseq/F210I_M323K_paper/RNA_Maps/data/noheader/M323K_adult_brain_se_intron_skipped.bed"
# exonList <- "/SAN/vyplab/IoN_RNAseq/F210I_M323K_paper/RNA_Maps/data/noheader/M323K_adult_brain_flank0_se_skipped.bed"

outFolder <- "/Users/Jack/google_drive/TDP_paper/intron_annotation/"

#annotation_code <- "/SAN/vyplab/HuRNASeq/leafcutter/leafcutter/data/gencode_mm10"
annotation_code <- paste0(outFolder,"/gencode_mm10")

pieList <- list()

annotateJunctions <- function( intronList, exonList, code){
  resultsFolder <- paste0(outFolder, code)
  if( ! dir.exists(resultsFolder)){
    dir.create(resultsFolder)
  }
  introns <- read.table(intronList, header=FALSE, stringsAsFactors=FALSE)
  exons <- read.table(exonList, header = FALSE , stringsAsFactors=FALSE)
  
  names(exons) <- c("chr","start","end","clusterID","score","strand")
  names(introns) <- c("chr","start","end","clusterID","score","strand")
  
  # remove duplicate exons - if there are multiple introns then they will be picked up.
  exons <- exons[ !duplicated(paste0(exons$chr, exons$start, exons$end) ), ]

  # F210I skipped is weird - need to match introns with exons
  # for each line in exons find a matching intron
  intronList <- lapply( 1:nrow(exons), FUN = function(i){
    ex <- exons[i,]
    int <- introns[ introns$chr == ex$chr &
                      introns$start < ex$start & introns$end > ex$end &
                      introns$strand == ex$strand, ]
    # add extra rows to ex if multiple introns are found - treat as separate junctions
    for( i in 1:nrow(int)){
      ex[i,] <- ex[1,]
    }
    return(list(ex, int) )
    
  })
  
  introns <- do.call( rbind, lapply(intronList, FUN = function(x){ x[[2]] }) )
  exons <- do.call( rbind, lapply(intronList, FUN = function(x){ x[[1]] }))
  # fix row names
  row.names(introns) <- 1:nrow(introns)
  row.names(exons) <- 1:nrow(exons)
  
  exons$clusterID <- paste0(exons$clusterID, "_", row.names(exons))
  introns$clusterID <- paste0(introns$clusterID, "_", row.names(introns))
  
  exon_fiveprime <- exons[,1:6]
  exon_threeprime <- exons[,1:6]
  
  # convert exon coordinates into list of junctions
  for( i in 1:nrow(introns)){
    exon_fiveprime$start[i] <- introns$start[i]
    exon_fiveprime$end[i] <- exons$start[i]
    exon_threeprime$start[i] <- exons$end[i]
    exon_threeprime$end[i] <- introns$end[i]
  }
  
  all <- rbind( introns, exon_fiveprime, exon_threeprime)
  
  names(all) <- c("chr","start","end","clusterID","score","strand")
  
  
  all.fiveprime <- data.frame( chr = all$chr,
                               start = all$start,
                               end = as.numeric( as.character(all$start) ) + 1,
                               clusterID = all$clusterID)
  all.threeprime <- data.frame( chr = all$chr,
                                start = all$end,
                                end = as.numeric( as.character(all$end) ) + 1,
                                clusterID = all$clusterID)
  
  all.bed <- select(all, chr, start, end, clusterID)
  
  all.file <- paste0(resultsFolder, "/all.bed")
  all.fiveprime.file <- paste0(resultsFolder, "/all.fiveprime.bed")
  all.threeprime.file <- paste0(resultsFolder, "/all.threeprime.bed")
  
  write.table( all.bed, all.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  write.table( all.threeprime, all.threeprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  write.table( all.fiveprime, all.fiveprime.file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t" )
  
  print( "BedTools intersect junctions with list of known splice sites")
  
  all.cmd <- paste0("bedtools intersect -a ", all.file, " -b ", annotation_code,"_all_introns.bed.gz", " -wa -wb -loj -f 1" )
  
  all_intersect <- fread(all.cmd)
  
  # intersect with bedtools to find the annotations of each splice site
  threeprime.cmd <- paste0( "bedtools intersect -a ", all.threeprime.file, " -b ", annotation_code,"_threeprime.bed.gz", " -wa -wb -loj -f 1" )
  
  threeprime_intersect <- fread(threeprime.cmd)
  
  fiveprime.cmd <- paste0( "bedtools intersect -a ", all.fiveprime.file, " -b ", annotation_code,"_fiveprime.bed.gz", " -wa -wb -loj -f 1" )
  
  fiveprime_intersect <- fread(fiveprime.cmd)
  
  # remove temporary files
  rm.cmd <- paste("rm ", all.fiveprime.file, all.threeprime.file, all.file) 
  system(rm.cmd)
  
  # now I have two lists of splice site annotation
  # for testing
  #cluster <- all[ all$clusterID == "clu_4879" , ]
  
  print("Annotating junctions")
  
  verdict.list <- list()
  coord.list <- list()
  gene.list <- list()
  ensemblID.list <- list()
  transcripts.list <- list()
  constitutive.list <- list()
  classification.list <- list()
  
  clusters <- unique( all$clusterID ) 
  for( clu in clusters ){
    # for each intron in the cluster, check for coverage of both
    # output a vector of string descriptions 
    cluster <- all[ all$clusterID == clu , ]
    #print(nrow(cluster))
    # for each intron in the cluster:
    #   create vector of overlapping splice sites, indexed by the row of the intersect
    # five prime splice sites
    fprime <- apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      start <- which( names(cluster) == "start" )
      fiveprime_intersect[   
        fiveprime_intersect$V1 == x[chr] & 
          fiveprime_intersect$V2 == as.numeric( x[start] )&
          fiveprime_intersect$V4 == clu,]
    } )
    # three prime splice sites
    tprime <- apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      end <- which( names(cluster) == "end" )
      threeprime_intersect[   
        threeprime_intersect$V1 == x[chr] & 
          threeprime_intersect$V2 == as.numeric( x[end] ) &
          threeprime_intersect$V4 == clu,]
    } )
    
    # both splice sites 
    bothSS <-  apply( cluster, MAR = 1, FUN = function(x) {
      chr <- which( names(cluster) == "chr" )
      start <- which(names(cluster) == "start")
      end <- which( names(cluster) == "end" )
      all_intersect[   
        all_intersect$V1 == x[chr] &
          all_intersect$V2 == x[start] &
          all_intersect$V3 == as.numeric( x[end] ) &
          all_intersect$V4 == clu,]
    } )
    
    # find gene and ensemblID by the most represented gene among all the splice sites
    cluster_genes <- names(sort(table(do.call( what = rbind, tprime )$V8), decreasing = TRUE ))
    
    cluster_gene <- cluster_genes[ cluster_genes != "." ][1]
    # if no cluster gene found then leave as "."
    if( length(cluster_gene) == 0){
      cluster_gene == "."
    }
    # do the same for EnsemblID
    cluster_ensemblIDs <- names(sort(table(do.call( what = rbind, tprime )$V9), decreasing = TRUE ))
    cluster_ensemblID <- cluster_ensemblIDs[ cluster_ensemblIDs != "." ][1]
    if( length( cluster_ensemblID ) == 0 ){
      cluster_ensemblID == "."
    }
    
    verdict <- c()
    coord <- c()
    gene <- c()
    ensemblID <- c()
    transcripts <- list() 
    
    for( intron in 1:nrow(cluster) ){
      coord[intron] <- paste(cluster[intron,]$chr,cluster[intron,]$start, cluster[intron,]$end )
      # record all transcripts that use the splice sites for each intron
      
      # tgene <- names(sort(table( tprime[[intron]]$V8 ), decreasing = TRUE)[1])
      # fgene <- names(sort(table( fprime[[intron]]$V8 ), decreasing = TRUE)[1])
      
      # tensemblID <- names(sort(table( tprime[[intron]]$V8 ), decreasing = TRUE)[1])
      # fensemblID <- names(sort(table( fprime[[intron]]$V8 ), decreasing = TRUE)[1])
      
      
      # gene[intron] <- ifelse( tgene == ".",  no = fgene, yes = tgene )
      # ensemblID[intron]<- ifelse( tensemblID == ".", no = fensemblID, yes = tensemblID )
      
      gene[intron] <- cluster_gene
      ensemblID[intron] <- cluster_ensemblID
      
      # for each intron create vector of all transcripts that contain both splice sites
      transcripts[[intron]] <- unique( bothSS[[intron]]$V10 ) 
      
      verdict[intron] <- "error"
      if(
        all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 == "." )
      ){ verdict[intron] <- "cryptic_unanchored"
      }
      if(
        all( tprime[[intron]]$V5 == ".") & all( fprime[[intron]]$V5 != "." )
      ){ verdict[intron] <- "cryptic_threeprime"
      }
      if(
        all( tprime[[intron]]$V5 != ".") & all( fprime[[intron]]$V5 == "." )
      ){ verdict[intron] <- "cryptic_fiveprime"
      }
      if(
        all( tprime[[intron]]$V5 != "." ) & all( fprime[[intron]]$V5 != "." )
      ){ 
        # test if the splice sites are paired in a known intron
        #tp <- paste( tprime[[intron]]$V9, tprime[[intron]]$V10 )
        #fp <- paste( fprime[[intron]]$V9, fprime[[intron]]$V10 )
        if( all(bothSS[[intron]]$V5 != ".") ){
          verdict[intron] <- "annotated"
        }else{
          verdict[intron] <- "skiptic"
        }
        
        # if( length( intersect(fp,tp) ) > 0 ){
        #   verdict[intron] <- "annotated"
        # }else{
        #   verdict[intron] <- "skiptic"
        # }
      }
      verdict.list[[clu]] <- verdict
      coord.list[[clu]] <- coord
      gene.list[[clu]] <- gene
      ensemblID.list[[clu]] <- ensemblID
      #transcripts.list[[clu]] <- transcripts
      
      # once all the transcripts for all the introns are found, go back and work out how many constitutive each junction is. Does the junction appear in every transcript? 
      
      if( intron == nrow(cluster)){ # only on final intron
        all_transcripts <- unique( unlist( transcripts ) )
        # remove "." - non-existent transcripts
        all_transcripts <- all_transcripts[ all_transcripts != "." ]
        
        constitutive <- lapply( transcripts, FUN = function(x) {
          # for each intron how many transcripts is it seen in?
          x <- x[ x != "." ]
          length(x) / length( all_transcripts)
          
        })
        
        constitutive.list[[clu]] <- constitutive
        
        # collapse all transcripts for each intron into a single string
        transcripts.list[[clu]] <- lapply(transcripts, FUN = function(x) paste( x, collapse = "+" ) )
        
      }
      
    }
    
    # predicting the event type from the shape of the junctions
    # easy start - cassette exons
    #print(clu)
    
    if( nrow(cluster) != 3){ 
      classification.list[[clu]] <- "." 
      next
    }else{
      classification.list[[clu]] <- "."
      
      tab <- select(cluster, start, end)
      tab$verdict <- verdict.list[[clu]]
      
      # the junctions are sorted by start and end coordinates
      tab <- arrange(tab, start, end)
      
      
      # check for the presence of a junction that spans the entire length of the cluster
      if( !any(  which( tab$start == min(tab$start) ) %in% which( tab$end == max(tab$end) )  ) ){
        classification.list[[clu]] <- "."
        next
      }
      
      # therefore for a cassette exon arrangement the longest junction always comes second 
      if( which( tab$start ==  min(tab$start) & tab$end == max(tab$end ) ) != 2 ){
        classification.list[[clu]] <- "." 
        next 
      }
      
      # now we know that junction 2 is the parent, junction 1 is the left most child and junction 3 is the right most
      # check that the end of junction 1 comes before the start of junction 3
      
      if( tab[1,"end"] > tab[3,"start"] ){
        classification.list[[clu]] <- "."
        next
      }
      
      # double check the starts and ends
      if( tab[1, "start"] != tab[2,"start"] | tab[3,"end"] != tab[2,"end"] ){
        classification.list[[clu]] <- "."
        next
      }
      
      # work out direction of change
      # if( cluster[1, "deltapsi"] > 0 & cluster[3, "deltapsi"] > 0 & cluster[2,"deltapsi"] < 0){
      #   classification.list[[clu]] <- "cassette exon - increased"
      # }
      # if( cluster[1, "deltapsi"] < 0 & cluster[3, "deltapsi"] < 0 & cluster[2,"deltapsi"] > 0){
      #   classification.list[[clu]] <- "cassette exon - decreased"
      # }
      
      classification.list[[clu]] <- "cassette"
      
      # work out annotation status
      if( all( tab$verdict == "annotated") ){
        classification.list[[clu]] <- paste0( classification.list[[clu]], " - annotated")
      }
      
      if( tab$verdict[2] == "annotated" & tab$verdict[1] != "annotated" & tab$verdict[3] != "annotated"  ){
        classification.list[[clu]] <- paste0( classification.list[[clu]], " - cryptic")
      }
      
      if( tab$verdict[2] != "annotated" & tab$verdict[1] == "annotated" & tab$verdict[3] == "annotated"  ){
        classification.list[[clu]] <- paste0( classification.list[[clu]], " - skiptic")
      }
      
      if( tab$verdict[2] == "annotated" & 
          ( tab$verdict[1] != "annotated" & tab$verdict[3] == "annotated" ) |
          ( tab$verdict[1] == "annotated" & tab$verdict[3] != "annotated" ) 
      ){
        classification.list[[clu]] <- "cryptic extension"
      }
      # anything weird 
      if( classification.list[[clu]] == "cassette" ){
        classification.list[[clu]] <- "complex annotation"
      }
      # print(n_parent)
      # print(children_list)
      
    }
    
    # print(clu)
    # # print(n_parent)
    # # print(children_list)
    # print(classification.list[[clu]])
    
  }
  print("Preparing results")
  
  # match all the lists together
  all$verdict <- unlist(verdict.list)[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]
  
  all$gene <- unlist(gene.list)[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]
  
  all$ensemblID <- unlist(ensemblID.list)[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]
  
  all$transcripts <- unlist( transcripts.list )[ match( paste( all$chr, all$start, all$end ), unlist(coord.list)) ]
  
  all$prediction <-  unlist( classification.list )[ match( all$clusterID, names(classification.list) )]
  
  annotations <- data.frame(
    geneID = names(classification.list),
    class = unlist(classification.list) )
  
  annotations$geneID <- str_split_fixed(annotations$geneID, "_", 2)[,1]
  
  annotations.out <- paste0(resultsFolder, "/cluster_annotations.tab")
  write.table( annotations, annotations.out, col.names=TRUE, row.names =FALSE, quote = FALSE, sep = "\t")
  
  all.out <- paste0(resultsFolder, "/junction_annotations.tab")
  write.table( all, all.out, col.names=TRUE, row.names =FALSE, quote = FALSE, sep = "\t" )
  
  # make a pie
  pieLabels <- paste0( names(table(annotations$class)), "\n(", table(annotations$class), ")" )
  
  pdf( paste0(resultsFolder, "/", code, "_pie.pdf"))
  junctionPie <- pie( table( annotations$class ), labels = pieLabels, main = gsub("_", " ", code), cex = 0.75) 
  print(junctionPie)
  dev.off()
  
  annotations$is.annotated <- ifelse( annotations$class == "cassette - annotated",
                                      "annotated", "novel")
  annotations$group <- code
  
  annotations$is.annotated <- ifelse( annotations$class == "cassette - annotated",
                                      "annotated", "novel")
  pieLabels <- paste0( names(table(annotations$is.annotated)), "\n(", table(annotations$is.annotated), ")" )
  print(pieLabels)
  
  pdf( paste0(resultsFolder, "/", code, "_clean_pie.pdf"))
  cleanPie <- pie( table( annotations$is.annotated ), labels = pieLabels, main = gsub("_", " ", code), cex = 0.75)
  print(cleanPie)
  dev.off()
  
  return(annotations)
}


makeAllPie <- function(annotations){
  pieLabels <- paste0( names(table(annotations$class)), "\n(", table(annotations$class), ")" )
  junctionPie <- pie( table( annotations$class ), labels = pieLabels, main = gsub("_", " ", code), cex = 0.75) 
  print(junctionPie)
}

makeCleanPie <- function(annotations){
  annotations$is.annotated <- ifelse( annotations$class == "cassette - annotated",
                                      "annotated", "novel")
  pieLabels <- paste0( names(table(annotations$is.annotated)), "\n(", table(annotations$is.annotated), ")" )
  junctionPie <- pie( table( annotations$is.annotated ), labels = pieLabels, main = gsub("_", " ", code), cex = 0.75)
  return(junctionPie)
}

#UNCOMMENT THIS
intronList <- paste0(outFolder, "M323K_adult_brain_se_intron_skipped.bed" )
exonList <- paste0(outFolder, "M323K_adult_brain_flank0_se_skipped.bed" )
code <- "M323K_skipped"

pieList[[1]] <- annotateJunctions(intronList = intronList, exonList = exonList, code = code)

intronList <- paste0(outFolder, "M323K_adult_brain_se_intron_included.bed" )
exonList <- paste0(outFolder, "M323K_adult_brain_flank0_se_included.bed" )
code <- "M323K_included"

pieList[[2]] <- annotateJunctions(intronList = intronList, exonList = exonList, code = code)

intronList <- paste0(outFolder, "F210I_embryonic_brain_se_intron_included.bed" )
exonList <- paste0(outFolder, "F210I_embryonic_brain_flank0_se_included.bed" )
code <- "F210I_included"

pieList[[3]] <- annotateJunctions(intronList = intronList, exonList = exonList, code = code)

intronList <- paste0(outFolder, "F210I_embryonic_brain_se_intron_skipped.bed" )
exonList <- paste0(outFolder, "F210I_embryonic_brain_flank0_se_skipped.bed" )
code <- "F210I_skipped"

pieList[[4]] <- annotateJunctions(intronList = intronList, exonList = exonList, code = code)

# F210I cryptic exons
intronList <- "/Users/Jack/google_drive/TDP_paper/RNA_maps/cryptic_skiptic_only/F210I_cryptic_introns.bed"
exonList <- "/Users/Jack/google_drive/TDP_paper/cryptics.bed"
code <- "F210I_cryptics"
pieList[[5]] <- annotateJunctions(intronList = intronList, exonList = exonList, code = code)


intronList <- "/Users/Jack/google_drive/TDP_paper/RNA_maps/cryptic_skiptic_only/M323K_skiptic_introns.bed"
exonList <- "/Users/Jack/google_drive/TDP_paper/skiptics.bed"
code <- "M323K_skiptics"
pieList[[6]] <- annotateJunctions(intronList = intronList, exonList = exonList, code = code)



