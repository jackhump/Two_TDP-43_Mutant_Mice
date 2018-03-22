# NMD prediction 
# Jack Humphrey 2017
# takes as input lists of exon as Bed files

# exons are either skipped or included

# use the splice junction lists from STAR to find the connecting junctions

# for skipped exons use the control junctions

# for included exons use the case junctions

# use annotation to find the outer boundaries that make up the connecting exons 
# and the reading frame of the leading exon

# translate to FASTA sequence

# compare inclusion vs exclusion sequences, look for frameshifts and PTCs



library(optparse)

library(data.table,quietly=T)
library(dplyr)
library(stringr)
library(GenomicRanges,quietly=T)
library(ggplot2)
library(Biostrings)

options(echo = TRUE)

N_PERMUTE=100
zero_base_fuckery=TRUE

################
### FUNCTIONS
#################



files.exist <- function(files.list){
  for(my.file in files.list){
    if(!file.exists(my.file)){
      stop(paste(my.file,"doesn't exist!"))
    }
  }
}
files.are.empty <- function(files.list){
  if(length(files.list) == 0){
    stop(paste(files.list,"is empty!"))
  }
  
}

# This function reads in all the SJ.tab files given in the list and merges them all together to give a list of unique SJs with total numbers of occurences across all the datasets in the list.
merge.SJ.files <- function(SJ.tab.list){
  # First read in all the files in the list
  SJ.file.list <- list()
  print( "reading in junctions")
  for(i in 1:length(SJ.tab.list)){
          SJ <- as.data.frame(fread(SJ.tab.list[i]))
          SJ.file.list[[i]] <- SJ
        }
  # rbind all together
  SJ.merge <- rbindlist(SJ.file.list)
  SJ.merge <- as.data.frame(SJ.merge)
  #Use dplyr to group the total list of splice junctions by unique SJ and then sum the total of the cases.
  print( "merging junctions with dplyr")
  by_coord <- group_by(SJ.merge, paste(V1,V2,V3))
  SJ.summary <- summarise(by_coord,
          count.unique = sum(V7),
          count.multi = sum(V8),
          strand = mean(V4),
      intron.motif = mean(V5)
          )
  names(SJ.summary)[1] <- "coord"
  #split the coordinate reference into columns
  SJ.summary$chr <- str_split_fixed(SJ.summary$coord, " ", 3)[,1]
  SJ.summary$start <- str_split_fixed(SJ.summary$coord, " ", 3)[,2]
  SJ.summary$end <- str_split_fixed(SJ.summary$coord, " ", 3)[,3]
  #this is where columns are dropped?
  SJ.summary <- SJ.summary[,c(6,7,8,2,3,4,5)]
  SJ.summary <- as.data.frame(SJ.summary)
  SJ.summary$start <- as.numeric(SJ.summary$start)
  SJ.summary$end <- as.numeric(SJ.summary$end)
  # make strand readable
  SJ.summary[SJ.summary$strand == 1,]$strand <- "+"
  SJ.summary[SJ.summary$strand == 2,]$strand <- "-"
    SJ.summary[SJ.summary$strand == 0,]$strand <- "*"
    #make intron motif readable
  motif.list <- c("non-canonical","GT/AG","CT/AC","GC/AG","CT/GC","AT/AC","GT/AT")
  for(i in 1:7){
    if(length(SJ.summary[SJ.summary$intron.motif == (i - 1),]$intron.motif) > 0){
      SJ.summary[SJ.summary$intron.motif == (i - 1),]$intron.motif <- motif.list[i]
    } 
  }
  #unique SJs will become "score" field in GRanges object
  names(SJ.summary)[4] <- "score"
  return(SJ.summary)
}

upstream_junction_query <- function(CE.chr,CE.start,CE.end, CE.ID, SJ.GRange){
    junction <- SJ.GRange[ 
      seqnames(SJ.GRange) == CE.chr & 
      #start(SJ.GRange) == as.numeric(canonical.start) & 
      end(SJ.GRange) >= as.numeric(CE.start) - 1 &  # used to be >= start - 1
      end(SJ.GRange) < as.numeric(CE.end) 
    ]

    junction <- head(junction[order(score(junction),decreasing=T)],1)
    if( length(junction) == 0 ){
      # return null junction
      junction <-  SJ.GRange[1]
      start(junction) <- 1
      end(junction) <- 2
      CE.ID <- "."
    }
    elementMetadata(junction)$ID <- as.character(CE.ID)
    return(junction)
}

downstream_junction_query <- function(CE.chr,CE.start,CE.end, CE.ID, SJ.GRange){
    junction <- SJ.GRange[
      seqnames(SJ.GRange) == CE.chr & 
      start(SJ.GRange) > as.numeric(CE.start) & 
      start(SJ.GRange) <= as.numeric(CE.end) + 1 
      #& start(SJ.GRange) >= as.numeric(CE.end) - 10# used to be <=
      ]
    junction <- head(junction[order(score(junction),decreasing=T)],1)
    if( length(junction) == 0 ){
      # return null junction
      junction <-  SJ.GRange[1]
      start(junction) <- 1
      end(junction) <- 2
      CE.ID <- "."
    }
    elementMetadata(junction)$ID <- CE.ID
    return(junction)
}

#This function assumes that the results file has been appended with the canonical start and end coordinates at positions ??? and ??? respectively.
bridging_junction_finder <- function(SJ.summary, results.df, query.type){

  GRanges_object <-  makeGRangesFromDataFrame(SJ.summary,keep.extra.columns=T)
  # create null chromosome for returning null results
  #seqlevels(GRanges_object) <- c( seqlevels(GRanges_object), "chrNull")
  if(query.type == "downstream"){
    junctions.list <- apply(results.df, MAR=1,FUN=function(x) downstream_junction_query(x[1], x[2], x[3], x[7],GRanges_object))
  }
  if(query.type == "upstream"){
    junctions.list <- apply(results.df, MAR=1,FUN=function(x) upstream_junction_query(x[1],x[2],x[3], x[7],GRanges_object))
  }
  #output is a list of GRange objects - unuseable.
  junctions.list <- unlist(GRangesList(junctions.list))
  #convert into a dataframe, extracting the relevent information from the GRanges object.
  #names(GRanges) is a vector of rownames, confusingly.
  if(query.type == "downstream"){
    bridging_junctions.df <- data.frame(
      #row.names=names(junctions.list),
      chr=seqnames(junctions.list),
      central.end=start(junctions.list),
      downstream.start=end(junctions.list),
      downstream.unique.count = score(junctions.list),
      downstream.strand = strand(junctions.list),
      intron.motif = mcols(junctions.list)[3],
      ID = elementMetadata(junctions.list)$ID)
  }
  if(query.type == "upstream"){
    bridging_junctions.df <- data.frame(
      #row.names=names(junctions.list),
      chr=seqnames(junctions.list),
      upstream.end=start(junctions.list),
      central.start=end(junctions.list),
      upstream.unique.count = score(junctions.list),
      upstream.strand = strand(junctions.list),
      intron.motif = mcols(junctions.list)[3],
    ID = elementMetadata(junctions.list)$ID)
  }
  return(bridging_junctions.df)
}

resolve_junctions <- function(junction.chr, junction.start, junction.end, junction.ID, cds.grange, mode){
  if( mode == "upstream"){
    exon <- cds.grange[ 
      seqnames(cds.grange) == as.character(junction.chr) & 
      end(cds.grange) == (as.numeric(junction.start) - 1)
    ]
  }
  if(mode == "downstream"){
    exon <- cds.grange[ 
      seqnames(cds.grange) == junction.chr & 
      start(cds.grange) == ( as.numeric(junction.end) + 1 )
    ]
  }
  if( length(exon) == 0 ){
    exon <- cds.grange[1]
    start(exon) <- 1
    end(exon) <- 2
    junction.ID <- "."
  }
  exon <- head(exon,1)
  elementMetadata(exon)$ID <- junction.ID

return(exon)
}

# rewrite translate toxic function
#testing
    # x <- clean_exons[12,]
    # upstream_fasta <- as.character( x$upstream_seq)
    # central_fasta <- as.character( x$central_exon)
    # downstream_fasta <-  as.character(x$downstream_seq)
    # upstream_phase <-  as.numeric(as.character(x$upstream_phase))
    # downstream_phase <- as.numeric(as.character(x$downstream_phase))
    # strand <-  as.character(x$strand)

translate_toxic <- function(upstream_fasta, central_fasta, downstream_fasta, upstream_phase, downstream_phase,strand, permute=FALSE){
    upstream_fasta <- as.character(upstream_fasta)
    central_fasta <- as.character(central_fasta)
    downstream_fasta <- as.character(downstream_fasta)
    upstream_phase <- as.numeric(as.character(upstream_phase) )
    downstream_phase <- as.numeric(as.character(downstream_phase) )
    strand <- as.character(strand)
    print(central_fasta)
    if( strand == "+"){
        a <- DNAString( upstream_fasta )
        a_b <- DNAString( paste0(upstream_fasta, central_fasta) )
        a_c <- DNAString( paste0(upstream_fasta,downstream_fasta) )
        a_b_c <- DNAString( paste0(upstream_fasta,central_fasta,downstream_fasta) )
        c <- reverseComplement( DNAString( downstream_fasta ) ) 
        phase <- upstream_phase + 1
      #c.phase <- downstream_phase
    }
    if( strand == "-"){
        a <- reverseComplement( DNAString( downstream_fasta ) ) 
        a_b <- reverseComplement( DNAString( paste0(central_fasta,downstream_fasta) ) )
        a_c <- reverseComplement( DNAString( paste0(upstream_fasta,downstream_fasta) ) )
        a_b_c <- reverseComplement( DNAString( paste0(upstream_fasta,central_fasta,downstream_fasta) ) )
        c <- reverseComplement( DNAString( upstream_fasta ) ) 
        phase <- downstream_phase + 1
      #c.phase <- upstream_phase
    }

    translate_orf <- function( seq, codon_phase ){
        protein <- suppressWarnings(translate( subseq(seq, start=codon_phase) ))
        return(protein)
    }
    count_PTC <- function( protein ){
      section_length <- str_length(protein)
      if( section_length > 17 ){
        section_length <- section_length - 17
        return( countPattern("*", subseq(protein, end = section_length ) ) )
      }else{
        # if exon is 17 amino acids or smaller then by definition any stop codons should escape NMD
        return(0)
      }
    }

    compare_c_sections <- function( a_b, a_b_c, c, phase){
        c.length <- floor( str_length(c) / 3)
        c_with_b <- translate_orf(a_b_c, phase)
        # cut out translated C section and compare
        c_with_b <- subseq( c_with_b, start =  ( str_length(c_with_b) - c.length + 2 )  )
        c_without_b <- translate_orf(a_c, phase)
        c_without_b <- subseq(c_without_b, start = ( str_length(c_without_b) - c.length + 2 ) )
        # true if the two are different
        c.conserved <- ( c_with_b != c_without_b )

        return( c.conserved )
    }
    # translate each segment
    protein_seq <- lapply( list(a, a_b, a_c, a_b_c ), FUN = function(seq) translate_orf(seq, phase) )
    
    functionality <- lapply( protein_seq, FUN = function(seq) count_PTC(seq) )

    c.conserved <- compare_c_sections( a_b, a_b_c, c, phase )

    functionality <- c(functionality, c.conserved)
    # invert logic - no PTCs (score 0) equals a functional transcript
    # 5 should be TRUE if C's sequence is maintained
    functionality <- !as.logical(functionality)
    # catch sneaky NMD escaping stop codons in C
    if( all(functionality == c(TRUE,TRUE,TRUE,TRUE,FALSE) ) & countPattern("*", protein_seq[[4]]) > 0 ){
      functionality <- c(F,F,F,F,F) # changes it to mean that there is a frameshift causing an NMD-scaping stop codon
    }


    return(functionality) 
  
}



# # F210I
# control_SJ_list <- c(
#   "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_CTL_1_norm/F210I_CTL_1_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_CTL_2_norm/F210I_CTL_2_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_CTL_3_norm/F210I_CTL_3_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_CTL_4_norm/F210I_CTL_4_normSJ.out.tab"
# )

# case_SJ_list <- c(
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_HOM_1_norm/F210I_HOM_1_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_HOM_2_norm/F210I_HOM_2_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_HOM_3_norm/F210I_HOM_3_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_HOM_4_norm/F210I_HOM_4_normSJ.out.tab"
#   )

# code <- "F210I"

# # exon lists!
# included_exon_list <- "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/F210I_se_included.bed"
# skipped_exon_list <- "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/F210I_se_skipped.bed"


# #############
# # M323K
# ##############

# control_SJ_list <- c(
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_WT_1/M323K_WT_1SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_WT_2/M323K_WT_2SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_WT_3/M323K_WT_3SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_WT_4/M323K_WT_4SJ.out.tab"
# )

# case_SJ_list <- c(
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_HOM_1/M323K_HOM_1SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_HOM_2/M323K_HOM_2SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_HOM_3/M323K_HOM_3SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_HOM_4/M323K_HOM_4SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_HOM_5/M323K_HOM_5SJ.out.tab"
# )

# included_exon_list <- "/SAN/vyplab/IoN_RNAseq/Kitty/M323K/sgseq/M323K_adult_brain_se_included.bed"
# skipped_exon_list <- "/SAN/vyplab/IoN_RNAseq/Kitty/M323K/sgseq/M323K_adult_brain_se_skipped.bed" 

# code <- "M323K"


# ###########
# ## NULL EXONS
# ###########

# #all_exons_f210i <- "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/F210I_embryonic_brain_res_clean.tab"
# all_mouse_exons <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/all_mouse_exons_merged.bed"
# all <- as.data.frame(fread(all_mouse_exons) )
# n.sample = 100
# all_sample <- sample_n(all, size = n.sample, replace = FALSE)
# # fix exon starts 
# all_sample$V2 <- all_sample$V2-1
# # sort out columns
# all_sample$p.value <- "."
# all_sample$gene_id <- "."
# all_sample <- select( all_sample, V1, V2, V3,gene_id,p.value,V4)


# # use all control samples for null exons
# control_SJ_list <- c(  
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_CTL_1_norm/F210I_CTL_1_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_CTL_2_norm/F210I_CTL_2_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_CTL_3_norm/F210I_CTL_3_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/F210I/New_embryonic_brain/processed/F210I_CTL_4_norm/F210I_CTL_4_normSJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_WT_1/M323K_WT_1SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_WT_2/M323K_WT_2SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_WT_3/M323K_WT_3SJ.out.tab",
# "/SAN/vyplab/IoN_RNAseq/M323K/New_adult_brain/processed/M323K_WT_4/M323K_WT_4SJ.out.tab" )

# code <- "all_exons_sampled"
# exon_code <- "sample_all_exons"
# junction_list <- control_SJ_list


# outFolder="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/cryptic_skiptic/NMD_prediction"
# code="skiptic"
# species="mouse"
# support.frame="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/NMD_prediction/m323k_support.tab"

# null_exons="/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/cryptic_skiptic/skiptics_nohead.bed"



#########
## BEGIN



  
# set skipped exons as central_exons
#exon_code <- c("skipped_exons", "included_exons","sampled_null")
#junction_list <- c(control_SJ_list, case_SJ_list, control_SJ_list )
#exon_list <- c(skipped_exon_list, included_exon_list,all_sample)



# outFolder <- paste0( "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/F210I_M323K/NMD_prediction/", code)




# options(echo=T) 






# UNCOMMENT THIS!!
# Rscript $script --outFolder $outFolder \
#                 --code $code \
#                 --species $species \
#                 --support $support \
#                 --skipped $skipped_exons \
#                 --included $included_exons



option_list <- list(
    make_option(c('--support'), help=''),
    make_option(c('--code'), help=''),
    make_option(c('--skipped'), help=''),
    make_option(c('--included'), help=''),
    make_option(c('--null_exons'), help = ''),
    make_option(c('--outFolder'), help=''),
    make_option(c('--species'), help=''),
    make_option(c('--permute'), help = 'should you permute the central exons?') 
)

########################## read arguments
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

if(length(opt) > 1){
  support.frame <- opt$support
  code <- opt$code
  skipped_exon_list <- opt$skipped
  included_exon_list <-  opt$included
  outFolder <- opt$outFolder
  species <- opt$species
  null_exons <- opt$null_exons
  permute <- opt$permute
}

print(species)

  
data_list <- c()

if( !is.null(null_exons) ){
  mode <- c("null")
  data_list <- c(null_exons)
}

if( !is.null( skipped_exon_list ) & !is.null(included_exon_list) ){
  mode <- c(mode, "both")
  data_list <- c(data_list, skipped_exon_list, included_exon_list)
}else{
  if( !is.null( skipped_exon_list ) & is.null(included_exon_list) ){
    mode <- c(mode, "skipped")
    data_list <- c(data_list, skipped_exon_list )
  }

  if( is.null( skipped_exon_list ) & !is.null(included_exon_list) ){
    mode <- c(mode, "included")
    data_list <- c(data_list, included_exon_list)
  }
}

print("exons to look at:")
print(data_list)

print("mode is:")
print(mode)

#### read support frame

support <- read.table(support.frame, header = FALSE, stringsAsFactors = FALSE)

#print(support)

control_SJ_list <- support[ grepl("control|ctl", support$V2, ignore.case = TRUE), ]$V1
case_SJ_list <- support[ !grepl("control|ctl", support$V2, ignore.case = TRUE), ]$V1


print(control_SJ_list)
print(case_SJ_list)


if( ! dir.exists( outFolder) ){
  dir.create(outFolder)
}

outFolder <- paste0(outFolder, "/",code)

toxic.outFolder <- paste(outFolder,"/toxic_exons", sep = "/")
fasta.outFolder <- paste(outFolder, "/DNAfasta",sep="/")

for(folder in c(toxic.outFolder, fasta.outFolder, outFolder)){
  if (! dir.exists(folder)) dir.create(folder,recursive=T)
}

print(data_list)


for( file in c(case_SJ_list, control_SJ_list, unlist(data_list) ) ){
    if( ! file.exists(file) ){
      stop( paste0( file, " do not exist"))  
    }
}



##############################
# SPECIES-SPECIFIC ANNOTATION 
##############################

if( species == "mouse" ){
  exon.gtf <- "/SAN/vyplab/HuRNASeq/GENCODE/gencode.vM12.annotation.gtf"
  genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/mm10.fa"
  # back up genome
  genome.fa <- "/SAN/vyplab/HuRNASeq/reference_datasets/mm10.fa"
}

if( species == "human" ){
  exon.gtf <- "/SAN/vyplab/HuRNASeq/GENCODE/gencode.v25.annotation.gtf"
  genome.fa <- "/SAN/vyplab/HuRNASeq/reference_datasets/hg38.fa"
}

# read in exon lists

exon_list <- list()

for( i in 1:length(data_list) ){
  exon_list[[i]] <- as.data.frame(fread( data_list[i] ))
}

exon_code <- c()

n.sample <- 1000
# if null mode then random sample of the null exons
# Null exons must be in tab delimited format
#   chr, start, end, strand 
if( "null" %in% mode ){
  if( nrow( exon_list[[1]] ) > n.sample ){
  all_sample <- sample_n(exon_list[[1]], size = n.sample, replace = FALSE)
  }else{ 
    all_sample <- exon_list[[1]]
  }
  all_sample$p.value <- "."
  all_sample$gene_id <- "."
  all_sample <- select( all_sample, V1, V2, V3,gene_id,p.value,V4)
  exon_list[[1]] <- all_sample

  exon_code <- "null_exons"
}

if("both" %in% mode){
  exon_code <- c(exon_code, "skipped_exons", "included_exons")
}
if("skipped" %in% mode){
  exon_code <- c(exon_code, "skipped_exons")
}
if("included" %in% mode){
  exon_code <- c(exon_code, "included_exons")
}


print(exon_code)

##########
# BEGIN
##########


for( exon_type in 1:length(exon_list) ){

  central_exons <- exon_list[[ exon_type ]]

  names(central_exons) <- c("chr","central.start","central.end", "ensemblID","p-value","strand")
  
  central_exons$ID <- paste(central_exons$chr, central_exons$central.start, central_exons$central.end, sep = ":")

  # the start is fudged so the fasta sequence lines up

  #central_exons$central.start <- central_exons$central.start +1
 

  # control_SJs.out <- paste0(outFolder,"/",code,"_SJs_control.tab")
  # case_SJs.out <- paste0(outFolder,"/",code,"_SJs_case.tab")
  # do the same to case_SJs!

  # remove exons that have 0 values - these are extensions rather than cassettes
  print("total number of exons")
  print(nrow(central_exons))
  central_exons <- central_exons[ central_exons$central.start > 100 & central_exons$central.end > 100, ]
  print("number of exons after removing zero values")
  print(nrow(central_exons))

  if( exon_code[exon_type] == "skipped_exons" ){
    sj_merge_out <- paste0( outFolder,"/exclusion_splice_junctions.RData")
    if( file.exists(sj_merge_out) ){
      load(sj_merge_out)
    }else{
      control_total_SJ_counts <- merge.SJ.files(control_SJ_list)

      # remove any junctions that have counts less than the number of samples from the merge
      min_number <- length( control_SJ_list )
      control_total_SJ_counts <- control_total_SJ_counts[ control_total_SJ_counts$score >= min_number, ]

      SJ_merge <- control_total_SJ_counts
      save(SJ_merge, file = sj_merge_out)
    }
  }

  if( exon_code[exon_type] == "included_exons" ){
    sj_merge_out <- paste0( outFolder,"/inclusion_splice_junctions.RData")
    if( file.exists(sj_merge_out) ){
      load(sj_merge_out)
    }else{

      case_total_SJ_counts <- merge.SJ.files(case_SJ_list)

    # remove any junctions that have counts less than the number of samples from the merge
      min_number <- length( case_SJ_list )
      case_total_SJ_counts <- case_total_SJ_counts[ case_total_SJ_counts$score >= min_number, ]

      SJ_merge <- case_total_SJ_counts
      save(SJ_merge, file = sj_merge_out)
    }
  }

  if( exon_code[exon_type] == "null_exons" ){
      sj_merge_out <- paste0( outFolder,"/all_splice_junctions.RData")


    
    if( file.exists(sj_merge_out) ){
      load(sj_merge_out)
    }else{
      all_SJs <- merge.SJ.files( c(case_SJ_list, control_SJ_list) )
      #min_number <- length( c(case_SJ_list, control_SJ_list))
      min_number <- 5
      all_SJs <- all_SJs[ all_SJs$score >= min_number, ]
      SJ_merge <- all_SJs
      save(SJ_merge, file = sj_merge_out)
    }
  }

  # only accept exact matches!

  # upstream junctions in control
  upstream_results_control <- bridging_junction_finder(
                             SJ.summary = SJ_merge,
                             results.df = central_exons,
                             query.type = "upstream")

  # downstream junctions in control
  downstream_results_control <- bridging_junction_finder(
                               SJ.summary = SJ_merge, 
                               results.df = central_exons, 
                               query.type = "downstream")

  print("# how many exons can be matched to either an upstream or downstream junction?")
  print(table(upstream_results_control$upstream.end >1 | downstream_results_control$downstream.start > 2))
  # match the remaining exon coordinates and readling frame from a GTF file
  if( !exists("gtf")){
    gtf <- as.data.frame(fread(exon.gtf))
  }
  cds <- subset(gtf, V3 == "CDS")
  names(cds) = c("chr","origin","type","start","end","misc","strand","phase")
  cds <- select(cds, chr, start, end, strand, phase)

  # get UTRs
  UTR <- subset(gtf, V3 == "UTR")
  names(UTR) <- c("chr","origin","type","start","end","misc","strand","phase")
  UTR <- select(UTR, chr, start, end, strand, phase)
  
  # for each junction find the exon that is linked to it
  # output the exon coordinates and the codon phase
  
  cds.grange <- makeGRangesFromDataFrame(cds, keep.extra.columns = TRUE)
  UTR.grange <- makeGRangesFromDataFrame(UTR, keep.extra.columns = TRUE)

  # find upstream_exons
  upstream_exons <- apply( upstream_results_control, MAR = 1, FUN = function(x){
      resolve_junctions( as.character(x[1]),x[2],x[3],x[7], cds.grange, mode = "upstream")  
  })

  upstream_exons <- unlist(GRangesList(upstream_exons)) 
  upstream <- data.frame(chr = seqnames(upstream_exons),
      start = (start(upstream_exons) -1 ), # used to be -1
      end = end(upstream_exons),
      strand = strand(upstream_exons),
      phase =  elementMetadata(upstream_exons)$phase,
      ID = elementMetadata(upstream_exons)$ID
      )

  # find downstream exons
  downstream_exons <- apply( downstream_results_control, MAR = 1, FUN = function(x){
      resolve_junctions( as.character(x[1]),x[2],x[3],x[7], cds.grange, mode = "downstream")  
  })
  downstream_exons <- unlist(GRangesList(downstream_exons)) 
  downstream <- data.frame(chr = seqnames(downstream_exons),
      start = ( start(downstream_exons) -1   ), # used to be - 1
      end = end(downstream_exons),
      strand = strand(downstream_exons),
      phase =  elementMetadata(downstream_exons)$phase,
      ID = elementMetadata(downstream_exons)$ID
      )


  # find UTR exons
  upstream_UTR <- apply( upstream_results_control, MAR = 1, FUN = function(x){
      resolve_junctions( as.character(x[1]),x[2],x[3],x[7], UTR.grange, mode = "upstream")  
  })

  upstream_UTR <- unlist(GRangesList(upstream_UTR)) 
  upstream_UTR <- data.frame(chr = seqnames(upstream_UTR),
        start = (start(upstream_UTR) -1 ), # used to be -1
        end = end(upstream_UTR),
        strand = strand(upstream_UTR),
        phase =  elementMetadata(upstream_UTR)$phase,
        ID = elementMetadata(upstream_UTR)$ID
        )

    # find downstream UTR
  downstream_UTR <- apply( downstream_results_control, MAR = 1, FUN = function(x){
      resolve_junctions( as.character(x[1]),x[2],x[3],x[7], UTR.grange, mode = "downstream")  
  })
  downstream_UTR <- unlist(GRangesList(downstream_UTR)) 
  downstream_UTR <- data.frame(chr = seqnames(downstream_UTR),
      start = ( start(downstream_UTR) -1   ), # used to be - 1
      end = end(downstream_UTR),
      strand = strand(downstream_UTR),
      phase =  elementMetadata(downstream_UTR)$phase,
      ID = elementMetadata(downstream_UTR)$ID
      )

  # if an exon has no entries in upstream or downstream but does in upstream or downstream UTR -> exon (maybe) not in CDS) 
  # if  an exon has ( no upstream exon or upstream CDS exon ) or (no downstream exon or downstream UTR ) -> problematic, evades classification
  UTR.verdicts <- rep(NA, nrow(central_exons))
  
  for( i in 1:nrow(central_exons) ){
    if( 
      (upstream$start[i] == 0 & upstream_UTR$start[i] != 0) |
      ( downstream$start[i] == 0 & downstream_UTR$start[i] != 0) 
      ){
      UTR.verdicts[i] <- "exon in UTR"
    }
    if(
      (upstream$start[i] == 0 & upstream_UTR$start[i] == 0) |
      ( downstream$start[i] == 0 & downstream_UTR$start[i] == 0) 
      ){
      UTR.verdicts[i] <- "exon evades classification"
    }
  }


  # fudge central exon start by 1 due to change of 1-base to 0-base
  # so fucking annoying!
  central_exons$central.start <- central_exons$central.start -1

  # some central exons cannot be matched to CDS sequence
  # this can be due to their presence in ncRNA
  # or they are in UTRs
  # or they're not a cassette - an extension or something
  # does an inclusion or exclusion of this exon change the position of the start codon?


  # now I have exon coordinates for all three I can put them all together in one large table linked by the ID column
  table_list <- list(central = central_exons, upstream = upstream, downstream = downstream)
  fasta_list <- list()

  # first transcibe FASTA for each table
  for( i in 1:length(table_list) ){
    # write table
    table_name <- paste0( outFolder,"/", code, "_", names(table_list)[i], "_fasta.bed" )
    write.table(table_list[[i]], table_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE) 
    fasta_out <- sub(".bed", ".fa", table_name)
    print(table_name)
    fasta.cmd <- paste0("bedtools getfasta -name -fi ",genome.fa," -bed ",table_name ," -fo ",fasta_out)
    system(fasta.cmd)

    fasta <- fread(fasta_out)
    

    fasta <- data.frame(sequence =  fasta[seq(1,nrow(fasta),2),])
    names(fasta) <- names(table_list)[i]

    fasta.name <- paste0(names(table_list)[i], "_fasta")

    assign(fasta.name,fasta)

  #  fasta_list[[i]] <- fasta

  }

  included <- paste0(upstream_fasta$upstream, central_fasta$central, downstream_fasta$downstream)
  skipped <- paste0(upstream_fasta$upstream, downstream_fasta$downstream)

  # make one gigantic table
  central_exons$upstream_phase <- upstream$phase[ match( central_exons$ID, upstream$ID)]
  central_exons$downstream_phase <- downstream$phase[ match( central_exons$ID, downstream$ID)]

  central_exons$included_seq <- included
  central_exons$skipped_seq <- skipped

  central_exons$upstream_seq <- upstream_fasta
  central_exons$central_exon <- central_fasta
  central_exons$downstream_seq <- downstream_fasta


  # remove exons where there is no up and downstream exon found
  clean_exons <- central_exons[ central_exons$upstream_seq != "NN" & central_exons$downstream_seq != "NN" & central_exons$skipped_seq != "NN" & central_exons$strand %in% c("+","-"),]
  # for testing
  #clean_exons <- head(clean_exons,10)

  print("got this far!")
  # F210I skipped exons fail at this entry due to C comparison on a very short downstream exon
  #clean_exons <- central_exons[48,]

  print(head(central_exons))

  print(head(central_exons[,1:8]))
  
  print(head(clean_exons))


  # only need to do once
  res <- lapply(1:nrow(clean_exons), FUN = function(x) translate_toxic( 
                        clean_exons[x,12], clean_exons[x,13],clean_exons[x,14], clean_exons[x,8], clean_exons[x,9], clean_exons[x,6]) )
  results <- do.call( what = rbind, args = res)
  
  if(permute == TRUE){
    # reorder the central exon sequence and run again - 100 times
    set.seed(1234)

    permute_res <- list()
    for( i in 1:N_PERMUTE ){
      # this shuffles the central exons vector
      #permute_exons <- clean_exons$central_exon[sample(1:nrow(clean_exons), size = nrow(clean_exons), replace = FALSE ),]
      # this instead creates random sequences of the same length
      permute_exons <- unname( sapply( clean_exons$central_exon[[1]], FUN = function(y){
           paste( sample( c("A","C","G","T"), size = str_length(y), replace = TRUE ), collapse = "" )
        } ) )


      permute_res[[i]] <- lapply(1:nrow(clean_exons), FUN = function(x) translate_toxic( 
                                  clean_exons[x,12], permute_exons[x], clean_exons[x,14], clean_exons[x,8], clean_exons[x,9], clean_exons[x,6]) )
      # rbind each run and then rbind all runs
      }
    permute_results <- do.call( rbind, lapply( 1:length(permute_res), FUN = function(x) do.call( rbind, permute_res[[x]] ) ) )
  }



  illegal.combos <- list( 
    c(T,F,T,T,T),
    c(T,T,F,T,T),
    c(T,F,F,F,F),
    c(T,F,F,T,T),
    c(T,F,T,T,F),
    c(T,F,F,T,F),
    c(T,F,F,T,F),
    c(T,F,F,F,T)
    )

  possible.combos <- list(
    c(T,T,T,T,T), # both isoforms functional 
    c(T,T,T,T,F), # B benign, frameshifts C to benign new function
    c(T,T,F,F,F), # B benign, frameshifts C from poison to poison # I checked and this is probably due to noisy splicing
    c(T,T,T,F,F), # B benign, frameshfts C from functional to poison
    c(T,F,T,F,T), # B toxic, doesn't frameshift C
    c(T,T,F,F,T), # B benign, C is toxic, doesn't frameshift C
    c(T,T,F,T,F), # B benign, frameshifts C from toxic to functional ( B is probably constitutive )
    c(T,F,T,F,F),  # B toxic, frameshifts C (from functional to toxic)
    c(T,T,T,F,T), # B contains PTC that escapes NMD
    c(F,F,F,F,F) # B frameshifts C to create an NMD-escaping stop codon
    )
  combo.names <- c(
    "Both transcripts functional",
    "Both transcripts function with different downstream protein sequences",
    "Both transcripts non-functional",
    "Included exon causes downstream frameshift",
    "Included exon contains PTC, does not frameshift",
    "Exclusion transcript contains PTC, Included exon benign but doesn't alter PTC",
    "Exclusion transcript contains PTC, Inclusion transcript is functional",
    "Included exon contains PTC and frameshifts",
    "Included exon has PTC that may escape NMD",
    "Included exon frameshifts to create a downstream potentially NMD-escaping stop codon",
    "Evades classification"
    )

  group.names <- c(
    "Transcript switch is neutral",
    "Transcript switch is neutral",
    "Transcript switch is neutral",
    "Inclusion destabilises transcript",
    "Inclusion destabilises transcript",
    "Transcript switch is neutral",
    "Skipping destabilises transcript",
    "Inclusion destabilises transcript",
    "Inclusion destablises transcript (possibly non-NMD)",
    "Inclusion destablises transcript (possibly non-NMD)",
    "Evades classification"
    )

  frameshift.names <- c(
    FALSE,
    TRUE,
    NA,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    FALSE,
    TRUE,
    NA
    )




# for checking whether tht combination
#apply( results, MAR = 1, FUN = function(x) all( x = illegal) )
# does that combination appear?
  illegals <- lapply( illegal.combos, FUN = function(combo) any(apply( results, MAR = 1, FUN = function(x) all( x == combo) )) )

  verdicts <- lapply( possible.combos, FUN = function(combo){ 
    sum( apply( results, MAR = 1, FUN = function(x) all( x == combo) ) )
    })

  # tack on illegal combos
  verdicts[[11]] <-  sum(unlist(illegals))

  verdicts <- data.frame( verdict = combo.names, count = as.numeric( do.call( what = rbind, verdicts) ), group = group.names, causes.frameshift = frameshift.names, stringsAsFactors = FALSE)

  #not_in_CDS <- c( "Exon not in CDS", nrow(central_exons) - nrow(clean_exons), "Exon not in coding sequence" )
  #verdicts <- rbind( not_in_CDS, verdicts)

  verdicts$count <- as.numeric(as.character(verdicts$count))
  verdicts$prop <- verdicts$count / sum(verdicts$count)

  verdicts.table <- paste0(outFolder,"/", code, "_", exon_code[exon_type], "_results.tab")

  print( paste0("writing ", verdicts.table) )

  write.table( verdicts, verdicts.table, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # match the exons with their verdicts
  #results$verdict <- 

  print("saving image for debugging")
  save.image( paste0( outFolder, "/all_objects_not_final.Rdata") )


  if( permute == TRUE ){
    permute_verdicts <- lapply( possible.combos, FUN = function(combo){ 
                          sum( apply( permute_results, MAR = 1, FUN = function(x) all( x == combo) ) )
    }) 

    permute_verdicts <- data.frame( verdict = combo.names, count = as.numeric( do.call( what = rbind, permute_verdicts) ), group = group.names, causes.frameshift = frameshift.names, stringsAsFactors = FALSE)

    #not_in_CDS <- c( "Exon not in CDS", nrow(central_exons) - nrow(clean_exons), "Exon not in coding sequence" )
    #verdicts <- rbind( not_in_CDS, verdicts)

    permute_verdicts$count <- as.numeric(as.character(permute_verdicts$count))
    permute_verdicts$prop <- permute_verdicts$count / sum( permute_verdicts$count)


    permute_verdicts.table <- paste0(outFolder,"/", code, "_", exon_code[exon_type], "_results_permuted",N_PERMUTE,"X.tab")

    print( paste0("writing ", permute_verdicts.table) )

    write.table( permute_verdicts, permute_verdicts.table, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  


  }

  # create results table for each exon

  results_list <- list()
  #apply(results_list, 1, FUN = function(x){
  
  verdict.list <- sapply( 1:nrow(results), FUN = function(x){
                   # check illegal combos first
                   illegals <- which(sapply(illegal.combos, FUN = function(y){
                      all(y == results[x,] ) 
                   }) )
                   if( length(illegals) > 0 ){
                    return(11)
                  }

                   which(sapply(possible.combos, FUN = function(y){
                      all(y == results[x,] ) 
                   }) )
                   # if illegal combo
  }) 
  
  # what about integer(0)s in the verdict.list?

  clean_exons$verdict <- combo.names[ unlist(verdict.list) ]

  # re-load the initial file without any filtering

  central_exons <- exon_list[[ exon_type ]]
  names(central_exons) <- c("chr","central.start","central.end", "ensemblID","p-value","strand")
  central_exons$ID <- paste(central_exons$chr, central_exons$central.start, central_exons$central.end, sep = ":")
  central_exons$verdict <- clean_exons$verdict[ match( paste0( central_exons$ensemblID, central_exons$p.value), paste0( clean_exons$ensemblID, clean_exons$p.value) )]

  central_exons[ is.na(central_exons$verdict) ,]$verdict <- "Exon not in CDS"

  central_exons$UTR.verdict <- UTR.verdicts

  central_exons$group <- verdicts$group[ match(central_exons$verdict, verdicts$verdict)]

  out_table <- paste0( outFolder, "/", code, "_", exon_code[exon_type], "_individual_exon_results.tab")

  print( "saving individual exon verdicts")
  write.table( central_exons, out_table, sep = "\t", row.names = FALSE)

  # if( permute == TRUE ){
  #   print("chi-squared test comparing the cryptic exon category frequencies with the random sequence percentages")

  #   chis <- chisq.test( verdicts$count, permute_verdicts$count )
  #   print(chis)

  #   print("saving image for debugging")
  #   save.image( paste0( outFolder, "/all_objects.Rdata") )

  # }

}
