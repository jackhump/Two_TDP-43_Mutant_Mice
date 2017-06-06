# Entire intron coverage created by whole_intron_coverage.py
# new method - using coverage across the entire intron
library("smoother")
library(ggplot2)
library(stringr)
library(optparse)
# Jack Humphrey 

smoothing <- 10
flank = 300
intron_length = 100
exon_length = 100
exon_height <- 0.1
# This needs to be done three times:
# Skipped exons
# Included exons
# Control exons
options(echo=TRUE)

opt <- parse_args(
  OptionParser(option_list=list(
    make_option( "--included", type="character", default=NULL, "coverage over the introns encompassing included exons"),
    make_option( "--skipped",type="character", default=NULL, "coverage over the introns encompassing skipped exons"),
    make_option( "--control",type="character", default=NULL, "coverage over a set of introns containing non-regulated exons"),
    make_option( "--skipped_exons",type="character", default=NULL, "a bed file of the skipped exons"),
    make_option( "--included_exons",type="character", default=NULL, "a bed file of the included exons"),
    make_option( "--control_exons",type="character", default=NULL, "a bed file of the control exons"),
    make_option( "--code",type="character", default=NULL, help = "the same dataset-specific code used throughout the pipeline"),
    make_option( "--outFolder",type="character", default=NULL, help = "where you want the plot to go")
  )
  ))


outFolder <- opt$outFolder
included <- opt$included
skipped <- opt$skipped
control <- opt$control
skipped_exons <- opt$skipped_exons
included_exons <- opt$included_exons
control_exons <- opt$control_exons
code <- opt$code

exists <- 0
for( file in c(skipped,included,control,skipped_exons,included_exons,control_exons)){
  if(!file.exists(file)){
    print(paste0(file," does not exist"))
    file <- file + 1
  }
}

stopifnot(exists == 0)


# FUNCTIONS

scale_exon <- function(x, len=100){
  get_clusters <- function( clusters ){
    # takes a vector produced by which(x == 1) 
    cluster_start <- c()
    cluster_length <- c()
    for( i in 1:length(clusters) ){
      if(i==1){
        n_cluster <- 1
        cluster_start[n_cluster] <- clusters[i]
        cluster_length[n_cluster] <- 1
        clu <- clusters[i]
      }
      if( clusters[i] == clu + 1){
        cluster_length[n_cluster] <- cluster_length[n_cluster] + 1
        clu <- clu + 1
        next
      }
      if( clusters[i] > clu + 1){
        n_cluster <- n_cluster + 1
        cluster_start[n_cluster] <- clusters[i]
        cluster_length[n_cluster] <- 1
        clu <- clusters[i]
      }
    }
    return( list(n_clusters = n_cluster, start = cluster_start, length = cluster_length))
  }
  # initialise new empty vector
  y <- rep(0, len)
  # x is a vector of 0 and 1 encoding the positions of clusters along a sequence
  # if x is all zero then return resized zero vector
  if( length( which(x==1) ) == 0 ){
    return(y)
  }
  L <- length(x)
  clusters <- get_clusters( which(x>0) )  
  #print(clusters)
  # scale each cluster
  
  for( i in 1:clusters$n_clusters){
    s <- clusters$start[i]
    e <- clusters$start[i] + ( clusters$length[i] - 1 )
    l <- clusters$length[i]
    pos <- (s-1) / L
    S <- (pos * len) + 1
    E <- ( e / L) * len
    y[S:E] <- 1
  }
  
  return(y)
}

create_intron_df <- function(coverage){
  # create data frame of introns in the coverage list
  introns <- unlist(lapply(coverage, FUN = function(x) return(x[[1]]) ))
  names(introns) <- NULL
  split <- str_split_fixed(introns, ":", 2)
  intron_df <- data.frame(
    chr = split[,1],
    start = as.numeric( str_split_fixed(split[,2], "-", 2 )[,1] ),
    end = as.numeric( str_split_fixed(split[,2], "-", 2)[,2] ),
    stringsAsFactors = FALSE
  )
  return(intron_df)
}

create_matched_exons <- function(exons, intron_df){
  # match exons with introns - pick one exon per intron
  matched_exons <- apply( intron_df, MAR = 1, FUN = function(x){
    exon <-  exons[exons$V1 == x[[1]] & exons$V2 > as.numeric( x[[2]]) & exons$V3 < as.numeric(x[[3]]) ,]
    if( nrow(exon) > 1){
      exon <- exon[1,]
    }
    if( nrow(exon) == 0 ){
      exon <- exons[1,]
      exon[,1:ncol(exon)] <- NA
    }
    return(exon)
  } )
  matched_exons <- do.call(rbind, args = matched_exons)
  names(matched_exons)[1:3] <- c("chr", "start", "end")
  
  return(matched_exons)
}

smooth_scaled_coverage <- function( coverage, intron_df, matched_exons, num.exons ){
  coverage <- lapply( coverage, FUN = function(x) as.numeric(x[2:length(x)]))
  scaled_coverage <- lapply(1:length(coverage), FUN = function(i) {
    print(i)
    cov <- coverage[[i]]
    exon <- matched_exons[i,]
    intron <- intron_df[i,]
    # check for sanity
    print(exon$start - intron$start)
    total_length <- exon_length + flank + intron_length + flank + exon_length + flank + intron_length + flank + exon_length 
    if( is.na(exon$chr) ){
      return( rep(0, total_length ))
    }
    stopifnot( length(cov) == intron$end - intron$start) # sanity check
    stopifnot( exon$chr == intron$chr & intron$start < exon$start & intron$end > exon$end) # double check the matching
    l <- length(cov)
    S <- intron$start
    E <- intron$end
    strand <- exon$V6
    #print(strand)
    if( strand == "+"){
      s <- exon$start - S
      e <- exon$end - S
    }else{
      s <- E - exon$end
      e <- E - exon$start
    }
    # extract raw coverage
    upstream_exon_with_flank <- cov[1:(exon_length+flank)]
    # danger of flank being larger than the distance between the start of the exon and the start of the intron
    if( s > flank){
      print("threeSS looks good")
      threeSS_flank <- cov[(s-flank):(s-1)]
    }else{
      threeSS_flank <- rep(0,flank)
      threeSS_flank[(flank-(s-2) ):flank] <- cov[1:(s-1)]
    }
    # exon 62 is partly NA for this
    if( e+flank < l){
      print("fiveSS looks good")
      fiveSS_flank <- cov[(e+1):(e+flank)]
    }else{
      fiveSS_flank <- rep(0,flank)
      fiveSS_flank
    }
    downstream_exon_with_flank <- cov[ (l - (flank + exon_length - 1) ):l]
    # extract scaled coverage
    # introns may not exist due to overlapping flanking areas so in that case add fake coverage
    if( (1+exon_length+flank) < (s-flank-1) ){
      upstream_intron <- scale_exon( cov[(1+exon_length+flank):(s-flank-1)], intron_length)
    }else{
      upstream_intron <- rep(0,intron_length)
      # what to put in here?
    }
    if( (e+flank) < (l-flank)){
      downstream_intron <- scale_exon( cov[(e+flank):(l-flank)], intron_length)
    }else{
      downstream_intron <- rep(0,intron_length)
    }
    exon <- scale_exon( cov[s:e], exon_length)
    
    # paste together all the parts
    total <- c(
      upstream_exon_with_flank,
      upstream_intron,
      threeSS_flank,
      exon,
      fiveSS_flank,
      downstream_intron,
      downstream_exon_with_flank
    )
    return(total)
  } )
  
  scaled <- as.data.frame( do.call( what = rbind, args = scaled_coverage) )
  scaled.count <- colSums(scaled,na.rm = TRUE)
  scaled.norm <- scaled.count / num.exons
  scaled.smooth <- smth(scaled.norm, window = smoothing, method = "gaussian", tails = TRUE)
  
  return(scaled.smooth)
}

# run loop three times 
create_coverage <- function(coverage_file, exons_bed){
  coverage <- readLines(coverage_file)
  coverage <- sapply(coverage, FUN = function(x) str_split(x, ","))
  
  intron_df <- create_intron_df(coverage)
  
  exons <- read.table(exons_bed, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  matched_exons <- create_matched_exons(exons, intron_df)
  num.exons <- nrow(exons)
  
  smooth_coverage <- smooth_scaled_coverage(coverage, intron_df, matched_exons, num.exons )
  return(smooth_coverage)
}

get_exon_number <- function(exons_bed){
  exons <- read.table(exons_bed, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  num_exons <- nrow(exons)
  return(num_exons)
}

# create numeric coverage vectors
# coverage_file <- "/Users/Jack/Google Drive/TDP_paper/RNA_maps/results/F210I_included_100_coverage.csv"
# exons_bed <- "/Users/Jack/SAN/IoN_RNAseq/RNA_Maps/data//noheader/F210I_embryonic_brain_se_intron_included.bed"
# test <- create_coverage( coverage_file, exons_bed)
# 

# coverage_file <- "/Users/Jack/Documents/Misc/intron_coverage.csv"
# exons_bed <- "/Users/Jack/SAN/IoN_RNAseq/RNA_Maps/data/noheader/F210I_embryonic_brain_flank0_se_included.bed"



# run each set through
cov_list <- list()
cov_list[[1]] <- create_coverage(included,included_exons)
cov_list[[2]] <- create_coverage(skipped,skipped_exons)
cov_list[[3]] <- create_coverage(control,control_exons)

# get exon numbers for plot
exon_num_list <- list()
exon_num_list[[1]] <- get_exon_number(included_exons)
exon_num_list[[2]] <- get_exon_number(skipped_exons)
exon_num_list[[3]] <- get_exon_number(control_exons)


results <- paste0(outFolder,"/coverage_data.Rdata")

save.image(results)

quit()

load("/Users/Jack/Google Drive/TDP_paper/RNA_maps/results/coverage_data.Rdata")


# plotting
my_ymax <- max(unlist(lapply(cov_list, max)))

x_breaks <- c( 1+exon_length,
               1+exon_length+flank,
               1+exon_length+flank+intron_length,
               1+exon_length+flank+intron_length+flank,
               1+exon_length+flank+intron_length+flank+exon_length,
               1+exon_length+flank+intron_length+flank+exon_length+flank,
               1+exon_length+flank+intron_length+flank+exon_length+flank+intron_length,
               1+exon_length+flank+intron_length+flank+exon_length+flank+intron_length+flank
)


x_labels <- c("5'",flank, -flank, "3'",  "5'", flank, -flank, "3'" )

centre_point <- 1+exon_length+flank+intron_length+flank+(exon_length/2)
# draw box on plot to represent exon

exon_df <- data.frame(
  x = c( 
    x_breaks[4],
    x_breaks[5],
    x_breaks[5],
    x_breaks[4],
    
    1,
    x_breaks[1],
    x_breaks[1],
    1,
    
    x_breaks[8], 
    x_breaks[8] + exon_length, 
    x_breaks[8] + exon_length, 
    x_breaks[8]
    
  ),
  y = rep( c( exon_height * my_ymax, exon_height * my_ymax, -exon_height * my_ymax, -exon_height * my_ymax), 3 ),
  group = c(
    1,1,1,1,
    2,2,2,2,
    3,3,3,3
  )
)
included_df <- data.frame(
  x = c(
    x_breaks[1],
    x_breaks[5] ), 
  xend = c(
    x_breaks[4],
    x_breaks[8]
  )
)
skipped_df <- data.frame(
  x = x_breaks[1],
  xend = x_breaks[8]
)

p <- ggplot() + theme_classic() +
  # included exons
  geom_area( aes(1:length(cov_list[[1]]), cov_list[[1]] ), color = NA, fill = "red", alpha = 0.9) + 
  # control exons
  geom_area( aes(1:length(cov_list[[3]]), cov_list[[3]] ), color = NA, fill = "gray", alpha = 0.75 ) +
  # skipped exons
  geom_area( aes(1:length(cov_list[[2]]), -1*cov_list[[2]] ), color = NA, fill = "blue", alpha = 0.9) + 
  # control exons
  geom_area( aes(1:length(cov_list[[3]]), -1*cov_list[[3]] ), color = NA, fill = "gray", alpha = 0.75 ) +
  
  scale_x_continuous(
    "",
    breaks = x_breaks,
    label = x_labels) + 
  
  scale_y_continuous(
    "Per-nucleotide normalised iCLIP cluster coverage",
    limits = c(-my_ymax, my_ymax ),
    labels = scales::percent 
  ) +
  labs(title = paste0(code, " iCLIP cluster RNA map")) +
  annotate("text", y = my_ymax, x = (exon_length/2),  label = paste0("Included exons: ", exon_num_list[[1]]), colour = "red" ) +
  annotate("text", y = my_ymax - (1/3*my_ymax), x = (exon_length/2),  label = paste0("Skipped exons: ", exon_num_list[[2]]), colour = "blue" ) +
  annotate("text", y = my_ymax - (2/3*my_ymax), x = (exon_length/2),  label = paste0("Control exons: ", exon_num_list[[3]]), colour = "black" ) +
  # add exon 
  geom_polygon(data=exon_df, aes(x,y, group = group), fill = NA, colour = "black", linetype="dashed", size = 0.75) + 
  
  # add intron lines
  geom_curve( data=included_df, aes(x=x,xend=xend, y=0,yend=0), curvature = -0.3, size = 0.75 ) +
  geom_curve( data=skipped_df, aes(x=x, xend=xend, y=0,yend=0), curvature = 0.25, size = 0.75) +
  # add dashed lines for locations
  geom_vline( data = data.frame(xintercept = x_breaks, slope=0), aes(xintercept=xintercept),colour = "black", linetype="dashed",alpha=0.5 )



ggsave(p,
       units = "in",
       width = 15, 
       height = 10,
       filename = paste0(outFolder, "/", gsub(" ", "_", code), "_RNAMap_whole_intron_scaled.pdf")
)


