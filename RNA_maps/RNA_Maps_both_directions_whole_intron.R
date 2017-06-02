#########################
# WHOLE INTRON RNA MAP
#####################
library(gplots)
library("smoother")
library(ggplot2)
library(stringr)
# Jack Humphrey 
# based on a script by Nejc Haberman

# the visualisation for 5' and 3' were done separately for a single set of exons
# I want to combine both splice sites and plot included and skipped exons on opposite ends of the y axis
# this should be done for F210I and M323K
flank <- 300
exon_flank <- 10
smoothing <- 30 # default is 10
spacer <- 100 # to separate the 3' and 5' exon coverage
exon_height <- 0.1 # proportion of highest peak in plot
fake_exon_size <- 100 # how large should we guess the flanking exons are for coverage estimates?
bar_width <- 66 # how wide should the exon bars be?
##################
#   all the data #
##################

# F210I
included_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_exons/F210I_embryonic_brain_flank0_se_included.bed-3SS-flanked300-clusters-map_positions-merged.csv"
included_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_exons/F210I_embryonic_brain_flank0_se_included.bed-5SS-flanked300-clusters-map_positions-merged.csv"
control_included_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_exons/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv"
control_included_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_exons/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv"
skipped_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_skipped_exons/F210I_embryonic_brain_flank0_se_skipped.bed-3SS-flanked300-clusters-map_positions-merged.csv"
skipped_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_skipped_exons/F210I_embryonic_brain_flank0_se_skipped.bed-5SS-flanked300-clusters-map_positions-merged.csv"
control_skipped_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_skipped_exons/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv"
control_skipped_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_skipped_exons/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv"
included_intron_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_introns/F210I_embryonic_brain_se_intron_included.bed-3SS-flanked300-clusters-map_positions-merged.csv"
included_intron_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_introns/F210I_embryonic_brain_se_intron_included.bed-5SS-flanked300-clusters-map_positions-merged.csv"
control_included_intron_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_introns/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv"
control_included_intron_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_introns/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv"
skipped_intron_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_skipped_introns/F210I_embryonic_brain_se_intron_skipped.bed-3SS-flanked300-clusters-map_positions-merged.csv"
skipped_intron_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_skipped_introns/F210I_embryonic_brain_se_intron_skipped.bed-5SS-flanked300-clusters-map_positions-merged.csv"
control_skipped_intron_3SS <-"/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_skipped_introns/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv"
control_skipped_intron_5SS <-"/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_skipped_introns/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv"
code <- "F210I whole intron"
iCLIP.code <- "F210I WT + M323K WT + M323K MUT + TDP-43 embryonic brain + TDP-43 E18 brain merge"


# M323K
included_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_included_exons/M323K_adult_brain_flank0_se_included.bed-3SS-flanked300-clusters-map_positions-merged.csv"
included_5SS  <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_included_exons/M323K_adult_brain_flank0_se_included.bed-5SS-flanked300-clusters-map_positions-merged.csv"
control_included_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_included_exons/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv"
control_included_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_included_exons/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv"
skipped_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_skipped_exons/M323K_adult_brain_flank0_se_skipped.bed-3SS-flanked300-clusters-map_positions-merged.csv"
skipped_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_skipped_exons/M323K_adult_brain_flank0_se_skipped.bed-5SS-flanked300-clusters-map_positions-merged.csv"
control_skipped_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_skipped_exons/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv"
control_skipped_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_skipped_exons/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv"
included_intron_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_included_introns/M323K_adult_brain_se_intron_included.bed-3SS-flanked300-clusters-map_positions-merged.csv"
included_intron_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_included_introns/M323K_adult_brain_se_intron_included.bed-5SS-flanked300-clusters-map_positions-merged.csv"
control_included_intron_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_included_introns/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv"
control_included_intron_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_included_introns/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv"
skipped_intron_3SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_skipped_introns/M323K_adult_brain_se_intron_skipped.bed-3SS-flanked300-clusters-map_positions-merged.csv"
skipped_intron_5SS <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_skipped_introns/M323K_adult_brain_se_intron_skipped.bed-5SS-flanked300-clusters-map_positions-merged.csv"
control_skipped_intron_3SS <-"/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_skipped_introns/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv"
control_skipped_intron_5SS <-"/SAN/vyplab/IoN_RNAseq/RNA_Maps/M323K/M323K_skipped_introns/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv"
code <- "M323K whole intron"
iCLIP.code <- "F210I WT + M323K WT + M323K MUT + TDP-43 embryonic brain + TDP-43 E18 brain merge"
###################
# regulated exons #
###################

# threeSS.csv="/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/F210I_embryonic_brain_se_sig_nooffset.bed-3SS-flanked300-clusters-map_positions-merged.csv" 
# fiveSS.csv="/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/F210I_embryonic_brain_se_sig_nooffset.bed-5SS-flanked300-clusters-map_positions-merged.csv"
# control.threeSS="/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/random_exons.bed-3SS-flanked300-clusters-map_positions-merged.csv" 
# control.fiveSS="/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/random_exons.bed-5SS-flanked300-clusters-map_positions-merged.csv" 

outFolder="/Users/Jack/Google Drive/TDP_paper/RNA_maps/"

# dir <- "/SAN/vyplab/IoN_RNAseq/RNA_Maps/F210I/F210I_included_exons/"
# exon_stem <- "F210I_embryonic_brain_flank0_se_included.bed"
# control_stem <- "random_exons.bed"
# 
# threeSS.csv <- paste0(dir, exon_stem, "-3SS-flanked300-clusters-map_positions-merged.csv" )
# fiveSS.csv <- paste0(dir, exon_stem, "-5SS-flanked300-clusters-map_positions-merged.csv" )
# control.threeSS <- paste0(dir, control_stem, "-3SS-flanked300-clusters-map_positions-merged.csv")
# control.fiveSS <- paste0( dir, control_stem, "-5SS-flanked300-clusters-map_positions-merged.csv")
# 
files_to_process <- c(
                      included_3SS,
                      included_5SS,
                      control_included_3SS,
                      control_included_5SS,
                      skipped_3SS,
                      skipped_5SS,
                      control_skipped_3SS,
                      control_skipped_5SS
)
intron_files_to_process <- c(
                      included_intron_3SS,
                      included_intron_5SS,
                      control_included_intron_3SS,
                      control_included_intron_5SS,
                      skipped_intron_3SS,
                      skipped_intron_5SS,
                      control_skipped_intron_3SS,
                      control_skipped_intron_5SS
                      )
# 
# 
# output1=paste0(output, "coverage_3ss.pdf") 
# output2=paste0(output, "coverage_5ss.pdf") 
# output3=paste0(output, "HeatMap.pdf") 

# if doing this through RStudio then switch the paths 
# if( .Platform$GUI == "RStudio"){
#   for( file in c("threeSS.csv","fiveSS.csv","control.threeSS","control.fiveSS")){
#     var <- get(file)
#     var <- gsub( "/SAN/vyplab", "/Users/Jack/SAN/", var )
#     print(file.exists(var))
#     assign(file,var)
#   }
# }


# set the boundaries 
fiveSS_boundary <- c( ( 7 + (3 * flank) - exon_flank ) : ( 7 + (4 * flank)  )   )
threeSS_boundary <- c( 4:(flank + exon_flank + 4) )
# intron boundaries
fiveSS_intron_boundary <- c( (5 + (2*flank) ):( 5+(3*flank)+exon_flank )   )
threeSS_intron_boundary <- c( (4 + flank - exon_flank): (4 + (2*flank)) )

#################
# functions     #
#################

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

#get_clusters( c(2,3,4,5,8,9,10,11,100,110,111,112))



scale_exons <- function(x, len=100){
  # initialise new empty vector
  y <- rep(0, len)
  # x is a vector of 0 and 1 encoding the positions of clusters along a sequence
  # if x is all zero then return resized zero vector
  if( length( which(x==1) ) == 0 ){
    return(y)
  }
  L <- length(x)
  clusters <- get_clusters( which(x==1) )  
  print(clusters)
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



smooth_coverage <- function( exon_3SS, exon_5SS, spacer ){
  # check if files exist
  if( all( !file.exists( c(exon_3SS, exon_5SS) ))){
    exon_3SS <- gsub( "/SAN/vyplab", "/Users/Jack/SAN/", exon_3SS )
    exon_5SS <- gsub( "/SAN/vyplab", "/Users/Jack/SAN/", exon_5SS )
  }
  # try substituting the mount name in 
  for( file in c(exon_3SS, exon_5SS)){
    if( !file.exists(file)){
      print(paste0(file," does not exist"))
    }
  }
  
  # import and merge of 3SS with 5SS
  exon3SSinput <- read.table(exon_3SS, header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
  exon5SSinput <- read.table(exon_5SS, header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
  exon3SS.5SSA <- merge(exon3SSinput, exon5SSinput, by="V1", all.x=TRUE)
  exon3SS.5SSA <- exon3SS.5SSA[order(exon3SS.5SSA$V2.x),]
  
  num.exons <- nrow(exon3SS.5SSA)
  
  # removal of transcripts with clusters outside of heat map region
  exon3SS.5SS <- exon3SS.5SSA[complete.cases(exon3SS.5SSA),]

  # get exon length 
  exon_length <- data.frame(str_split_fixed(str_split_fixed(exon3SSinput$V1, ":",2)[,2], ":", 2), stringsAsFactors = FALSE)
  exon_length <- as.numeric(exon_length$X2) - as.numeric(exon_length$X1)
  
  # extract coverage information across the exon
  exon_cov_list <- list()
  for( exon in 1:length(exon_length)){
    #print(exon)
    if( exon_length[exon] > flank ){
      exon_length[exon] <- flank
    }
    exon_cov_list[[exon]] <- exon3SSinput[ exon, (4 + flank):(4+flank+exon_length[exon] ) ]
  }
  
  # two options: are there iCLIP clusters in the exon?
  proportion_exons <- as.logical( lapply( exon_cov_list, FUN = max ) ) 
  # create proportion
  proportion_exons <- sum(proportion_exons) / length(proportion_exons)
  
 
  
  
  # keep both upstream and downstream intron in the same object
  input <- exon3SS.5SS[, c(threeSS_boundary, fiveSS_boundary)]
  
  # set all nonzero counts to one - otherwise exons with many iCLIP clusters will bias the plot
  input[ input[,] < 0 ]  <- 0
  input[ input[,] > 0 ]  <- 1
  
  # normalise and smooth the sum of each position
  input.sum <- colSums(input, na.rm = FALSE, dims = 1)
  input.norm <- input.sum / num.exons
  input.smooth <- smth(input.norm, window = smoothing, method = "gaussian", tails = TRUE) # tails option keeps the tail ends of the data which are probably out a bit
  
  #print(length(input.smooth) )
  
  # add spacer of 0 for plotting
  input.spaced <- rep(0, length(input.smooth) )
  input.spaced[1:(1+flank+exon_flank)] <- input.smooth[1:(1+flank+exon_flank)]
  input.spaced[ (1+flank+exon_flank+spacer):( (2*flank)+(2*exon_flank) + spacer + 2 )] <- input.smooth[ (1+flank+exon_flank):( (2*flank)+(2*exon_flank) + 2) ]

  
  # get ymax
  ymax.exons <- max(input.smooth, na.rm = TRUE)

  # return smoothed coverage and ymax to list
  return( list( input.spaced, ymax.exons, num.exons, proportion_exons))
}

###############


smooth_intron_coverage <- function(exon_3SS, exon_5SS, spacer ){
  # a way of hacking together intron coverage by treating the introns as large exons
  # Nejc's script thus finds coverage around the start and end of the intron
  # the intron start is here called exon_3SS and the intron end is called exon_5SS - actually opposite way round
  if( all( !file.exists( c(exon_3SS, exon_5SS) ))){
    exon_3SS <- gsub( "/SAN/vyplab", "/Users/Jack/SAN/", exon_3SS )
    exon_5SS <- gsub( "/SAN/vyplab", "/Users/Jack/SAN/", exon_5SS )
  }
  # try substituting the mount name in 
  for( file in c(exon_3SS, exon_5SS)){
    if( !file.exists(file)){
      print(paste0(file," does not exist"))
    }
  }
  
  # import and merge of 3SS with 5SS
  exon3SSinput <- read.table(exon_3SS, header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
  exon5SSinput <- read.table(exon_5SS, header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
  exon3SS.5SSA <- merge(exon3SSinput, exon5SSinput, by="V1", all.x=TRUE)
  exon3SS.5SSA <- exon3SS.5SSA[order(exon3SS.5SSA$V2.x),]
  
  num.exons <- nrow(exon3SS.5SSA)
  
  # removal of transcripts with clusters outside of heat map region
  exon3SS.5SS <- exon3SS.5SSA[complete.cases(exon3SS.5SSA),]
  
  #keep both upstream and downstream intron in the same object
  input <- exon3SS.5SS[, c(threeSS_intron_boundary, fiveSS_intron_boundary)]
  
  # set all nonzero counts to one - otherwise exons with many iCLIP clusters will bias the plot
  input[ input[,] < 0 ]  <- 0
  input[ input[,] > 0 ]  <- 1
  
  # normalise and smooth the sum of each position
  input.sum <- colSums(input, na.rm = FALSE, dims = 1)
  input.norm <- input.sum / num.exons
  input.smooth <- smth(input.norm, window = smoothing, method = "gaussian", tails = TRUE) # tails option keeps the tail ends of the data which are probably out a bit
  
  #print(length(input.smooth) )
  
  # add spacer of 0s for plotting - introns need a lot of space!
  intron_spacer <- (5*spacer) + (4*flank) + (4*exon_flank) + 4
  input.spaced <- rep(0, intron_spacer )
  input.spaced[(1+spacer):(1+flank+exon_flank+spacer)] <- input.smooth[1:(1+flank+exon_flank)] # 5' end of the intron
  input.spaced[ ( (4*spacer) + (3*flank) + (3*exon_flank) + 4  ):( (4*spacer) + (4*flank) +  (4*exon_flank) + 4 )] <- input.smooth[ (2+flank+exon_flank):( (2*flank)+(2*exon_flank) + 2) ]
  
  # get iCLIP coverage within flanking exons
  upstream_exon_cov_list <- list()
  downstream_exon_cov_list <- list()
  for( exon in 1:num.exons){
    #print(exon)
    upstream_exon_cov_list[[exon]] <- exon3SSinput[ exon, (4 + flank - fake_exon_size):(4+flank) ]
    downstream_exon_cov_list[[exon]] <- exon5SSinput[ exon, (4 + flank):(4 + flank + fake_exon_size) ] 
  }
  
  # two options: are there iCLIP clusters in the exon?
  upstream_exon_prop <- as.logical( lapply( upstream_exon_cov_list, FUN = max))
  downstream_exon_prop <- as.logical( lapply( downstream_exon_cov_list, FUN = max ) ) 
  # create proportion
  upstream_exon_prop <- sum(upstream_exon_prop) / length(upstream_exon_prop)
  downstream_exon_prop <- sum(downstream_exon_prop) / length(downstream_exon_prop)
  
  # get ymax
  ymax.exons <- max(input.smooth, na.rm = TRUE)
  
  # return smoothed coverage and ymax to list
  return( list( input.spaced, ymax.exons, num.exons, upstream_exon_prop, downstream_exon_prop))
  
}


###################
# Do Calculations #
####################


# loop through every pair of files
cov_list <- list()
ymax_list <- list()
n_exon_list <- list()
prop_exon_list <- list()
i <- 1
for( pair in seq(1, length(files_to_process), 2) ){
  print(i)
  print(paste(files_to_process[pair], files_to_process[pair+1]) )
  data <- smooth_coverage( files_to_process[pair], files_to_process[pair+1], spacer )
  cov_list[[i]] <- data[[1]]
  ymax_list[[i]] <- data[[2]]
  n_exon_list[[i]] <- data[[3]]
  prop_exon_list[[i]] <- data[[4]]
  i <- i + 1
}


# do the same for the introns
intron_cov_list <- list()
intron_ymax_list <- list()
upstream_prop_list <- list()
downstream_prop_list <- list()
i <- 1
for( pair in seq(1, length(intron_files_to_process), 2) ){
  print(i)
  print(paste(intron_files_to_process[pair], intron_files_to_process[pair+1]) )
  data <- smooth_intron_coverage( intron_files_to_process[pair], intron_files_to_process[pair+1], spacer )
  intron_cov_list[[i]] <- data[[1]]
  intron_ymax_list[[i]] <- data[[2]]
  upstream_prop_list[[i]] <- data[[4]]
  downstream_prop_list[[i]] <- data[[5]]
  i <- i + 1
}



merge_exons_into_intron <- function( exon_cov, intron_cov ){
  exon_start <- (2*spacer) + exon_flank + flank + 1 
  exon_end <- (3*spacer) + (3*flank) + (3*exon_flank) + 2
  intron_cov[exon_start:exon_end] <- exon_cov
  return(intron_cov)
  }

merged_list <- list()
for( i in 1:length(cov_list)){
  merged_list[[i]] <- merge_exons_into_intron( cov_list[[i]], intron_cov_list[[i]] )
}

my_ymax <- max(c( unlist(ymax_list),unlist(intron_ymax_list) ) )
# now have a list of coverage for:
## 1. included exons 2. control 3. skipped exons 4. control
# create numerical scale
#x_scale <- c(-(flank):exon_flank, -exon_flank:flank)

x_breaks <- c( 1+spacer+exon_flank,
              flank + 1 + spacer + exon_flank , 
              (2*spacer) + flank + exon_flank + 1,
              (2*spacer) + (2*flank) + exon_flank + 1,
              (3*spacer) + (2*flank) + (3*exon_flank) + 3,
              (3*spacer) + (3*flank) + (3*exon_flank) + 3,
              (4*spacer) + (3*flank) + (3*exon_flank) + 3,
              (4*spacer) + (4*flank) + (3*exon_flank) + 4
            )
x_labels <- c("5'",flank, -flank, "3'",  "5'", flank, -flank, "3'" )

centre_point <- 2 + (2*flank) + (2*spacer) + (2*exon_flank) + (spacer/2) 
# draw box on plot to represent exon

exon_df <- data.frame(
    x = c( 
      x_breaks[4],
      x_breaks[5],
      x_breaks[5],
      x_breaks[4],
      
      1,x_breaks[1],x_breaks[1],1,
      
      x_breaks[8], x_breaks[8] + spacer, x_breaks[8] + spacer, x_breaks[8]
      
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
  geom_area( aes(1:length(merged_list[[1]]), merged_list[[1]] ), color = NA, fill = "red", alpha = 0.9) + 
  # control exons
  geom_area( aes(1:length(merged_list[[2]]), merged_list[[2]] ), color = NA, fill = "gray", alpha = 0.75 ) +
  # skipped exons
  geom_area( aes(1:length(merged_list[[3]]), -1*merged_list[[3]] ), color = NA, fill = "blue", alpha = 0.9) + 
  # control exons
  geom_area( aes(1:length(merged_list[[4]]), -1*merged_list[[4]] ), color = NA, fill = "gray", alpha = 0.75 ) +
  scale_x_continuous(
    "",
    breaks = x_breaks,
    label = x_labels) + 
  
  scale_y_continuous(
    "Per-nucleotide normalised iCLIP cluster coverage",
    limits = c(-max( ymax_list[[3]], intron_ymax_list[[3]] ), max( ymax_list[[1]], intron_ymax_list[[1]] ) ),
    labels = scales::percent,
    sec.axis = sec_axis(~./ymax_list[[1]],
                        name = "Per-exon iCLIP cluster coverage", 
                        labels = scales::percent)) +
  labs(title = paste0(code, " iCLIP cluster RNA map"), caption = paste0("iCLIP peaks from ",iCLIP.code )) +
  annotate("text", y = max( ymax_list[[1]], intron_ymax_list[[1]] ), x = (spacer/2),  label = paste0("Included exons: ", n_exon_list[[1]]), colour = "red" ) +
  annotate("text", y = -max( ymax_list[[3]], intron_ymax_list[[3]] ), x = (spacer/2),  label = paste0("Skipped exons: ", n_exon_list[[3]]), colour = "blue" ) +

  # add exon coverage proportion
  geom_bar( aes( x = centre_point, y = prop_exon_list[[1]] * ymax_list[[1]]), stat = "identity" , fill = "red", alpha = 0.9, width = bar_width, colour = "black") +
  geom_bar( aes( x = centre_point, y = prop_exon_list[[2]] * ymax_list[[1]]), stat = "identity" , fill = "gray", alpha = 0.75, width = bar_width, colour = "black") +
  geom_bar( aes( x = centre_point, y = -prop_exon_list[[3]] * ymax_list[[1]]), stat = "identity" , fill = "blue", alpha = 0.9, width = bar_width, colour = "black") +
  geom_bar( aes( x = centre_point, y = -prop_exon_list[[4]] * ymax_list[[1]]), stat = "identity" , fill = "gray", alpha = 0.75, width = bar_width, colour = "black") +

  # flanking exons
  geom_bar( aes( x = (spacer/2), y = upstream_prop_list[[1]] * ymax_list[[1]]), stat = "identity" , fill = "red", alpha = 0.9, width = bar_width, colour = "black") +
  geom_bar( aes( x = (spacer/2), y = upstream_prop_list[[2]] * ymax_list[[1]]), stat = "identity" , fill = "gray", alpha = 0.75, width = bar_width, colour = "black") +
  geom_bar( aes( x =  (spacer/2), y = -upstream_prop_list[[3]] * ymax_list[[1]]), stat = "identity" , fill = "blue", alpha = 0.9, width = bar_width, colour = "black") +
  geom_bar( aes( x =  (spacer/2), y = -upstream_prop_list[[4]] * ymax_list[[1]]), stat = "identity" , fill = "gray", alpha = 0.75, width = bar_width, colour = "black") +
    
  geom_bar( aes( x = x_breaks[8] + (spacer/2), y = downstream_prop_list[[1]] * ymax_list[[1]]), stat = "identity" , fill = "red", alpha = 0.9, width = bar_width, colour = "black") +
  geom_bar( aes( x = x_breaks[8] + (spacer/2), y = downstream_prop_list[[2]] * ymax_list[[1]]), stat = "identity" , fill = "gray", alpha = 0.75, width = bar_width, colour = "black") +
  geom_bar( aes( x = x_breaks[8] + (spacer/2), y = -downstream_prop_list[[3]] * ymax_list[[1]]), stat = "identity" , fill = "blue", alpha = 0.9, width = bar_width, colour = "black") +
  geom_bar( aes( x = x_breaks[8] + (spacer/2), y = -downstream_prop_list[[4]] * ymax_list[[1]]), stat = "identity" , fill = "gray", alpha = 0.75, width = bar_width, colour = "black") + 
  
  # add exon 
  geom_polygon(data=exon_df, aes(x,y, group = group), fill = NA, colour = "black", linetype="dashed", size = 0.75) + 
  
  # add intron lines
  geom_curve( data=included_df, aes(x=x,xend=xend, y=0,yend=0), curvature = -0.3, size = 0.75 ) +
  geom_curve( data=skipped_df, aes(x=x, xend=xend, y=0,yend=0), curvature = 0.25, size = 0.75) 
  


ggsave(p,
       units = "in",
       width = 15, 
       height = 10,
       filename = paste0(outFolder, "/", gsub(" ", "_", code), "_",gsub(" ","_",iCLIP.code), "_RNAMap_whole_intron.pdf")
       )

# including exon coverage
# 
# 
# p2 <- ggplot() + theme_classic() +
#   # included exons
#   geom_area( aes(1:length(cov_list[[1]]), cov_list[[1]] ), color = NA, fill = "red", alpha = 0.9) + 
#   # control exons
#   geom_area( aes(1:length(cov_list[[2]]), cov_list[[2]] ), color = NA, fill = "gray", alpha = 0.75 ) +
#   # skipped exons
#   geom_area( aes(1:length(cov_list[[3]]), -1*cov_list[[3]] ), color = NA, fill = "blue", alpha = 0.9) + 
#   # control exons
#   geom_area( aes(1:length(cov_list[[4]]), -1*cov_list[[4]] ), color = NA, fill = "gray", alpha = 0.75 ) +
#   scale_x_continuous(
#     "",
#     breaks = x_breaks,
#     label = x_labels) + 
#   
#   scale_y_continuous(
#     "Per-nucleotide normalised iCLIP cluster coverage",
#     limits = c(-ymax_list[[3]], ymax_list[[1]]),
#     labels = scales::percent,
#     sec.axis = sec_axis(~./ymax_list[[1]],
#                name = "Per-exon iCLIP cluster coverage", 
#                labels = scales::percent)) +
#   # add exon coverage bars
#   geom_bar( aes( x = centre_point, y = prop_exon_list[[1]] * ymax_list[[1]]), stat = "identity" , fill = "red", alpha = 0.9, width = 30, colour = "black") +
#   geom_bar( aes( x = centre_point, y = prop_exon_list[[2]] * ymax_list[[1]]), stat = "identity" , fill = "gray", alpha = 0.75, width = 30, colour = "black") +
#   geom_bar( aes( x = centre_point, y = -prop_exon_list[[3]] * ymax_list[[1]]), stat = "identity" , fill = "blue", alpha = 0.9, width = 30, colour = "black") +
#   geom_bar( aes( x = centre_point, y = -prop_exon_list[[4]] * ymax_list[[1]]), stat = "identity" , fill = "gray", alpha = 0.75, width = 30, colour = "black") +
#   # add exon 
#   geom_polygon(data=exon_df, aes(x,y), fill = NA, colour = "black", linetype="dashed", size = 0.75) + 
#   # add intron lines
#   geom_curve( data=included_df, aes(x=x,xend=xend, y=0,yend=0), curvature = -0.3, size = 0.75 ) +
#   geom_curve( data=skipped_df, aes(x=x, xend=xend, y=0,yend=0), curvature = 0.25, size = 0.75) +
#   labs(title = paste0(code, " iCLIP cluster RNA map"), caption = paste0("iCLIP peaks from ",iCLIP.code )) +
#   annotate("text", y = ymax_list[[1]], x = 100,  label = paste0("Included exons: ", n_exon_list[[1]]), colour = "red" ) +
#   annotate("text", y = -ymax_list[[3]], x = 100,  label = paste0("Skipped exons: ", n_exon_list[[3]]), colour = "blue" )
# 
# 
# 
# 
# 
# 


# 
# 
# # ggplot cluster density
# 
# 
# gg.3SS <- ggplot() + theme_bw() + 
#   geom_area(aes(c(-300:50), sum.3SSA.norm.smooth), color = "grey", fill = "red", alpha = 0.3) + 
#   geom_area(aes(c(-300:50), sum.control.3SSA.norm.smooth), color="grey", alpha=0.3) + 
#   geom_vline(xintercept = 0, alpha = 0.5, linetype="longdash") +
#   ggtitle("Cluster coverage") + 
#   xlab("position relative to 3'SS") + 
#   ylab("normalized density of XL clusters") + 
#   theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
#   scale_x_continuous(limits = c(-300, 50)) +
#   scale_y_continuous(limits = c(0, ymax))
# 
# ggsave(output1)
# 
# gg.5SS <- ggplot() + theme_bw() + 
#   geom_area(aes(c(-50:300), sum.control.5SSA.norm.smooth), color="grey", alpha=0.3) + 
#   geom_line(aes(c(-50:300), sum.5SSA.norm.smooth)) + 
#   geom_vline(xintercept = 0, alpha = 0.5, linetype="longdash") +
#   ggtitle("Cluster coverage") + 
#   xlab("position relative to 5'SS") + 
#   ylab("normalized density of XL clusters") + 
#   theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
#   scale_x_continuous(limits = c(-50, 300)) +
#   scale_y_continuous(limits = c(0, ymax))
# 
# ggsave(output2) 
# 
# # separate 3SS from 5SS and convert it to matrix
# rnames <- data3SS.5SS[,2]
# data3SS <- data3SS.5SS[,4:354]
# data3SS[,300] <- -1  #exon start
# data5SS <- data3SS.5SS[,857:1206]
# data5SS[,50] <- -1  #exon start
# 
# data3SS <- as.data.frame(data3SS)
# data_matrix <- data.matrix(data3SS)
# 
# 
# pdf(paste(output3, sep=""), width = 8.5, height = 11)
# 
# heatmap.2(
#   data_matrix+1,
#   dendrogram = "none",
#   scale      = "none",
#   trace      = "none",
#   Rowv = FALSE,
#   Colv = FALSE,
#   key        = FALSE,
#   #labRow     = FALSE,
#   labRow     = data3SS.5SS$V1,
#   #labCol     = c(-300:50),
#   labCol = FALSE,
#   col    = c("#FFFFFF", "#FFBDC0", "#F27F89", "#A12931", "#270004"),
#   main="",
#   key.xlab="pentamer density",
#   key.ylab="",
#   key.title="",
#   cexRow=0.2,
#   cexCol=0.2
# )
# 
# #gg.3SS
# 
# # 5SS HeatMap
# data5SS <- as.data.frame(data5SS)
# data_matrix <- data.matrix(data5SS)
# 
# heatmap.2(
#   data_matrix+1,
#   dendrogram = "none",
#   scale      = "none",
#   trace      = "none",
#   Rowv = FALSE,
#   Colv = FALSE,
#   key        = FALSE,
#   labRow     = data3SS.5SS$V1,
#   #labCol     = c(-50:300),
#   labCol = FALSE,
#   col    = c("#FFFFFF", "#FFBDC0", "#F27F89", "#A12931", "#270004"),
#   main="",
#   key.xlab="pentamer denstiy",
#   key.ylab="",
#   key.title="",
#   cexRow=0.2,
#   cexCol=0.2 
# )
# 
# #gg.5SS
# 
# # distance between regulated and controls
# dist3SS.reg.vs.control <- norm.3SSA.coverage - norm.control.3SSA.coverage
# dist5SS.reg.vs.control <- norm.5SSA.coverage - norm.control.5SSA.coverage
# dist.exon.reg.vs.control <- norm.exon.coverage - norm.control.exon.coverage
# 
# # enrichment between regulated and controls
# enr3SS.reg.vs.control <- norm.3SSA.coverage / norm.control.3SSA.coverage
# enr5SS.reg.vs.control <- norm.5SSA.coverage / norm.control.5SSA.coverage
# enr.exon.reg.vs.control <- norm.exon.coverage / norm.control.exon.coverage
# 
# # Table 1
# reg.exons.cov <- c(norm.3SSA.coverage, norm.exon.coverage, norm.5SSA.coverage)
# control.exons.cov <- c(norm.control.3SSA.coverage, norm.control.exon.coverage, norm.control.5SSA.coverage)
# enrichment <- c(norm.3SSA.coverage/norm.control.3SSA.coverage, norm.exon.coverage/norm.control.exon.coverage, norm.5SSA.coverage/norm.control.5SSA.coverage)
# distance <- c(dist3SS.reg.vs.control, dist.exon.reg.vs.control, dist5SS.reg.vs.control)
# table1 = data.frame(reg.exons.cov, control.exons.cov, enrichment, distance) 
# rownames(table1) <- c("3ss","exon","5ss")
# colnames(table1) <- c("reg.exons.cov","control.exons.cov","enrichment","distance")
# 
# # Table 2
# ratio.regulated <- c(norm.3SSA.coverage/norm.5SSA.coverage, norm.exon.coverage/norm.3SSA.coverage, norm.exon.coverage/norm.5SSA.coverage, norm.5SSA.coverage/norm.3SSA.coverage)
# ratio.contols <- c(norm.control.3SSA.coverage/norm.control.5SSA.coverage, norm.control.exon.coverage/norm.control.3SSA.coverage, norm.control.exon.coverage/norm.control.5SSA.coverage, norm.control.5SSA.coverage/norm.control.3SSA.coverage)
# table2 = data.frame(ratio.regulated, ratio.contols) 
# rownames(table2) <- c("3ss vs. 5ss","exon vs. 3ss", "exon vs. 5ss","5ss vs. 3ss")
# colnames(table2) <- c("regulated","control")
# 
# library(gridExtra)
# gt1 <- tableGrob(round(table1,2))
# gt2 <- tableGrob(round(table2,2))
# grid.arrange(gt1, gt2)
# dev.off() 
# 





x <- c(0,1,1,0)
y <- scale_exons( x, 1000)
