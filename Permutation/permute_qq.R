library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(readr)

# find files and order by permutation number
result_files <- list.files(pattern="_res_clean_novel.tab", full.names = TRUE, recursive=TRUE)
result_ID <- as.numeric(str_split_fixed(basename(result_files), pattern = "_", n = 4)[,3])

result_files <- result_files[order(result_ID)]

support_files <- list.files(pattern="_support.tab.permute", full.names=TRUE,recursive=TRUE)
support_ID <- as.numeric(str_split_fixed(basename(support_files), pattern = "\\.", n = 5)[,4])

support_files <- support_files[order(support_ID)]


# for each result,
# read in files
# extract p value vector

results <- purrr::map( result_files, ~{
  print(.x)
  d <- fread(.x, data.table = FALSE)
  return(d)
})


# do same for support file
supports <- purrr::map( support_files, ~{
	s <- fread(.x, data.table = FALSE)
	return(s)
})



create_quantiles <- function(results){
  pvector <- results$pvalue
  observed <- -log10(sort(pvector, decreasing = FALSE))
  return(observed)
}

get_sample_list <- function(support){
  samples <- support[support[,3] == "CONTROL",]$sample_name
  samples <- samples[order(samples)]
}

sample_orders <- purrr::map( supports, get_sample_list)

# look for order of wildtype samples
support_df <- tibble(
  replicate = support_ID,
  order = map_chr(sample_orders, ~paste( .x, collapse = "+") )
  )

magic_order <- c(
  "M323K_WT_1+M323K_WT_2+M323K_WT_3+M323K_WT_4",
  "M323K_HOM_1+M323K_HOM_2+M323K_HOM_3+M323K_HOM_4"
  )

support_df$correct.ordering <- ifelse( support_df$order %in% magic_order,
  TRUE, FALSE)


quantiles <- map(results, create_quantiles)
names(quantiles) <- result_ID

# bind all quantiles together and add permutation order info for plotting

df <- 
  bind_rows(quantiles) %>%
  mutate( expected = -log10(ppoints(nrow(.)))) %>%
  select( expected, everything() ) %>%
  tidyr::gather( "replicate", "observed", -expected) %>%
  left_join( support_df, by = "replicate")

df_distinct <- 
  df %>%
  select( -replicate) %>%
  distinct()

save( df_distinct, file = "all_permutations.Rdata")

# plot big QQ - doesn't work on cluster! bullshit
# library(ggplot2)
# png("test.png")

# ggplot( df, aes( x = expected, y = observed)) +#, group = replicate )) + 
#   geom_point( aes(colour = correct.ordering)) +
#   geom_abline( slope = 1, intercept = 0, linetype = 3) +
#   xlab("Expected quantiles") +
#   ylab("Observed quantiles") +
#   theme_bw()

# dev.off()

##----------------- skiptic counting ---------------

# iterate through
# support and results - must match!

# for each results, count the number of skiptic exons


# d - an SGSeq results table
# groupIDs - those events of interest
# sample_list - ?
createVarTable <- function(d, groupIDs,sample_list){
  # conversion table for SGSeq categories to human readable names
  variantType_df <- data.frame(
    SGSeq = c("SE", "S2E", "RI", "A3SS", "A5SS", "ALE", "AFE", "AE", "AS", "MXE"),
    human = c("cassette_exon", "multi-cassette_exon", "retained_intron", "alt_3_splice_site", "alt_5_splice_site", "alt_last_exon", "alt_first_exon", "multi-alt_last_exon", "multi_alt_first_exon", "mut_exclusive_exons"),
    stringsAsFactors = FALSE
  )
  varTable <- lapply(groupIDs, FUN = function(i){
    event <- d[ d$groupID == i,]
    varType <- event$variantType
    varTypeSplit <- str_split_fixed(varType, ":",2)[,1]
    varTypeHuman <- NA
    # simple variants first:
    if( !any(grepl("\\+", varType)) ){
      ref <- "unassigned"
      if( nrow(event) == 1){
        ref <- 1
      }else{
        # for each event, rule which of the two is the reference
        # if any other type (ALE, AFE, MXE etc) or a mixture then just pick the most represented (meanBase)
        if( all("SE" %in% varTypeSplit) ){
          ref <- which( varType == "SE:I")
        }
        if( all("S2E" %in% varTypeSplit) ){
          ref <- which( varType == "S2E:I")
        }
        if( all("RI" %in% varTypeSplit) ){
          ref <- which( varType == "RI:R")
        }
        if( all("A5SS" %in% varTypeSplit) ){
          ref <- which( varType == "A5SS:P")
        }
        if( all("A3SS" %in% varTypeSplit) ){
          ref <- which( varType == "A3SS:P")
        }
      if( length(ref) > 1 | any(ref == "unassigned") ){
          ref <- which( event$exonBaseMean == max(event$exonBaseMean))
        }
      }
    }else{
      # for complex variants 
      #set whichever row has the highest meanBase as ref
      ref <- head(which( event$exonBaseMean == max(event$exonBaseMean)), 1)
      #return( list( i, varType, ref, "complex"))
      varTypeHuman <- "complex"
    }
    # create dPSI from the reference row
    control_psi <- event[ref, paste0(sample_list[[1]], "_psi") ]
    case_psi <- event[ref, paste0(sample_list[[2]], "_psi") ]
    dPSI <- mean(as.numeric(case_psi), na.rm = TRUE) - mean(as.numeric(control_psi), na.rm = TRUE)

    # translate varType to understandable
    if(is.na(varTypeHuman) ){
      if( !any(grepl("\\+", varType)) ){
        varTypeHuman <- variantType_df$human[ match( unique(varTypeSplit), variantType_df$SGSeq)]
        if( all(is.na(varTypeHuman) ) | length(varTypeHuman) > 1 ){
          varTypeHuman <- "complex"
        }
      }else{
        varTypeHuman <- "complex"
      }
    }
    # pick the largest coordinates to encompass all variants
    coords <- do.call(rbind,strsplit(event$coords, split = ":|-"))
    coordSizes <- apply(coords, 1, FUN = function(x) as.numeric(x[3]) - as.numeric(x[2]) )
    biggestCoords <- head(event$coords[order(coordSizes)],1)

    # annotation - are both the reference and alternate forms annotated?
    ref_anno <- event$txName[ref]
    alt_anno <- event$txName[-ref]
    if( all(ref_anno == "") ){ ref_anno <- "novel"}else{ ref_anno <- "annotated"}
    if( all(alt_anno == "") ){ alt_anno <- "novel"}else{ alt_anno <- "annotated"}
    # stranding! This should have been created in step2b oh well
    strand <- unique( str_split_fixed(c(event$to,event$from), ":", 4)[,4] )
    if( length(strand) > 1){
      strand <- "*"
    }
    varTable_row <- NULL
    tryCatch(
# assemble into a table
    varTable_row <- data.frame(
      groupID = i,
      gene = event$geneName[ref], # in case of weirdness
      EnsemblID = event$ensemblName[ref],
      coords = biggestCoords,
      strand = strand,
      variant_type = varTypeHuman,
      control_PSI = signif( mean(as.numeric(control_psi), na.rm = TRUE), 3),
      case_PSI = signif( mean(as.numeric(case_psi), na.rm = TRUE), 3),
      dPSI = signif(dPSI,3),
      FDR = signif(min(event$padj, na.rm= TRUE),3),
      ref_anno,
      alt_anno,
      stringsAsFactors=FALSE
    ), error = function(e){
      print(paste0("event: ", i, " in ", event$geneName[ref], " with ref = ", ref ))
     }
    )
    return(varTable_row)
  })
  varTable <- do.call(rbind,varTable)
  return(varTable)
}


#support <- supports[[1]]
#res <- results[[1]]

countEvents <- function(support, res, sample_order){
  sample_list <- split(support, support$condition_hom) %>% map("sample_name")


  sigGroupIDs <- unique( dplyr::filter(res, padj < 0.05 & geneName != "")$groupID)
  print("Number of unique events:")
  print(length(sigGroupIDs))

  # create variant table for significant events
  sigVarTable <- createVarTable(res, sigGroupIDs, sample_list)

  # apply cryptic and skiptic criteria

  sig <- nrow(sigVarTable)

  cryptic <- filter(sigVarTable, variant_type == "cassette_exon", control_PSI <= 0.05,   abs(dPSI) >= 0.05 ) %>% nrow()
  skiptic <- filter(sigVarTable, variant_type == "cassette_exon", control_PSI >= 0.95,   abs(dPSI) >= 0.05 ) %>% nrow()

  output <- tibble(order = sample_order, significant = sig, cryptic = cryptic, skiptic = skiptic)
  return(list(count_table = output, event_table = sigVarTable ) )
}

#magic <- "M323K_WT_1+M323K_WT_2+M323K_WT_3+M323K_WT_4"

#magic.order <- which(ordering == magic)[1]

#test <- countEvents( supports[[magic.order]], results[[magic.order]], ordering[[magic.order]])

ordering <- map_chr(sample_orders, ~paste( .x, collapse = "+") )

all_results  <- 
  pmap( list(supports, results, ordering), countEvents )

count_results <- 
  map_df( all_results, "count_table") %>%
  arrange(desc(significant)) %>%
  distinct()

write_csv(count_results, path = "permutation_event_counts.csv")

# save a bunch of var tables
var_tables <- map(all_results, "event_table")

save(list = var_tables, supports, ordering, file = "all_permutation_variant_tables.Rdata" )





# create_fake_p_vectors <- function(
#   difference,
#   n.repeats = 10,
#   n.genes = 20000,
#   n.samples = 4
#   ){
#   pvectors <- map(1:n.repeats, 
#     ~{
#       map_dbl( 1:n.genes, 
#       ~{
#         t.test( rnorm(n.samples, mean = 0, sd =1), rnorm(n.samples, 1 + difference) )$p.value
#       })
#   })
#   return(pvectors)
# }

# #fake_quantiles <- map(1:n.replicates, ~create_fake_p_vectors)

# fake_quantiles <- create_fake_p_vectors(difference = 0.2)



# names(fake_quantiles) <- as.character(1:n.replicates)


# qq_df <- 
#   map_df(fake_quantiles, create_quantiles ) %>%
#   mutate(expected = -log10(ppoints(nrow(.)))) %>%
#   select( expected, everything() ) %>%
#   tidyr::gather( "replicate", "observed", -expected) %>%
#   filter( expected >= 3)






map( fake_quantiles, )

