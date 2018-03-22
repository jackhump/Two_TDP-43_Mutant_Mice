# make MA plots from SGSeq result tables
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)

# M323K adult brain

m323k_file <- "../data/SGSeq/M323K_New_CONTROL_HOM_res_clean_novel.tab"

m323k <- read_tsv( m323k_file )


mean_reads <- 
  m323k %>% 
  select( M323K_WT_1, M323K_WT_2, M323K_WT_3, M323K_WT_4,  M323K_HOM_1, M323K_HOM_2, M323K_HOM_3, M323K_HOM_4 ) %>%
  rowMeans()

m323k_MA <- 
  m323k %>%
  mutate( mean_reads = mean_reads) %>%
  select( groupID, padj, log2FC = log2fold_HOM_CONTROL, mean_reads ) #%>%

p <- 
  ggplot( m323k_MA, aes( x = log10(mean_reads), y = log2FC )) +
  geom_point( data = m323k_MA[ m323k_MA$padj >= 0.05,], colour = "black", size = 0.5 ) +
  geom_point( data = m323k_MA[ m323k_MA$padj < 0.05,], colour = "red", size = 0.5 ) +
  theme_classic() +
  xlab("log10(mean reads per junction)") +
  ylab("log2(fold change)") +
  ggtitle( "LCDmut Adult SGSeq - all splicing") +
  ylim( -5,5)

ggsave(filename = "m323k_adult_sgseq_MA_plot.png", device = "png", width = 7, height = 7 )

# F210I Embryonic brain

f210i_file <- "../data/SGSeq/F210I_embryonic_brain_res_clean_novel.tab"


f210i <- read_tsv(f210i_file)

f210i_mean_reads <- f210i %>% 
  select( F210I_CTL_1_norm, F210I_CTL_2_norm, F210I_CTL_3_norm, F210I_CTL_4_norm,  F210I_HOM_1_norm, F210I_HOM_2_norm, F210I_HOM_3_norm, F210I_HOM_4_norm ) %>%
  rowMeans()

f210i_MA <- f210i %>%
  mutate( mean_reads = f210i_mean_reads) %>%
  select( groupID, padj, log2FC = log2fold_HOM_CTL, mean_reads ) #%>%

p <- ggplot( f210i_MA, aes( x = log10(mean_reads), y = log2FC )) +
  geom_point( data = f210i_MA[ f210i_MA$padj >= 0.05,], colour = "black", size = 0.5 ) +
  geom_point( data = f210i_MA[ f210i_MA$padj < 0.05,], colour = "red", size = 0.5 ) +
  theme_classic() +
  xlab("log10(mean reads per junction)") +
  ylab("log2(fold change)") +
  ggtitle( "RRM2mut Embryo SGSeq - all splicing") +
  ylim( -5,5)

ggsave(filename = "f210i_adult_sgseq_MA_plot.png", device = "png", width = 7, height = 7 )

