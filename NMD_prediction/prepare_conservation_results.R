library(dplyr)
library(data.table)
conservation_table <- "/Users/Jack/google_drive/TDP_paper/NMD_prediction/conservation/extreme_splicing_conservation_table.tab.gz"

d <- fread(paste0("zless ",conservation_table))
res <- group_by(d, exon.type) %>% 
        summarise( 
          n = n() , 
          mean = mean(phyloP.score, na.rm= TRUE), 
          sd = sd(phyloP.score, na.rm = TRUE), 
          sem = sd / sqrt(n) ) %>% 
        mutate(pretty = paste0(signif(mean,3),"±",signif(sem,3)) )


toPlot <- d
bed.names <- res$exon.type
bed.labels <- paste0(res$exon.type, "\n(",res$n, ")")

p <- ggplot(toPlot, aes(x = gsub(".","\n",exon.type,fixed =T), y = phyloP.score, fill = exon.type, colour = exon.type)) + 
  geom_violin(colour = "black", fill = "white") +
  geom_boxplot(colour = "black", fill = "skyblue", notch = T, width = 0.1 ) +
  xlab("") + 
  ylab("per exon mean phyloP score") + 
  ylim(c(-1,6)) +
  #scale_fill_manual(values = c("skyblue","firebrick3","skyblue","skyblue")) +
  scale_x_discrete(name=" ",
          limits=bed.names,
                    labels=bed.labels) +
  # ggtitle(paste0(graph_title,"\naverage phyloP conservation score per exon")) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
          axis.line = element_line(colour = "black"),
          legend.position="none"
          )

  other <- as.data.frame(fread("/Users/Jack/project/F210I_M323K/cryptic_skiptic/conservation/extreme_splicing_exon_conservation.tab"))
