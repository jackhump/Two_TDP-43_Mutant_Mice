#Make pie chart of the event types in the new datasets 

library(SGSeq) 
library(dplyr) 
library(ggplot2) 
library(data.table)

# f210i.dir <- "/SAN/vyplab/IoN_RNAseq/Kitty/F210I/sgseq/f210i/"
# m323k.dir <- "/SAN/vyplab/IoN_RNAseq/Kitty/M323K/sgseq/" 

f210i.dir <- "/Users/Jack/google_drive/TDP_paper/Sgseq/F210I_embryonic_brain/"
m323k.dir <- "/Users/Jack/google_drive/TDP_paper/Sgseq/M323K_adult_brain/"

f210i.res <- paste0(f210i.dir, "F210I_embryonic_brain_res_clean_novel.RData") 
m323k.res <- paste0(m323k.dir, "M323K_adult_brain_res_clean_novel.RData")

load(f210i.res) 

f210i.res.clean <- res.clean 

load(m323k.res) 

m323k.res.clean <- res.clean 

#f210i.res.clean <- read.table(f210i.res,header=TRUE)
#m323k.res.clean <- read.table(m323k.res,header=TRUE)

makePieChart <- function(sgseqRes, title, FDRlimit, outFolder){
# filter by FDR < 0.05
res.sig <- dplyr::filter(sgseqRes, FDR < FDRlimit) %>% select(one_of(c("groupID", "variantType", "FDR")) ) 
# for each groupID take one event
res.sig.by.group <- res.sig %>% group_by(groupID) %>% do(head(.,1))  
num.events <- nrow(res.sig.by.group)
print(paste0("total number of events: ", num.events ))
print(table(res.sig.by.group$variantType))
# treat each class of event as separate
res.events <- table(unlist(strsplit(res.sig.by.group$variantType, "+", fixed= TRUE) ) ) 
print(res.events)
print( sum(res.events) )
res.events.plot <- c() 
res.events.plot["Cassette exons"] <- res.events["S2E:I"] + res.events["S2E:S"] + 
    res.events["SE:I"] + res.events["SE:S"]
res.events.plot["Retained introns"] <- res.events["RI:E"] + res.events["RI:R"] 
res.events.plot["Alternative first exon"] <- res.events["AS"] + res.events["AFE"] 
res.events.plot["Alternative last exon"] <- res.events["AE"] + res.events["ALE"] 
res.events.plot["Alternative 3' site"] <- res.events["A3SS:P"] + res.events["A3SS:D"] 
res.events.plot["Alternative 5' site"] <- res.events["A5SS:P"] + res.events["A5SS:D"] 
res.events.plot["Mutually exclusive exons"] <- res.events["MXE"]
res.events.plot <- as.data.frame(res.events.plot) 
names(res.events.plot) <- "varcounts" 
res.events.plot$variantType <- row.names(res.events.plot)
res.events.plot <- res.events.plot[ order(res.events.plot$varcounts,decreasing=FALSE),]  
res.events.plot$variantType <- factor(res.events.plot$variantType, levels = rev(res.events.plot$variantType) )
res.events.plot <- dplyr::mutate(res.events.plot, pos = cumsum(varcounts) - 0.5*varcounts) 
res.events.plot$prop <- signif( (res.events.plot$varcounts / sum(res.events.plot$varcounts) ) * 100, 3) 
res.events.plot$prop <- paste0(res.events.plot$prop, "%")

pie <- ggplot(res.events.plot, aes(x="", y = varcounts, fill = variantType)) + 
       geom_bar(width = 1, stat="identity")  + 
       geom_text(aes(x=1.6, y = pos, label = prop), size = 3 ) +
       coord_polar(theta="y") + 
       scale_fill_brewer("",palette="Dark2", direction = 1) + 
       theme_void() +
       ggtitle(title) +
       annotate("text", x = 1.8, y = sum(res.events.plot$varcounts) / 2, 
        label = paste0("total number of events at adjusted p < ",FDRlimit,": ", num.events ) )

print(pie)
ggsave(paste0(outFolder,"/", title, "_sgseq_pie_chart.pdf") )

write.table( res.events.plot, paste0(outFolder, "/", title, "_sgseq_variant_type_table.tab"), col.names = TRUE, sep = "\t")
} 

makePieChart(f210i.res.clean, "RRM2mut", 0.05, f210i.dir)
makePieChart(m323k.res.clean, "LCDmut", 0.05, m323k.dir) 

# m323k.sig <- dplyr::filter(m323k.res.clean, FDR < 0.05) %>% select(one_of(c("groupID", "variantType", "FDR")) ) 
# m323k.sig.by.group <- m323k.sig %>% group_by(groupID) %>% do(sample_n(.,1))  
# m323k.events <- table(unlist(strsplit(m323k.sig.by.group$variantType, "+", fixed= TRUE) ) ) 
# m323k.events.plot <- c() 
# m323k.events.plot["Cassette exons"] <- m323k.events["S2E:I"] + m323k.events["S2E:S"] + m323k.events["SE:I"] + m323k.events["SE:S"] 
# m323k.events.plot["Retained introns"] <- m323k.events["RI:E"] + m323k.events["RI:R"] 
# m323k.events.plot["Alternative first exon"] <- m323k.events["AS"] + m323k.events["AFE"] 
# m323k.events.plot["Alternative last exon"] <- m323k.events["AE"] + m323k.events["ALE"] 
# m323k.events.plot["Alternative 3' site"] <- m323k.events["A3SS:P"] + m323k.events["A3SS:D"] 
# m323k.events.plot["Alternative 5' site"] <- m323k.events["A5SS:P"] + m323k.events["A5SS:D"] 
# m323k.events.plot <- as.data.frame(m323k.events.plot) 
# names(m323k.events.plot) <- "varcounts"  
# m323k.events.plot$variantType <- factor(row.names(m323k.events.plot), level = row.names(m323k.events.plot) ) 
# m323k.events.plot <- dplyr::mutate(m323k.events.plot, pos = cumsum(varcounts) - 0.5*varcounts) 

# bar <- ggplot(m323k.events.plot, aes(x="", y = varcounts, fill = variantType)) + geom_bar(width = 1, stat="identity") + 
#        geom_text(aes(x="", y = pos, label = varcounts) ) 
# pie <- bar + coord_polar(theta="y") + scale_fill_brewer(palette="Dark2") + theme_void()  
# ggsave("m323k_new_emb_pie.pdf") 
 
