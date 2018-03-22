#Make pie chart of the event types in the new datasets 

#library(SGSeq) 
library(dplyr) 
library(ggplot2) 
library(data.table)

#--- data

f210i.res <- "../data/SGSeq/F210I_embryonic_brain_res_clean_novel.tab"
m323k.res <- "../data/SGSeq/M323K_New_CONTROL_HOM_res_clean_novel.tab"

f210i.res.clean <- read.table(f210i.res,header=TRUE,stringsAsFactors=FALSE)
m323k.res.clean <- read.table(m323k.res,header=TRUE, stringsAsFactors=FALSE)


#------ functions
 
makePieChart <- function(sgseqRes, title, FDRlimit, outFolder){
	# filter by FDR < 0.01!
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

	pie <- 
		ggplot(res.events.plot, aes(x="", y = varcounts, fill = variantType)) + 
       		geom_bar(width = 1, stat="identity")  + 
       		geom_text(aes(x=1.6, y = pos, label = prop), size = 3 ) +
       		coord_polar(theta="y") + 
       		scale_fill_brewer("",palette="Dark2", direction = 1) + 
       		theme_void() +
       		ggtitle(title) +
       		annotate("text", x = 1.8, y = sum(res.events.plot$varcounts) / 2, 
        	label = paste0("total number of events at adjusted p < ",FDRlimit,": ", num.events ) )

	#print(pie)
	ggsave(paste0(outFolder, title, "_sgseq_pie_chart.pdf") )

	write.table( res.events.plot, paste0(outFolder, title, "_sgseq_variant_type_table.tab"), col.names = TRUE, sep = "\t")
} 

#---- run script

makePieChart(f210i.res.clean, "RRM2mut", 0.01, "")
makePieChart(m323k.res.clean, "LCDmut", 0.01, "") 

