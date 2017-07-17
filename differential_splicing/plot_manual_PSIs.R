# a set of known TDP-regulated exons
# splice junction counts come from the old embryonic data
library(dplyr)
library(ggplot2)
iFolder <- "/Users/Jack/google_drive/TDP_paper/"
manual_results <- paste0(iFolder, "/differential_splicing/manually_curated_psi_values.tab")

data <- read.table(manual_results,header=TRUE)

data$IID <- paste0(data$Number, data$GeneID)

data$my.B <- pmax(data$InclusionA, data$InclusionB)

data$cc <- 2- as.numeric(data$Sample)

data <- data[complete.cases(data),]

data.F210I <- filter(data, Dataset == "F210I")
data.M323K <- filter(data,Dataset=="M323K")

testPSI <- function(data, code){

  data.filtered <- filter(data, Dataset == code)

  mod <- glm(data = data.filtered, family = binomial(link = "logit"), formula = 'cbind(my.B, Exclusion) ~ cc*GeneID - cc - 1')
  print(summary(mod))
  data.filtered.sum <- dplyr::mutate(data.filtered, fitted.psi.log = boot::logit(fitted(mod))) %>%
    dplyr::select(GeneID, fitted.psi.log, cc) %>% 
    dplyr::distinct()

  my.std <- coef(summary(mod))[,"Std. Error"]
  my.pvalue <- coef(summary(mod))[,"Pr(>|z|)"]
  
  data.filtered.sum$std <- my.std[ paste0('cc:GeneID', data.filtered.sum$GeneID) ]
  data.filtered.sum$pvalue <- my.pvalue[ paste0('cc:GeneID', data.filtered.sum$GeneID) ]
  
  data.filtered.sum <- dplyr::mutate(data.filtered.sum,
                                  low.CI.log = fitted.psi.log - 1.96*std,
                                  upp.CI.log = fitted.psi.log + 1.96*std,
                                  low.CI = boot::inv.logit(low.CI.log),
                                  upp.CI = boot::inv.logit(upp.CI.log),
                                  fitted.psi = boot::inv.logit(fitted.psi.log))
  # split into case and control to calculate dPSI
  cases <- dplyr::filter(data.filtered.sum, cc == 1)
  controls <- dplyr::filter(data.filtered.sum, cc == 0)
  dPSI <- dplyr::mutate(cases, 
                        dPSI = cases$fitted.psi - controls$fitted.psi,
                        low.CI = cases$low.CI - controls$fitted.psi,
                        upp.CI = cases$upp.CI - controls$fitted.psi)
  
  
  dPSI$signif.code <- ifelse( dPSI$pvalue > 0.05, yes = "ns",
                              no = ifelse( dPSI$pvalue > 0.01, yes = "*",
                                           no = ifelse( dPSI$pvalue > 0.001, yes = "**",
                                                        no = "***" )))
  dPSI$dataset <- code
  return(dPSI)
}
F210I <- testPSI(data, "F210I")
M323K <- testPSI(data, "M323K")

dPSI <- rbind(F210I, M323K)

sigInOne <- group_by(dPSI, GeneID) %>%
            summarise( minP = min(pvalue)) %>%
            filter(minP < 0.05) %>%
            pull(GeneID)

dPSI_sigInOne <- filter(dPSI, GeneID %in% sigInOne)

dPSI_sigInOne$dataset <- gsub("F210I", "RRM2mut", dPSI_sigInOne$dataset)
dPSI_sigInOne$dataset <- gsub("M323K", "LCDmut", dPSI_sigInOne$dataset)

direction_in_KD <- read.table("google_drive/TDP_paper/differential_splicing/target_gene_direction_in_KD.txt", header=FALSE, sep = "\t")
dPSI_sigInOne$direction <- direction_in_KD$V2[match(dPSI_sigInOne$GeneID, direction_in_KD$V1)]  

p <- ggplot( dPSI_sigInOne, aes(x = rev(as.factor(dataset)), y = dPSI, label = signif.code, fill = dataset)) + 
    geom_col() + 
    geom_errorbar(aes(ymin = low.CI, ymax = upp.CI), width = 0.2) + 
    facet_wrap(~direction+GeneID) + 
    geom_text( y = 0.6) +
    ylim( -0.6, 0.65) +
    geom_hline(yintercept=0, linetype=2) + 
    #scale_x_discrete(levels= rev(levels(as.factor(dPSI_sigInOne$dataset)))) +
    theme_bw() + 
    xlab("") +
    ylab("delta PSI")

#data.F210I$fitted.psi <- fitted(mod)
#print(summary(mod))
#print(confint(mod))


pdf(paste0("/Users/Jack/google_drive/TDP_paper/differential_splicing/manual_dPSI_chosen_genes_embryo_heads.pdf"))
print(p)
dev.off()
