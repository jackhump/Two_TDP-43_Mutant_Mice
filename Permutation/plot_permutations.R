library(tidyverse)

data <- "all_permutations.Rdata"

if( !exists("df_distinct")){
load(data)
}
d <- df_distinct
# for testing

#d <- filter(d, expected > 3)

# standard scatter plot - colour by order

p1 <- ggplot( d, aes(x = expected, y = observed, group = order, colour = correct.ordering) ) +
	geom_point() +
	geom_abline( slope = 1, intercept = 0, linetype = 3 ) +
	scale_colour_manual("true sample order", values = c("grey", "purple") ) +
	guides(colour=FALSE) +
	theme_classic() +
	ylab("observed quantiles") +
	xlab("expected quantiles")

# scatter plus line

p2 <- ggplot( d, aes(x = expected, y = observed, group = order, colour = correct.ordering) ) +
	geom_line() +
	geom_point() +
	geom_abline( slope = 1, intercept = 0, linetype = 3 ) +
	scale_colour_manual("true sample order", values = c("grey", "purple") ) +
	theme_classic() +
        ylab("observed quantiles") +
        xlab("expected quantiles")

# box plots for permutations and line/points for true order

p3 <- ggplot( d, aes(x = expected, y = observed, group = expected, colour = correct.ordering) ) +
	geom_line(data = d[d$correct.ordering,], aes(group = correct.ordering), colour = "purple" ) +
	geom_point(data = d[d$correct.ordering,], colour = "purple"  )  +
	geom_boxplot(data = d[!d$correct.ordering,], width = 0.1, colour = "black", fill = "grey") +
	geom_abline( slope = 1, intercept = 0, linetype = 3 ) +
	theme_classic() +
        ylab("observed quantiles") +
        xlab("expected quantiles")

# geom line + geom_ribbon for mean and sd

d_mean <- 
	filter(d, correct.ordering == FALSE ) %>%
	group_by(expected) %>%
	summarise( mean_o = mean(observed), sd_o = sd(observed), sem_o = sd(observed) / sqrt(n() ) )

p4 <- ggplot( ) +
        #geom_line(data = d[d$correct.ordering,], aes(x = expected, y = observed, group = correct.ordering), colour = "purple" ) +
        geom_point(data = d[d$correct.ordering,], aes(x = expected, y = observed), colour = "purple"  )  +
        geom_ribbon( data = d_mean, aes( x = expected, ymin = mean_o - sd_o, ymax = mean_o + sd_o), colour = NA, fill = "lightgrey" ) +
	geom_line(data = d_mean, aes(x = expected, y = mean_o), colour = "black" ) +
	
        geom_abline( slope = 1, intercept = 0, linetype = 3 ) +
        theme_classic() +
        ylab("observed quantiles") +
        xlab("expected quantiles")


# save plots

#png("permutation_scatter.png")
#print(p1)
#dev.off()

#png("permutation_scatter_line.png")
#print(p2)
#dev.off()

#png("permutation_boxplot.png")
#print(p3)
#dev.off()

svg("permutation_ribbon.svg")
print(p4)
dev.off()


