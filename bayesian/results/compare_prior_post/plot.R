library(ggplot2)
library(reshape2)

df = read.csv("prior_post_p14.csv")
names(df) = c("p14_a1_post", "p14_a10_post", "p14_a1_prior", "p14_a10_prior")

molten.df = melt(df)
molten.df$alpha = "mu = 1"
molten.df$alpha[(molten.df$variable=="p14_a10_post") | (molten.df$variable=="p14_a10_prior")]="mu = 10"
levels(molten.df$variable) = c("Posterior",   "Posterior",  "Prior",  "Prior")
labs = c("Posterior","Prior")

cbPalette <- c(## "#000000",
                "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p = ggplot(molten.df) + geom_density(aes(value,colour=variable,fill=variable),alpha=0.5) + theme(axis.title.x=element_blank(),axis.title.y = element_blank(),legend.title = element_blank()) + facet_wrap(~alpha) + theme_bw(base_size = 18) + theme(legend.position = "bottom", legend.title = element_blank()) + xlab("P(Incubation >= 14 days)") + ylab("Density") + scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette)
p

ggsave(filename = "../../Figures/compare_prior_post.pdf", p, width = 13, height = 7)

df$day = (0:29)+0.5
