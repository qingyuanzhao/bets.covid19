
library(ggplot2)

df = read.csv("posterior_w.csv",header = F,sep="")
names(df) = c("pos_mean","q025","q975")
df$day = (0:29)+0.5

df$Model <- "Nonparametric"

x = (0:(30*1000))/1000
y = dgamma(x,shape=1.8632188,rate=0.3347591)
df2 = data.frame(x,y)
df2$Model <- "Parametric"

## cbPalette <- c(## "#000000",
##                 "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p = ggplot(df,aes(x=day,y=pos_mean, color = Model)) + geom_point() + geom_errorbar(aes(ymin=q025,ymax=q975)) + geom_line(data=df2,aes(x=x,y=y)) + scale_color_manual(values = c("#56B4E9", "#E69F00"))
p = p + xlab("Incubation period (days)") + ylab("Density") + theme_bw(base_size = 15) + theme(legend.position = "bottom")

p

ggsave(filename = "../../Figures/posterior_of_w.pdf", p)
