
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

ggsave(filename = "../../../Figures/posterior_of_w.pdf", p)


## Gender
df = read.csv("posterior_w_gender.csv",header = T,sep=",")
df$day = (0:29)+0.5

library(reshape2)
df_long <- melt(df, id.vars = "day")
df_long$Gender <- rep(c("Men", "Women"), each = nrow(df_long) / 2)
df_long$variable <- rep(rep(c("mean", "lower", "upper"), each = nrow(df_long) / 6), 2)
df <- dcast(df_long, Gender + day ~ variable)

library(ggplot2)
p <- ggplot(df, aes(x = day,y = mean, ymin = lower, ymax = upper, fill = Gender, linetype = Gender)) + geom_line(aes(color = Gender), size = 1.5) + geom_ribbon(alpha = 0.5) + scale_color_manual(values = c("#56B4E9", "#E69F00")) + scale_fill_manual(values = c("#56B4E9", "#E69F00")) + theme_bw(base_size = 18) + xlab("Incubation period") + ylab("Density (days)") + theme(legend.position = "bottom")

ggsave(filename = "../../../Figures/posterior_of_w_gender.pdf", p)

## Age
df = read.csv("posterior_w_age.csv",header = T,sep=",")
df$day = (0:29)+0.5

library(reshape2)
df_long <- melt(df, id.vars = "day")
df_long$Age <- rep(c(">=50", "<50"), each = nrow(df_long) / 2)
df_long$variable <- rep(rep(c("mean", "lower", "upper"), each = nrow(df_long) / 6), 2)
df <- dcast(df_long, Age + day ~ variable)

library(ggplot2)
p <- ggplot(df, aes(x = day,y = mean, ymin = lower, ymax = upper, fill = Age, linetype = Age)) + geom_line(aes(color = Age), size = 1.5) + geom_ribbon(alpha = 0.5) + scale_color_manual(values = c("#56B4E9", "#E69F00")) + scale_fill_manual(values = c("#56B4E9", "#E69F00")) + theme_bw(base_size = 18) + xlab("Incubation period") + ylab("Density (days)") + theme(legend.position = "bottom")

ggsave(filename = "../../../Figures/posterior_of_w_age.pdf", p)
