source("functions.R")

data <- read.table("Feb5 - Outside China.tsv", sep = "\t", header = TRUE)

data$Confirmed <- date.process(data$Confirmed)
data$Arrived <- date.process(data$Arrived)
data$Symptom <- date.process(data$Symptom)
data$Initial <- date.process(data$Initial)
data$Hospital <- date.process(data$Hospital)

data <- parse.infected(data)

## Only focus on the following countries/regions: Japan, Singapore, Taiwan, HongKong, Macau, Korea
data$Country_or_Region <- do.call(rbind, strsplit(as.character(data$Case), "-"))[, 1]
data <- subset(data, Country_or_Region %in% c("Japan", "Singapore", "Taiwan", "Korea", "HongKong", "Macau"))

## Only consider cases who come from Wuhan on or before January 23
data <- subset(data, ! (Outside %in% c("Y", "L")) & Arrived <= 23+31 ## & Confirmed <= 28+31
               )

table(data$Country_or_Region)


#### Data visualisation

#' Sample the trajectory
sample.trajectory <- function(data, m = 100, ...) {
    library(mice)

    data_imputed <- mice(data[, c("Arrived", "Symptom", "Initial", "Hospital", "Confirmed")], m, print = FALSE)
    data_output <- list()
    for (i in 1:m) {
        infected_imputed <- impute.infected(complete(data_imputed, i)$Symptom,
                                            data$Infected_first, data$Infected_last, ...)
        data_output[[i]] <- data.frame(Case = data$Case)
        data_output[[i]]$iter <- i
        data_output[[i]]$Infected <- infected_imputed
        data_output[[i]] <- cbind(data_output[[i]], complete(data_imputed, i))

    }

    do.call(rbind, data_output)
}

set.seed(20200206)
m <- 1000
data_imputed <- sample.trajectory(data, m)

data_imputed_long <- reshape2::melt(data_imputed, id.vars = c("Case", "iter"))
data_imputed_long$value <- data_imputed_long$value - 1 + as.Date("2019-12-1")

data_long <- reshape2::melt(data[, setdiff(colnames(data_imputed), c("iter", "Infected"))], id.vars = c("Case"))
data_long$value <- data_long$value - 1 + as.Date("2019-12-1")

simplify.case <- function(case) {
    tmp <- strsplit(as.character(case), "-")[[1]]
    paste0(switch(tmp[1],
                  Japan = "JP",
                  Korea = "KR",
                  Macau = "MO",
                  HongKong = "HK",
                  Taiwan = "TW",
                  Singapore = "SG"),
           "-", tmp[2])
}
data_imputed_long$Case <- sapply(data_imputed_long$Case, simplify.case)
data_long$Case <- sapply(data_long$Case, simplify.case)

library(ggplot2)
library(ggthemes)
library(ggstance)
pl <- ggplot() +
    stat_bin(binwidth = 1,
             aes(x = value, y = ..ncount.., fill = variable),
             data = subset(data_imputed_long, variable == "Infected"),
             position = position_dodge(width = 2),
             alpha = 0.5) +
    scale_fill_manual(values = c("grey25")) +
    facet_grid(Case ~ ., scale = "free")
pl <- pl + geom_point(aes(x = value, y = 0.5, color = variable),
                      data = data_long, alpha = 0.7, stroke = 1.5, size = 1.5, shape = 4,
                      position = position_dodgev(height = 0.6))
pl <- pl + theme_void() + theme(axis.text.x = element_text(),
                                legend.title=element_blank(),
                                legend.position = "bottom")

ggsave(filename = "Feb4_data.pdf", pl, height = 12, width = 7)

#### Imputed infection time matrix
as.count <- function(infected, last_date = 23+31) {
    table(factor(infected, levels = 1:last_date))
}
N <- 23 + 31

library(reshape2)
library(ggplot2)
library(ggthemes)
library(data.table)

OI_imputed <- matrix(0, N, m)
for (i in 1:m) {
    OI_imputed[, i] <- as.count(subset(data_imputed, iter == i)$Infected)
}

df <- melt(OI_imputed[, 1:20])
names(df) <- c("date", "impute", "count")
df$date <- df$date - 1 + as.Date("2019-12-01")
df_summary <- data.frame(date = 1:N - 1 + as.Date("2019-12-01"),
                         count_mean = apply(OI_imputed, 1, mean),
                         count_se = apply(OI_imputed, 1, sd) / sqrt(m))

df$lm_data <- df$date > as.Date("2020-1-1") & df$date < as.Date("2020-1-16")
df_summary$lm_data <- df_summary$date > as.Date("2020-1-1") & df_summary$date < as.Date("2020-1-16")

compute.offset <- function(t, or = c(1, 2, 3)) {

    if (or == 1) {
        or <- rep(1, N)
    } else if (or == 2) {
        or <- c(rep(1, N - 4), rep(2, 4))
    } else if (or == 3) {
        or <- c(rep(0, N - 1), N)
    }

    log(rev(cumsum(rev(or)))[t])

}

df_summary$t <- as.numeric(df_summary$date - as.Date("2019-12-01") + 1)
df_summary$offset1 <- compute.offset(df_summary$t, 1)
df_summary$offset2 <- compute.offset(df_summary$t, 2)
df_summary$offset3 <- compute.offset(df_summary$t, 3)


pl <- ggplot() + aes(x = date) +
  geom_point(data = df, aes(y = count), position = position_jitter(width = 0.1), alpha = 0.2) +
  geom_bar(stat = "identity", data = df_summary, aes(y = count_mean), alpha = 0.5, col = "white") +
    geom_errorbar(data = df_summary, aes(ymin = count_mean - 1 * count_se, ymax = count_mean + 1 * count_se, col = !lm_data)) + guides(col = FALSE) +
    theme_tufte() + xlab("") + ylab("count")

lancet.r <- log(2) / 6.4 ## https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)30260-9/fulltext#

pl1 <- ggplot(data = df_summary) +
    aes(x = date,
        y = log(count_mean) - offset1,
        ymin = log(count_mean - 1 * count_se) - offset1,
        ymax = log(count_mean + 1 * count_se) - offset1) +
    geom_point(aes(col = !lm_data)) + geom_errorbar(aes(col = !lm_data)) +
    theme_tufte() + xlab("") + ylab("log(mean count) - offset1")
fit1 <- lm(log(count_mean) - offset1 ~ date, subset(df_summary, lm_data))
pl1 <- pl1 + geom_abline(intercept = fit1$coef[1], slope = fit1$coef[2], col = "red")
pl1 <- pl1 + geom_text(x = as.Date("2020-01-05"), y = 0, label = paste("r =", round(fit1$coef[2], 2)), col = "red") + guides(col = FALSE)
fit1.lancet <- lm(log(count_mean) - offset1 - I(lancet.r * as.numeric(date)) ~ 1, subset(df_summary, lm_data))
pl1 <- pl1 + geom_abline(intercept = fit1.lancet$coef[1], slope = lancet.r, col = "blue")
pl1 <- pl1 + geom_text(x = as.Date("2020-01-18"), y = -3, label = paste("r =", round(lancet.r, 2)), col = "blue") + guides(col = FALSE)

pl2 <- ggplot(data = df_summary) +
    aes(x = date,
        y = log(count_mean) - offset2,
        ymin = log(count_mean - 1 * count_se) - offset2,
        ymax = log(count_mean + 1 * count_se) - offset2) +
    geom_point(aes(col = !lm_data)) + geom_errorbar(aes(col = !lm_data)) +
    theme_tufte() + xlab("") + ylab("log(mean count) - offset2")
fit2 <- lm(log(count_mean) - offset2 ~ date, subset(df_summary, lm_data))
pl2 <- pl2 + geom_abline(intercept = fit2$coef[1], slope = fit2$coef[2], col = "red")
pl2 <- pl2 + geom_text(x = as.Date("2020-01-05"), y = -0.5, label = paste("r =", round(fit2$coef[2], 2)), col = "red") + guides(col = FALSE)
fit2.lancet <- lm(log(count_mean) - offset2 - I(lancet.r * as.numeric(date)) ~ 1, subset(df_summary, lm_data))
pl2 <- pl2 + geom_abline(intercept = fit2.lancet$coef[1], slope = lancet.r, col = "blue")
pl2 <- pl2 + geom_text(x = as.Date("2020-01-18"), y = -4, label = paste("r =", round(lancet.r, 2)), col = "blue") + guides(col = FALSE)

pl3 <- ggplot(data = df_summary) +
    aes(x = date,
        y = log(count_mean) - offset3,
        ymin = log(count_mean - 1 * count_se) - offset3,
        ymax = log(count_mean + 1 * count_se) - offset3) +
    geom_point(aes(col = !lm_data)) + geom_errorbar(aes(col = !lm_data)) +
    theme_tufte() + xlab("") + ylab("log(mean count)")
fit3 <- lm(log(count_mean) - offset3 ~ date, subset(df_summary, lm_data))
pl3 <- pl3 + geom_abline(intercept = fit3$coef[1], slope = fit3$coef[2], col = "red")
pl3 <- pl3 + geom_text(x = as.Date("2020-01-05"), y = -2.5, label = paste("r =", round(fit3$coef[2], 2)), col = "red") + guides(col = FALSE)
fit3.lancet <- lm(log(count_mean) - offset3 - I(lancet.r * as.numeric(date)) ~ 1, subset(df_summary, lm_data))
pl3 <- pl3 + geom_abline(intercept = fit3.lancet$coef[1], slope = lancet.r, col = "blue")
pl3 <- pl3 + geom_text(x = as.Date("2020-01-18"), y = -5, label = paste("r =", round(lancet.r, 2)), col = "blue") + guides(col = FALSE)

g <- ggplotGrob(pl + ggtitle("A"))
g1 <- ggplotGrob(pl1 + ggtitle("B"))
g2 <- ggplotGrob(pl2 + ggtitle("C"))
g3 <- ggplotGrob(pl3 + ggtitle("D"))

library(gridExtra)
pdf("Feb4_fit.pdf", height = 6, width = 9)
grid.arrange(grobs = list(g, g1, g2, g3), widths = c(1, 1, 1), heights = c(1,1), layout_matrix = rbind(c(1,1,1),c(2,3,4)))
dev.off()
