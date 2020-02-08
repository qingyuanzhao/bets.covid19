source("functions.R")

data <- read.table("Feb5 - In China.tsv", sep = "\t", header = TRUE)

data$Confirmed <- date.process(data$Confirmed) 
data$Arrived <- date.process(data$Arrived) 
data$Symptom <- date.process(data$Symptom)
data$Initial <- date.process(data$Initial) 
data$Hospital <- date.process(data$Hospital)

data <- data[-72, ] # Don't know how to parse this infected date yet
data <- parse.infected(data)

## Only consider cases with known symptom onset, arrived on or before
##January 23 
data <- subset(data, !is.na(Symptom))
data <- subset(data, Arrived <= 23+31) 
data <- subset(data, !(is.na(Arrived) & Infected_first == 1 & Infected_last == Symptom)) # remove cases with no information

#' Compute the likelihood
#'
#' GT is a discretized distribution for the generation time
#'
infection.likelihood <- function(symptom, infected_first, infected_last, GT) {
    loglike <- 0 
    for (i in 1:nrow(data)) {
        min.incub <- symptom[i] - infected_last[i] 
        max.incub <- symptom[i] - infected_first[i]
        loglike <- loglike + log(sum(GT$GT[1 +(min.incub):(max.incub)])) 
    } 
    loglike 
}

myfun <- function(par) { 
    GT <- R0::generation.time("gamma", par,truncate = 100); 
    infection.likelihood(data$Symptom, data$Infected_first, data$Infected_last, GT) 
}

fit <- optim(c(7.5, 3.4), myfun, control = list(fnscale = -1))

pars <- expand.grid(mean = seq(5, 15, 0.1), sd = seq(3, 10, 0.1))

pars$loglike <- apply(pars, 1, myfun)

pars$in.CR <- (pars$loglike > fit$value - qchisq(0.95, 1) / 2)

library(ggplot2) 
ggplot(pars) + aes(x = mean, y = sd, fill = in.CR) +geom_tile()

## MLE for incubation period: mean = 6.8, sd = 4.2 CI for mean
## incubation period: about 5.5 to 9.0
