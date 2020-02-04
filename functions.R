#' Processing age to print its distribution
#'
#' @param age a vector of age, each entry is either a number (like 34) or age group (like 30s)
#'
#' @return each age is either repeated 10 times or expanded to 10 numbers (for example, 30s is expanded to 30, 31, ..., 39).
age.process <- function(age) {
    if (length(age) > 1) { # a vector
        return(as.vector(sapply(age, age.process)))
    }
    age <- as.character(age)
    if (substr(age, nchar(age), nchar(age)) == "s") {
        group <- as.numeric(substr(age, 1, nchar(age) - 1))
        return(group + (0:9))
    } else {
        return(rep(as.numeric(age), 10))
    }
}

#' Transform date to numeric
#'
#' For example, "23-Jan" to 23+31 = 52 (start date is set to December 1st, 2019).
#'
date.process <- function(date) {

    tmp <- as.numeric(as.Date(date, "%d-%b") - as.Date("2019-12-01")) + 1

    for (i in which(!is.na(tmp))) {
        if (tmp[i] > 365) {
            tmp[i] <- tmp[i] - 365
        }
    }

    tmp

}

#' Parse the infected date
parse.infected <- function(data) {
    parse.one.infected <- function(infected) {
        tmp <- strsplit(as.character(infected), "to")[[1]]
        if (length(tmp) == 0) {
            return(c(-Inf, Inf))
        } else if (length(tmp) == 1) {
            return(rep(date.process(tmp), 2))
        } else {
            return(date.process(tmp))
        }
    }

    infected_interval <- t(sapply(data$Infected, parse.one.infected))

    data$Infected_first <- pmax(infected_interval[, 1], 1)
    data$Infected_last <- pmin(infected_interval[, 2],
                               data$Arrived, data$Symptom,
                               data$Initial, data$Confirmed, na.rm = TRUE)

    data
}

#' Simple imputation
#'
simple.impute.onset <- function(data) {

    na.ind <- is.na(data$Symptom)
    print(paste(sum(na.ind), "of", nrow(data), "symptom onset date(s) are missing."))
    avg.symptom.to.initial <- round(mean(data$Initial - data$Symptom, na.rm = TRUE))
    print(paste("Imputing by initial medical visit date minus", avg.symptom.to.initial, "days..."))
    data$Symptom[na.ind] <- data$Initial[na.ind] - avg.symptom.to.initial

    na.ind <- is.na(data$Symptom)
    if (sum(na.ind) > 0) {
        print(paste(sum(na.ind), "of", nrow(data), "symptom onset date(s) are still missing."))
        avg.symptom.to.confirm <- round(mean(data$Confirmed - data$Symptom, na.rm = TRUE))
        print(paste("Imputing by confirmation date minus", avg.symptom.to.confirm, "days..."))
        data$Symptom[na.ind] <- data$Confirmed[na.ind] - avg.symptom.to.confirm
    }

    na.ind <- is.na(data$Symptom)
    print(paste(sum(na.ind), "sympton onset date(s) are missing after imputation."))

    data$Symptom
}

#' Multiple imputation
#'
multiple.impute.onset <- function(data, m = 50) {
    library(mice)
    data.imputed <- mice(data[, c("Arrived", "Symptom", "Initial", "Hospital", "Confirmed")], m, print = FALSE)

    sapply(1:m, function(i) complete(data.imputed, i)$Symptom)
}

#' Impute infected using symptom onset and the distribution of the incubation period, respecting information about the infected time
#'
#' @param symptom symptom onset date
#' @param infected_first beginning of possible infection time
#' @param infected_end end of the possible infection time
#' @param incubation_alpha alpha parameter in the gamma distribution of incubation
#' @param incubation_beta beta parameter in the gamma distribution of incubation
#'
impute.infected <- function(symptom, infected_first, infected_last,
                            incubation_alpha = 1.92, incubation_beta = 0.37 ## not exactly the same as NEJM paper
                            ) {

    stopifnot(length(symptom) == length(infected_first))
    stopifnot(length(symptom) == length(infected_last))
    if (sum(is.na(symptom)) > 0) {
        stop("Symptom onset date cannot be missing. Use simple.impute.onset or multiple.imput.onset to impute the missing values!")
    }

    infected <- rep(NA, length(symptom))
    for (i in 1:length(symptom)) {
        infected[i] <- round(symptom[i] - rgamma(1, incubation_alpha, incubation_beta))
        while (infected[i] < infected_first[i] ||
               infected[i] > infected_last[i]) {
            infected[i] <- ceiling(symptom[i] - rgamma(1, incubation_alpha, incubation_beta))
        }
    }

    infected
}
