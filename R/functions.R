#' Processing age to print its distribution
#'
#' @param age a vector of age, each entry is either a number (like 34) or age by decade (like 30s)
#'
#' @return each age is either repeated 10 times or expanded to 10 numbers (for example, 30s is expanded to 30, 31, ..., 39).
#'
#' @export
#'
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
#' @param date a vector of dates of the form "DD-MMM" (for example, 23-Jan).
#'
#' @return a vector of days since December 1st, 2019 (or example, 23-Jan is converted to 23+31 = 52).
#'
#' @export
#'
date.process <- function(date) {

    options(tz="UTC")
    Sys.setenv(TZ = "UTC")

    tmp <- as.numeric(as.Date(date, "%d-%b") - as.Date("2019-12-01")) + 1

    for (i in which(!is.na(tmp))) {
        if (tmp[i] > 365) {
            tmp[i] <- tmp[i] - 365
        }
    }

    tmp

}


#' Simple imputation of symptom onset date
#'
#' @param data a data frame with the following columns: Symptom, Initial, Confirmed.
#'
#' @return a vector of symptom onset dates with all the missing values imputed using initial medical visit date or confirmation date.
#'
#' @export
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

#' Multiple imputation of symptom onset date
#'
#' @param data a data frame with the following columns: Arrived, Symptom, Initial, Hospital, Confirmed.
#'
#' @return a matrix of symptom onset dates with missing values imputed using the mice package.
#'
#' @import mice
#' @export
#'
multiple.impute.onset <- function(data, m = 50) {
    data.imputed <- mice(data[, c("Arrived", "Symptom", "Initial", "Hospital", "Confirmed")], m, print = FALSE)

    sapply(1:m, function(i) complete(data.imputed, i)$Symptom)
}

#' Simulated infection time using symptom onset
#'
#' Uses distribution of the incubation period and respects information about the infected time
#'
#' @param symptom a vector of symptom onset dates
#' @param infected_first a vector of the first possible infection dates
#' @param infected_end a vector of the last possible infection dates
#' @param incubation_alpha alpha parameter in the gamma distribution of incubation period
#' @param incubation_beta beta parameter in the gamma distribution of incubation period
#'
#' @return a vector of simulated infection dates.
#'
#' @importFrom R0 generation.time
#' @export
#'
impute.infected <- function(symptom, infected_first, infected_last,
                            incubation_mean = 6.5, incubation_sd = 2.6) {
    ## default incubation period is gamma-distributed

    stopifnot(length(symptom) == length(infected_first))
    stopifnot(length(symptom) == length(infected_last))
    if (sum(is.na(symptom)) > 0) {
        stop("Symptom onset date cannot be missing. Use simple.impute.onset or multiple.imput.onset to impute the missing values!")
    }

    mGT <- generation.time("gamma", c(incubation_mean, incubation_sd))

    infected <- rep(NA, length(symptom))
    for (i in 1:length(symptom)) {
        infected[i] <- symptom[i] - sample(1:length(mGT$time) - 1, 1, prob = mGT$GT)
        boundary.keep <- sample(c(TRUE, FALSE), 1)
        while (infected[i] < infected_first[i] ||
               infected[i] > infected_last[i] ||
               (infected[i] == infected_first[i] & !boundary.keep) ||
               (infected[i] == infected_last[i] & !boundary.keep)) {
                   infected[i] <- symptom[i] - sample(1:length(mGT$time) - 1, 1, prob = mGT$GT)
                   boundary.keep <- sample(c(TRUE, FALSE), 1)
        }
    }

    infected
}
