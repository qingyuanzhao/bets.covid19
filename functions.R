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

#' Transform date to numeric ("23-Jan" to 23)
date.process <- function(date) {
    as.numeric(as.Date(date, "%d-%b") - as.Date("2020-01-01")) + 1
}

#' Simple imputation
#'
simple.impute <- function(data) {

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
multiple.impute <- function(data, m = 50) {
    library(mice)
    data.imputed <- mice(data[, c("Arrive", "Symptom", "Initial", "Hospital", "Confirmed")], m, print = FALSE)

    sapply(1:m, function(i) complete(data.imputed, i)$Symptom)
}
