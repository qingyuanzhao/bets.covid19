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

    tmp <- as.numeric(as.Date(date, "%d-%b") - as.Date("2019-12-01")) + 1

    for (i in which(!is.na(tmp))) {
        if (tmp[i] > 365) {
            tmp[i] <- tmp[i] - 365
        }
    }

    tmp[date == "NO"] <- Inf
    tmp

}
