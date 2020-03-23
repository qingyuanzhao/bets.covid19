#' Parse the infected date
#'
#' Deprecated. Used in the previous analysis, now replaced by \code{preprocess.data}.
#'
#' @param data a data frame with the following columns: Infected, Arrived, Symptom, Initial, Confirmed.
#'
#' @return the data frame with two new columns, Infected_first and Infected_last
#'
#' @export
#'
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

    data$Infected_first <- pmax(infected_interval[, 1], 1, na.rm = TRUE)
    data$Infected_last <- pmin(infected_interval[, 2],
                               data$Arrived, data$Symptom,
                               data$Initial, data$Confirmed, na.rm = TRUE)

    data
}


#' Parse infected date (basic)
#'
#' @param infected A string of the form "DATE1" or "DATE1 to DATE2".
#' @return A vector of length 2 for the infection window
#'
parse.one.infected <- function(infected) {
    tmp <- strsplit(as.character(infected), "&")[[1]]
    if (length(tmp) > 1) {
        res <- sapply(tmp, parse.one.infected)
        return(c(res[1, 1], res[nrow(res), 2]))
    }
    tmp <- strsplit(as.character(infected), "to")[[1]]
    if (length(tmp) == 0) { ## if no information
        return(c(-Inf, Inf))
    } else if (length(tmp) == 1) { ## if format is date1
        return(rep(date.process(tmp), 2))
    } else {## if format is date1 to date2
        tmp <- c(tmp[1], tmp[length(tmp)]) ## For an edge case like Hefei-72: 15-Jan to 16-Jan & 19-Jan to 21-Jan
        return(date.process(tmp))
    }
}

#' Prepare data frame for analysis
#'
#' @param data A data frame
#' @param infected.in Either "Wuhan" or "Outside"
#'
#' @return A data frame
#' @details here is a summary of the procedures
#' \enumerate{
#'   \item Convert all dates to number of days since 1-Dec-2019.
#'   \item Restrict to cases with a known symtom onset date.
#'   \item Separates data into those returned from Wuhan and those infected outside of wuhan.
#'   \item Parse column 'Infected' into two columns: Infected_first and Infected_last.
#'   \item For all cases, set Infected_first to 1 if it is missing.
#'   \item For outside cases, set Infected_last to be no later than symptom onset.
#'   \item For Wuhan-exported cases, set Infected_last to no later than symptom onset and arrival.
#' }
#'
#' @author Nianqiao Ju <nju@g.harvard.edu>
#' @examples
#'
#' data <- cases.in.china
#' head(preprocess.data(data))
#'
#' @export
#'
preprocess.data <- function(data, infected_in = c("Wuhan", "Outside")){

    infected_in <- match.arg(infected_in, c("Wuhan", "Outside"))

    ## Step 1: Convert all dates to 'number of days since 1-Dec-2019'
    for (column in c("Confirmed", "Arrived", "Symptom", "Initial", "Hospital")) {
        data[[column]] <- date.process(data[[column]])
    }
    ## For visualization, remove some columns
    data$Source <- data$Cluster <- data$Death <- data$Discharged <- data$Verified <- data$Note <- NULL
    ## Initialize columns
    data$Infected_last  <-  data$Infected_first <- NA

    ## Step 2: Considers only cases with a known symtom onset date
    print(paste("Removing", sum(is.na(data$Symptom)), "cases with unknown symptom onset dates:", sum(!is.na(data$Symptom)), "cases left."))
    data <- subset(data, !is.na(Symptom))
    data$Infected <- as.character(data$Infected)

    ## Step 3: separates data into two parts: Wuhan-exported and local transmitted
    data <- switch(infected_in,
                   Wuhan = subset(data, Outside == '' ),
                   Outside = subset(data, Outside != '' ))
    print(paste0("Only keeping cases who were infected in ", infected_in, ": ", nrow(data), " cases left."))

    ## Step 4: Parse column 'Infected' into two columns: Infected_first and Infected_last.
    infected_interval <- t(sapply(data$Infected, parse.one.infected))

    ## Step 5: For all cases, set Infected_first to 1 if it is missing.
    data$Infected_first <- pmax(infected_interval[, 1], 1 , na.rm = TRUE)

    ## Step 6: For outside cases, set Infected_last to be no greater than Symptom onset date.
    if (infected_in == "Outside") {
        data$Infected_last <- pmin(infected_interval[ ,2], data$Symptom, na.rm = TRUE)
    }

    ## Step 7: For Wuhan-exported cases, set Infected_last to be no greater than min(Symptom, Arrived)
    if (infected_in == "Wuhan") {
        data$Infected_last <- pmin(infected_interval[ ,2], data$Symptom, data$Arrived, na.rm = TRUE)
    }

    return(data)
}
