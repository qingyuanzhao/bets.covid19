#' Parse the infected date
#'
#' Deprecated. Used in the previous analysis, now replaced by \code{preprocess.data}.
#'
#' @param data a data frame with the following columns: Infected, Arrived, Symptom, Initial, Confirmed.
#'
#' @return the data frame with two new columns, Infected_first and Infected_last
#'
#' @keywords internal
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
#' @keywords internal
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
#' @param infected_in Either "Wuhan" or "Outside"
#' @param symptom_impute Whether to use initial medical visit and confirmation to impute missing symptom onset.
#'
#' @return A data frame
#' @details A summary of the procedures:
#' \enumerate{
#'   \item Convert all dates to number of days since 1-Dec-2019.
#'   \item Separates data into those returned from Wuhan and those infected outside of wuhan.
#'   \item Restrict to cases with a known symptom onset date.
#'   \item Parse column 'Infected' into two columns: Infected_first and Infected_last.
#'   \item For all cases, set Infected_first to 1 if it is missing.
#'   \item For outside cases, set Infected_last to be no later than symptom onset.
#'   \item For Wuhan-exported cases, set Infected_last to no later than symptom onset and end of Wuhan stay.
#' }
#'
#' @author Nianqiao Ju <nju@g.harvard.edu>, Qingyuan Zhao <qyzhao@statslab.cam.ac.uk>
#' @examples
#'
#' data(covid19_data)
#' head(data <- preprocess.data(covid19_data))
#'
#' \donttest{ ## This is how the wuhan_exported data frame is created
#' data <- subset(data, Symptom < Inf)
#' data <- subset(data, Arrived <= 54)
#' data$Location <- do.call(rbind, strsplit(as.character(data$Case), "-"))[, 1]
#' wuhan_exported <- data.frame(Location = data$Location,
#'                              B = data$Begin_Wuhan,
#'                              E = data$End_Wuhan,
#'                              S = data$Symptom)
#' ## devtools::use_data(wuhan_exported)
#' }
#'
#' @export
#'
preprocess.data <- function(data, infected_in = c("Wuhan", "Outside"), symptom_impute = FALSE){

    infected_in <- match.arg(infected_in, c("Wuhan", "Outside"))

    ## Step 1: Convert all dates to 'number of days since 1-Dec-2019'
    for (column in c("Confirmed", "Begin_Wuhan", "End_Wuhan", "Arrived", "Symptom", "Initial", "Hospital")) {
        data[[column]] <- date.process(data[[column]])
    }
    ## For visualization, remove some columns
    data$Source <- data$Death <- data$Discharged <- data$Verified <- data$Note <- NULL
    ## Initialize columns
    data$Infected_last  <-  data$Infected_first <- NA

    ## Step 2: separates data into two parts: Wuhan-exported and local transmitted
    data <- switch(infected_in,
                   Wuhan = data[data$Outside == '',],
                   Outside = data[data$Outside != '',])
    message(paste0("Only keeping cases who were infected in ", infected_in, ": ", nrow(data), " cases left."))

    ## Step 3: Considers only cases with a known symtom onset date
    if (!symptom_impute) {
        message(paste("Removing", sum(is.na(data$Symptom)), "cases with unknown symptom onset dates:", sum(!is.na(data$Symptom)), "cases left."))
        data <- data[!is.na(data$Symptom), ]
    } else {
        message(paste("Imputing", sum(is.na(data$Symptom)), "cases with unknown symptom onset dates."))
        data$Symptom[is.na(data$Symptom)] <- data$Initial[is.na(data$Symptom)] - 2
        data$Symptom[is.na(data$Symptom)] <- data$Confirmed[is.na(data$Symptom)] - 7
    }

    ## Step 4: Parse column 'Infected' into two columns: Infected_first and Infected_last.
    data$Infected <- as.character(data$Infected)
    infected_interval <- t(sapply(data$Infected, parse.one.infected))

    ## Step 5: For all cases, set Infected_first to 1 if it is missing.
    data$Infected_first <- pmax(infected_interval[, 1], 1 , na.rm = TRUE)

    ## Step 6: For outside cases, set Infected_last to be no greater than Symptom onset date.
    if (infected_in == "Outside") {
        data$Infected_last <- pmin(infected_interval[ ,2], data$Symptom, na.rm = TRUE)
    }

    ## Step 7: For Wuhan-exported cases, set Infected_last to be no greater than min(Symptom, End_Wuhan)
    if (infected_in == "Wuhan") {
        data$Begin_Wuhan[is.na(data$Begin_Wuhan)] <- 0
        data$End_Wuhan[is.na(data$End_Wuhan)] <- 54
        data$Arrived[is.na(data$End_Wuhan)] <- 54
        data$Infected_first <- pmax(infected_interval[ ,1], data$Begin_Wuhan, na.rm = TRUE)
        data$Infected_last <- pmin(infected_interval[ ,2], data$Symptom, data$End_Wuhan, na.rm = TRUE)
    }

    return(data)
}
