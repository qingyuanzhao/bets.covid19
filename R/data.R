#' Confirmed cases of COVID-19
#'
#' A dataset containing the trajectory of cases of COVID-19.
#'
#' @format A data frame with 1091 rows and 20 variables:
#' \describe{
#'   \item{Case}{Label of the case, in the format of Country-Case number.}
#'   \item{Nationality/Residence}{Nationality or residence of the patient.}
#'   \item{Gender}{Male (M) or Female (F).}
#'   \item{Age}{Age of the patient, either an integer or age by decade (for example, 40s).}
#'   \item{Cluster}{Other confirmed cases that this patient had contacts with.}
#'   \item{Known Contact}{Whether the case has contact with earlier confirmed cases or visited Hubei province.}
#'   \item{Outside}{Was the patient infected outside Wuhan? Yes (Y), Likely (L), or No (empty string and the default).}
#'   \item{Begin_Wuhan}{Begin of stay in Wuhan.}
#'   \item{End_Wuhan}{End of Stay in Wuhan.}
#'   \item{Infected}{When was the patient infected? Can be an interval or multiple dates.}
#'   \item{Arrived}{When did the patient arrive in the country where he/she was confirmed a 2019-nCoV case?}
#'   \item{Symptom}{When did the patient first show symptoms of 2019-nCoV (cough, fever, fatigue, etc.)?}
#'   \item{Initial}{After developing symptoms, when was the patient first went to (or taken to) a medical institution?}
#'   \item{Hospital}{If the patient was not admitted to or isolated in a hospital after the initial medical visit, when was the patient finally admitted or isolated?}
#'   \item{Confirmed}{When was the patient confirmed as a case of 2019-nCoV?}
#'   \item{Discharged}{When was the patient discharged from hospital?}
#'   \item{Death}{When did the patient die?}
#'   \item{Verified}{Has this information been verified by another data collector?}
#'   \item{Source}{URLs to the information recorded (usually government websites or news reports).}
#' }
"covid19_data"

#' COVID-19 exported from Wuhan
#'
#' Constructed from \code{covid19_data}, see \code{example(preprocess.data)}.
#'
#' @format A data frame with 378 rows and 4 variables:
#' \describe{
#'   \item{Location}{Where the case is confirmed.}
#'   \item{Gender}{Gender of the patient.}
#'   \item{Age}{Age of the patient.}
#'   \item{B}{Beginning of stay in Wuhan.}
#'   \item{E}{End of stay in Wuhan.}
#'   \item{S}{Symptom onset.}
#' }
"wuhan_exported"
