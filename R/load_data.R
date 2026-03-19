#' Example of input data
#'
#' This is an example of the input data used in the \code{seroreconstruct} function. This data frame illustrates the format of the input data.
#' @docType data
#' @usage data(inputdata)
#' @format A data frame with 9 variables, where each row represents an individual:
#' \describe{
#'   \item{age_group}{0: children, 1: adults, 2: older adults}
#'   \item{start_time}{start of follow-up}
#'   \item{end_time}{end of follow-up}
#'   \item{time1}{date of first serum collection}
#'   \item{time2}{date of second serum collection}
#'   \item{time3}{date of third serum collection}
#'   \item{HAI_titer_1}{HAI titer for first serum collection}
#'   \item{HAI_titer_2}{HAI titer for second serum collection}
#'   \item{HAI_titer3}{HAI titer for third serum collection}
#' }
#' @family inputdata
"inputdata"

#' Example of flu activity data
#'
#' This is an example of the flu activity data used in the \code{seroreconstruct} function. This data frame specifies the format of the flu activity data.
#' @docType data
#' @usage data(flu_activity)
#' @format A data frame with 1 variable, where each row represents a date, and it should match the date in the input data:
#' \describe{
#'   \item{h1.activity}{This is the influenza activity from surveillance data. It can be on a relative scale, as the model includes a scale parameter to estimate infection probability.}
#'}
#' @family example_data
"flu_activity"

#' Example of parameter vector for the main model
#'
#' This is an example of the parameter vector for the main model used in the \code{seroreconstruct} function. This data frame specifies the format of the parameter vector for the main model.
#' @docType data
#' @usage data(para1)
#' @format A vector with 10 elements, where each of them is a model parameter:
#' \describe{
#'   \item{element 1}{the parameter value of the random measurement error}
#'   \item{element 2}{the parameter value of the 2-fold error}
#'   \item{element 3}{the boosting in HAI titer after infection for children (in log2 unit)}
#'   \item{element 4}{the waning in HAI titer for children (in log2 unit)}
#'   \item{element 5}{the boosting in HAI titer after infection for adults (in log2 unit)}
#'   \item{element 6}{the waning in HAI titer for adults (in log2 unit)}
#'   \item{element 7}{the scale parameter for children}
#'   \item{element 8}{the scale parameter for adults}
#'   \item{element 9}{the scale parameter for older adults}
#'   \item{element 10}{the log of risk ratio of 2-fold increase in baseline HAI titer}
#'}
#' @family example_data
"para1"

#' Example of parameter vector for the baseline HAI titer for the main model
#'
#' This is an example of the parameter vector for the baseline HAI titer for the main model used in the \code{seroreconstruct} function. This data frame specifies the format of the parameter vector for the baseline HAI titer for the main model.
#' @docType data
#' @usage data(para2)
#' @format A vector with 10 elements, where each of them is a model parameter:
#' \describe{
#'   \item{element 1-10}{The probability that the HAI titer is equal to 0-9 for children}
#'   \item{element 11-20}{The probability that the HAI titer is equal to 0-9 for adults}
#'}
#' @family example_data
"para2"
