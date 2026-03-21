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
#' @format A numeric vector with 10 elements (for a single-season model, \code{S = 1}).
#'   The general length is \code{6 + 4*S} where \code{S} is the number of seasons.
#' \describe{
#'   \item{Elements 1--6 (shared)}{1) random measurement error, 2) 2-fold error,
#'     3) boosting for children (log2), 4) waning for children (log2),
#'     5) boosting for adults (log2), 6) waning for adults (log2).}
#'   \item{Elements 7--9 (per-season)}{infection risk scale parameters for
#'     children, adults, and older adults (3 per season).}
#'   \item{Element 10 (per-season)}{log risk ratio of 2-fold increase in
#'     baseline HAI titer (1 per season).}
#'}
#' @family example_data
"para1"

#' Example of parameter vector for the baseline HAI titer for the main model
#'
#' This is an example of the parameter vector for the baseline HAI titer for the main model used in the \code{seroreconstruct} function. This data frame specifies the format of the parameter vector for the baseline HAI titer for the main model.
#' @docType data
#' @usage data(para2)
#' @format A numeric vector with 20 elements (for a single-season model, \code{S = 1}).
#'   The general length is \code{20 * S} where \code{S} is the number of seasons.
#' \describe{
#'   \item{Elements 1--10}{probability that the HAI titer is 0--9 for children}
#'   \item{Elements 11--20}{probability that the HAI titer is 0--9 for adults}
#' }
#' For multi-season models, this pattern repeats for each season.
#' @family example_data
"para2"
