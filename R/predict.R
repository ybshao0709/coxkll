#' Predict method for objects of class "coxkl"
#'
#' Computes linear predictors (and, if desired, the log partial
#' likelihood) for new observations given a fitted Cox–KL model.
#'
#' @param object  An object of class \code{"coxkl"} returned by
#'                \code{\link{coxkl}}.
#' @param newz    Numeric matrix (n_new × p) of covariates for which
#'                predictions are required.  It must have the same
#'                number and ordering of columns used in training.
#' @param delta   Numeric event indicator vector (1 = event, 0 = censored)
#'                for the new data.  Required if \code{likelihood = TRUE}.
#' @param time    Numeric vector of event or censoring times for the
#'                new data.  Required if \code{likelihood = TRUE}.
#' @param likelihood Logical; if \code{TRUE} (default) the log partial
#'                likelihood is evaluated on the new data with
#'                \code{\link{pl_cal_theta}}.  Set to \code{FALSE} if
#'                you only need the linear predictors.
#' @param ...     Currently ignored; included for S3 compatibility.
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{LP_list}}{List of length \code{length(object$eta_list)}
#'           containing the linear predictors for each \eqn{\eta}.}
#'     \item{\code{beta_list}}{The coefficient vectors carried over from
#'           the fitted object.}
#'     \item{\code{eta_list}}{The \eqn{\eta} grid used in fitting.}
#'     \item{\code{likelihood}}{Numeric vector of log partial likelihood
#'           values on the new data (returned only if
#'           \code{likelihood = TRUE}).}
#'   }
#'
#' @examples
#' 
#' @method predict coxkl
#' @export
#' 
#' data(ExampleData)
#' library(survival)
#' fit <- coxkl(z = ExampleData$z, delta = ExampleData$status,
#'              time = ExampleData$time, beta = ExampleData$beta_external,
#'              eta_list = seq(0, 5, 1))
#' 
#' pred <- predict.coxkl(fit,
#'                       newz  = ExampleData$z,
#'                       delta = ExampleData$status,
#'                       time  = ExampleData$time)
predict.coxkl <- function(object,
                          newz,
                          delta     = NULL,
                          time      = NULL,
                          likelihood = TRUE,
                          ...) {

  if (!inherits(object, "coxkl"))
    stop("`object` must be of class \"coxkl\" (output of coxkl()).")

  if (!is.matrix(newz))
    newz <- as.matrix(newz)

  p_train <- nrow(object$beta_list[[1L]])
  if (ncol(newz) != p_train)
    stop(sprintf("`newz` has %d columns but the model was trained with %d.",
                 ncol(newz), p_train))

  if (likelihood && (is.null(delta) || is.null(time)))
    stop("`delta` and `time` must be supplied when likelihood = TRUE.")

  LP_list <- lapply(object$beta_list,
                    function(b) newz %*% b)

  if (likelihood) {
    lik_vec <- vapply(LP_list,
                      pl_cal_theta,
                      FUN.VALUE = numeric(1L),
                      delta = delta,
                      time  = time)
  } else {
    lik_vec <- NULL
  }

  out <- list(LP_list     = LP_list,
              beta_list   = object$beta_list,
              eta_list    = object$eta_list,
              likelihood  = lik_vec)

  class(out) <- "coxklpred"

  invisible(out)
}
