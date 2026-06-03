#' Constrained SVM objective function for README examples
#'
#' Internal/example helper exported for convenience.
#'
#' @keywords internal
#' @export
#' @noRd
svm_objective_constrained <- function(theta, X, y, folds,
                                      max_support_fraction = 0.2,
                                      max_training_time = 0.01,
                                      penalty_outside = 1e6) {
  
  out <- svm_objective(theta, X, y, folds)
  
  ## Force output to a numeric vector
  out <- as.numeric(out)
  
  ## Enforce objective names
  names(out) <- c(
    "CV_error",
    "CV_logloss",
    "Support_fraction",
    "Training_time"
  )
  
  ## Objective-space constraints:
  ## 1. Penalize overly complex SVM fits.
  ## 2. Penalize tuning choices that are too computationally expensive.
  if (out["Support_fraction"] > max_support_fraction ||
      out["Training_time"] > max_training_time) {
    
    out["CV_error"] <- penalty_outside
    out["CV_logloss"] <- penalty_outside
    out["Support_fraction"] <- penalty_outside
    out["Training_time"] <- penalty_outside
  }
  
  return(out)
}