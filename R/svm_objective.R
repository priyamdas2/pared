svm_objective <- function(theta, X, y, folds) {
  
  cost  <- 10^theta[1]
  gamma <- 10^theta[2]
  
  t0 <- proc.time()[3]
  
  cv_error     <- numeric(length(folds))
  cv_logloss   <- numeric(length(folds))
  support_frac <- numeric(length(folds))
  
  for (k in seq_along(folds)) {
    
    test_id  <- folds[[k]]
    train_id <- setdiff(seq_len(nrow(X)), test_id)
    
    fit <- e1071::svm(
      x = X[train_id, , drop = FALSE],
      y = y[train_id],
      kernel = "radial",
      cost = cost,
      gamma = gamma,
      probability = TRUE
    )
    
    pred <- predict(
      fit,
      X[test_id, , drop = FALSE],
      probability = TRUE
    )
    
    prob <- attr(pred, "probabilities")
    
    cv_error[k] <- mean(pred != y[test_id])
    
    true_class <- as.character(y[test_id])
    
    ## Defensive probability lookup
    if (is.null(prob) || is.null(colnames(prob))) {
      cv_logloss[k] <- 1e6
    } else {
      
      class_id <- match(true_class, colnames(prob))
      
      if (any(is.na(class_id))) {
        cv_logloss[k] <- 1e6
      } else {
        prob_true <- prob[cbind(seq_along(true_class), class_id)]
        prob_true <- pmax(prob_true, 1e-8)
        cv_logloss[k] <- -mean(log(prob_true))
      }
    }
    
    support_frac[k] <- length(fit$index) / length(train_id)
  }
  
  training_time <- proc.time()[3] - t0
  
  c(
    CV_error = mean(cv_error),
    CV_logloss = mean(cv_logloss),
    Support_fraction = mean(support_frac),
    Training_time = training_time
  )
}
