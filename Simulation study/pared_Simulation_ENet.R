rm(list = ls())
setwd("U:/pared_2026/Simulation study")

library(pared)

## -------------------------------------------------------------------------
## Small ENet simulation design
## -------------------------------------------------------------------------

make_block_cov <- function(p, block_size, rho) {
  
  Sigma <- diag(p)
  n_blocks <- p / block_size
  
  for (b in seq_len(n_blocks)) {
    ind <- ((b - 1) * block_size + 1):(b * block_size)
    Sigma[ind, ind] <- rho
    diag(Sigma)[ind] <- 1
  }
  
  Sigma
}

generate_enet_data <- function(n_train,
                               n_test,
                               p,
                               block_size,
                               rho,
                               snr,
                               seed) {
  
  set.seed(seed)
  
  Sigma <- make_block_cov(p = p, block_size = block_size, rho = rho)
  
  beta_true <- rep(0, p)
  active_set <- c(1, 2, 6, 7, 11, 16)
  beta_true[active_set] <- c(2.0, 1.5, -2.0, -1.5, 1.75, -1.75)
  
  X_train_raw <- MASS::mvrnorm(n_train, mu = rep(0, p), Sigma = Sigma)
  X_test_raw  <- MASS::mvrnorm(n_test,  mu = rep(0, p), Sigma = Sigma)
  
  x_center <- colMeans(X_train_raw)
  x_scale <- apply(X_train_raw, 2, sd)
  x_scale[x_scale == 0] <- 1
  
  X_train <- scale(X_train_raw, center = x_center, scale = x_scale)
  X_test  <- scale(X_test_raw,  center = x_center, scale = x_scale)
  
  signal_train <- as.numeric(X_train %*% beta_true)
  signal_test  <- as.numeric(X_test %*% beta_true)
  
  sigma_eps <- sqrt(var(signal_train) / snr)
  
  y_train <- signal_train + rnorm(n_train, sd = sigma_eps)
  y_test  <- signal_test  + rnorm(n_test,  sd = sigma_eps)
  
  list(
    X_train = X_train,
    y_train = as.numeric(y_train),
    X_test = X_test,
    y_test = as.numeric(y_test),
    beta_true = beta_true,
    active_set = active_set,
    sigma_eps = sigma_eps
  )
}

## -------------------------------------------------------------------------
## Model evaluation
## -------------------------------------------------------------------------

evaluate_enet_model <- function(X_train,
                                y_train,
                                X_test,
                                y_test,
                                beta_true,
                                alpha,
                                lambda,
                                method_name) {
  
  fit <- glmnet::glmnet(
    x = X_train,
    y = y_train,
    alpha = alpha,
    lambda = lambda,
    standardize = FALSE,
    intercept = TRUE
  )
  
  beta_hat <- as.numeric(fit$beta)
  
  pred_train <- as.numeric(predict(fit, newx = X_train, s = lambda))
  pred_test  <- as.numeric(predict(fit, newx = X_test,  s = lambda))
  
  selected <- which(abs(beta_hat) > 1e-6)
  true_active <- which(abs(beta_true) > 0)
  
  tp <- length(intersect(selected, true_active))
  fp <- length(setdiff(selected, true_active))
  fn <- length(setdiff(true_active, selected))
  tn <- length(setdiff(seq_along(beta_true), union(selected, true_active)))
  
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
  fpr <- ifelse(fp + tn == 0, 0, fp / (fp + tn))
  
  data.frame(
    Method = method_name,
    alpha = alpha,
    lambda = lambda,
    Test_MSE = mean((y_test - pred_test)^2),
    Train_MSE = mean((y_train - pred_train)^2),
    Deviance = deviance(fit)[1],
    CV_MSE = NA_real_,
    Model_size = length(selected),
    L2_norm = sqrt(sum(beta_hat^2)),
    TP = tp,
    FP = fp,
    FN = fn,
    FDR = fdr,
    TPR = recall,
    FPR = fpr,
    F1 = f1,
    stringsAsFactors = FALSE
  )
}

## -------------------------------------------------------------------------
## Fixed-candidate CV MSE
## -------------------------------------------------------------------------

compute_enet_cv_mse <- function(X_train,
                                y_train,
                                alpha,
                                lambda,
                                foldid) {
  
  nfolds <- length(unique(foldid))
  fold_mse <- numeric(nfolds)
  
  for (k in seq_len(nfolds)) {
    
    test_id <- which(foldid == k)
    train_id <- setdiff(seq_len(nrow(X_train)), test_id)
    
    fit <- glmnet::glmnet(
      x = X_train[train_id, , drop = FALSE],
      y = y_train[train_id],
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      intercept = TRUE
    )
    
    pred <- as.numeric(
      predict(fit, newx = X_train[test_id, , drop = FALSE], s = lambda)
    )
    
    fold_mse[k] <- mean((y_train[test_id] - pred)^2)
  }
  
  mean(fold_mse)
}

## -------------------------------------------------------------------------
## Conventional ENet competitors
## -------------------------------------------------------------------------

fit_enet_competitors <- function(X_train,
                                 y_train,
                                 X_test,
                                 y_test,
                                 beta_true,
                                 alpha_grid,
                                 nfolds,
                                 ebic_gamma,
                                 seed) {
  
  set.seed(seed)
  
  n <- nrow(X_train)
  p <- ncol(X_train)
  
  foldid <- sample(rep(seq_len(nfolds), length.out = n))
  
  all_candidates <- list()
  row_id <- 1
  
  for (alpha in alpha_grid) {
    
    cvfit <- glmnet::cv.glmnet(
      x = X_train,
      y = y_train,
      alpha = alpha,
      foldid = foldid,
      type.measure = "mse",
      standardize = FALSE,
      intercept = TRUE
    )
    
    lambda_seq <- cvfit$lambda
    
    fit_path <- glmnet::glmnet(
      x = X_train,
      y = y_train,
      alpha = alpha,
      lambda = lambda_seq,
      standardize = FALSE,
      intercept = TRUE
    )
    
    beta_mat <- as.matrix(fit_path$beta)
    pred_train <- predict(fit_path, newx = X_train, s = lambda_seq)
    
    for (j in seq_along(lambda_seq)) {
      
      beta_j <- beta_mat[, j]
      df_j <- sum(abs(beta_j) > 1e-6)
      rss_j <- sum((y_train - pred_train[, j])^2)
      rss_j <- max(rss_j, .Machine$double.eps)
      
      bic_j <- n * log(rss_j / n) + log(n) * df_j
      
      ebic_j <- bic_j
      if (df_j > 0 && df_j < p) {
        ebic_j <- ebic_j + 2 * ebic_gamma * lchoose(p, df_j)
      }
      
      all_candidates[[row_id]] <- data.frame(
        alpha = alpha,
        lambda = lambda_seq[j],
        CV_MSE = cvfit$cvm[j],
        CV_SE = cvfit$cvsd[j],
        df = df_j,
        EBIC = ebic_j
      )
      
      row_id <- row_id + 1
    }
  }
  
  cand <- do.call(rbind, all_candidates)
  
  idx_cv_min <- which.min(cand$CV_MSE)
  cv_min <- cand[idx_cv_min, ]
  
  cv_threshold <- cand$CV_MSE[idx_cv_min] + cand$CV_SE[idx_cv_min]
  cand_1se <- cand[cand$CV_MSE <= cv_threshold, ]
  cand_1se <- cand_1se[order(cand_1se$df, cand_1se$CV_MSE), ]
  cv_1se <- cand_1se[1, ]
  
  ebic_min <- cand[which.min(cand$EBIC), ]
  
  performance <- rbind(
    evaluate_enet_model(
      X_train, y_train, X_test, y_test, beta_true,
      alpha = cv_min$alpha,
      lambda = cv_min$lambda,
      method_name = "CV-min"
    ),
    evaluate_enet_model(
      X_train, y_train, X_test, y_test, beta_true,
      alpha = cv_1se$alpha,
      lambda = cv_1se$lambda,
      method_name = "CV-1SE"
    ),
    evaluate_enet_model(
      X_train, y_train, X_test, y_test, beta_true,
      alpha = ebic_min$alpha,
      lambda = ebic_min$lambda,
      method_name = "EBIC-min"
    )
  )
  
  performance$CV_MSE <- c(
    cv_min$CV_MSE,
    cv_1se$CV_MSE,
    ebic_min$CV_MSE
  )
  
  list(
    performance = performance,
    candidates = cand,
    foldid = foldid
  )
}

## -------------------------------------------------------------------------
## pared ENet Pareto set
## -------------------------------------------------------------------------

fit_pared_enet_sim <- function(X_train,
                               y_train,
                               X_test,
                               y_test,
                               beta_true,
                               Pareto_budget,
                               lower,
                               upper,
                               nfolds,
                               seed,
                               cv_mse_tolerance,
                               deviance_tolerance) {
  
  set.seed(seed)
  
  foldid <- sample(rep(seq_len(nfolds), length.out = nrow(X_train)))
  
  objective_fun <- function(alpha_loglambda) {
    
    alpha <- alpha_loglambda[1]
    lambda <- 10^alpha_loglambda[2]
    
    fit <- glmnet::glmnet(
      x = X_train,
      y = y_train,
      alpha = alpha,
      lambda = lambda,
      standardize = FALSE,
      intercept = TRUE
    )
    
    beta_hat <- as.numeric(fit$beta)
    
    num_nonzero <- sum(abs(beta_hat) > 1e-6)
    l2_norm <- sqrt(sum(beta_hat^2))
    dev <- deviance(fit)[1]
    
    c(num_nonzero, l2_norm, dev)
  }
  
  res <- pared_optimize(
    objective_fun = objective_fun,
    lower = lower,
    upper = upper,
    budget = Pareto_budget,
    parameter_names = c("alpha", "log10_lambda"),
    objective_names = c("Num_nonzero", "L2_norm", "Deviance"),
    objective_directions = c("min", "min", "min"),
    control = list(),
    round_digits = 6,
    verbose = FALSE
  )
  
  pareto_par <- res$par
  pareto_eval <- vector("list", nrow(pareto_par))
  
  for (i in seq_len(nrow(pareto_par))) {
    
    alpha_i <- pareto_par[i, 1]
    lambda_i <- 10^pareto_par[i, 2]
    
    cur_eval <- evaluate_enet_model(
      X_train, y_train, X_test, y_test, beta_true,
      alpha = alpha_i,
      lambda = lambda_i,
      method_name = "pared-candidate"
    )
    
    cur_eval$CV_MSE <- compute_enet_cv_mse(
      X_train = X_train,
      y_train = y_train,
      alpha = alpha_i,
      lambda = lambda_i,
      foldid = foldid
    )
    
    pareto_eval[[i]] <- cur_eval
  }
  
  pareto_metrics <- do.call(rbind, pareto_eval)
  
  min_cv <- min(pareto_metrics$CV_MSE, na.rm = TRUE)
  cv_cutoff <- min_cv * (1 + cv_mse_tolerance)
  
  eligible_cv <- pareto_metrics$CV_MSE <= cv_cutoff
  eligible_cv_metrics <- pareto_metrics[eligible_cv, , drop = FALSE]
  eligible_cv_metrics <- eligible_cv_metrics[
    order(
      eligible_cv_metrics$Model_size,
      eligible_cv_metrics$CV_MSE,
      eligible_cv_metrics$Deviance
    ),
    ,
    drop = FALSE
  ]
  
  sparse_cv <- eligible_cv_metrics[1, , drop = FALSE]
  sparse_cv$Method <- "pared sparse within 10% CV MSE"
  
  min_dev <- min(pareto_metrics$Deviance, na.rm = TRUE)
  dev_cutoff <- min_dev * (1 + deviance_tolerance)
  
  eligible_dev <- pareto_metrics$Deviance <= dev_cutoff
  eligible_dev_metrics <- pareto_metrics[eligible_dev, , drop = FALSE]
  eligible_dev_metrics <- eligible_dev_metrics[
    order(
      eligible_dev_metrics$Model_size,
      eligible_dev_metrics$Deviance,
      eligible_dev_metrics$CV_MSE
    ),
    ,
    drop = FALSE
  ]
  
  sparse_dev <- eligible_dev_metrics[1, , drop = FALSE]
  sparse_dev$Method <- "pared sparse within 10% deviance"
  
  idx_best_test <- which.min(pareto_metrics$Test_MSE)
  best_test <- pareto_metrics[idx_best_test, , drop = FALSE]
  best_test$Method <- "pared oracle best test MSE on front"
  
  idx_best_f1 <- which.max(pareto_metrics$F1)
  best_f1 <- pareto_metrics[idx_best_f1, , drop = FALSE]
  best_f1$Method <- "pared oracle best F1 on front"
  
  set_level <- data.frame(
    Pareto_points = nrow(pareto_metrics),
    Min_CV_MSE = min(pareto_metrics$CV_MSE, na.rm = TRUE),
    Median_CV_MSE = median(pareto_metrics$CV_MSE, na.rm = TRUE),
    Min_test_MSE = min(pareto_metrics$Test_MSE, na.rm = TRUE),
    Median_test_MSE = median(pareto_metrics$Test_MSE, na.rm = TRUE),
    Max_F1 = max(pareto_metrics$F1, na.rm = TRUE),
    Min_model_size = min(pareto_metrics$Model_size, na.rm = TRUE),
    Median_model_size = median(pareto_metrics$Model_size, na.rm = TRUE),
    Max_model_size = max(pareto_metrics$Model_size, na.rm = TRUE)
  )
  
  list(
    performance = rbind(sparse_cv, sparse_dev, best_test, best_f1),
    pareto_all = pareto_metrics,
    set_level = set_level,
    pared_result = res
  )
}

## -------------------------------------------------------------------------
## One replicate
## -------------------------------------------------------------------------

run_one_replicate <- function(rep_id,
                              n_train,
                              n_test,
                              p,
                              block_size,
                              rho,
                              snr,
                              Pareto_budget,
                              alpha_grid,
                              nfolds,
                              ebic_gamma,
                              lower,
                              upper,
                              cv_mse_tolerance,
                              deviance_tolerance) {
  
  dat <- generate_enet_data(
    n_train = n_train,
    n_test = n_test,
    p = p,
    block_size = block_size,
    rho = rho,
    snr = snr,
    seed = 1000 + rep_id
  )
  
  comp <- fit_enet_competitors(
    X_train = dat$X_train,
    y_train = dat$y_train,
    X_test = dat$X_test,
    y_test = dat$y_test,
    beta_true = dat$beta_true,
    alpha_grid = alpha_grid,
    nfolds = nfolds,
    ebic_gamma = ebic_gamma,
    seed = 2000 + rep_id
  )
  
  pared <- fit_pared_enet_sim(
    X_train = dat$X_train,
    y_train = dat$y_train,
    X_test = dat$X_test,
    y_test = dat$y_test,
    beta_true = dat$beta_true,
    Pareto_budget = Pareto_budget,
    lower = lower,
    upper = upper,
    nfolds = nfolds,
    seed = 3000 + rep_id,
    cv_mse_tolerance = cv_mse_tolerance,
    deviance_tolerance = deviance_tolerance
  )
  
  perf <- rbind(comp$performance, pared$performance)
  perf$Replicate <- rep_id
  
  set_level <- pared$set_level
  set_level$Replicate <- rep_id
  
  list(
    performance = perf,
    set_level = set_level
  )
}

## -------------------------------------------------------------------------
## Run simulation
## -------------------------------------------------------------------------

B <- 10
Pareto_budget <- 50

n_train <- 200
n_test <- 200
p <- 20
block_size <- 5
rho <- 0.65
snr <- 2.5

alpha_grid <- seq(0, 1, by = 0.1)
nfolds <- 5
ebic_gamma <- 0.5

lower <- c(0, -3)
upper <- c(1, 1)

cv_mse_tolerance <- 0.10
deviance_tolerance <- 0.10

all_perf <- vector("list", B)
all_set <- vector("list", B)

for (b in seq_len(B)) {
  
  cat("\nRunning replicate", b, "of", B, "\n")
  
  ans_b <- run_one_replicate(
    rep_id = b,
    n_train = n_train,
    n_test = n_test,
    p = p,
    block_size = block_size,
    rho = rho,
    snr = snr,
    Pareto_budget = Pareto_budget,
    alpha_grid = alpha_grid,
    nfolds = nfolds,
    ebic_gamma = ebic_gamma,
    lower = lower,
    upper = upper,
    cv_mse_tolerance = cv_mse_tolerance,
    deviance_tolerance = deviance_tolerance
  )
  
  all_perf[[b]] <- ans_b$performance
  all_set[[b]] <- ans_b$set_level
}

sim_perf <- do.call(rbind, all_perf)
sim_set <- do.call(rbind, all_set)

## -------------------------------------------------------------------------
## Manuscript-style summary tables
## -------------------------------------------------------------------------

mean_sd <- function(x, digits) {
  paste0(
    round(mean(x, na.rm = TRUE), digits),
    " (",
    round(sd(x, na.rm = TRUE), digits),
    ")"
  )
}

summarize_by_method <- function(dat, digits) {
  
  method_order <- c(
    "CV-min",
    "CV-1SE",
    "EBIC-min",
    "pared sparse within 10% CV MSE",
    "pared sparse within 10% deviance",
    "pared oracle best test MSE on front",
    "pared oracle best F1 on front"
  )
  
  out <- lapply(method_order, function(m) {
    
    tmp <- dat[dat$Method == m, , drop = FALSE]
    
    data.frame(
      Method = m,
      Test_MSE = mean_sd(tmp$Test_MSE, digits),
      CV_MSE = mean_sd(tmp$CV_MSE, digits),
      Model_size = mean_sd(tmp$Model_size, digits),
      TP = mean_sd(tmp$TP, digits),
      FP = mean_sd(tmp$FP, digits),
      FDR = mean_sd(tmp$FDR, digits),
      TPR = mean_sd(tmp$TPR, digits),
      F1 = mean_sd(tmp$F1, digits),
      L2_norm = mean_sd(tmp$L2_norm, digits),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

summarize_pareto_set <- function(dat, digits) {
  
  data.frame(
    Quantity = c(
      "Number of Pareto solutions",
      "Minimum CV MSE on Pareto set",
      "Median CV MSE on Pareto set",
      "Minimum test MSE on Pareto set",
      "Median test MSE on Pareto set",
      "Maximum F1 on Pareto set",
      "Minimum model size on Pareto set",
      "Median model size on Pareto set",
      "Maximum model size on Pareto set"
    ),
    Mean_SD = c(
      mean_sd(dat$Pareto_points, digits),
      mean_sd(dat$Min_CV_MSE, digits),
      mean_sd(dat$Median_CV_MSE, digits),
      mean_sd(dat$Min_test_MSE, digits),
      mean_sd(dat$Median_test_MSE, digits),
      mean_sd(dat$Max_F1, digits),
      mean_sd(dat$Min_model_size, digits),
      mean_sd(dat$Median_model_size, digits),
      mean_sd(dat$Max_model_size, digits)
    ),
    stringsAsFactors = FALSE
  )
}

digits <- 3

Table_ENet_Performance <- summarize_by_method(sim_perf, digits)
Table_ENet_Pareto_Set <- summarize_pareto_set(sim_set, digits)

print(Table_ENet_Performance)
print(Table_ENet_Pareto_Set)

write.csv(
  Table_ENet_Performance,
  file = "Table_ENet_simulation_performance.csv",
  row.names = FALSE
)

# write.csv(
#   Table_ENet_Pareto_Set,
#   file = "Table_ENet_pareto_set_summary.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   sim_perf,
#   file = "Raw_ENet_simulation_performance.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   sim_set,
#   file = "Raw_ENet_pareto_set_summary.csv",
#   row.names = FALSE
# )