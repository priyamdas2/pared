rm(list=ls())
setwd("U:/pared_2026/Simulation study")
# pak::pak("priyamdas2/pared")
library(pared)

## -------------------------------------------------------------------------
## Generate GGL data with known truth
## -------------------------------------------------------------------------
set.seed(1)

generate_sample_with_truth <- function(sample_sizes = c(100, 200, 50, 150),
                                       rand_seed = 1) {
  
  set.seed(rand_seed)
  
  K <- 4
  p <- 20
  
  A1 <- toeplitz(c(1, 0.5, 0.4, rep(0, p - 3)))
  
  n_edges <- (sum(A1 != 0) - p) / 2
  n_possible <- p * (p - 1) / 2
  
  pos_A1 <- which(upper.tri(A1) & A1 != 0, arr.ind = TRUE)
  zero_A1 <- which(upper.tri(A1) & A1 == 0, arr.ind = TRUE)
  
  pos_inds <- sample(seq_len(nrow(pos_A1)), n_edges)
  zero_inds <- sample(seq_len(nrow(zero_A1)), n_possible - n_edges)
  
  A2 <- A1
  
  for (j in 1:5) {
    
    sign_j <- ifelse(runif(1) > 0.5, 1, -1)
    val_j <- runif(1, 0.4, 0.6) * sign_j
    
    r_add <- zero_A1[zero_inds[j], 1]
    c_add <- zero_A1[zero_inds[j], 2]
    
    A2[r_add, c_add] <- val_j
    A2[c_add, r_add] <- val_j
    
    r_del <- pos_A1[pos_inds[j], 1]
    c_del <- pos_A1[pos_inds[j], 2]
    
    A2[r_del, c_del] <- 0
    A2[c_del, r_del] <- 0
  }
  
  A3 <- A2
  
  for (j in 6:15) {
    
    sign_j <- ifelse(runif(1) > 0.5, 1, -1)
    val_j <- runif(1, 0.4, 0.6) * sign_j
    
    r_add <- zero_A1[zero_inds[j], 1]
    c_add <- zero_A1[zero_inds[j], 2]
    
    A3[r_add, c_add] <- val_j
    A3[c_add, r_add] <- val_j
    
    r_del <- pos_A1[pos_inds[j], 1]
    c_del <- pos_A1[pos_inds[j], 2]
    
    A3[r_del, c_del] <- 0
    A3[c_del, r_del] <- 0
  }
  
  A4 <- A3
  
  for (j in 16:20) {
    
    sign_j <- ifelse(runif(1) > 0.5, 1, -1)
    val_j <- runif(1, 0.4, 0.6) * sign_j
    
    r_add <- zero_A1[zero_inds[j], 1]
    c_add <- zero_A1[zero_inds[j], 2]
    
    A4[r_add, c_add] <- val_j
    A4[c_add, r_add] <- val_j
    
    r_del <- pos_A1[pos_inds[j], 1]
    c_del <- pos_A1[pos_inds[j], 2]
    
    A4[r_del, c_del] <- 0
    A4[c_del, r_del] <- 0
  }
  
  A2 <- fix_matrix(A2, 1)
  A3 <- fix_matrix(A3, 1)
  A4 <- fix_matrix(A4, 1)
  
  Precision_true <- array(0, dim = c(p, p, K))
  
  Precision_true[, , 1] <- A1
  Precision_true[, , 2] <- A2
  Precision_true[, , 3] <- A3
  Precision_true[, , 4] <- A4
  
  sample_list <- vector("list", K)
  
  for (k in seq_len(K)) {
    sample_list[[k]] <- mvtnorm::rmvnorm(
      n = sample_sizes[k],
      mean = rep(0, p),
      sigma = solve(Precision_true[, , k])
    )
  }
  
  list(
    sample_list = sample_list,
    Precision_true = Precision_true
  )
}

## -------------------------------------------------------------------------
## Metrics for GGL
## -------------------------------------------------------------------------

upper_edge_vec <- function(A, threshold = 1e-4) {
  abs(A[upper.tri(A)]) > threshold
}

count_total_edges <- function(theta_list, threshold = 1e-4) {
  
  total_edges <- 0
  
  for (k in seq_along(theta_list)) {
    total_edges <- total_edges + sum(upper_edge_vec(theta_list[[k]], threshold))
  }
  
  total_edges
}

count_shared_edges <- function(theta_list, threshold = 1e-4) {
  
  K <- length(theta_list)
  edge_mat <- NULL
  
  for (k in seq_len(K)) {
    edge_k <- upper_edge_vec(theta_list[[k]], threshold)
    
    if (k == 1) {
      edge_mat <- matrix(edge_k, nrow = 1)
    } else {
      edge_mat <- rbind(edge_mat, edge_k)
    }
  }
  
  sum(colSums(edge_mat) == K)
}

theta_list_to_edge_vector <- function(theta_list, threshold = 1e-4) {
  
  out <- NULL
  
  for (k in seq_along(theta_list)) {
    out <- c(out, upper_edge_vec(theta_list[[k]], threshold))
  }
  
  out
}

precision_array_to_edge_vector <- function(Precision_true, threshold = 1e-12) {
  
  K <- dim(Precision_true)[3]
  out <- NULL
  
  for (k in seq_len(K)) {
    out <- c(out, upper_edge_vec(Precision_true[, , k], threshold))
  }
  
  out
}

compute_edge_recovery <- function(theta_list,
                                  Precision_true,
                                  threshold_est = 1e-4,
                                  threshold_true = 1e-12) {
  
  est_edges <- theta_list_to_edge_vector(theta_list, threshold_est)
  true_edges <- precision_array_to_edge_vector(Precision_true, threshold_true)
  
  tp <- sum(est_edges & true_edges)
  fp <- sum(est_edges & !true_edges)
  fn <- sum(!est_edges & true_edges)
  tn <- sum(!est_edges & !true_edges)
  
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  fdr <- ifelse(tp + fp == 0, 0, fp / (tp + fp))
  fpr <- ifelse(fp + tn == 0, 0, fp / (fp + tn))
  
  list(
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    FDR = fdr,
    TPR = recall,
    FPR = fpr,
    F1 = f1
  )
}

compute_jgl_aic <- function(theta_list, sample_list, threshold = 1e-4) {
  
  K <- length(sample_list)
  AIC_sum <- 0
  
  for (k in seq_len(K)) {
    
    n_k <- nrow(sample_list[[k]])
    theta_k <- theta_list[[k]]
    
    det_k <- determinant(theta_k, logarithm = TRUE)
    if (det_k$sign <= 0) return(Inf)
    
    logdet_k <- as.numeric(det_k$modulus)
    S_k <- crossprod(sample_list[[k]]) / n_k
    
    trace_term <- sum(diag(S_k %*% theta_k))
    num_edges_k <- sum(upper_edge_vec(theta_k, threshold))
    
    neg2loglik_k <- n_k * (trace_term - logdet_k)
    AIC_k <- neg2loglik_k + 2 * num_edges_k
    
    AIC_sum <- AIC_sum + AIC_k
  }
  
  AIC_sum
}

fit_group_jgl_once <- function(sample_list,
                               log10_lambda1,
                               log10_lambda2,
                               jgl_control = list()) {
  
  lambda1 <- 10^log10_lambda1
  lambda2 <- 10^log10_lambda2
  
  args <- c(
    list(
      Y = sample_list,
      penalty = "group",
      lambda1 = lambda1,
      lambda2 = lambda2,
      return.whole.theta = TRUE
    ),
    jgl_control
  )
  
  fit <- do.call(JGL, args)
  
  fit$theta
}

evaluate_ggl_model <- function(theta_list,
                               sample_list,
                               Precision_true,
                               method_name,
                               log10_lambda1 = NA,
                               log10_lambda2 = NA,
                               runtime_sec = NA) {
  
  aic <- compute_jgl_aic(theta_list, sample_list)
  
  total_edges <- count_total_edges(theta_list)
  shared_edges <- count_shared_edges(theta_list)
  
  rec <- compute_edge_recovery(theta_list, Precision_true)
  
  data.frame(
    Method = method_name,
    log10_lambda1 = log10_lambda1,
    log10_lambda2 = log10_lambda2,
    lambda1 = ifelse(is.na(log10_lambda1), NA, 10^log10_lambda1),
    lambda2 = ifelse(is.na(log10_lambda2), NA, 10^log10_lambda2),
    AIC = aic,
    Total_edges = total_edges,
    Shared_edges = shared_edges,
    TP = rec$TP,
    FP = rec$FP,
    FN = rec$FN,
    FDR = rec$FDR,
    TPR = rec$TPR,
    FPR = rec$FPR,
    F1 = rec$F1,
    Runtime_sec = runtime_sec,
    stringsAsFactors = FALSE
  )
}

## -------------------------------------------------------------------------
## AIC-selected GGL baseline by grid search
## -------------------------------------------------------------------------

fit_aic_selected_ggl <- function(sample_list,
                                 Precision_true,
                                 log10_lambda1_grid = seq(-2, 0, length.out = 6),
                                 log10_lambda2_grid = seq(-2, 0, length.out = 6),
                                 jgl_control = list()) {
  
  grid <- expand.grid(
    log10_lambda1 = log10_lambda1_grid,
    log10_lambda2 = log10_lambda2_grid
  )
  
  grid_results <- vector("list", nrow(grid))
  
  start_time <- Sys.time()
  
  for (i in seq_len(nrow(grid))) {
    
    log1 <- grid$log10_lambda1[i]
    log2 <- grid$log10_lambda2[i]
    
    theta_i <- tryCatch(
      {
        fit_group_jgl_once(
          sample_list = sample_list,
          log10_lambda1 = log1,
          log10_lambda2 = log2,
          jgl_control = jgl_control
        )
      },
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(theta_i)) {
      
      grid_results[[i]] <- data.frame(
        log10_lambda1 = log1,
        log10_lambda2 = log2,
        AIC = Inf,
        Total_edges = NA,
        Shared_edges = NA,
        F1 = NA
      )
      
    } else {
      
      eval_i <- evaluate_ggl_model(
        theta_list = theta_i,
        sample_list = sample_list,
        Precision_true = Precision_true,
        method_name = "grid candidate",
        log10_lambda1 = log1,
        log10_lambda2 = log2
      )
      
      grid_results[[i]] <- eval_i[, c(
        "log10_lambda1",
        "log10_lambda2",
        "AIC",
        "Total_edges",
        "Shared_edges",
        "F1"
      )]
    }
  }
  
  grid_table <- do.call(rbind, grid_results)
  
  best_idx <- which.min(grid_table$AIC)
  
  best_log1 <- grid_table$log10_lambda1[best_idx]
  best_log2 <- grid_table$log10_lambda2[best_idx]
  
  theta_best <- fit_group_jgl_once(
    sample_list = sample_list,
    log10_lambda1 = best_log1,
    log10_lambda2 = best_log2,
    jgl_control = jgl_control
  )
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  best_eval <- evaluate_ggl_model(
    theta_list = theta_best,
    sample_list = sample_list,
    Precision_true = Precision_true,
    method_name = "AIC-selected GGL",
    log10_lambda1 = best_log1,
    log10_lambda2 = best_log2,
    runtime_sec = total_time
  )
  
  list(
    performance = best_eval,
    grid_table = grid_table
  )
}

## -------------------------------------------------------------------------
## pared GGL Pareto candidates and representatives
## -------------------------------------------------------------------------

fit_pared_ggl_with_metrics <- function(sample_list,
                                       Precision_true,
                                       Pareto_budget = 100,
                                       lb = c(-2, -2),
                                       ub = c(0, 0),
                                       aic_tolerance = 0.1,
                                       jgl_control = list()) {
  
  start_time <- Sys.time()
  
  pareto_res <- pared_JGL(
    sample_list = sample_list,
    method = "group",
    lb = lb,
    ub = ub,
    Pareto_budget = Pareto_budget,
    jgl_control = jgl_control,
    draw_projection = 0
  )
  
  pared_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  raw_par <- pareto_res$pared_result$par
  
  pareto_metrics <- vector("list", nrow(raw_par))
  
  for (i in seq_len(nrow(raw_par))) {
    
    log1 <- raw_par[i, 1]
    log2 <- raw_par[i, 2]
    
    theta_i <- fit_group_jgl_once(
      sample_list = sample_list,
      log10_lambda1 = log1,
      log10_lambda2 = log2,
      jgl_control = jgl_control
    )
    
    pareto_metrics[[i]] <- evaluate_ggl_model(
      theta_list = theta_i,
      sample_list = sample_list,
      Precision_true = Precision_true,
      method_name = "pared GGL candidate",
      log10_lambda1 = log1,
      log10_lambda2 = log2,
      runtime_sec = pared_time
    )
  }
  
  pareto_table <- do.call(rbind, pareto_metrics)
  
  min_aic <- min(pareto_table$AIC, na.rm = TRUE)
  aic_cutoff <- min_aic * (1 + aic_tolerance)
  
  eligible <- pareto_table$AIC <= aic_cutoff
  
  if (sum(eligible) == 0) {
    eligible <- rep(TRUE, nrow(pareto_table))
  }
  
  eligible_table <- pareto_table[eligible, , drop = FALSE]
  
  sparse_rep <- eligible_table[which.min(eligible_table$Total_edges), , drop = FALSE]
  sparse_rep$Method <- "pared sparse GGL within 10% AIC"
  
  shared_rep <- eligible_table[which.max(eligible_table$Shared_edges), , drop = FALSE]
  shared_rep$Method <- "pared shared-edge GGL within 10% AIC"
  
  best_f1_rep <- pareto_table[which.max(pareto_table$F1), , drop = FALSE]
  best_f1_rep$Method <- "pared best-F1 GGL on front"
  
  performance <- rbind(
    sparse_rep,
    shared_rep,
    best_f1_rep
  )
  
  set_summary <- data.frame(
    Pareto_points = nrow(pareto_table),
    Min_AIC = min(pareto_table$AIC, na.rm = TRUE),
    Median_AIC = median(pareto_table$AIC, na.rm = TRUE),
    Min_total_edges = min(pareto_table$Total_edges, na.rm = TRUE),
    Median_total_edges = median(pareto_table$Total_edges, na.rm = TRUE),
    Max_shared_edges = max(pareto_table$Shared_edges, na.rm = TRUE),
    Median_shared_edges = median(pareto_table$Shared_edges, na.rm = TRUE),
    Max_F1 = max(pareto_table$F1, na.rm = TRUE),
    Runtime_sec = pared_time
  )
  
  list(
    performance = performance,
    pareto_table = pareto_table,
    set_summary = set_summary,
    pared_result = pareto_res
  )
}

## -------------------------------------------------------------------------
## Run one in-silico GGL study
## -------------------------------------------------------------------------

set.seed(1)

dat <- generate_sample_with_truth()

sample_data <- dat$sample_list
Precision_true <- dat$Precision_true

## AIC-selected grid baseline
ggl_aic <- fit_aic_selected_ggl(
  sample_list = sample_data,
  Precision_true = Precision_true,
  log10_lambda1_grid = seq(-2, 0, length.out = 6),
  log10_lambda2_grid = seq(-2, 0, length.out = 6)
)

## pared GGL
ggl_pared <- fit_pared_ggl_with_metrics(
  sample_list = sample_data,
  Precision_true = Precision_true,
  Pareto_budget = 100,
  lb = c(-2, -2),
  ub = c(0, 0),
  aic_tolerance = 0.1
)

Table_GGL_Performance <- rbind(
  ggl_aic$performance,
  ggl_pared$performance
)

## Round for manuscript-style display
Table_GGL_Performance_Display <- Table_GGL_Performance

num_cols <- sapply(Table_GGL_Performance_Display, is.numeric)

Table_GGL_Performance_Display[num_cols] <-
  lapply(Table_GGL_Performance_Display[num_cols], function(x) round(x, 3))

Table_GGL_Pareto_Set <- ggl_pared$set_summary

num_cols2 <- sapply(Table_GGL_Pareto_Set, is.numeric)

Table_GGL_Pareto_Set[num_cols2] <-
  lapply(Table_GGL_Pareto_Set[num_cols2], function(x) round(x, 3))

print(Table_GGL_Performance_Display)
print(Table_GGL_Pareto_Set)

write.csv(
  Table_GGL_Performance_Display,
  file = "Table_GGL_simulation_performance.csv",
  row.names = FALSE
)

# write.csv(
#   Table_GGL_Pareto_Set,
#   file = "Table_GGL_pareto_set_summary.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   ggl_pared$pareto_table,
#   file = "Raw_GGL_pareto_candidates.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   ggl_aic$grid_table,
#   file = "Raw_GGL_AIC_grid_candidates.csv",
#   row.names = FALSE
# )

## Optional: view the pared GGL figure
ggl_pared$pared_result$figure