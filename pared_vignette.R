rm(list=ls())

# remove.packages("pared")
# remove.packages("pared2026v1")

# pak::pak("priyamdas2/pared")

library(pared)

help(package = "pared")

?pared_ENet
?pared_FLasso
?pared_JGL
?pared_optimize
?plot_pared_2d
?plot_pared_3d

################################################################################
# 1. Toy example
################################################################################

toy_objective <- function(z) c(Mean = mean(z),
                               Harmonic_mean = length(z) / sum(1 / z),
                               MAD_from_1.5 = mean(abs(z - 1.5)))

set.seed(1)

toy_res <- pared_optimize(objective_fun = toy_objective,
                          lower = c(1, 1, 1),
                          upper = c(2, 2, 2),
                          budget = 30,
                          parameter_names = c("z1", "z2", "z3"),
                          objective_names = c("Mean", "Harmonic mean",
                                              "MAD around 1.5"),
                          objective_directions = c("max", "min", "min"),
                          verbose = FALSE)

print(toy_res$summary_table)

plot_pared_2d(toy_res,
              x_objective = "Mean",
              y_objective = "MAD around 1.5",
              plot_title = "Toy Pareto front")

plot_pared_3d(toy_res,
              x_objective = "Mean",
              y_objective = "Harmonic mean",
              z_objective = "MAD around 1.5",
              plot_title = "Toy Pareto front",
              plot_marker_color = "blue")

################################################################################
# 2. SVM example
################################################################################

library(e1071)
data(iris)

X <- scale(as.matrix(iris[, 1:4]))
y <- as.factor(iris$Species)

set.seed(1)
K <- 5
fold_id <- rep(NA_integer_, length(y))
for (cl in levels(y)) fold_id[y == cl] <- sample(rep(seq_len(K), length.out = sum(y == cl)))
folds <- lapply(seq_len(K), function(k) which(fold_id == k))

svm_objective(c(0.5, -1.5), X, y, folds)

## Run 4-objective Pareto optimization
set.seed(123)

res_svm <- pared_optimize(
  objective_fun = function(theta) svm_objective(theta, X, y, folds),
  lower = c(-2, -4), upper = c(3, 1), budget = 50,
  objective_directions = c("min", "min", "min", "min"),
  objective_names = c("CV_error", "CV_logloss", "Support_fraction", "Training_time"),
  parameter_names = c("log10_cost", "log10_gamma"),
  control = list(), verbose = TRUE
)

## Inspect Pareto-optimal solutions
res_svm$summary_table

## 2D projection of the full 4-objective Pareto result
fig2d <- plot_pared_2d(
  res_svm,
  x_objective = "Support_fraction", y_objective = "Training_time",
  plot_title = "SVM Pareto front:<br>Complexity versus computation",
  xaxis_title = "Supp. vec. frac.", yaxis_title = "Training time",
  plot_marker_color = "rgba(46, 139, 87, 0.7)", plot_marker_symbol = "diamond"
)

fig2d

## 3D projections of the full 4-objective Pareto result
fig3d_1 <- plot_pared_3d(
  res_svm,
  x_objective = "CV_logloss", y_objective = "Support_fraction", z_objective = "Training_time",
  plot_title = "SVM Pareto front:<br>Log-loss, complexity, and computation",
  xaxis_title = "CV log-loss", yaxis_title = "Supp. vec. frac.", zaxis_title = "Training time",
  plot_marker_color = "rgba(46, 139, 87, 0.7)", plot_marker_symbol = "square",
  draw_projection = 0
)

fig3d_1

fig3d_2 <- plot_pared_3d(
  res_svm,
  x_objective = "CV_error", y_objective = "Support_fraction", z_objective = "Training_time",
  plot_title = "SVM Pareto front:<br>Accuracy, complexity, and computation",
  xaxis_title = "CV error", yaxis_title = "Supp. vec. frac.", zaxis_title = "Training time",
  plot_marker_color = "rgba(30, 144, 255, 0.7)",
  plot_marker_symbol = "diamond", plot_marker_size = 10
)

fig3d_2

## Run Pareto optimization with objective-space constraints
set.seed(123)

res_svm_constrained <- pared_optimize(
  objective_fun = function(theta) svm_objective_constrained(
    theta, X, y, folds, max_support_fraction = 0.2, max_training_time = 0.02
  ),
  lower = c(-2, -4), upper = c(3, 1), budget = 50,
  objective_directions = c("min", "min", "min", "min"),
  objective_names = c("CV_error", "CV_logloss", "Support_fraction", "Training_time"),
  parameter_names = c("log10_cost", "log10_gamma"),
  control = list(), verbose = TRUE
)

res_svm_constrained$summary_table

## 3D projection of the constrained 4-objective Pareto result
fig3d_2_constrained <- plot_pared_3d(
  res_svm_constrained,
  x_objective = "CV_error", y_objective = "Support_fraction", z_objective = "Training_time",
  plot_title = "Constrained SVM Pareto front:<br>Accuracy, complexity, and computation",
  xaxis_title = "CV error", yaxis_title = "Supp. vec. frac.", zaxis_title = "Training time",
  plot_marker_color = "rgba(30, 144, 255, 0.7)",
  plot_marker_symbol = "diamond", plot_marker_size = 10
)

fig3d_2
fig3d_2_constrained


################################################################################
# 3. JGL
################################################################################

set.seed(1)

sample_data <- generate_sample(sample_sizes = c(30, 50, 40, 70), rand_seed = 1)

## Group graphical LASSO
jgl_group_res <- pared_JGL(sample_list = sample_data, method = "group", Pareto_budget = 50)
jgl_group_res$summary_table
jgl_group_res$figure

## Fused graphical LASSO
jgl_fused_res <- pared_JGL(sample_data, method = "fused", Pareto_budget = 50,
                           plot_marker_color = "blue", plot_marker_symbol = "diamond")
jgl_fused_res$summary_table
jgl_fused_res$figure

################################################################################
# 4. Elastic Net
################################################################################

set.seed(1)

n <- 100; p <- 20
X <- matrix(rnorm(n * p), nrow = n)
beta_true <- c(2, -1.5, 1, rep(0, p - 3))
y <- as.numeric(X %*% beta_true + rnorm(n))

enet_res <- pared_ENet(X, y, Pareto_budget = 40)

enet_res$summary_table
enet_res$figure


################################################################################
# 5. Fused Lasso
################################################################################


set.seed(123)

n <- 100; p <- 10
X <- matrix(rnorm(n * p), nrow = n)
beta_true <- c(0, 0, 1, 1, 1, 0, 0, 0, 2, 2)
y <- as.numeric(X %*% beta_true + rnorm(n))

flasso_res <- pared_FLasso(X, y, Pareto_budget = 20)

flasso_res$summary_table
flasso_res$figure


