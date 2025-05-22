#remove.packages("pared")

# Install packages
devtools::install_github("priyamdas2/pared", force = TRUE)
# load library
library(pared)


################################################################################
### JGL ########################################################################
################################################################################

## DEMO: Group JGL

?pared_JGL

# Generate sample
sample_data <- generate_sample(sample_sizes = c(30, 50, 40, 70), rand_seed = 123)
# Finding optimals on Pareto-frontier for 'group' JGL
result <- pared_JGL(sample_list = sample_data, method = "group", Pareto_budget = 50)
# Extract list of optimal tuning parameters
result$summary_table
# Extract interactive figure showing optimal points on Pareto-frontier
result$figure

### Demo: Fused JGL 

# Generate sample
sample_data2 <- generate_sample()
resultFused <- pared_JGL(sample_list = sample_data2, method = "fused", 
                         plot_marker_symbol = 'diamond', plot_marker_color = 'blue', 
                         plot_marker_size = 7, Pareto_budget = 40)

resultFused$summary_table
resultFused$figure


################################################################################
### ENet #######################################################################
################################################################################

# See help for pared_ENet
?pared_ENet


set.seed(1)
p <- 5
X <- matrix(rnorm(100), ncol = p)
n <- dim(X)[1]
beta.true <- matrix(c(1,2), ncol = 1)  # only first few coordinates are non-zero
y <- X[, 1:2] %*% beta.true + rnorm(n)

A <- pared_ENet(X, y, Pareto_budget = 50)
A$summary_table
A$fig


################################################################################
### FLasso #####################################################################
################################################################################

# See help for pared_FLasso
?pared_FLasso

set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta_true <- c(0, 0, 1, 1, 1, 0, 0, 0, 2, 2)
y <- X %*% beta_true + rnorm(n)


result <- pared_FLasso(X, y, Pareto_budget = 80, plot_marker_symbol = 'square', plot_marker_size = 7)
result$summary_table
result$figure


# Fitting fused lasso with desired lambda_1 and lambda_2 values

log_lam1 <- log10(0.146)
log_lam2 <- log10(0.002)

g <- function(beta) FLasso_objective(X, y, beta, log_lam1, log_lam2)

beta_init <- rep(0, p)
opt_result <- optim(par = beta_init, fn = g,  method = "BFGS", control = list(maxit = 1000))
beta_opt <- opt_result$par
print(beta_opt)




