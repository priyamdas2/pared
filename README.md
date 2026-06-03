# `pared`: Pareto-based multi-objective tuning for statistical models

`pared` is an R package for model selection and hyperparameter tuning using multi-objective optimization. Instead of reducing model selection to a single criterion, `pared` identifies Pareto-optimal solutions that represent trade-offs among competing objectives, such as prediction error, sparsity, model complexity, coefficient magnitude, graph structure, and goodness of fit.

## 📚 Citation

> Das P, Robinson S, Peterson C (2026).  
> *pared: Model selection using multi-objective optimization*.  
> **ArXiv**: [https://arxiv.org/abs/2505.21730](https://arxiv.org/abs/2505.21730)

---

## Table of Contents
- [Main functions](#main-functions)
- [Installation](#installation)
- [Basic Pareto Front Example](#toy-example)
- [General Pareto Front Example: SVM](#general-pareto-front-example-svm)
- [Joint Graphical LASSO](#joint-graphical-lasso)
- [Elastic-Net](#elastic-net)
- [Fused LASSO](#fused-lasso)
- [Case-study: Fitting JGL to cancer proteomics dataset](#case-study-fitting-jgl-to-cancer-proteomics-dataset)

---

## Main functions

The package is built around six main user-facing functions.

### `pared_optimize()`

`pared_optimize()` is the general-purpose Pareto optimization engine. Users provide a custom objective function, lower and upper bounds for the tuning/search parameters, objective names, and optimization directions. Each objective can be minimized or maximized. The function returns a `pared_result` object containing the Pareto-optimal parameter values, corresponding objective values, and a summary table.

This function can be used for any user-defined tuning problem, including models or objectives not included among the built-in wrappers.

### `pared_ENet()`

`pared_ENet()` applies Pareto-based tuning to Elastic Net regression. The tuning parameters are the Elastic Net mixing parameter `alpha` and the regularization parameter `lambda`, represented internally through `log10(lambda)`. The default objectives are:

- number of non-zero coefficients,
- L2 norm of the coefficient vector,
- deviance.

The function returns an interactive 3D Pareto-front plot, a summary table of Pareto-optimal tuning values, and the full `pared_result` object.

### `pared_FLasso()`

`pared_FLasso()` applies Pareto-based tuning to Fused LASSO regression. The tuning parameters are the LASSO penalty `lambda1` and the fusion penalty `lambda2`, both represented internally on the `log10` scale. The default objectives are:

- number of non-zero coefficients,
- residual sum of squares,
- coefficient roughness.

This wrapper is useful when one wants to balance sparsity, model fit, and smoothness of adjacent regression coefficients.

### `pared_JGL()`

`pared_JGL()` applies Pareto-based tuning to Joint Graphical LASSO models for multiple related groups. The tuning parameters are `lambda1` and `lambda2`, both represented internally on the `log10` scale. The function supports both:

- `method = "group"` for group graphical LASSO,
- `method = "fused"` for fused graphical LASSO.

For group graphical LASSO, the default objectives are total number of edges, shared edges, and AIC. For fused graphical LASSO, the default objectives are total number of edges, precision-matrix MAD, and AIC.

### `plot_pared_2d()`

`plot_pared_2d()` creates an interactive two-dimensional Plotly scatter plot from a `pared_result` object. Users select two objective columns to display on the x- and y-axes. Hover labels show the corresponding tuning-parameter values. This function is useful for visualizing pairwise projections of Pareto fronts, especially when more than two objectives are considered.

### `plot_pared_3d()`

`plot_pared_3d()` creates an interactive three-dimensional Plotly scatter plot from a `pared_result` object. Users select three objective columns to display on the x-, y-, and z-axes. Optional projection lines can be drawn from each Pareto point to the xy-plane. This function is useful for visualizing three-objective Pareto fronts or three-dimensional projections of higher-dimensional Pareto results.

---

## Installation

R package `pared` can be directly installed from GitHub as follows.

```r
# Install packages
devtools::install_github("priyamdas2/pared2026v1", force = TRUE)
# load library
library(pared2026v1)
```
To view the package documentation and function-level help pages, use:

```r
help(package = "pared2026v1")

?pared_ENet
?pared_FLasso
?pared_JGL
?pared_optimize
?plot_pared_2d
?plot_pared_3d
```
---

## Basic Pareto Front Example

## 📈 Basic Pareto Front Example

The main function in `pared` is `pared_optimize()`. It is designed for situations where a user wants to tune one or more parameters, but the quality of a solution is measured using more than one criterion. Instead of combining all criteria into a single number, `pared_optimize()` searches for Pareto-optimal solutions, meaning solutions for which no objective can be improved without making at least one other objective worse.

To use `pared_optimize()`, the user needs to provide an objective function. This objective function should take a numeric vector of tuning/search parameters as input and return a numeric vector containing multiple objective values. Each objective can be set to either `"min"` or `"max"` depending on whether smaller or larger values are preferred.

In this toy example, we optimize a three-dimensional vector

\[
z = (z_1,z_2,z_3),
\]

where each coordinate is constrained to lie between 1 and 2. Therefore, the search space is

\[
1 \leq z_j \leq 2, \qquad j=1,2,3.
\]

We define three objectives:

1. the arithmetic mean of \(z_1,z_2,z_3\), which we want to maximize;
2. the harmonic mean of \(z_1,z_2,z_3\), which we want to minimize;
3. the mean absolute deviation from 1.5, which we want to minimize.

Although this example is not based on a statistical model, it illustrates the required input structure for `pared_optimize()`.

```r
toy_objective <- function(z) c(Mean = mean(z),
                              Harmonic_mean = length(z) / sum(1 / z),
                              MAD_from_1.5 = mean(abs(z - 1.5)))
```

The function `toy_objective()` takes a candidate vector `z` and returns three objective values. The names inside the returned vector are optional, but they make the code easier to read. The important point is that the order of the returned objectives must match the order used in `objective_names` and `objective_directions` below.

Next, we run the Pareto optimization. The arguments `lower` and `upper` define the search region for the three components of `z`. Since `lower = c(1, 1, 1)` and `upper = c(2, 2, 2)`, the optimizer searches over vectors with all three coordinates between 1 and 2.

```r
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
```

Here, `budget = 30` specifies the number of objective-function evaluations used in the Gaussian-process-based Pareto search. The argument `parameter_names` labels the search parameters, while `objective_names` labels the three objectives returned by `toy_objective()`.

The argument `objective_directions` is especially important. In this example,

```r
objective_directions = c("max", "min", "min")
```

means that we want to maximize the first objective, minimize the second objective, and minimize the third objective. Thus, `pared_optimize()` can handle mixed optimization goals across objectives.

The Pareto-optimal solutions can be printed as a summary table.

```r
print(toy_res$summary_table)
```

Each row of `toy_res$summary_table` corresponds to a Pareto-optimal candidate. The table includes both the selected parameter values, `z1`, `z2`, and `z3`, and the corresponding objective values. These are the trade-off solutions identified by the optimizer.

For a pairwise view of the Pareto front, we can plot two objectives at a time. The following plot displays the trade-off between maximizing the mean and minimizing the deviation around 1.5.

```r
plot_pared_2d(toy_res,
              x_objective = "Mean",
              y_objective = "MAD around 1.5",
              plot_title = "Toy Pareto front")
```

The 2D plot is useful when the user wants to inspect a specific pair of objectives. However, since the full example contains three objectives, we can also visualize all three objectives simultaneously using `plot_pared_3d()`.

```r
plot_pared_3d(toy_res,
              x_objective = "Mean",
              y_objective = "Harmonic mean",
              z_objective = "MAD around 1.5",
              plot_title = "Toy Pareto front",
              plot_marker_color = "blue")
```

This 3D plot shows the Pareto-optimal trade-offs among all three objectives. Points on the plot correspond to candidate vectors \(z\) that lie on the Pareto front. Hovering over a point displays the corresponding parameter values, allowing the user to inspect which choices of \(z_1,z_2,z_3\) produce each trade-off.

This basic example illustrates the general workflow of `pared_optimize()`:

1. define a search space using `lower` and `upper`;
2. write an objective function that returns multiple criteria;
3. specify whether each criterion should be minimized or maximized;
4. run `pared_optimize()`;
5. inspect the Pareto-optimal solutions numerically and graphically.

---


## General Pareto Front Example: SVM

## Joint Graphical LASSO

## Elastic-Net

## Fused LASSO

## Case-study: Fitting JGL to cancer proteomics dataset
