# `pared`: Pareto-based multi-objective tuning for statistical models

`pared` is an R package for model selection and hyperparameter tuning using multi-objective optimization. Instead of reducing model selection to a single criterion, `pared` identifies Pareto-optimal solutions that represent trade-offs among competing objectives, such as prediction error, sparsity, model complexity, coefficient magnitude, graph structure, and goodness of fit.

The package is built around six main user-facing functions.

---
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



## Main functions

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

## Installation

## Basic Pareto Front Example

## General Pareto Front Example: SVM

## Joint Graphical LASSO

## Elastic-Net

## Fused LASSO

## Case-study: Fitting JGL to cancer proteomics dataset
