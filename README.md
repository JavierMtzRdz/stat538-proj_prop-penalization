
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Final Project for STAT 538A (*Generalized Linear Models*)

<!-- badges: start -->
<!-- badges: end -->

This project aims to measure the performance of Sorted L-One Penalized
Estimation (SLOPE) for variable selection in the context of
high-dimensional count data modeled with Poisson regression. SLOPE
extends the LASSO penalty by introducing a sequence of non-increasing
regularization parameters, resulting in the penalty
$\sum_{i=1}^{p}\lambda_{i}|\hat{\beta}_{(i)}|$, where the coefficients
$|\hat{\beta}|$ are sorted in descending order and
$\lambda_1 \ge \dots \ge \lambda_p \ge 0$ (Bogdan et al., 2015). While
SLOPE, inspired by the Benjamini-Hochberg (BH) procedure, has shown
theoretical False Discovery Rate (FDR) control in Gaussian linear
models, its behaviour in non-Gaussian settings, such as Poisson
regression, remains less explored. For this reason, this project aims to
determine how the performance of SLOPE regarding variable selection
accuracy (FDR and Power) compares to LASSO and Adaptive LASSO when
applied to high-dimensional Poisson regression.

This project uses simulations to compare the variable selection accuracy
(measured by FDR and Power) of SLOPE against standard and adaptive LASSO
across scenarios varying in predictor dimensionality ($p/n$ ratio),
inter-predictor correlation ($\rho$), sparsity ($k$), and signal
strength, using a target FDR ($q=0.1$) for SLOPE and cross-validation
for LASSO variants.

Main document in the repository:

- Document
  - [JN source code](/docs/SLOPE-proj_Hotz-Lau-Martinez.ipynb)
  - [HTML doc](/docs/SLOPE-proj_Hotz-Lau-Martinez.html)
- Presentation
  - [Quarto source code](index.qmd)
  - [HTML](https://javiermtzrdz.github.io/stat538-proj_prop-penalization/#/title-slide)
- [References](/refs/refs.bib)
