# Setup ----
## Packages to use ----
if (!require("pacman")) install.packages("pacman")
# if (!require("mytidyfunctions")) remotes::install_github("JavierMtzRdz/mytidyfunctions")

pacman::p_load(tidyverse, janitor, 
               SLOPE, glmnet, MASS,
               # mytidyfunctions,
               progress,
               patchwork, here)

## Load fonts ----
extrafont::loadfonts(quiet = TRUE)

## Set theme ------
#mytidyfunctions::set_mytheme(text = element_text(family = "Lato"))


n <- 1000         
# p_values <- c(500, 1000, 2000)
# rho_values <- c(0, 0.5, 0.8)
k_values <- c(10, 20, 50, 100) # Non-zero betas
signal_strengths <- list( 
  weak = list(beta_min = 0.1, beta_max = 0.5)
  #strong = list(beta_min = 0.5, beta_max = 1.5)
)
R <- 15
q_fdr <- 0.05       # q parameter for SLOPE
adapt_lasso_gamma <- 1 # ALasso weights
beta0 <- 0.1
alpha <- 1

generate_data <- function(n, p, rho, k, signal_info, beta0) {
  
  beta_true <- numeric(p)
  
  if (k > 0) {
    non_zero_indices <- sample(1:p, k)
    magnitudes <- runif(k, min = signal_info$beta_min, max = signal_info$beta_max)
    signs <- sample(c(-1, 1), k, replace = TRUE)
    beta_true[non_zero_indices] <- magnitudes * signs
  }
  true_support <- which(beta_true != 0)
  
  # Generate X
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  X <- scale(X)
  
  # Count response 
  lambda <- exp(beta0 + X %*% beta_true)
  # lambda <- pmin(lambda, some_large_value)
  y <- rpois(n, lambda)
  
  return(list(X = X, y = y, beta_true = beta_true, true_support = true_support, beta0 = beta0, lambda = lambda))
}

calculate_metrics <- function(selected_indices, true_indices, p) {
  true_positives <- length(intersect(selected_indices, true_indices))
  false_positives <- length(setdiff(selected_indices, true_indices))
  # Power
  power <- ifelse(length(true_indices) == 0, NA, true_positives / length(true_indices))
  # FDR
  fdr <- ifelse((true_positives + false_positives) == 0, 0, false_positives / (true_positives + false_positives))
  return(list(FDR = fdr, Power = power, TP = true_positives, FP = false_positives, Selected_Count = length(selected_indices)))
}

# Loop ---
results_list <- list()
total_runs <- length(p_values) * length(rho_values) * length(k_values) * length(signal_strengths) * R
cli::cli_progress_bar("Cleaning data", total = total_runs)
iter <- 0


# Loop ---
results_list <- list()
total_runs <- length(p_values) * length(rho_values) * length(k_values) * length(signal_strengths) * R
cli::cli_progress_bar("Cleaning data", total = total_runs)
iter <- 0

for (p in p_values) {
  for (rho in rho_values) {
    set.seed(538)
    for (k in k_values) {
      for (signal_name in names(signal_strengths)) {
        signal_info <- signal_strengths[[signal_name]]
        for (rep in 1:R) {
          cli::cli_progress_update()
          iter <- iter + 1
          
          # Generate Data
          sim_data <- generate_data(n = n, p = p, rho = rho, k = k,
                                    signal_info = signal_info, beta0 = beta0)
          X <- sim_data$X
          y <- sim_data$y
          lambda_val <- sim_data$lambda
          true_support <- sim_data$true_support # Indices of non-zero elements in beta_true
          # Models
          
          # SLOPE
          selected_slope <- integer(0) # Initialize as empty
          fdr_slope <- NA
          power_slope <- NA
          slope_error <- FALSE
          slope_fit <- SLOPE::SLOPE(X, y, family = "poisson", q = q_fdr, lambda = "gaussian", alpha = sqrt(exp(beta0)/n)+0.02,
                                    max_passes = 5000)
          slope_coeffs <- coef(slope_fit)
          
          selected_slope <- which(abs(slope_coeffs[-1]) > 1e-6)
          metrics_slope <- calculate_metrics(selected_slope, true_support, p)
          fdr_slope <- metrics_slope$FDR
          power_slope <- metrics_slope$Power
          
          cli::cli_inform("SLOPE finished")
          
          # SLOPE CV
          #selected_slope_cv <- integer(0) # Initialize as empty
          #fdr_slope_cv <- NA
          #power_slope_cv <- NA
          #slope_cv_error <- FALSE
          #slope_cv_fit <- SLOPE::SLOPE(X, y, family = "poisson", q = q_fdr, lambda = "gaussian", alpha = sqrt(exp(beta0)/n)+0.02,
          #max_passes = 5000)
          #slope_cv_coeffs <- coef(slope_cv_fit)
          #selected_slope_cv <- which(abs(slope_cv_coeffs[-1]) > 1e-6)
          #metrics_slope <- calculate_metrics(selected_slope_cv, true_support, p)
          #fdr_slope_cv <- metrics_slope$FDR
          #power_slope_cv <- metrics_slope$Power
          
          #cli::cli_inform("SLOPE CV finished")
          
          # LASSO 
          selected_lasso <- integer(0)
          fdr_lasso <- NA
          power_lasso <- NA
          lasso_error <- FALSE
          cv_lasso_fit <- cv.glmnet(X, y, family = "poisson", alpha = 1, standardize = FALSE)
          lasso_coeffs <- coef(cv_lasso_fit, s = "lambda.min")
          selected_lasso <- which(lasso_coeffs[-1] != 0) 
          metrics_lasso <- calculate_metrics(selected_lasso, true_support, p)
          fdr_lasso <- metrics_lasso$FDR
          power_lasso <- metrics_lasso$Power
          cli::cli_inform("LASSO finished")
          
          # Adaptive LASSO
          selected_adapt <- integer(0)
          fdr_adapt <- NA
          power_adapt <- NA
          adapt_error <- FALSE
          
          cv_ridge_fit <- cv.glmnet(X, y, family = "poisson", alpha = 0, standardize = FALSE) 
          ridge_coeffs <- coef(cv_ridge_fit, s = "lambda.min")[-1] 
          
          # Calculate weights 
          weights <- 1 / (abs(ridge_coeffs) + .Machine$double.eps)^adapt_lasso_gamma
          weights <- pmin(weights, 1e10) 
          
          # Fit weighted LASSO using CV
          cv_adapt_fit <- cv.glmnet(X, y, family = "poisson", alpha = 1,
                                    penalty.factor = weights, standardize = FALSE)
          adapt_coeffs <- coef(cv_adapt_fit, s = "lambda.min")
          selected_adapt <- which(adapt_coeffs[-1] != 0) 
          
          metrics_adapt <- calculate_metrics(selected_adapt, true_support, p)
          fdr_adapt <- metrics_adapt$FDR
          power_adapt <- metrics_adapt$Power
          
          cli::cli_inform("ALASSO finished")
          
          # Save results
          results_list[[iter]] <- data.frame(
            n = n, p = p, rho = rho, k = k, signal = signal_name, replication = rep,
            Method = c("SLOPE", "LASSO", "AdaptiveLASSO"),
            FDR = c(fdr_slope, fdr_lasso, fdr_adapt),
            Power = c(power_slope, power_lasso, power_adapt),
            SelectedCount = c(length(selected_slope), length(selected_lasso), length(selected_adapt),
                              metrics = c(metrics_slope, metrics_lasso, metrics_adapt))
          )
          write_rds(results_list[[iter]], here("data-prcssd", paste0("results-", iter, "_p-",p,
                                                                     "_rho-", rho, ".rds")))
          cli::cli_inform("Iteration {iter}/{total_runs} finished\n\n {results_list[[iter]]}")
          
        }}}}} 
cli::cli_progress_done()
