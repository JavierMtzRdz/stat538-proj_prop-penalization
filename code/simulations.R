#' To install mytidyfunctions, you need 
if (!require("pacman")) install.packages("pacman")
if (!require("mytidyfunctions")) remotes::install_github("JavierMtzRdz/mytidyfunctions")
pacman::p_load(tidyverse, janitor, scales,
               mytidyfunctions, here, extrafont,
               logbin, brm)

## Set theme and options ------
set_mytheme(text = element_text(family = "Times New Roman"))


# Functions

sigmoid <- rbrm:::sigmoid

generate_data2 <- function(pa, pb, n, n_test, alpha, beta, gamma, misspec_nuisance = FALSE, misspec_propensity = FALSE) {
  
  # Simulate training data
  v.train <- matrix(stats::runif(n * pa, min = -2, max = 2), nrow = n, ncol = pa)
  
  # Generate propensity scores for training data
  if (misspec_propensity) {
    # Misspecified propensity model: use only the first half of covariates
    gamma_true <- c(0, gamma[1:(pa/2)], rep(0, pa/2))
  } else {
    # Correctly specified propensity model: use all covariates
    gamma_true <- c(0, gamma)
  }
  pscore.true <- sigmoid(cbind(rep(1, n), v.train) %*% gamma_true)
  x.train <- stats::rbinom(n, 1, pscore.true)
  
  # Generate outcome probabilities for training data
  if (misspec_nuisance) {
    # Misspecified nuisance model: use only the first half of covariates
    alpha_true <- c(alpha[1:(pa/2)], rep(0, pa/2))
    beta_true <- c(beta[1:(pa/2)], rep(0, pa/2))
  } else {
    # Correctly specified nuisance model: use all covariates
    alpha_true <- alpha
    beta_true <- beta
  }
  p0p1.true <- brm::getProbRR(v.train %*% alpha_true, v.train %*% beta_true)
  pA.true <- p0p1.true[, 1] # P(Y=1|X=0)
  pA.true[x.train == 1] <- p0p1.true[x.train == 1, 2] # P(Y=1|X=1)
  y.train <- stats::rbinom(n, 1, pA.true) # P(Y=1|X)
  
  # Simulate test data
  v.test <- matrix(stats::runif(n_test * pa, min = -2, max = 2), nrow = n_test, ncol = pa)
  
  # Generate propensity scores for test data
  pscore.true.test <- sigmoid(cbind(rep(1, n_test), v.test) %*% gamma_true)
  x.test <- stats::rbinom(n_test, 1, pscore.true.test)
  
  # Generate outcome probabilities for test data
  p0p1.true.test <- brm::getProbRR(v.test %*% alpha_true, v.test %*% beta_true)
  pA.true.test <- p0p1.true.test[, 1] # P(Y=1|X=0)
  pA.true.test[x.test == 1] <- p0p1.true.test[x.test == 1, 2] # P(Y=1|X=1)
  y.test <- stats::rbinom(n_test, 1, pA.true.test) # P(Y=1|X)
  
  # Return training and test data
  return(list(v.train = v.train, x.train = x.train, y.train = y.train, 
              p0.train = p0p1.true[, 1], p1.train = p0p1.true[, 2],
              v.test = v.test, x.test = x.test, y.test = y.test,
              p0.test = p0p1.true.test[, 1], p1.test = p0p1.true.test[, 2]))
}


set.seed(123)
results <- tibble()

for (i in 1:50) {
# Simple setting
samplesize <- c(500)
dimensions <- c(2)

true_alphas <- c(1, 
                 -1)

true_betas <- c(-0.5, 
                1)


true_gammas <- c(0.1, 
                 -0.5)

dat <- generate_data2(
  pa = dimensions,
  pb = dimensions, n = samplesize, n_test = 100,
  alpha = true_alphas,
  beta = true_betas,
  gamma = true_gammas
)

Y <- dat$y.train
A <- dat$x.train
V <- dat$v.train
rr <-  mean(dat$p1.train)/mean(dat$p0.train)
data <- tibble(Y, A, V1 = V[,1], V2 = V[,2] )

poisson_model <- glm(Y ~ A*V, family = "poisson")
poisson_rr <- exp(coef(poisson_model)["A"])

logbin_model <- logbin(Y ~ .,
                       data = data,
                       method = "em")
logbin_rr <- exp(coef(logbin_model)["A"])

mle_model <- brm::brm(Y, A, V, V, param = "RR",
                  est.method = "MLE")

mle_p <- colMeans(predict(mle_model))

mle_rr <-  mle_p[2]/mle_p[1]

dr_model <- brm::brm(Y, A, V, V, param = "RR",
                      est.method = "DR")

dr_p <- colMeans(predict(dr_model))

dr_rr <-  dr_p[2]/dr_p[1]

set_1 <- tibble(sim = i, setting = "Correct",
                rr = rr, 
                poisson_rr, logbin_rr , mle_rr, dr_rr)



# Missespacification setting
samplesize <- c(500)
dimensions <- c(3)

true_alphas <- c(1, 
                 -1, 
                 0)

true_betas <- c(-0.5, 
                1, 
                0)


true_gammas <- c(0.1, 
                 -0.5, 
                 0)


dat <- generate_data2(
  pa = dimensions,
  pb = dimensions, n = samplesize, n_test = 100,
  alpha = true_alphas,
  beta = true_betas,
  gamma = true_gammas,
)

Y <- dat$y.train
A <- dat$x.train
V <- dat$v.train#[,c(1, 3)]

rr <-  mean(dat$p1.train)/mean(dat$p0.train)
data <- tibble(Y, A, V1 = V[,1], V2 = V[,2], V3 = V[,2] )

poisson_model <- glm(Y ~ A*V, family = "poisson")
poisson_rr <- exp(coef(poisson_model)["A"])

logbin_model <- logbin(Y ~ .,
                       data = data,
                       method = "em")
logbin_rr <- exp(coef(logbin_model)["A"])

mle_model <- brm::brm(Y, A, V, V, param = "RR",
                      est.method = "MLE")

mle_p <- colMeans(predict(mle_model))

mle_rr <-  mle_p[2]/mle_p[1]

dr_model <- brm::brm(Y, A, V, V, param = "RR",
                     est.method = "DR")

dr_p <- colMeans(predict(dr_model))

dr_rr <-  dr_p[2]/dr_p[1]

set_2 <- tibble(sim = i, setting = "Misspecification",
                rr = rr, 
                poisson_rr, logbin_rr, mle_rr, dr_rr)

# In the analizis the second variable is removed. 

# Setting with extreme values 

samplesize <- c(500)
dimensions <- c(2)

true_alphas <- c(3, 
                 3)

true_betas <- c(5, 
                5)


true_gammas <- c(-3, 
                 3)

dat <- generate_data2(
  pa = dimensions,
  pb = dimensions, n = samplesize, n_test = 100,
  alpha = true_alphas,
  beta = true_betas,
  gamma = true_gammas
)

mean(dat$p1.train)
mean(dat$p0.train)

# In this setting p_1 = 0.09 and p_0 = 0.32

Y <- dat$y.train
A <- dat$x.train
V <- dat$v.train
rr <-  mean(dat$p1.train)/mean(dat$p0.train)
data <- tibble(Y, A, V1 = V[,1], V2 = V[,2] )

poisson_model <- glm(Y ~ A*V, family = "poisson")
poisson_rr <- exp(coef(poisson_model)["A"])

logbin_model <- logbin(Y ~ .,
                       data = data,
                       method = "em")

logbin_rr <- exp(coef(logbin_model)["A"])

mle_model <- brm::brm(Y, A, V, V, param = "RR",
                      est.method = "MLE")

mle_p <- colMeans(predict(mle_model))

mle_rr <-  mle_p[2]/mle_p[1]

dr_model <- brm::brm(Y, A, V, V, param = "RR",
                     est.method = "DR")

dr_p <- colMeans(predict(dr_model))

dr_rr <-  dr_p[2]/dr_p[1]

set_3 <- tibble(sim = i, setting = "Extreme",
                rr = rr, 
                poisson_rr, logbin_rr, mle_rr, dr_rr)

results <- bind_rows(results, 
                     set_1, set_2, set_3)

cli::cli_alert_success("Iteration {i}")

}
results2 <-results %>% 
  mutate(setting = case_match(setting, "Correct" ~ "Correct\nspecification",
                                        "Misspecification" ~ "Misspecification\nsetting",
                                        "Extreme" ~ "Extreme\nProb. Setting",
                              T))

# results <- results %>% 
#   mutate(setting = rep(c("Correct\nspecification", "Misspecification\nsetting",
#                          "Extreme\nProb. Setting"), 50))

saveRDS(results2, "final-proj/data-prcssd/results.rds")

results_rr <- results2 %>% 
  group_by(setting) %>% 
  summarise(value = mean(rr))

results2 %>% 
  pivot_longer(poisson_rr:dr_rr) %>% 
  mutate(name = case_match(name, 
                           "dr_rr" ~ "DR",
                           "mle_rr" ~ "MLE",
                           "logbin_rr" ~ "Log-binomial",
                           "poisson_rr" ~ "Poisson")) %>% 
  ggplot(aes(x = name, y = value, color = name)) +
  geom_boxplot(fill = "transparent",
               show.legend = F) +
  geom_hline(data = results_rr, 
             aes(yintercept = value)) +
  facet_wrap(~setting) 
  lab()
