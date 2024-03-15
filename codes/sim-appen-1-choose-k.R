library(nonprobsvy)
library(doSNOW)
library(progress)
library(tidyverse)

set.seed(2024)
N <- 1e5
n <- 500
KK2 <- 5

sigma <- diag(1, nrow = 5)
sigma[upper.tri(sigma)] <- runif(n = (5^2 - 5) / 2, max = .5, min = -.5)
sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
x1 <- MASS::mvrnorm(n = N / 5, mu = rep(1, 5), Sigma = sigma) |> as.vector()
x2 <- rexp(n = N, rate = 1)
sigma <- diag(2, nrow = 5)
sigma[upper.tri(sigma)] <- runif(n = (5^2 - 5) / 2, max = 1, min = -.7)
sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
epsilon <- MASS::mvrnorm(n = N / 5, mu = rep(0, 5), Sigma = sigma) |> as.vector()

p1 <- exp(x2)/(1+exp(x2))
p2 <- exp(x1)/(1+exp(x1))
population <- data.frame(
  x1,
  x2,
  y1 = 1 + x1 * .5 + x2 * .35 + epsilon,
  y2 = -1.2 + (x1 - 0.5) ^ 2 + atan(x2) ^ (3 + sin(x1 + x2)) + sin(x1) * cos(x2) + epsilon, # weird
  p1 = p1,
  base_w_srs = N/n
)


# simulation --------------------------------------------------------------
<<<<<<< HEAD
cores <- 8
sims <- 500
=======

cores <- 8
>>>>>>> 356b1ccc6ae55d379c4409fab3cf8130b3aacf22
cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy"),
               .options.snow = opts) %dopar% {
  flag_srs <- rbinom(n = N, size = 1, prob = n / N)
  
  flag_bd1 <- pmin(
    rbinom(n = 1:N, size = 1, prob = p1),
    epsilon > quantile(epsilon, .8) |
      quantile(epsilon, .2) > epsilon,
    rbinom(n = 1:N, size = 1, prob = p2)
  )
  
  base_w_bd <- N/sum(flag_bd1)
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = subset(population, flag_srs == 1))
  
  ## y1
  glm1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  REP <- TRUE
  KK <- 0
  pmm1a_dyn <- NULL
  while (REP) {
    se_prev <- if (is.null(pmm1a_dyn)) Inf else pmm1a_dyn$output$SE
    pmm1a_dyn <- nonprob(
      outcome = y1 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "pmm",
      pop_size = N,
      family_outcome = "gaussian",
      control_outcome = controlOut(k = KK + 1, predictive_match = 1)
    )
    REP <- pmm1a_dyn$output$SE < se_prev
    KK <- KK + 1
  }
  KK <- KK - 1
  ## pmm A (yhat-y) with dynamic kk
  pmm1a_dyn <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1)
  )
  
  
  REP <- TRUE
  KK <- 0
  pmm1b_dyn <- NULL
  while (REP) {
    se_prev <- if (is.null(pmm1b_dyn)) Inf else pmm1b_dyn$output$SE
    pmm1b_dyn <- nonprob(
      outcome = y1 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "pmm",
      pop_size = N,
      family_outcome = "gaussian",
      control_outcome = controlOut(k = KK + 1, predictive_match = 2),
      control_inference = controlInf(pmm_exact_se = TRUE)
    )
    REP <- pmm1b_dyn$output$SE < se_prev
    KK <- KK + 1
  }
  KK <- KK - 1
  
  ## pmm B (yhat-yhat) with dynamic kk
  pmm1b_dyn <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  ## pmm A (yhat-y) with fixed kk=4
  pmm1a_fixed <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK2, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## pmm B (yhat-yhat) with fixed kk=4
  pmm1b_fixed <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK2, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## y2
  glm2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  REP <- TRUE
  KK <- 0
  pmm2a_dyn <- NULL
  while (REP) {
    se_prev <- if (is.null(pmm2a_dyn)) Inf else pmm2a_dyn$output$SE
    pmm2a_dyn <- nonprob(
      outcome = y2 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "pmm",
      pop_size = N,
      family_outcome = "gaussian",
      control_outcome = controlOut(k = KK + 1, predictive_match = 1),
      control_inference = controlInf(pmm_exact_se = TRUE)
    )
    REP <- pmm2a_dyn$output$SE < se_prev
    KK <- KK + 1
  }
  KK <- KK - 1
  ## pmm A (yhat-y) with dynamic k
  pmm2a_dyn <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  REP <- TRUE
  KK <- 0
  pmm2b_dyn <- NULL
  while (REP) {
    se_prev <- if (is.null(pmm2b_dyn)) Inf else pmm2b_dyn$output$SE
    pmm2b_dyn <- nonprob(
      outcome = y2 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "pmm",
      pop_size = N,
      family_outcome = "gaussian",
      control_outcome = controlOut(k = KK + 1, predictive_match = 2),
      control_inference = controlInf(pmm_exact_se = TRUE)
    )
    REP <- pmm2b_dyn$output$SE < se_prev
    KK <- KK + 1
  }
  KK <- KK - 1
  ## pmm B (yhat-yhat) with dynamic k
  pmm2b_dyn <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## pmm A (yhat-y) with fixed k=4
  pmm2a_fixed <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK2, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## pmm B (yhat-yhat) with fixed k=4
  pmm2b_fixed <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK2, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  # we need some extractors for this
  # quite unsightly that call
  
  data.frame(
    k = k,
    y = c("y1", "y2"),
    trues = c(mean(population$y1), mean(population$y2)),
    glm = c(glm1$output$mean, glm2$output$mean),
    pmmA_dyn = c(pmm1a_dyn$output$mean, pmm2a_dyn$output$mean),
    pmmA_fix = c(pmm1a_fixed$output$mean, pmm2a_fixed$output$mean),
    pmmB_dyn = c(pmm1b_dyn$output$mean, pmm2b_dyn$output$mean),
    pmmB_fix = c(pmm1b_fixed$output$mean, pmm2b_fixed$output$mean),
    glm_ci = c(
      glm1$confidence_interval[, 1] < mean(population$y1) & 
      mean(population$y1) < glm1$confidence_interval[, 2],
      glm2$confidence_interval[, 1] < mean(population$y2) & 
      mean(population$y2) < glm2$confidence_interval[, 2]
    ),
    pmmA_dyn_ci = c(
      pmm1a_dyn$confidence_interval[, 1] < mean(population$y1) & 
      mean(population$y1) < pmm1a_dyn$confidence_interval[, 2],
      pmm2a_dyn$confidence_interval[, 1] < mean(population$y2) & 
      mean(population$y2) < pmm2a_dyn$confidence_interval[, 2]
    ),
    pmmA_fix_ci = c(
      pmm1a_fixed$confidence_interval[, 1] < mean(population$y1) & 
      mean(population$y1) < pmm1a_fixed$confidence_interval[, 2],
      pmm2a_fixed$confidence_interval[, 1] < mean(population$y2) & 
      mean(population$y2) < pmm2a_fixed$confidence_interval[, 2]
    ),
    pmmB_dyn_ci =  c(
      pmm1b_dyn$confidence_interval[, 1] < mean(population$y1) & 
      mean(population$y1) < pmm1b_dyn$confidence_interval[, 2], 
      pmm2b_dyn$confidence_interval[, 1] < mean(population$y2) & 
      mean(population$y2) < pmm2b_dyn$confidence_interval[, 2]
    ),
    pmmB_fix_ci =  c(
      pmm1b_fixed$confidence_interval[, 1] < mean(population$y1) & 
      mean(population$y1) < pmm1b_fixed$confidence_interval[, 2],
      pmm2b_fixed$confidence_interval[, 1] < mean(population$y2) & 
      mean(population$y2) < pmm2b_fixed$confidence_interval[, 2]
    )
  )
}

stopCluster(cl)

## processing results
setDT(res)

results_simulation1_process <- res |> melt(id.vars = 1:3)
results_simulation1_process[, c("est", "sample", "ci"):=tstrsplit(variable, "_")]
results_simulation1_process[sample == "ci", ci := "ci"]
results_simulation1_process[sample == "ci", sample := NA]
saveRDS(results_simulation1_process, file = "results/sim-appen1-choose-k-results.RDS")

## stats
tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues))*100, 
                                                 se = sd(value)*100, 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))*100), 
                                    keyby=.(y, est,sample)] 

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(y, est,sample)] 

tab1[tab2, on = c("y", "est", "sample")] |>
  xtable(digits = 4) |>
  print.xtable(include.rownames = FALSE)

