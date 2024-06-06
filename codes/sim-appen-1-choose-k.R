set.seed(seed_for_sim)
N <- 1e5
n <- 500
KK_pmm2 <- 5
KK_nn2 <- 5

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
population <- data.table(
  x1,
  x2,
  y1 = 1 + x1 * .5 + x2 * .35 + epsilon,
  y2 = -1.2 + (x1 - 0.5) ^ 2 + atan(x2) ^ (3 + sin(x1 + x2)) + sin(x1) * cos(x2) + epsilon, # weird
  p1 = p1,
  p2 = p2,
  base_w_srs = N/n
)


# simulation --------------------------------------------------------------
cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy"),
               .options.snow = opts) %dopar% {
  
  flag_srs <- sample(1:N, n)
  
  flag_bd1 <- pmin(
    rbinom(n = 1:N, size = 1, prob = p1),
    epsilon > quantile(epsilon, .8) |
      quantile(epsilon, .2) > epsilon,
    rbinom(n = 1:N, size = 1, prob = p2)
  )
  
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = population[flag_srs, ])
  

   # variable: y1 ------------------------------------------------------------

  glm_y1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  # NN 
  REP_nn <- TRUE
  KK_nn <- 0
  nn_y1_dyn <- NULL
  while (REP_nn) {
    se_prev <- if (is.null(nn_y1_dyn)) Inf else nn_y1_dyn$output$SE
    nn_y1_dyn <- nonprob(
      outcome = y1 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "nn",
      pop_size = N,
      control_outcome = controlOut(k = KK_nn + 1)
    )
    REP_nn <- nn_y1_dyn$output$SE < se_prev
    KK_nn <- KK_nn + 1
  }

  KK_nn_y1_dyn <- KK_nn - 1
  ## nn with dynamic KK_nn
  nn_y1_dyn <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "nn",
    pop_size = N,
    control_outcome = controlOut(k = KK_nn_y1_dyn)
  )
  
  ## nn with fixed KK_nn2
  nn_y1_fixed <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "nn",
    pop_size = N,
    control_outcome = controlOut(k = KK_nn2)
  )
  
  # PMM_A -- yhat -yhat
  REP_pmm <- TRUE
  KK_pmm <- 0
  pmm_a_y1_dyn <- NULL
  while (REP_pmm) {
    se_prev <- if (is.null(pmm_a_y1_dyn)) Inf else pmm_a_y1_dyn$output$SE
    pmm_a_y1_dyn <- nonprob(
      outcome = y1 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "pmm",
      pop_size = N,
      family_outcome = "gaussian",
      control_outcome = controlOut(k = KK_pmm + 1, predictive_match = 2),
      control_inference = controlInf(pmm_exact_se = TRUE)
    )
    REP_pmm <- pmm_a_y1_dyn$output$SE < se_prev
    KK_pmm <- KK_pmm + 1
  }
  
  KK_pmm_a_y1_dyn <- KK_pmm - 1
  
  ## pmm a (yhat-yhat) with dynamic KK_pmm
  pmm_a_y1_dyn <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK_pmm_a_y1_dyn, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## pmm A (yhat-y) with fixed KK_pmm=5
  pmm_a_y1_fixed <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK_pmm2, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  # PMM_B -- yhat -y
  REP_pmm <- TRUE
  KK_pmm <- 0
  pmm_b_y1_dyn <- NULL
  
  while (REP_pmm) {
    se_prev <- if (is.null(pmm_b_y1_dyn)) Inf else pmm_b_y1_dyn$output$SE
    pmm_b_y1_dyn <- nonprob(
      outcome = y1 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "pmm",
      pop_size = N,
      family_outcome = "gaussian",
      control_outcome = controlOut(k = KK_pmm + 1, predictive_match = 1),
      control_inference = controlInf(pmm_exact_se = TRUE)
    )
    REP_pmm <- pmm_b_y1_dyn$output$SE < se_prev
    KK_pmm <- KK_pmm + 1
  }
  
  KK_pmm_b_y1_dyn <- KK_pmm - 1
  
  ## pmm B (yhat-y) with dynamic KK_pmm
  pmm_b_y1_dyn <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK_pmm_b_y1_dyn, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )


  ## pmm b (yhat-yhat) with fixed KK_pmm=5
  pmm_b_y1_fixed <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK_pmm2, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

   # variable: y2 ------------------------------------------------------------

  glm_y2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  # NN 
  REP_nn <- TRUE
  KK_nn <- 0
  nn_y2_dyn <- NULL
  while (REP_nn) {
    se_prev <- if (is.null(nn_y2_dyn)) Inf else nn_y2_dyn$output$SE
    nn_y2_dyn <- nonprob(
      outcome = y2 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "nn",
      pop_size = N,
      control_outcome = controlOut(k = KK_nn + 1)
    )
    REP_nn <- nn_y2_dyn$output$SE < se_prev
    KK_nn <- KK_nn + 1
  }
  
  KK_nn_y2_dyn <- KK_nn - 1
  ## nn with dynamic KK_nn
  nn_y2_dyn <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "nn",
    pop_size = N,
    control_outcome = controlOut(k = KK_nn_y2_dyn)
  )
  
  ## nn with fixed KK_nn2
  nn_y2_fixed <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "nn",
    pop_size = N,
    control_outcome = controlOut(k = KK_nn2)
  )
  
  # PMM_A -- yhat -yhat
  REP_pmm <- TRUE
  KK_pmm <- 0
  pmm_a_y2_dyn <- NULL
  while (REP_pmm) {
    se_prev <- if (is.null(pmm_a_y2_dyn)) Inf else pmm_a_y2_dyn$output$SE
    pmm_a_y2_dyn <- nonprob(
      outcome = y2 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "pmm",
      pop_size = N,
      family_outcome = "gaussian",
      control_outcome = controlOut(k = KK_pmm + 1, predictive_match = 2),
      control_inference = controlInf(pmm_exact_se = TRUE)
    )
    REP_pmm <- pmm_a_y2_dyn$output$SE < se_prev
    KK_pmm <- KK_pmm + 1
  }
  KK_pmm_a_y2_dyn <- KK_pmm - 1
  
  ## pmm b (yhat-yhat) with dynamic KK_pmm
  pmm_a_y2_dyn <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK_pmm_a_y2_dyn, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## pmm a (yhat-yhat) with fixed KK_pmm=5
  pmm_a_y2_fixed <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK_pmm2, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  # PMM_B -- yhat -y
  REP_pmm <- TRUE
  KK_pmm <- 0
  pmm_b_y2_dyn <- NULL
  
  while (REP_pmm) {
    se_prev <- if (is.null(pmm_b_y2_dyn)) Inf else pmm_b_y2_dyn$output$SE
    pmm_b_y2_dyn <- nonprob(
      outcome = y2 ~ x1 + x2,
      data = population[flag_bd1 == 1, , drop = FALSE],
      svydesign = sample_prob,
      method_outcome = "pmm",
      pop_size = N,
      family_outcome = "gaussian",
      control_outcome = controlOut(k = KK_pmm + 1, predictive_match = 1),
      control_inference = controlInf(pmm_exact_se = TRUE)
    )
    REP_pmm <- pmm_b_y2_dyn$output$SE < se_prev
    KK_pmm <- KK_pmm + 1
  }
  
  KK_pmm_b_y2_dyn <- KK_pmm - 1
  
  ## pmm B (yhat-y) with dynamic KK_pmm
  pmm_b_y2_dyn <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK_pmm_b_y2_dyn, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## pmm b (yhat-y) with fixed KK_pmm=5
  pmm_b_y2_fixed <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK_pmm2, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  # save results ------------------------------------------------------------

  # we need some extractors for this
  # quite unsightly that call
  
  data.frame(
    k = k,
    y = c("y1", "y2"),
    trues      = c(mean(population$y1), mean(population$y2)),
    
    glm        = c(glm_y1$output$mean, glm_y2$output$mean),
    
    nn_dyn     = c(nn_y1_dyn$output$mean, nn_y2_dyn$output$mean),
    nn_dyn_k   = c(KK_nn_y1_dyn, KK_nn_y2_dyn),
    nn_fix     = c(nn_y1_fixed$output$mean, nn_y2_fixed$output$mean),
    
    pmmA_dyn   = c(pmm_a_y1_dyn$output$mean, pmm_a_y2_dyn$output$mean),
    pmmA_dyn_k = c(KK_pmm_a_y1_dyn, KK_pmm_a_y2_dyn),
    pmmA_fix   = c(pmm_a_y1_fixed$output$mean, pmm_a_y2_fixed$output$mean),
    
    pmmB_dyn   = c(pmm_b_y1_dyn$output$mean, pmm_b_y2_dyn$output$mean),
    pmmB_dyn_k = c(KK_pmm_b_y1_dyn, KK_pmm_b_y2_dyn),
    pmmB_fix   = c(pmm_b_y1_fixed$output$mean, pmm_b_y2_fixed$output$mean),
    
    glm_ci = c(
      glm_y1$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < glm_y1$confidence_interval[, 2],
      glm_y2$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < glm_y2$confidence_interval[, 2]
    ),
    
    nn_dyn_ci = c(
      nn_y1_dyn$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < nn_y1_dyn$confidence_interval[, 2],
      nn_y2_dyn$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < nn_y2_dyn$confidence_interval[, 2]
    ),
    
    nn_fix_ci = c(
      pmm_a_y1_fixed$confidence_interval[, 1] < mean(population$y1) &  mean(population$y1) < pmm_a_y1_fixed$confidence_interval[, 2],
      pmm_a_y2_fixed$confidence_interval[, 1] < mean(population$y2) &  mean(population$y2) < pmm_a_y2_fixed$confidence_interval[, 2]
    ),
    
    pmmA_dyn_ci = c(
      pmm_a_y1_dyn$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmm_a_y1_dyn$confidence_interval[, 2],
      pmm_a_y2_dyn$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmm_a_y2_dyn$confidence_interval[, 2]
    ),
    
    pmmA_fix_ci = c(
      pmm_a_y1_fixed$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmm_a_y1_fixed$confidence_interval[, 2],
      pmm_a_y2_fixed$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmm_a_y2_fixed$confidence_interval[, 2]
    ),
    
    pmmB_dyn_ci =  c(
      pmm_b_y1_dyn$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmm_b_y1_dyn$confidence_interval[, 2],
      pmm_b_y2_dyn$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmm_b_y2_dyn$confidence_interval[, 2]
    ),
    pmmB_fix_ci =  c(
      pmm_b_y1_fixed$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmm_b_y1_fixed$confidence_interval[, 2],
      pmm_b_y2_fixed$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmm_b_y2_fixed$confidence_interval[, 2]
    )
  )
}

stopCluster(cl)

## processing results
setDT(res)

results_simulation1_process <- res |> melt(id.vars = 1:3)
results_simulation1_process[, c("est", "sample", "ci"):=tstrsplit(variable, "_")]
results_simulation1_process[sample == "ci", ci := "ci"]
results_simulation1_process[ci == "k", k_sel := "k"]
results_simulation1_process[sample == "ci", sample := NA]
results_simulation1_process[ci == "k", ci := NA]

saveRDS(results_simulation1_process, file = "results/sim-appen1-choose-k-results.RDS")

