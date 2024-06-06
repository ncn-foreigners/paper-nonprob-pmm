set.seed(seed_for_sim)

N <- 1e5
n <- 500
KK <- 5

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
  y2 = -1.2 + (x1 - 0.5) ^ 2 + atan(x2) ^ (3 + sin(x1 + x2)) + sin(x1) * cos(x2) + epsilon, 
  y3 = x1 * x2 * epsilon, # multiplicative
  p1 = p1,
  p2 = p2,
  base_w_srs = N/n
)


cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy"),
               .errorhandling = "remove",
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
  
  ## Y1
  glm1_y1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
  
  nn1_y1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "nn",
    family_outcome = "gaussian",
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  # 
  
  pmmB_y1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmmA_y1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## 
  pmmB_y1_l <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmmA_y1_l <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## Y2
  glm1_y2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
  
  nn1_y2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "nn",
    family_outcome = "gaussian",
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmmB_y2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmmA_y2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  ##
  pmmB_y2_l <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmmA_y2_l <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## Y3
  glm1_y3 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
  
  nn1_y3 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "nn",
    family_outcome = "gaussian",
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmmB_y3 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmmA_y3 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  ##
  pmmB_y3_l <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmmA_y3_l <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  data.frame(
    k = k,
    y = c("y1", "y2", "y3"),
    trues =  c(mean(population$y1), mean(population$y2), mean(population$y3)),
    glm =    c(glm1_y1$output$mean, glm1_y2$output$mean, glm1_y2$output$mean),
    nn =    c(nn1_y1$output$mean, nn1_y2$output$mean, nn1_y2$output$mean),
    pmmA =   c(pmmA_y1$output$mean, pmmA_y2$output$mean, pmmA_y3$output$mean),
    pmmA_l = c(pmmA_y1_l$output$mean, pmmA_y2_l$output$mean, pmmA_y3_l$output$mean),
    pmmB =   c(pmmB_y1$output$mean, pmmB_y2$output$mean, pmmB_y3$output$mean),
    pmmB_l =  c(pmmB_y1_l$output$mean, pmmB_y2_l$output$mean, pmmB_y3_l$output$mean),
    
    glm_ci = c(glm1_y1$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < glm1_y1$confidence_interval[, 2],
               glm1_y2$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < glm1_y2$confidence_interval[, 2],
               glm1_y3$confidence_interval[, 1] < mean(population$y3) & mean(population$y3) < glm1_y3$confidence_interval[, 2]),
    
    nn_ci = c(nn1_y1$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < nn1_y1$confidence_interval[, 2],
               nn1_y2$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < nn1_y2$confidence_interval[, 2],
               nn1_y3$confidence_interval[, 1] < mean(population$y3) & mean(population$y3) < nn1_y3$confidence_interval[, 2]),
    
    pmmA_ci = c(pmmA_y1$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmmA_y1$confidence_interval[, 2],
                pmmA_y2$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmmA_y2$confidence_interval[, 2],
                pmmA_y3$confidence_interval[, 1] < mean(population$y3) & mean(population$y3) < pmmA_y3$confidence_interval[, 2]),

    pmmA_l_ci = c(pmmA_y1_l$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmmA_y1_l$confidence_interval[, 2],
                  pmmA_y2_l$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmmA_y2_l$confidence_interval[, 2],
                  pmmA_y3_l$confidence_interval[, 1] < mean(population$y3) & mean(population$y3) < pmmA_y3_l$confidence_interval[, 2]),

    pmmB_ci =  c(pmmB_y1$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmmB_y1$confidence_interval[, 2],
                 pmmB_y2$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmmB_y2$confidence_interval[, 2],
                 pmmB_y3$confidence_interval[, 1] < mean(population$y3) & mean(population$y3) < pmmB_y3$confidence_interval[, 2]),

    pmmB_l_ci = c(pmmB_y1_l$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmmB_y1_l$confidence_interval[, 2],
                  pmmB_y2_l$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmmB_y2_l$confidence_interval[, 2],
                  pmmB_y3_l$confidence_interval[, 1] < mean(population$y3) & mean(population$y3) < pmmB_y3_l$confidence_interval[, 2])
  )
}

stopCluster(cl)


## processing results
setDT(res)

results_simulation1_process <- res |> melt(id.vars = 1:3)
results_simulation1_process[, c("est", "var", "ci"):=tstrsplit(variable, "_")]
results_simulation1_process[var == "ci", ci :="ci"]
results_simulation1_process[var == "ci", var:= NA]

saveRDS(results_simulation1_process, file = "results/sim-appen3-nonparam.RDS")


