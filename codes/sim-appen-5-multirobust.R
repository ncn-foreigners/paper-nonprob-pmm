library(nonprobsvy)
library(sampling)
library(doSNOW)
library(progress)
library(foreach)
library(data.table)
library(xtable)

set.seed(2024)

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
  y1 = 1 + x1 * .2 + x2 * 5 + epsilon,
  y2 = -2 + 5 * (x1 - 0.5)^5 + x2 ^ 3 + epsilon,
  p1 = p1,
  p2 = p2,
  base_w_srs = N/n
)

sims <- 500
cores <- 8

cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy", "dbscan"),
               .options.snow = opts,
               .errorhandling = "remove") %dopar% {
  flag_bd1 <- pmin( # planned size ~~ 9K
    rbinom(n = 1:N, size = 1, prob = p1),
    rbinom(n = 1:N, size = 1, prob = p2)
  )
  base_w_bd <- N/n
  
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = population[sample(1:N, n),])

  glm1 <- nonprob(
    outcome = y1 + y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  nn1 <- nonprob(
    outcome = y1 + y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "nn",
    control_outcome = controlOut(k = KK),
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  pmmA <- nonprob(
    outcome = y1 + y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(
      pmm_exact_se = TRUE
    )
  )

  pmmB <- nonprob(
    outcome = y1 + y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(
      pmm_exact_se = TRUE
    )
  )
  
  mm1 <- lm(cbind(y1, y2) ~ x1 + x2, 
            population[flag_bd1 == 1, , drop = FALSE])
  
  mm2 <- lm(cbind(y1, y2) ~ I((x1 - 0.5)^5) + I(x2^3), 
            population[flag_bd1 == 1, , drop = FALSE])
  
  ddf <- cbind(population[flag_bd1 == 1, , drop = FALSE],
               "preds_lin" = predict(mm1),
               "preds_nlin" = predict(mm2))
  
  sample_prob$variables <- cbind(sample_prob$variables, 
                                 "preds_lin" = predict(mm1, sample_prob$variables),
                                 "preds_nlin" = predict(mm2, sample_prob$variables))
  
  pmmA_y1_mv <- nonprob(
    outcome = y1 ~ preds_lin.y1 + preds_nlin.y1 - 1,
    data = ddf,
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(
      pmm_exact_se = TRUE
    )
  )

  pmmA_y2_mv <- nonprob(
    outcome = y2 ~ preds_lin.y2 + preds_nlin.y2 - 1,
    data = ddf,
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(
      pmm_exact_se = TRUE
    )
  )
  
  pmmB_y1_mv <- nonprob(
    outcome = y1 ~ preds_lin.y1 + preds_nlin.y1 - 1,
    data = ddf,
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(
      pmm_exact_se = TRUE
    )
  )

  pmmB_y2_mv <- nonprob(
    outcome = y2 ~ preds_lin.y2 + preds_nlin.y2 - 1,
    data = ddf,
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(
      pmm_exact_se = TRUE
    )
  )
  
  ## save results
  
  data.frame(
    k = k,
    y = c("y1", "y2"),
    trues =  c(mean(population$y1), mean(population$y2)),
    glm =    glm1$output$mean,
    nn =   nn1$output$mean,
    pmmA =   pmmA$output$mean,
    pmmA_mv = c(pmmA_y1_mv$output$mean, pmmA_y2_mv$output$mean),
    pmmB =   pmmB$output$mean,
    pmmB_mv =  c(pmmB_y1_mv$output$mean, pmmB_y2_mv$output$mean),
    
    glm_ci = c(glm1$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < glm1$confidence_interval[1, 2],
               glm1$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < glm1$confidence_interval[2, 2]),
    
    nn_ci = c(nn1$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < nn1$confidence_interval[1, 2],
               nn1$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < nn1$confidence_interval[2, 2]),
    
    pmmA_ci = c(pmmA$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < pmmA$confidence_interval[1, 2],
                pmmA$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < pmmA$confidence_interval[2, 2]),

    pmmA_mv_ci = c(pmmA_y1_mv$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmmA_y1_mv$confidence_interval[, 2],
                   pmmA_y2_mv$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmmA_y2_mv$confidence_interval[, 2]),

    pmmB_ci =  c(pmmB$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < pmmB$confidence_interval[1, 2],
                 pmmB$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < pmmB$confidence_interval[2, 2]),

    pmmB_mv_ci = c(pmmB_y1_mv$confidence_interval[, 1] < mean(population$y1) & mean(population$y1) < pmmB_y1_mv$confidence_interval[, 2],
                   pmmB_y2_mv$confidence_interval[, 1] < mean(population$y2) & mean(population$y2) < pmmB_y2_mv$confidence_interval[, 2])
  )
}

stopCluster(cl)

## processing results
setDT(res)

results_simulation1_process <- res |> melt(id.vars = 1:3)
results_simulation1_process[, c("est", "var", "ci"):=tstrsplit(variable, "_")]
results_simulation1_process[var == "ci", ci :="ci"]
results_simulation1_process[var == "ci", var:= NA]

saveRDS(results_simulation1_process, file = "results/sim-appen5-multirobust.RDS")

## just for knn added
# results_simulation1_process_old <- readRDS( "results/sim-appen5-multirobust.RDS")
# results_simulation1_process <- rbind(results_simulation1_process_old, results_simulation1_process)
# saveRDS(results_simulation1_process, "results/sim-appen5-multirobust.RDS")


tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues)), 
                                                 se = sd(value), 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))), 
                                    keyby=.(y, est, var)]

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(y, est, var)]


tab1[tab2, on = c("y", "est", "var")] |>
  xtable(digits = 2) |>
  print.xtable(include.rownames = F)

