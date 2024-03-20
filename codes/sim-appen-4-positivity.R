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
x3 <- rnbinom(n = N, mu = 10, size = 4)+18
x3 <- as.factor(ifelse(x3 > 40, 40, x3))

sigma <- diag(2, nrow = 5)
sigma[upper.tri(sigma)] <- runif(n = (5^2 - 5) / 2, max = 1, min = -.7)
sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
epsilon <- MASS::mvrnorm(n = N / 5, mu = rep(0, 5), Sigma = sigma) |> as.vector()

x3_effect <- runif(n = NROW(unique(x3)), min = -6, max = 10)

p1 <- exp(x2 - x1 - 2) / (1 + exp(x2 - x1 - 2))
p2 <- exp(x1 * .6 - x2 - 2) / (1 + exp(x1 * .6 - x2 - 2))

x3_dat <- model.matrix(~ x3 - 1, data.frame(x3 = x3), 
                       contrasts.arg = list("x3" = contrasts(x3, FALSE)))

population <- data.frame(
  x1,
  x2,
  x3,
  y1 = 6 * x1 - 5 * x2 - 7 + x3_dat %*% x3_effect + 15 * epsilon,
  y2 = -2 + .37 * (x1 - 0.5) ^ 2 + x2 ^ 2 + x3_dat %*% x3_effect + 5 * epsilon,
  p1 = p1,
  p2 = p2,
  base_w_srs = N/n
)

# "samplable" population expected size ~~ 600_000

xx1 <- rbinom(prob = plogis(
  x3_dat %*% seq(from = .75, to = .1,
                 length.out = NROW(unique(x3)))
), size = 1, n = N) |>
  as.logical()

# "samplable" population of size 75% of original population

xx2 <- x3_dat %*% seq(
  from = .75, to = .1, 
  length.out = NROW(unique(x3))
) + model.matrix(~x1 + x2, population) %*% runif(n = 3)
xx2 <- as.vector(xx2)
xx2 <- xx2 > quantile(xx2, .25)


# stochastic
population1 <- population[xx1, ]
# deterministic
population2 <- population[xx2, ]

sims <- 500
cores <- 8
cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy"),
               .options.snow = opts) %dopar% {
  # planned size ~~ .9% of size "samplable" population (population1)
  flag_bd1 <- pmax(
    rbinom(n = 1:NROW(population1), size = 1, prob = population1$p1),
    rbinom(n = 1:NROW(population1), size = 1, prob = population1$p2)
  )
  # planned size ~~ .8% of size "samplable" population (population2)
  flag_bd2 <- pmax(
    rbinom(n = 1:NROW(population2), size = 1, prob = population2$p1),
    rbinom(n = 1:NROW(population2), size = 1, prob = population2$p2)
  )
  base_w_bd <- N/n
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = population[sample(1:N, n),])
  
  ## stoch
  glm_s <- nonprob(
    outcome = y1 + y2 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  pmmA_s <- nonprob(
    outcome = y1 + y2 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmmB_s <- nonprob(
    outcome = y1 + y2 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## determ
  glm_d <- nonprob(
    outcome = y1 + y2 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  pmmA_d <- nonprob(
    outcome = y1 + y2 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmmB_d <- nonprob(
    outcome = y1 + y2 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  data.frame(
    k = k,
    y = c("y1", "y2", "y1", "y2"),
    type = rep(c("s", "d"), each = 2),
    trues = c(mean(population$y1), mean(population$y2), mean(population$y1), mean(population$y2)),
    naive = c(mean(population1$y1[flag_bd1 == 1]), mean(population1$y2[flag_bd1 == 1]), 
              mean(population2$y1[flag_bd2 == 1]), mean(population2$y2[flag_bd2 == 1])),
    glm =   c(glm_s$output$mean, glm_d$output$mean),
    pmmA =  c(pmmA_s$output$mean, pmmA_d$output$mean),
    pmmB =  c(pmmB_s$output$mean, pmmB_d$output$mean),
    
    glm_ci = c(glm_s$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < glm_s$confidence_interval[1, 2],
               glm_s$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < glm_s$confidence_interval[2, 2],
               glm_d$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < glm_d$confidence_interval[1, 2],
               glm_d$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < glm_d$confidence_interval[2, 2]),
    
    pmmA_ci = c(pmmA_s$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < pmmA_s$confidence_interval[1, 2],
                pmmA_s$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < pmmA_s$confidence_interval[2, 2],
                pmmA_d$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < pmmA_d$confidence_interval[1, 2],
                pmmA_d$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < pmmA_d$confidence_interval[2, 2]),
    
    pmmB_ci = c(pmmB_s$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < pmmB_s$confidence_interval[1, 2],
                pmmB_s$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < pmmB_s$confidence_interval[2, 2],
                pmmB_d$confidence_interval[1, 1] < mean(population$y1) & mean(population$y1) < pmmB_d$confidence_interval[1, 2],
                pmmB_d$confidence_interval[2, 1] < mean(population$y2) & mean(population$y2) < pmmB_d$confidence_interval[2, 2])
  )
}

stopCluster(cl)

## processing results
setDT(res)

results_simulation1_process <- res |> melt(id.vars = 1:4)
results_simulation1_process[, c("est", "ci"):=tstrsplit(variable, "_")]
results_simulation1_process[ci == "ci", ci :="ci"]

saveRDS(results_simulation1_process, file = "results/sim-appen4-positivity.RDS")

tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues))*100, 
                                                 se = sd(value)*100, 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))*100), 
                                    keyby=.(y, type, est)]

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(y, type, est)]


tab1[tab2, on = c("y", "type", "est")] |>
  xtable(digits = 2) |>
  print.xtable(include.rownames = F)

