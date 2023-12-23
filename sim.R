# instalacja z tego branchu
#remotes::install_github("ncn-foreigners/nonprobsvy@dev_2")
library(nonprobsvy)
library(doParallel)
library(foreach)

set.seed(123)

cores <- detectCores()[1] - 1

sims <- 7 * 100
N <- 1e5
n <- 100
KK <- 2

x1 <- rnorm(n = N, mean = 1, sd = 1)
x2 <- rexp(n = N, rate = 1)
epsilon <- rnorm(n = N, sd = 2.5) # rnorm(N)
p1 <- exp(x2)/(1+exp(x2))
population <- data.frame(
  x1,
  x2,
  y1 = 1 + x1 * .2 + x2 * .1 + epsilon,
  y2 = -2 + 0.1 * (x1 - 0.5)^2 + x2 * .5 + epsilon,
  p1 = p1,
  base_w_srs = N/n
)

cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoParallel(cl)

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy"),
               .verbose = TRUE) %dopar% {
  flag_srs <- rbinom(n = N, size = 1, prob = n / N)
  flag_bd1 <- rbinom(n = N, size = 1, prob = p1)
  base_w_bd <- N/sum(flag_bd1)
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = subset(population, flag_srs == 1))

  pmm1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1)
  )

  pmm1.1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2)
  )

  pmm1.2 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmm1.3 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmm2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1)
  )

  pmm2.1 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2)
  )

  pmm2.2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmm2.3 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  # we need some extractors for this
  # quite unsightly that call

  cbind(
    pmm1$output$mean, pmm1.1$output$mean,
    pmm1.2$output$mean, pmm1.3$output$mean,
    pmm1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1$confidence_interval[, 2],
    pmm1.1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.1$confidence_interval[, 2],
    pmm1.2$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.2$confidence_interval[, 2],
    pmm1.3$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.3$confidence_interval[, 2],
    mean(population$y1),
    pmm2$output$mean, pmm2.1$output$mean,
    pmm2.2$output$mean, pmm2.3$output$mean,
    pmm2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2$confidence_interval[, 2],
    pmm2.1$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.1$confidence_interval[, 2],
    pmm2.2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.2$confidence_interval[, 2],
    pmm2.3$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.3$confidence_interval[, 2],
    mean(population$y2),
    pmm1$output$SE, pmm1.1$output$SE,
    pmm1.2$output$SE, pmm1.3$output$SE,
    pmm2$output$SE, pmm2.1$output$SE,
    pmm2.2$output$SE, pmm2.3$output$SE
  )
}

stopCluster(cl)
colnames(res) <- c(
  "sample 1 - yhat - y match",#1
  "sample 1 - yhat - yhat match",#2
  "sample 1 - yhat - y match - third term included",#3
  "sample 1 - yhat - yhat match - third term included",#4
  "sample 1 - yhat - y match - coverage",#5
  "sample 1 - yhat - yhat match - coverage",#6
  "sample 1 - yhat - y match - coverage - third term included",#7
  "sample 1 - yhat - yhat match - coverage - third term included",#8
  "sample 1 - mean of y",#9
  "sample 2 - yhat - y match",#10
  "sample 2 - yhat - yhat match",#11
  "sample 2 - yhat - y match - third term included",#12
  "sample 2 - yhat - yhat match - third term included",#13
  "sample 2 - yhat - y match - coverage",#14
  "sample 2 - yhat - yhat match - coverage",#15
  "sample 2 - yhat - y match - coverage - third term included",#16
  "sample 2 - yhat - yhat match - coverage - third term included",#17
  "sample 2 - mean of y",#18
  "sample 1 - est se yhat - y",#19
  "sample 1 - est se yhat - yhat",#20
  "sample 1 - est se yhat - y - third term included",#21
  "sample 1 - est se yhat - yhat - third term included",#22
  "sample 2 - est se yhat - y",#23
  "sample 2 - est se yhat - yhat",#24
  "sample 2 - est se yhat - y - third term included",#25
  "sample 2 - est se yhat - yhat - third term included"#26
)

df <- data.frame(
  bias = c(
    mean(res[,  1] - res[,  9]),
    mean(res[,  2] - res[,  9]),
    mean(res[,  3] - res[,  9]),
    mean(res[,  4] - res[,  9]),
    mean(res[, 10] - res[, 18]),
    mean(res[, 11] - res[, 18]),
    mean(res[, 12] - res[, 18]),
    mean(res[, 13] - res[, 18])
  ),
  mse = c(
    mean((res[,  1] - res[,  9]) ^ 2),
    mean((res[,  2] - res[,  9]) ^ 2),
    mean((res[,  3] - res[,  9]) ^ 2),
    mean((res[,  4] - res[,  9]) ^ 2),
    mean((res[, 10] - res[, 18]) ^ 2),
    mean((res[, 11] - res[, 18]) ^ 2),
    mean((res[, 12] - res[, 18]) ^ 2),
    mean((res[, 13] - res[, 18]) ^ 2)
  ),
  mae = c(
    mean(abs(res[,  1] - res[,  9])),
    mean(abs(res[,  2] - res[,  9])),
    mean(abs(res[,  3] - res[,  9])),
    mean(abs(res[,  4] - res[,  9])),
    mean(abs(res[, 10] - res[, 18])),
    mean(abs(res[, 11] - res[, 18])),
    mean(abs(res[, 12] - res[, 18])),
    mean(abs(res[, 13] - res[, 18]))
  ),
  variance = c(
    var(res[,  1] - res[,  9]),
    var(res[,  2] - res[,  9]),
    var(res[,  3] - res[,  9]),
    var(res[,  4] - res[,  9]),
    var(res[, 10] - res[, 18]),
    var(res[, 11] - res[, 18]),
    var(res[, 12] - res[, 18]),
    var(res[, 13] - res[, 18])
  ),
  se = c(
    sd(res[,  1] - res[,  9]),
    sd(res[,  2] - res[,  9]),
    sd(res[,  3] - res[,  9]),
    sd(res[,  4] - res[,  9]),
    sd(res[, 10] - res[, 18]),
    sd(res[, 11] - res[, 18]),
    sd(res[, 12] - res[, 18]),
    sd(res[, 13] - res[, 18])
  ),
  mean_se_est = c(
    mean(res[, 19]),
    mean(res[, 20]),
    mean(res[, 21]),
    mean(res[, 22]),
    mean(res[, 23]),
    mean(res[, 24]),
    mean(res[, 25]),
    mean(res[, 26])
  ),
  coverage = c(
    mean(res[,  5]),
    mean(res[,  6]),
    mean(res[,  7]),
    mean(res[,  8]),
    mean(res[, 14]),
    mean(res[, 15]),
    mean(res[, 16]),
    mean(res[, 17])
  ),
  row.names = c(
    "sample 1 - yhat - y match",
    "sample 1 - yhat - yhat match",
    "sample 1 - yhat - y match - third term included",
    "sample 1 - yhat - yhat match - third term included",
    "sample 2 - yhat - y match",
    "sample 2 - yhat - yhat match",
    "sample 2 - yhat - y match - third term included",
    "sample 2 - yhat - yhat match - third term included"
  )
)

df
