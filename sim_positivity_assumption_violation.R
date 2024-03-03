# instalacja z tego branchu
#remotes::install_github("https://github.com/ncn-foreigners/nonprobsvy/tree/dev")
library(nonprobsvy)
library(doSNOW)
library(progress)
library(tidyverse)

set.seed(2051)

cores <- 8

sims <- 5 * 100
N <- 1e6
n <- 200
KK <- 3
x3_break_num <- 30

sigma <- diag(1, nrow = 5)
sigma[upper.tri(sigma)] <- runif(n = (5^2 - 5) / 2, max = .5, min = -.5)
sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
x1 <- MASS::mvrnorm(n = N / 5, mu = rep(1, 5), Sigma = sigma) |> as.vector()
x2 <- rexp(n = N, rate = 1)
x3 <- rnbinom(n = N, mu = 40, size = 6)
x3 <- cut(x3, breaks = quantile(x3, seq(0, 1, length.out = x3_break_num + 1)), include.lowest = TRUE)
sigma <- diag(2, nrow = 5)
sigma[upper.tri(sigma)] <- runif(n = (5^2 - 5) / 2, max = 1, min = -.7)
sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
epsilon <- MASS::mvrnorm(n = N / 5, mu = rep(0, 5), Sigma = sigma) |> as.vector()

x3_effect <- runif(n = x3_break_num, min = -6, max = 10)

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
  base_w_srs = N/n
)

# inclusion probabilities

xx1 <- rbinom(prob = plogis(
                x3_dat %*% seq(from = .75, to = .1,
                               length.out = x3_break_num)
              ), size = 1, n = N) |>
  as.logical()

xx2 <- x3_dat %*% seq(
  from = .75, to = .1, 
  length.out = x3_break_num
) + model.matrix(~x1 + x2, population) %*% runif(n = 3)
xx2 <- as.vector(xx2)
xx2 <- xx2 > quantile(xx2, .25)

# "samplable" population expected size ~~ 600_000
population1 <- population[xx1, ]
# "samplable" population of size 75% of original population
population2 <- population[xx2, ]

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
    rbinom(n = 1:NROW(population1), size = 1, prob = p1[xx1]),
    rbinom(n = 1:NROW(population1), size = 1, prob = p2[xx1])
  )
  # planned size ~~ .8% of size "samplable" population (population2)
  flag_bd2 <- pmax(
    rbinom(n = 1:NROW(population2), size = 1, prob = p1[xx2]),
    rbinom(n = 1:NROW(population2), size = 1, prob = p2[xx2])
  )
  base_w_bd <- N/n
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = population[sample(1:N, n),])
  glm1_stochastic <- nonprob(
    outcome = y1 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  pmm1_stochastic <- nonprob(
    outcome = y1 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmm1.1_stochastic <- nonprob(
    outcome = y1 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  glm2_stochastic <- nonprob(
    outcome = y2 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )

  pmm2_stochastic <- nonprob(
    outcome = y2 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmm2.1_stochastic <- nonprob(
    outcome = y2 ~ x1 + x2 + x3,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  glm1_deterministic <- nonprob(
    outcome = y1 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  pmm1_deterministic <- nonprob(
    outcome = y1 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm1.1_deterministic <- nonprob(
    outcome = y1 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  glm2_deterministic <- nonprob(
    outcome = y2 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  pmm2_deterministic <- nonprob(
    outcome = y2 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm2.1_deterministic <- nonprob(
    outcome = y2 ~ x1 + x2 + x3,
    data = population2[flag_bd2 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  # we need some extractors for this
  # quite unsightly that call

  data.frame(
    "Stochastic undercoverage pmm yhat-y est y1"       = pmm1_stochastic$output$mean, 
    "Stochastic undercoverage pmm yhat-yhat est y1"    = pmm1.1_stochastic$output$mean,
    "Stochastic undercoverage pmm glm est y1"          = glm1_stochastic$output$mean,
    "Deterministic undercoverage pmm yhat-y est y1"    = pmm1_deterministic$output$mean, 
    "Deterministic undercoverage pmm yhat-yhat est y1" = pmm1.1_deterministic$output$mean,
    "Deterministic undercoverage glm est y1"           = glm1_deterministic$output$mean,
    "Stochastic undercoverage pmm yhat-y coverage y1" = 
      pmm1_stochastic$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1_stochastic$confidence_interval[, 2],
    "Stochastic undercoverage pmm yhat-yhat coverage y1" = 
      pmm1.1_stochastic$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.1_stochastic$confidence_interval[, 2],
    "Stochastic undercoverage glm coverage y1" = 
      glm1_stochastic$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < glm1_stochastic$confidence_interval[, 2],
    "Deterministic undercoverage pmm yhat-y coverage y1" = 
      pmm1_deterministic$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1_deterministic$confidence_interval[, 2],
    "Deterministic undercoverage pmm yhat-yhat coverage y1" = 
      pmm1.1_deterministic$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.1_deterministic$confidence_interval[, 2],
    "Deterministic undercoverage glm coverage y1" = 
      glm1_deterministic$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < glm1_deterministic$confidence_interval[, 2],
    "Population mean y1" = mean(population$y1),
    "Stochastic undercoverage pmm yhat-y est y2"       = pmm2_stochastic$output$mean, 
    "Stochastic undercoverage pmm yhat-yhat est y2"    = pmm2.1_stochastic$output$mean,
    "Stochastic undercoverage pmm glm est y2"          = glm2_stochastic$output$mean,
    "Deterministic undercoverage pmm yhat-y est y2"    = pmm2_deterministic$output$mean, 
    "Deterministic undercoverage pmm yhat-yhat est y2" = pmm2.1_deterministic$output$mean,
    "Deterministic undercoverage glm est y2"           = glm2_deterministic$output$mean,
    "Stochastic undercoverage pmm yhat-y coverage y2" = 
      pmm2_stochastic$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2_stochastic$confidence_interval[, 2],
    "Stochastic undercoverage pmm yhat-yhat coverage y2" = 
      pmm2.1_stochastic$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.1_stochastic$confidence_interval[, 2],
    "Stochastic undercoverage glm coverage y2" = 
      glm2_stochastic$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < glm2_stochastic$confidence_interval[, 2],
    "Deterministic undercoverage pmm yhat-y coverage y2" = 
      pmm2_deterministic$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2_deterministic$confidence_interval[, 2],
    "Deterministic undercoverage pmm yhat-yhat coverage y2" = 
      pmm2.1_deterministic$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.1_deterministic$confidence_interval[, 2],
    "Deterministic undercoverage glm coverage y2" = 
      glm2_deterministic$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < glm2_deterministic$confidence_interval[, 2],
    "Population mean y2" = mean(population$y2),
    "Stochastic undercoverage pmm yhat-y se_est y1"       = pmm1_stochastic$output$SE, 
    "Stochastic undercoverage pmm yhat-yhat se_est y1"    = pmm1.1_stochastic$output$SE,
    "Stochastic undercoverage pmm glm se_est y1"          = glm1_stochastic$output$SE,
    "Deterministic undercoverage pmm yhat-y se_est y1"    = pmm1_deterministic$output$SE, 
    "Deterministic undercoverage pmm yhat-yhat se_est y1" = pmm1.1_deterministic$output$SE,
    "Deterministic undercoverage glm se_est y1"           = glm1_deterministic$output$SE,
    "Stochastic undercoverage pmm yhat-y se_est y2"       = pmm2_stochastic$output$SE, 
    "Stochastic undercoverage pmm yhat-yhat se_est y2"    = pmm2.1_stochastic$output$SE,
    "Stochastic undercoverage pmm glm se_est y2"          = glm2_stochastic$output$SE,
    "Deterministic undercoverage pmm yhat-y se_est y2"    = pmm2_deterministic$output$SE, 
    "Deterministic undercoverage pmm yhat-yhat se_est y2" = pmm2.1_deterministic$output$SE,
    "Deterministic undercoverage glm se_est y2"           = glm2_deterministic$output$SE
  )
}

stopCluster(cl)

df <- data.frame(
  bias = c(
    mean(res[,  1] - res[, 13]),
    mean(res[,  2] - res[, 13]),
    mean(res[,  3] - res[, 13]),
    mean(res[,  4] - res[, 13]),
    mean(res[,  5] - res[, 13]),
    mean(res[,  6] - res[, 13]),
    mean(res[, 14] - res[, 26]),
    mean(res[, 15] - res[, 26]),
    mean(res[, 16] - res[, 26]),
    mean(res[, 17] - res[, 26]),
    mean(res[, 18] - res[, 26]),
    mean(res[, 19] - res[, 26])
  ),
  rel_bias_prec = c(
    100 * mean(res[,  1] - res[, 13]) / mean(res[, 13]),
    100 * mean(res[,  2] - res[, 13]) / mean(res[, 13]),
    100 * mean(res[,  3] - res[, 13]) / mean(res[, 13]),
    100 * mean(res[,  4] - res[, 13]) / mean(res[, 13]),
    100 * mean(res[,  5] - res[, 13]) / mean(res[, 13]),
    100 * mean(res[,  6] - res[, 13]) / mean(res[, 13]),
    100 * mean(res[, 14] - res[, 26]) / mean(res[, 26]),
    100 * mean(res[, 15] - res[, 26]) / mean(res[, 26]),
    100 * mean(res[, 16] - res[, 26]) / mean(res[, 26]),
    100 * mean(res[, 17] - res[, 26]) / mean(res[, 26]),
    100 * mean(res[, 18] - res[, 26]) / mean(res[, 26]),
    100 * mean(res[, 19] - res[, 26]) / mean(res[, 26])
  ),
  mse = c(
    mean((res[,  1] - res[, 13]) ^ 2),
    mean((res[,  2] - res[, 13]) ^ 2),
    mean((res[,  3] - res[, 13]) ^ 2),
    mean((res[,  4] - res[, 13]) ^ 2),
    mean((res[,  5] - res[, 13]) ^ 2),
    mean((res[,  6] - res[, 13]) ^ 2),
    mean((res[, 14] - res[, 26]) ^ 2),
    mean((res[, 15] - res[, 26]) ^ 2),
    mean((res[, 16] - res[, 26]) ^ 2),
    mean((res[, 17] - res[, 26]) ^ 2),
    mean((res[, 18] - res[, 26]) ^ 2),
    mean((res[, 19] - res[, 26]) ^ 2)
  ),
  mae = c(
    mean(abs(res[,  1] - res[, 13])),
    mean(abs(res[,  2] - res[, 13])),
    mean(abs(res[,  3] - res[, 13])),
    mean(abs(res[,  4] - res[, 13])),
    mean(abs(res[,  5] - res[, 13])),
    mean(abs(res[,  6] - res[, 13])),
    mean(abs(res[, 14] - res[, 26])),
    mean(abs(res[, 15] - res[, 26])),
    mean(abs(res[, 16] - res[, 26])),
    mean(abs(res[, 17] - res[, 26])),
    mean(abs(res[, 18] - res[, 26])),
    mean(abs(res[, 19] - res[, 26]))
  ),
  error_variance = c(
    var(res[,  1] - res[, 13]),
    var(res[,  2] - res[, 13]),
    var(res[,  3] - res[, 13]),
    var(res[,  4] - res[, 13]),
    var(res[,  5] - res[, 13]),
    var(res[,  6] - res[, 13]),
    var(res[, 14] - res[, 26]),
    var(res[, 15] - res[, 26]),
    var(res[, 16] - res[, 26]),
    var(res[, 17] - res[, 26]),
    var(res[, 18] - res[, 26]),
    var(res[, 19] - res[, 26])
  ),
  error_se = c(
    sd(res[,  1] - res[, 13]),
    sd(res[,  2] - res[, 13]),
    sd(res[,  3] - res[, 13]),
    sd(res[,  4] - res[, 13]),
    sd(res[,  5] - res[, 13]),
    sd(res[,  6] - res[, 13]),
    sd(res[, 14] - res[, 26]),
    sd(res[, 15] - res[, 26]),
    sd(res[, 16] - res[, 26]),
    sd(res[, 17] - res[, 26]),
    sd(res[, 18] - res[, 26]),
    sd(res[, 19] - res[, 26])
  ),
  mean_se_est = c(
    mean(res[, 27]),
    mean(res[, 28]),
    mean(res[, 29]),
    mean(res[, 30]),
    mean(res[, 31]),
    mean(res[, 32]),
    mean(res[, 33]),
    mean(res[, 34]),
    mean(res[, 35]),
    mean(res[, 36]),
    mean(res[, 37]),
    mean(res[, 38])
  ),
  coverage = c(
    mean(res[,  7]),
    mean(res[,  8]),
    mean(res[,  9]),
    mean(res[, 10]),
    mean(res[, 11]),
    mean(res[, 12]),
    mean(res[, 20]),
    mean(res[, 21]),
    mean(res[, 22]),
    mean(res[, 23]),
    mean(res[, 24]),
    mean(res[, 25])
  ),
  row.names = c(
    "Stochastic undercoverage sample 1 - yhat - y match",
    "Stochastic undercoverage sample 1 - yhat - yhat match",
    "Stochastic undercoverage sample 1 - glm",
    "Deterministic undercoverage sample 1 - yhat - y match",
    "Deterministic undercoverage sample 1 - yhat - yhat match",
    "Deterministic undercoverage sample 1 - glm",
    "Stochastic undercoverage sample 2 - yhat - y match",
    "Stochastic undercoverage sample 2 - yhat - yhat match",
    "Stochastic undercoverage sample 2 - yhat - glm",
    "Deterministic undercoverage sample 2 - yhat - y match",
    "Deterministic undercoverage sample 2 - yhat - yhat match",
    "Deterministic undercoverage sample 2 - yhat - glm"
  )
)

saveRDS(res, file = "results/custom-pmm-500-sims-check-positivity.rds")
res <- readRDS("results/custom-pmm-500-sims-check-positivity.rds")

df <- rbind(
  as.matrix(res[,c(c(1, 7, 27) +  0, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  1, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  2, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  3, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  4, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  5, 13)]),
  # ---------------------------------------------
  as.matrix(res[,c(c(14, 20, 33) +  0, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  1, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  2, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  3, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  4, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  5, 26)])
) |>
  as_tibble()

colnames(df) <- c("est", "coverage", "se", "true")

df[,"y_name"] <- paste0(rep(c("Stochastic undercoverage linear dependence", 
                              "Stochastic undercoverage non-linear dependence",
                              "Deterministic undercoverage linear dependence", 
                              "Deterministic undercoverage non-linear dependence"), 
                            each = NROW(df) / 4))
df[,"est_name"] <- rep(rep(c("yhat-y", "yhat-yhat", "glm"), each = 500), 4)

df <- df |>
  mutate(diff = est - true)

pp <- ggplot(data = df, aes(x = est_name, y = diff)) + 
  geom_violin(alpha = 0.8, draw_quantiles = 1:9 / 10, scale = "width") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point") +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  facet_wrap(~ y_name, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1)) +
  xlab("Estimator name") +
  ylab("Estimate error")

pp2 <- df |> 
  group_by(y_name, est_name) |>
  group_modify(.f = function(x, y) {
    xx <- binom.test(c(sum(x$coverage), sum(1 - x$coverage)), p = .95, n = sims)
    res <- data.frame(
      xx$conf.int[1],
      xx$conf.int[2],
      xx$estimate
    )
    colnames(res) <- c("lower", "upper", "mean")
    res
  }) |>
  #mutate(est_name = paste0(est_name, " - ", y_name)) |> 
  ggplot(aes(y = est_name, x = mean)) +
  geom_point(col = "blue", size = 5) +
  geom_errorbar(aes(xmin = lower, xmax = upper)) +
  geom_vline(aes(xintercept = .95), color = "red", linetype = "dashed") +
  facet_wrap(~ y_name, ncol = 2, scales = "free_x") +
  theme_bw() +
  xlab("Coverage") +
  ylab("Estimator and design")

ggsave("results/custom-pmm-500-sims-check-positivity-plot-errors.png", pp)
ggsave("results/custom-pmm-500-sims-check-positivity-plot-coverage.png", pp2)

res <- readRDS("results/custom-pmm-500-sims-check-positivity-really-large-sample.rds")

df <- rbind(
  as.matrix(res[,c(c(1, 7, 27) +  0, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  1, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  2, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  3, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  4, 13)]),
  as.matrix(res[,c(c(1, 7, 27) +  5, 13)]),
  # ---------------------------------------------
  as.matrix(res[,c(c(14, 20, 33) +  0, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  1, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  2, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  3, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  4, 26)]),
  as.matrix(res[,c(c(14, 20, 33) +  5, 26)])
) |>
  as_tibble()

colnames(df) <- c("est", "coverage", "se", "true")

df[,"y_name"] <- paste0(rep(c("Stochastic undercoverage linear dependence", 
                              "Stochastic undercoverage non-linear dependence",
                              "Deterministic undercoverage linear dependence", 
                              "Deterministic undercoverage non-linear dependence"), 
                            each = NROW(df) / 4))
df[,"est_name"] <- rep(rep(c("yhat-y", "yhat-yhat", "glm"), each = 500), 4)

df <- df |>
  mutate(diff = est - true)

pp <- ggplot(data = df, aes(x = est_name, y = diff)) + 
  geom_violin(alpha = 0.8, draw_quantiles = 1:9 / 10, scale = "width") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point") +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  facet_wrap(~ y_name, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1)) +
  xlab("Estimator name") +
  ylab("Estimate error")

pp2 <- df |> 
  group_by(y_name, est_name) |>
  group_modify(.f = function(x, y) {
    xx <- binom.test(c(sum(x$coverage), sum(1 - x$coverage)), p = .95, n = sims)
    res <- data.frame(
      xx$conf.int[1],
      xx$conf.int[2],
      xx$estimate
    )
    colnames(res) <- c("lower", "upper", "mean")
    res
  }) |>
  #mutate(est_name = paste0(est_name, " - ", y_name)) |> 
  ggplot(aes(y = est_name, x = mean)) +
  geom_point(col = "blue", size = 5) +
  geom_errorbar(aes(xmin = lower, xmax = upper)) +
  geom_vline(aes(xintercept = .95), color = "red", linetype = "dashed") +
  facet_wrap(~ y_name, ncol = 2, scales = "free_x") +
  theme_bw() +
  xlab("Coverage") +
  ylab("Estimator and design")

df |> 
  group_by(y_name, est_name) |>
  group_modify(.f = function(x, y) {
    xx <- binom.test(c(sum(x$coverage), sum(1 - x$coverage)), p = .95, n = sims)
    res <- data.frame(
      xx$conf.int[1],
      xx$conf.int[2],
      xx$estimate,
      mean(x$diff),
      mean((x$diff) ^ 2),
      mean(abs(x$diff)),
      sd(x$est)
    )
    res <- round(res, digits = 3)
    colnames(res) <- c("lower", "upper", "mean", "bias", "mse", "mae", "sd")
    res
  }) |>
  kableExtra::kable(format = "latex")

ggsave("results/custom-pmm-500-sims-check-positivity-really-large-sample-plot-errors.png", pp)
ggsave("results/custom-pmm-500-sims-check-positivity-really-large-sample-plot-coverage.png", pp2)
