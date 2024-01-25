library(nonprobsvy)
library(doSNOW)
library(progress)
library(tidyverse)

set.seed(stringr::str_split(lubridate::today(), "-") |> unlist() |> as.integer() |> sum())

cores <- 5

sims <- 5 * 100
N <- 1e5
n <- 200
KK <- 3

sigma <- diag(1, nrow = 5)
sigma[upper.tri(sigma)] <- runif(n = (5^2 - 5) / 2, max = .5, min = -.5)
sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
x1 <- MASS::mvrnorm(n = N / 5, mu = rep(1, 5), Sigma = sigma) |> as.vector()
x2 <- rexp(n = N, rate = 1)
sigma <- diag(2, nrow = 5)
sigma[upper.tri(sigma)] <- runif(n = (5^2 - 5) / 2, max = 1, min = -.7)
sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
epsilon <- MASS::mvrnorm(n = N / 5, mu = rep(0, 5), Sigma = sigma) |> as.vector()

p1 <- exp(x2 / 3 - x1 / 5 - 1.5) / (1 + exp(x2 / 3 - x1 / 5 - 1.5))
p2 <- exp(x1 / 3 - x2 / 5 - 1.5) / (1 + exp(x1 / 3 - x2 / 5 - 1.5))
population <- data.frame(
  x1,
  x2,
  y1 = 1 + x1 * .2 + x2 * .1 + epsilon, # linear
  y2 = -1.2 + (x1 - 0.5) ^ 2 + atan(x2) ^ (3 + sin(x1 + x2)) + sin(x1) * cos(x2) + epsilon, # weird
  y3 = x1 * x2 * epsilon, # multiplicative
  p1 = p1,
  base_w_srs = N/n
)

cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoSNOW(cl)

pb <- progress_bar$new(total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy", "dbscan", "gplm"),
               .options.snow = opts,
               .errorhandling = "stop") %dopar% {
  flag_srs <- rbinom(n = N, size = 1, prob = n / N)
  # flag_bd1 <- pmin(
  #   rbinom(n = 1:N, size = 1, prob = p1),
  #   epsilon > quantile(epsilon, .8) |
  #     quantile(epsilon, .2) > epsilon,
  #   rbinom(n = 1:N, size = 1, prob = p2)
  # )
  #flag_bd1 <- rbinom(n = 1:N, size = 1, prob = p2)
  # expected size ~~ 2k
  flag_bd1 <- pmin(rbinom(n = 1:N, size = 1, prob = p1 / 2),
                   rbinom(n = 1:N, size = 1, prob = p2))
  
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = subset(population, flag_srs == 1))
  
  glm1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
  
  pmm1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm1_loess <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2, 
                                 pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm1.1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm1.1_loess <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  glm2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
  
  pmm2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm2_loess <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm2.1 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm2.1_loess <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  glm3 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
  
  pmm3 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm3.1 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm3_loess <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  pmm3.1_loess <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1, pmm_reg_engine = "loess"),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  cbind(
    pmm1$output$mean, pmm1.1$output$mean,
    pmm1_loess$output$mean, pmm1.1_loess$output$mean,
    glm1$output$mean,
    pmm1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1$confidence_interval[, 2],
    pmm1.1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.1$confidence_interval[, 2],
    pmm1_loess$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1_loess$confidence_interval[, 2],
    pmm1.1_loess$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.1_loess$confidence_interval[, 2],
    glm1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < glm1$confidence_interval[, 2],
    mean(population$y1),
    pmm2$output$mean, pmm2.1$output$mean,
    pmm2_loess$output$mean, pmm2.1_loess$output$mean,
    glm2$output$mean,
    pmm2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2$confidence_interval[, 2],
    pmm2.1$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.1$confidence_interval[, 2],
    pmm2_loess$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2_loess$confidence_interval[, 2],
    pmm2.1_loess$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.1_loess$confidence_interval[, 2],
    glm2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < glm2$confidence_interval[, 2],
    mean(population$y2),
    pmm3$output$mean, pmm3.1$output$mean,
    pmm3_loess$output$mean, pmm3.1_loess$output$mean,
    glm3$output$mean,
    pmm3$confidence_interval[, 1] < mean(population$y3) &
      mean(population$y3) < pmm3$confidence_interval[, 2],
    pmm3.1$confidence_interval[, 1] < mean(population$y3) &
      mean(population$y3) < pmm3.1$confidence_interval[, 2],
    pmm3_loess$confidence_interval[, 1] < mean(population$y3) &
      mean(population$y3) < pmm3_loess$confidence_interval[, 2],
    pmm3.1_loess$confidence_interval[, 1] < mean(population$y3) &
      mean(population$y3) < pmm3.1_loess$confidence_interval[, 2],
    glm3$confidence_interval[, 1] < mean(population$y3) &
      mean(population$y3) < glm3$confidence_interval[, 2],
    mean(population$y3),
    pmm1$output$SE, pmm1.1$output$SE,
    pmm1_loess$output$SE, pmm1.1_loess$output$SE,
    glm1$output$SE, 
    pmm2$output$SE, pmm2.1$output$SE, 
    pmm2_loess$output$SE, pmm2.1_loess$output$SE, 
    glm2$output$SE, 
    pmm3$output$SE, pmm3.1$output$SE, 
    pmm3_loess$output$SE, pmm3.1_loess$output$SE, 
    glm3$output$SE
  )
}

stopCluster(cl)

df <- data.frame(
  bias = c(
    mean(res[,  1] - res[, 11]),
    mean(res[,  2] - res[, 11]),
    mean(res[,  3] - res[, 11]),
    mean(res[,  4] - res[, 11]),
    mean(res[,  5] - res[, 11]),
    mean(res[, 12] - res[, 22]),
    mean(res[, 13] - res[, 22]),
    mean(res[, 14] - res[, 22]),
    mean(res[, 15] - res[, 22]),
    mean(res[, 16] - res[, 22]),
    mean(res[, 23] - res[, 33]),
    mean(res[, 24] - res[, 33]),
    mean(res[, 25] - res[, 33]),
    mean(res[, 26] - res[, 33]),
    mean(res[, 27] - res[, 33])
  ),
  mse = c(
    mean((res[,  1] - res[, 11]) ^ 2),
    mean((res[,  2] - res[, 11]) ^ 2),
    mean((res[,  3] - res[, 11]) ^ 2),
    mean((res[,  4] - res[, 11]) ^ 2),
    mean((res[,  5] - res[, 11]) ^ 2),
    mean((res[, 12] - res[, 22]) ^ 2),
    mean((res[, 13] - res[, 22]) ^ 2),
    mean((res[, 14] - res[, 22]) ^ 2),
    mean((res[, 15] - res[, 22]) ^ 2),
    mean((res[, 16] - res[, 22]) ^ 2),
    mean((res[, 23] - res[, 33]) ^ 2),
    mean((res[, 24] - res[, 33]) ^ 2),
    mean((res[, 25] - res[, 33]) ^ 2),
    mean((res[, 26] - res[, 33]) ^ 2),
    mean((res[, 27] - res[, 33]) ^ 2)
  ),
  mae = c(
    mean(abs(res[,  1] - res[, 11])),
    mean(abs(res[,  2] - res[, 11])),
    mean(abs(res[,  3] - res[, 11])),
    mean(abs(res[,  4] - res[, 11])),
    mean(abs(res[,  5] - res[, 11])),
    mean(abs(res[, 12] - res[, 22])),
    mean(abs(res[, 13] - res[, 22])),
    mean(abs(res[, 14] - res[, 22])),
    mean(abs(res[, 15] - res[, 22])),
    mean(abs(res[, 16] - res[, 22])),
    mean(abs(res[, 23] - res[, 33])),
    mean(abs(res[, 24] - res[, 33])),
    mean(abs(res[, 25] - res[, 33])),
    mean(abs(res[, 26] - res[, 33])),
    mean(abs(res[, 27] - res[, 33]))
  ),
  error_variance = c(
    var(res[,  1] - res[, 11]),
    var(res[,  2] - res[, 11]),
    var(res[,  3] - res[, 11]),
    var(res[,  4] - res[, 11]),
    var(res[,  5] - res[, 11]),
    var(res[, 12] - res[, 22]),
    var(res[, 13] - res[, 22]),
    var(res[, 14] - res[, 22]),
    var(res[, 15] - res[, 22]),
    var(res[, 16] - res[, 22]),
    var(res[, 23] - res[, 33]),
    var(res[, 24] - res[, 33]),
    var(res[, 25] - res[, 33]),
    var(res[, 26] - res[, 33]),
    var(res[, 27] - res[, 33])
  ),
  error_se = c(
    sd(res[,  1] - res[, 11]),
    sd(res[,  2] - res[, 11]),
    sd(res[,  3] - res[, 11]),
    sd(res[,  4] - res[, 11]),
    sd(res[,  5] - res[, 11]),
    sd(res[, 12] - res[, 22]),
    sd(res[, 13] - res[, 22]),
    sd(res[, 14] - res[, 22]),
    sd(res[, 15] - res[, 22]),
    sd(res[, 16] - res[, 22]),
    sd(res[, 23] - res[, 33]),
    sd(res[, 24] - res[, 33]),
    sd(res[, 25] - res[, 33]),
    sd(res[, 26] - res[, 33]),
    sd(res[, 27] - res[, 33])
  ),
  mean_se_est = c(
    mean(res[, 34]),
    mean(res[, 35]),
    mean(res[, 36]),
    mean(res[, 38]),
    mean(res[, 38]),
    mean(res[, 39]),
    mean(res[, 40]),
    mean(res[, 41]),
    mean(res[, 42]),
    mean(res[, 43]),
    mean(res[, 44]),
    mean(res[, 45]),
    mean(res[, 46]),
    mean(res[, 47]),
    mean(res[, 48])
  ),
  coverage = c(
    mean(res[,  6]),
    mean(res[,  7]),
    mean(res[,  8]),
    mean(res[,  9]),
    mean(res[, 10]),
    mean(res[, 17]),
    mean(res[, 18]),
    mean(res[, 19]),
    mean(res[, 20]),
    mean(res[, 21]),
    mean(res[, 28]),
    mean(res[, 29]),
    mean(res[, 30]),
    mean(res[, 31]),
    mean(res[, 32])
  ),
  row.names = c(
    "y1 linear - yhat - yhat match",
    "y1 linear - yhat - y match",
    "y1 linear - yhat - yhat loess match",
    "y1 linear - yhat - y loess match",
    "y1 linear - glm",
    "y2 highly non linear - yhat - yhat match",
    "y2 highly non linear - yhat - y match",
    "y2 highly non linear - yhat - yhat loess match",
    "y2 highly non linear - yhat - y loess match",
    "y2 highly non linear - yhat - glm",
    "y3 multiplicative - yhat - yhat match",
    "y3 multiplicative - yhat - y match",
    "y3 multiplicative - yhat - yhat loess match",
    "y3 multiplicative - yhat - y loess match",
    "y3 multiplicative - yhat - glm"
  )
)

saveRDS(res, file = "results/custom-pmm-with-nonparametric-regression-500-sims.rds")


df <- rbind(
  as.matrix(res[,c(c(1, 6, 34) +  0, 11)]), # y1 linear - yhat - yhat match
  as.matrix(res[,c(c(1, 6, 34) +  1, 11)]), # y1 linear - yhat - y match
  as.matrix(res[,c(c(1, 6, 34) +  2, 11)]), # y1 linear - yhat - yhat match loess
  as.matrix(res[,c(c(1, 6, 34) +  3, 11)]), # y1 linear - yhat - y match loess
  as.matrix(res[,c(c(1, 6, 34) +  4, 11)]), # y1 linear - glm
  # ---------------------------------------------
  as.matrix(res[,c(c(12, 17, 39) +  0, 22)]), # y2 highly non linear - yhat - yhat match
  as.matrix(res[,c(c(12, 17, 39) +  1, 22)]), # y2 highly non linear - yhat - y match
  as.matrix(res[,c(c(12, 17, 39) +  2, 22)]), # y2 highly non linear - yhat - yhat match loess
  as.matrix(res[,c(c(12, 17, 39) +  3, 22)]), # y2 highly non linear - yhat - y match loess
  as.matrix(res[,c(c(12, 17, 39) +  4, 22)]), # y2 highly non linear - glm
  # ---------------------------------------------
  as.matrix(res[,c(c(23, 28, 44) +  0, 33)]), # y3 multiplicative - yhat - yhat match
  as.matrix(res[,c(c(23, 28, 44) +  1, 33)]), # y3 multiplicative - yhat - y match
  as.matrix(res[,c(c(23, 28, 44) +  2, 33)]), # y3 multiplicative - yhat - yhat match loess
  as.matrix(res[,c(c(23, 28, 44) +  3, 33)]), # y3 multiplicative - yhat - y match loess
  as.matrix(res[,c(c(23, 28, 44) +  4, 33)])  # y3 multiplicative - glm
) |>
  as_tibble()

colnames(df) <- c("est", "coverage", "se", "true")

df[,"y_name"] <- paste0(rep(c("linear", "highly non-linear", "multiplicative"), each = NROW(df) / 3))
df[,"est_name"] <- rep(rep(c("yhat-yhat", "yhat-y",
                             "yhat-yhat-loess", "yhat-y-loess",
                             "glm"), each = NROW(df) / 15), 3)

df <- df |>
  mutate(diff = est - true)

pp <- ggplot(data = df, aes(x = est_name, y = diff)) + 
  geom_violin(alpha = 0.8, draw_quantiles = 1:9 / 10, scale = "width") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point") +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  facet_wrap(~ y_name, ncol = 3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme_bw() +
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
  mutate(est_name = paste0(est_name, " - ", y_name)) |> 
  ggplot(aes(y = est_name, x = mean)) +
  geom_point(col = "blue", size = 5) +
  geom_errorbar(aes(xmin = lower, xmax = upper)) +
  geom_vline(aes(xintercept = .95), color = "red", linetype = "dashed") +
  theme_bw() +
  xlab("Coverage") +
  ylab("Estimator and design")

pp3 <- df |> 
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
  mutate(est_name = paste0(est_name, " - ", y_name)) |>
  filter(mean > .5) |> 
  ggplot(aes(y = est_name, x = mean)) +
  geom_point(col = "blue", size = 5) +
  geom_errorbar(aes(xmin = lower, xmax = upper)) +
  geom_vline(aes(xintercept = .95), color = "red", linetype = "dashed") +
  theme_bw() +
  xlab("Coverage") +
  ylab("Estimator and design")

ggsave("results/custom-pmm-with-nonparametric-regression-500-sims-plot-errors.png", pp)
ggsave("results/custom-pmm-with-nonparametric-regression-500-sims-plot-coverage.png", pp2)
ggsave("results/custom-pmm-with-nonparametric-regression-500-sims-plot-coverage-alt.png", pp3)
