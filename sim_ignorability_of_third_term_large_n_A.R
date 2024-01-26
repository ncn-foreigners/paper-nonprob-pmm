# instalacja z tego branchu
#remotes::install_github("https://github.com/ncn-foreigners/nonprobsvy/tree/dev")
library(nonprobsvy)
library(doSNOW)
library(progress)
library(tidyverse)

set.seed(stringr::str_split(lubridate::today(), "-") |> unlist() |> as.integer() |> sum())

cores <- 5

sims <- 5 * 100
N <- 1e6
n <- 500
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

p1 <- exp(x2)/(1+exp(x2))
p2 <- exp(x1)/(1+exp(x1))
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

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy"),
               .options.snow = opts) %dopar% {
  # flag_bd1 <- c(rbinom(n = 1:N, size = 1, prob = p1))
  #flag_bd1 <- c(rbinom(n = 1:(N - 10000), size = 1, prob = p1), rep(0, 10000))
  # flag_bd1 <- pmin(rbinom(n = 1:N, size = 1, prob = p1),
  #                  1 - rbinom(n = 1:N, size = 1, prob = p2))
  flag_bd1 <- pmax( # planned size ~~ 9K
    rbinom(n = 1:N, size = 1, prob = p1),
    rbinom(n = 1:N, size = 1, prob = p2)
  )
  base_w_bd <- N/n
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = population[sample(1:N, n),])

  glm1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
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
  
  glm2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
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
    glm1$output$mean,
    pmm1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1$confidence_interval[, 2],
    pmm1.1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.1$confidence_interval[, 2],
    pmm1.2$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.2$confidence_interval[, 2],
    pmm1.3$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.3$confidence_interval[, 2],
    glm1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < glm1$confidence_interval[, 2],
    mean(population$y1),
    pmm2$output$mean, pmm2.1$output$mean,
    pmm2.2$output$mean, pmm2.3$output$mean,
    glm2$output$mean,
    pmm2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2$confidence_interval[, 2],
    pmm2.1$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.1$confidence_interval[, 2],
    pmm2.2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.2$confidence_interval[, 2],
    pmm2.3$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.3$confidence_interval[, 2],
    glm2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < glm2$confidence_interval[, 2],
    mean(population$y2),
    pmm1$output$SE, pmm1.1$output$SE,
    pmm1.2$output$SE, pmm1.3$output$SE,
    glm1$output$SE, pmm2$output$SE, 
    pmm2.1$output$SE, pmm2.2$output$SE, 
    pmm2.3$output$SE, glm2$output$SE
  )
}

stopCluster(cl)
colnames(res) <- c(
  "sample 1 - yhat - y match",#1
  "sample 1 - yhat - yhat match",#2
  "sample 1 - yhat - y match - third term included",#3
  "sample 1 - yhat - yhat match - third term included",#4
  "sample 1 - glm",#5
  "sample 1 - yhat - y match - coverage",#6
  "sample 1 - yhat - yhat match - coverage",#7
  "sample 1 - yhat - y match - coverage - third term included",#8
  "sample 1 - yhat - yhat match - coverage - third term included",#9
  "sample 1 - glm - coverage",#10
  "sample 1 - mean of y",#11
  "sample 2 - yhat - y match",#12
  "sample 2 - yhat - yhat match",#13
  "sample 2 - yhat - y match - third term included",#14
  "sample 2 - yhat - yhat match - third term included",#15
  "sample 2 - glm",#16
  "sample 2 - yhat - y match - coverage",#17
  "sample 2 - yhat - yhat match - coverage",#18
  "sample 2 - yhat - y match - coverage - third term included",#19
  "sample 2 - yhat - yhat match - coverage - third term included",#20
  "sample 2 - glm - coverage",#21
  "sample 2 - mean of y",#22
  "sample 1 - est se yhat - y",#23
  "sample 1 - est se yhat - yhat",#24
  "sample 1 - est se yhat - y - third term included",#25
  "sample 1 - est se yhat - yhat - third term included",#26
  "sample 1 - est se glm",#27
  "sample 2 - est se yhat - y",#28
  "sample 2 - est se yhat - yhat",#29
  "sample 2 - est se yhat - y - third term included",#30
  "sample 2 - est se yhat - yhat - third term included",#31
  "sample 2 - est se glm"#32
)

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
    mean(res[, 16] - res[, 22])
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
    mean((res[, 16] - res[, 22]) ^ 2)
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
    mean(abs(res[, 16] - res[, 22]))
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
    var(res[, 16] - res[, 22])
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
    sd(res[, 16] - res[, 22])
  ),
  mean_se_est = c(
    mean(res[, 23]),
    mean(res[, 24]),
    mean(res[, 25]),
    mean(res[, 26]),
    mean(res[, 27]),
    mean(res[, 28]),
    mean(res[, 29]),
    mean(res[, 30]),
    mean(res[, 31]),
    mean(res[, 32])
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
    mean(res[, 21])
  ),
  row.names = c(
    "sample 1 - yhat - y match",
    "sample 1 - yhat - yhat match",
    "sample 1 - yhat - y match - third term included",
    "sample 1 - yhat - yhat match - third term included",
    "sample 1 - glm",
    "sample 2 - yhat - y match",
    "sample 2 - yhat - yhat match",
    "sample 2 - yhat - y match - third term included",
    "sample 2 - yhat - yhat match - third term included",
    "sample 2 - yhat - glm"
  )
)

saveRDS(res, file = "results/custom-pmm-500-sims-check-variance-large-nonprob.rds")
res <- readRDS("results/custom-pmm-500-sims-check-variance-large-nonprob.rds")

df <- rbind(
  as.matrix(res[,c(c(1, 6, 23) +  0, 11)]), # y1 linear - yhat - yhat match
  as.matrix(res[,c(c(1, 6, 23) +  1, 11)]), # y1 linear - yhat - y match
  as.matrix(res[,c(c(1, 6, 23) +  2, 11)]), # y1 linear - yhat - yhat match third term included
  as.matrix(res[,c(c(1, 6, 23) +  3, 11)]), # y1 linear - yhat - y match third term included
  as.matrix(res[,c(c(1, 6, 23) +  4, 11)]), # y1 linear - glm
  # ---------------------------------------------
  as.matrix(res[,c(c(12, 17, 28) +  0, 22)]), # y2 highly non linear - yhat - yhat match
  as.matrix(res[,c(c(12, 17, 28) +  1, 22)]), # y2 highly non linear - yhat - y match
  as.matrix(res[,c(c(12, 17, 28) +  2, 22)]), # y2 highly non linear - yhat - yhat match loess
  as.matrix(res[,c(c(12, 17, 28) +  3, 22)]), # y2 highly non linear - yhat - y match loess
  as.matrix(res[,c(c(12, 17, 28) +  4, 22)])  # y2 highly non linear - glm
) |>
  as_tibble()

colnames(df) <- c("est", "coverage", "se", "true")

df[,"y_name"] <- paste0(rep(c("linear", "non-linear"), each = NROW(df) / 2))
df[,"est_name"] <- rep(rep(c("yhat-y", "yhat-yhat",
                             "yhat-y-third term included", "yhat-yhat-third term included",
                             "glm"), each = NROW(df) / 10), 2)

df <- df |>
  mutate(diff = est - true)

pp <- ggplot(data = df, aes(x = est_name, y = diff)) + 
  geom_violin(alpha = 0.8, draw_quantiles = 1:9 / 10, scale = "width") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point") +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  facet_wrap(~ y_name, ncol = 3, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
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

ggsave("results/custom-pmm-500-sims-check-variance-large-nonprob-plot-errors.png", pp)
ggsave("results/custom-pmm-500-sims-check-variance-large-nonprob-plot-coverage.png", pp2)
