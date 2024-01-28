# instalacja z tego branchu
#remotes::install_github("https://github.com/ncn-foreigners/nonprobsvy/tree/dev")
library(nonprobsvy)
library(doSNOW)
library(progress)
library(tidyverse)

set.seed(2051)

cores <- 5

sims <- 5 * 100
N <- 1e6
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

p1 <- exp(x2 - x1 - 2)/(1 + exp(x2 - x1 - 2))
p2 <- exp(x1 * .6 - x2 - 2) / (1 + exp(x1 * .6 - x2 - 2))
population <- data.frame(
  x1,
  x2,
  y1 = 6 * x1 - 5 * x2 - 7 + 10 * epsilon,
  y2 = -2 + .37 * (x1 - 0.5) ^ 2 + x2 + epsilon,
  p1 = p1,
  base_w_srs = N/n
)

xx <- rbinom(prob = plogis(.5 * cos(x2) - .5 * sin(x1)), 
             size = 1, n = N) |> 
  as.logical()

# expected size about 500_000
population1 <- population[xx,]

cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy"),
               .options.snow = opts) %dopar% {
  flag_bd1 <- pmin( # planned size ~~ 4500
    rbinom(n = 1:NROW(population1), size = 1, prob = p1[xx]),
    rbinom(n = 1:NROW(population1), size = 1, prob = p2[xx])
  )
  base_w_bd <- N/n
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = population[sample(1:N, n),])

  glm1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  pmm1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmm1.1 <- nonprob(
    outcome = y1 ~ x1 + x2,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  glm2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )

  pmm2 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population1[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )

  pmm2.1 <- nonprob(
    outcome = y2 ~ x1 + x2,
    data = population1[flag_bd1 == 1, , drop = FALSE],
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
    glm1$output$mean,
    pmm1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1$confidence_interval[, 2],
    pmm1.1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < pmm1.1$confidence_interval[, 2],
    glm1$confidence_interval[, 1] < mean(population$y1) &
      mean(population$y1) < glm1$confidence_interval[, 2],
    mean(population$y1),
    pmm2$output$mean, pmm2.1$output$mean,
    glm2$output$mean,
    pmm2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2$confidence_interval[, 2],
    pmm2.1$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < pmm2.1$confidence_interval[, 2],
    glm2$confidence_interval[, 1] < mean(population$y2) &
      mean(population$y2) < glm2$confidence_interval[, 2],
    mean(population$y2),
    pmm1$output$SE, pmm1.1$output$SE,
    glm1$output$SE, pmm2$output$SE, 
    pmm2.1$output$SE, glm2$output$SE
  )
}

stopCluster(cl)
colnames(res) <- c(
  "sample 1 - yhat - y match",#1
  "sample 1 - yhat - yhat match",#2
  "sample 1 - glm",#3
  "sample 1 - yhat - y match - coverage",#4
  "sample 1 - yhat - yhat match - coverage",#5
  "sample 1 - glm - coverage",#6
  "sample 1 - mean of y",#7
  "sample 2 - yhat - y match",#8
  "sample 2 - yhat - yhat match",#9
  "sample 2 - glm",#10
  "sample 2 - yhat - y match - coverage",#11
  "sample 2 - yhat - yhat match - coverage",#12
  "sample 2 - glm - coverage",#13
  "sample 2 - mean of y",#14
  "sample 1 - est se yhat - y",#15
  "sample 1 - est se yhat - yhat",#16
  "sample 1 - est se glm",#17
  "sample 2 - est se yhat - y",#18
  "sample 2 - est se yhat - yhat",#19
  "sample 2 - est se glm"#20
)

df <- data.frame(
  bias = c(
    mean(res[,  1] - res[,  7]),
    mean(res[,  2] - res[,  7]),
    mean(res[,  3] - res[,  7]),
    mean(res[,  8] - res[, 14]),
    mean(res[,  9] - res[, 14]),
    mean(res[, 10] - res[, 14])
  ),
  mse = c(
    mean((res[,  1] - res[,  7]) ^ 2),
    mean((res[,  2] - res[,  7]) ^ 2),
    mean((res[,  3] - res[,  7]) ^ 2),
    mean((res[,  8] - res[, 14]) ^ 2),
    mean((res[,  9] - res[, 14]) ^ 2),
    mean((res[, 10] - res[, 14]) ^ 2)
  ),
  mae = c(
    mean(abs(res[,  1] - res[,  7])),
    mean(abs(res[,  2] - res[,  7])),
    mean(abs(res[,  3] - res[,  7])),
    mean(abs(res[,  8] - res[, 14])),
    mean(abs(res[,  9] - res[, 14])),
    mean(abs(res[, 10] - res[, 14]))
  ),
  error_variance = c(
    var(res[,  1] - res[,  7]),
    var(res[,  2] - res[,  7]),
    var(res[,  3] - res[,  7]),
    var(res[,  8] - res[, 14]),
    var(res[,  9] - res[, 14]),
    var(res[, 10] - res[, 14])
  ),
  error_se = c(
    sd(res[,  1] - res[,  7]),
    sd(res[,  2] - res[,  7]),
    sd(res[,  3] - res[,  7]),
    sd(res[,  8] - res[, 14]),
    sd(res[,  9] - res[, 14]),
    sd(res[, 10] - res[, 14])
  ),
  mean_se_est = c(
    mean(res[, 15]),
    mean(res[, 16]),
    mean(res[, 17]),
    mean(res[, 18]),
    mean(res[, 19]),
    mean(res[, 20])
  ),
  coverage = c(
    mean(res[,  4]),
    mean(res[,  5]),
    mean(res[,  6]),
    mean(res[, 11]),
    mean(res[, 12]),
    mean(res[, 13])
  ),
  row.names = c(
    "sample 1 - yhat - y match",
    "sample 1 - yhat - yhat match",
    "sample 1 - glm",
    "sample 2 - yhat - y match",
    "sample 2 - yhat - yhat match",
    "sample 2 - yhat - glm"
  )
)

saveRDS(res, file = "results/custom-pmm-500-sims-check-positivity.rds")
res <- readRDS("results/custom-pmm-500-sims-check-positivity.rds")

df <- rbind(
  as.matrix(res[,c(c(1, 4, 15) +  0, 7)]), # y1 linear - yhat - yhat match
  as.matrix(res[,c(c(1, 4, 15) +  1, 7)]), # y1 linear - yhat - y match
  as.matrix(res[,c(c(1, 4, 15) +  2, 7)]), # y1 linear - glm
  # ---------------------------------------------
  as.matrix(res[,c(c(8, 11, 18) +  0, 14)]), # y2 highly non linear - yhat - yhat match
  as.matrix(res[,c(c(8, 11, 18) +  1, 14)]), # y2 highly non linear - yhat - y match
  as.matrix(res[,c(c(8, 11, 18) +  2, 14)])  # y2 highly non linear - glm
) |>
  as_tibble()

colnames(df) <- c("est", "coverage", "se", "true")

df[,"y_name"] <- paste0(rep(c("linear", "non-linear"), each = NROW(df) / 2))
df[,"est_name"] <- rep(rep(c("yhat-y", "yhat-yhat", "glm"), each = 500), 2)

df <- df |>
  mutate(diff = est - true)

pp <- ggplot(data = df, aes(x = est_name, y = diff)) + 
  geom_violin(alpha = 0.8, draw_quantiles = 1:9 / 10, scale = "width") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point") +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  facet_wrap(~ y_name, ncol = 3, scales = "free_y") +
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
  mutate(est_name = paste0(est_name, " - ", y_name)) |> 
  ggplot(aes(y = est_name, x = mean)) +
  geom_point(col = "blue", size = 5) +
  geom_errorbar(aes(xmin = lower, xmax = upper)) +
  geom_vline(aes(xintercept = .95), color = "red", linetype = "dashed") +
  theme_bw() +
  xlab("Coverage") +
  ylab("Estimator and design")

ggsave("results/custom-pmm-500-sims-check-positivity-plot-errors.png", pp)
ggsave("results/custom-pmm-500-sims-check-positivity-plot-coverage.png", pp2)
