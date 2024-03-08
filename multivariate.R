library(nonprobsvy)
library(tidyverse)
library(dbscan)
library(doSNOW)
library(progress)

set.seed(2051)

sims <- 500
cores <- 6
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

p1 <- exp(x2 - x1 - 2)/(1+exp(x2 - x1 - 2))
p2 <- exp(x1 * .6 - x2 - 2)/(1+exp(x1 * .6 - x2 - 2))
population <- data.frame(
  x1,
  x2,
  y1 = 1 + x1 * .2 + x2 * 5 + epsilon,
  y2 = -2 + 5 * (x1 - 0.5)^5 + x2 ^ 3 + epsilon,
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
               .packages = c("survey", "nonprobsvy", "dbscan"),
               .options.snow = opts,
               .errorhandling = "stop") %dopar% {
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
    
    pmm1 <- nonprob(
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
    
    pmm2 <- nonprob(
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
    

    ddf <- cbind(population[flag_bd1 == 1, , drop = FALSE],
                 "preds" = predict(mm1))
    
    sample_prob$variables <- cbind(sample_prob$variables, "preds" = predict(mm1, sample_prob$variables))
    
    pmm1_mv <- nonprob(
      outcome = y1 + y2 ~ preds.y1 + preds.y2 - 1,
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
    
    pmm2_mv <- nonprob(
      outcome = y1 + y2 ~ preds.y1 + preds.y2 - 1,
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
    
    cbind(
      pmm1$output$mean |> as.vector() |> t(),
      pmm2$output$mean |> as.vector() |> t(),
      glm1$output$mean |> as.vector() |> t(),
      pmm1_mv$output$mean |> as.vector() |> t(),
      pmm2_mv$output$mean |> as.vector() |> t()
    )
}
stopCluster(cl)

saveRDS(res, "results/custom-pmm-500-sims-robust.rds")
#res <- readRDS("~/Desktop/nonprobsvy-predictive-mean-matching/results/custom-pmm-500-sims-robust.rds")

df <- data.frame(
  bias = c(
    mean(res[,  1] - mean(population$y1)),
    mean(res[,  2] - mean(population$y2)),
    mean(res[,  3] - mean(population$y1)),
    mean(res[,  4] - mean(population$y2)),
    mean(res[,  5] - mean(population$y1)),
    mean(res[,  6] - mean(population$y2)),
    mean(res[,  7] - mean(population$y1)),
    mean(res[,  8] - mean(population$y2)),
    mean(res[,  9] - mean(population$y1)),
    mean(res[,  10] - mean(population$y2))
  ),
  rel_bias_prec = c(
    100 * mean(res[,  1] - mean(population$y1)) / mean(population$y1),
    100 * mean(res[,  2] - mean(population$y2)) / mean(population$y2),
    100 * mean(res[,  3] - mean(population$y1)) / mean(population$y1),
    100 * mean(res[,  4] - mean(population$y2)) / mean(population$y2),
    100 * mean(res[,  5] - mean(population$y1)) / mean(population$y1),
    100 * mean(res[,  6] - mean(population$y2)) / mean(population$y2),
    100 * mean(res[,  7] - mean(population$y1)) / mean(population$y1),
    100 * mean(res[,  8] - mean(population$y2)) / mean(population$y2),
    100 * mean(res[,  9] - mean(population$y1)) / mean(population$y1),
    100 * mean(res[, 10] - mean(population$y2)) / mean(population$y2)
  ),
  mse = c(
    mean((res[,  1] - mean(population$y1)) ^ 2),
    mean((res[,  2] - mean(population$y2)) ^ 2),
    mean((res[,  3] - mean(population$y1)) ^ 2),
    mean((res[,  4] - mean(population$y2)) ^ 2),
    mean((res[,  5] - mean(population$y1)) ^ 2),
    mean((res[,  6] - mean(population$y2)) ^ 2),
    mean((res[,  7] - mean(population$y1)) ^ 2),
    mean((res[,  8] - mean(population$y2)) ^ 2),
    mean((res[,  9] - mean(population$y1)) ^ 2),
    mean((res[, 10] - mean(population$y2)) ^ 2)
  ),
  mae = c(
    mean(abs(res[,  1] - mean(population$y1))),
    mean(abs(res[,  2] - mean(population$y2))),
    mean(abs(res[,  3] - mean(population$y1))),
    mean(abs(res[,  4] - mean(population$y2))),
    mean(abs(res[,  5] - mean(population$y1))),
    mean(abs(res[,  6] - mean(population$y2))),
    mean(abs(res[,  7] - mean(population$y1))),
    mean(abs(res[,  8] - mean(population$y2))),
    mean(abs(res[,  9] - mean(population$y1))),
    mean(abs(res[, 10] - mean(population$y2)))
  ),
  row.names = c(
    "y1 yhat - y",
    "y2 yhat - y",
    "y1 yhat - yhat",
    "y2 yhat - yhat",
    "y1 glm",
    "y2 glm",
    "y1 yhat - y - robust",
    "y2 yhat - y - robust",
    "y1 yhat - yhat - robust",
    "y2 yhat - yhat - robust"
  )
)

ddf <- c(
  res[,  1] - mean(population$y1),
  res[,  2] - mean(population$y2),
  res[,  3] - mean(population$y1),
  res[,  4] - mean(population$y2),
  res[,  5] - mean(population$y1),
  res[,  6] - mean(population$y2),
  res[,  7] - mean(population$y1),
  res[,  8] - mean(population$y2),
  res[,  9] - mean(population$y1),
  res[, 10] - mean(population$y2)
) |> as_tibble()

colnames(ddf) <- "diff"

ddf[,"y_name"] <- rep(rep(c("y1", "y2"), each = 500), 5)
ddf[,"est_name"] <- rep(rep(c("yhat-y", "yhat-yhat", "glm", 
                              "yhat-y-robust", "yhat-yhat-robust"), 
                            each = 1000))

pp <- ggplot(data = ddf, aes(x = est_name, y = diff)) + 
  geom_violin(alpha = 0.8, draw_quantiles = 1:9 / 10, scale = "width") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point") +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed") +
  facet_wrap(~ y_name, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1)) +
  xlab("Estimator name") +
  ylab("Estimate error")

kableExtra::kable(df, format = "latex")
ggsave("results/custom-pmm-500-sims-robust-plot-errors.png", pp)
