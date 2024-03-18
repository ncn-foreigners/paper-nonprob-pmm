# coverage 
library(nonprobsvy)
library(sampling)
library(doSNOW)
library(progress)
library(foreach)
library(tidyverse)

seed_for_sim <- 2051
set.seed(seed_for_sim)

KK <- 1

N <- 100000
n_a <- 500
sims <- 500
n_b <- 1000
n_b1 <- 0.7 * n_b
n_b2 <- 0.3 * n_b

cores <- 5
cl <- makeCluster(cores)

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy", "sampling"),
               .options.snow = opts) %dopar% {
  ## generate sample
  x <- rnorm(N, 2, 1)
  e <- rnorm(N)

  y1 <- 1 + 2*x + e
  y2 <- 3 + x + 2*e
  y3 <- 2.5 + 0.5*x^2 + e

  strata <- x <= 2
  pop <- data.frame(x, y1, y2, y3, strata)
  sample_a <- pop[sample(1:N, n_a),]
  sample_a$w_a <- N/n_a
  svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_a)
  pop1 <- subset(pop, strata == TRUE)
  pop2 <- subset(pop, strata == FALSE)
  sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                    pop2[sample(1:nrow(pop2), n_b2), ])
  
  ## Mass imputation with pmm, yhat-yhat
  
  est_mi_pmm_yhat_yhat_y1 <- nonprob(
    outcome = as.formula(y1 ~ x),
    data = sample_b,
    svydesign = svy_a,
    pop_size = N,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_yhat_yhat_y2 <- nonprob(
    outcome = as.formula(y2 ~ x),
    data = sample_b,
    svydesign = svy_a,
    pop_size = N,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_yhat_yhat_y3 <- nonprob(
    outcome = as.formula(y3 ~ x),
    data = sample_b,
    svydesign = svy_a,
    pop_size = N,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  ## Mass imputation with pmm, yhat-y
  
  est_mi_pmm_yhat_y_y1 <- nonprob(
    outcome = as.formula(y1 ~ x),
    data = sample_b,
    svydesign = svy_a,
    pop_size = N,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_yhat_y_y2 <- nonprob(
    outcome = as.formula(y2 ~ x),
    data = sample_b,
    svydesign = svy_a,
    pop_size = N,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_yhat_y_y3 <- nonprob(
    outcome = as.formula(y3 ~ x),
    data = sample_b,
    svydesign = svy_a,
    pop_size = N,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  ## Mass imputation with glm
  
  est_mi_glm_y1 <- nonprob(
    outcome = as.formula(y1 ~ x),
    data = sample_b,
    pop_size = N,
    svydesign = svy_a,
    family_outcome = "gaussian"
  )
  
  est_mi_glm_y2 <- nonprob(
    outcome = as.formula(y2 ~ x),
    data = sample_b,
    pop_size = N,
    svydesign = svy_a,
    family_outcome = "gaussian"
  )
  
  est_mi_glm_y3 <- nonprob(
    outcome = as.formula(y3 ~ x),
    data = sample_b,
    pop_size = N,
    svydesign = svy_a,
    family_outcome = "gaussian"
  )
  
  cbind(
    # Y1
    est_mi_pmm_yhat_yhat_y1$output,
    est_mi_pmm_yhat_yhat_y1$confidence_interval,
    est_mi_pmm_yhat_y_y1$output,
    est_mi_pmm_yhat_y_y1$confidence_interval,
    est_mi_glm_y1$output,
    est_mi_glm_y1$confidence_interval,
    mean(pop$y1),
    # Y2
    est_mi_pmm_yhat_yhat_y2$output,
    est_mi_pmm_yhat_yhat_y2$confidence_interval,
    est_mi_pmm_yhat_y_y2$output,
    est_mi_pmm_yhat_y_y2$confidence_interval,
    est_mi_glm_y2$output,
    est_mi_glm_y2$confidence_interval,
    mean(pop$y2),
    # Y3
    est_mi_pmm_yhat_yhat_y3$output,
    est_mi_pmm_yhat_yhat_y3$confidence_interval,
    est_mi_pmm_yhat_y_y3$output,
    est_mi_pmm_yhat_y_y3$confidence_interval,
    est_mi_glm_y3$output,
    est_mi_glm_y3$confidence_interval,
    mean(pop$y3)
  )
}

stopCluster(cl)
res <- as.data.frame(res)

colnames(res) <- c(
  "Y_1-yhat-yhat-match-est",
  "Y_1-yhat-yhat-match-se",
  "Y_1-yhat-yhat-match-lower",
  "Y_1-yhat-yhat-match-upper",
  "Y_1-yhat-y-match-est",
  "Y_1-yhat-y-match-se",
  "Y_1-yhat-y-match-lower",
  "Y_1-yhat-y-match-upper",
  "Y_1-glm-est",
  "Y_1-glm-se",
  "Y_1-glm-lower",
  "Y_1-glm-upper",
  "Y_1_true_mean",
  "Y_2-yhat-yhat-match-est",
  "Y_2-yhat-yhat-match-se",
  "Y_2-yhat-yhat-match-lower",
  "Y_2-yhat-yhat-match-upper",
  "Y_2-yhat-y-match-est",
  "Y_2-yhat-y-match-se",
  "Y_2-yhat-y-match-lower",
  "Y_2-yhat-y-match-upper",
  "Y_2-glm-est",
  "Y_2-glm-se",
  "Y_2-glm-lower",
  "Y_2-glm-upper",
  "Y_2_true_mean",
  "Y_3-yhat-yhat-match-est",
  "Y_3-yhat-yhat-match-se",
  "Y_3-yhat-yhat-match-lower",
  "Y_3-yhat-yhat-match-upper",
  "Y_3-yhat-y-match-est",
  "Y_3-yhat-y-match-se",
  "Y_3-yhat-y-match-lower",
  "Y_3-yhat-y-match-upper",
  "Y_3-glm-est",
  "Y_3-glm-se",
  "Y_3-glm-lower",
  "Y_3-glm-upper",
  "Y_3_true_mean"
)

saveRDS(res, file = "results/kim2021-pmm-500-sims.rds")
res <- readRDS("results/kim2021-pmm-500-sims.rds")

df <- rbind(
  as.matrix(res[,c(1:4, 13)]), as.matrix(res[,c(1:4+4, 13)]), as.matrix(res[,c(1:4+8, 13)]), 
  as.matrix(res[,c(1:4+13, 26)]), as.matrix(res[,c(1:4+17, 26)]), as.matrix(res[,c(1:4+21, 26)]), 
  as.matrix(res[,c(1:4+26, 39)]), as.matrix(res[,c(1:4+30, 39)]), as.matrix(res[,c(1:4+34, 39)])
) |>
  as_tibble()

colnames(df) <- c("est", "se", "lower", "upper", "true")

df[,"y_name"] <- paste0("Y_", rep(1:3, each = 1500))
df[,"est_name"] <- rep(rep(c("yhat-yhat", "yhat-y", "glm"), each = 500), 3)

df <- df |>
  mutate(diff = est - true, covr = lower < true & true < upper)

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
    xx <- binom.test(c(sum(x$covr), sum(1 - x$covr)), p = .95, n = 500)
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

ggsave("results/kim2021-pmm-500-sims-plot-errors.png", pp)
ggsave("results/kim2021-pmm-500-sims-plot-coverage.png", pp2)
