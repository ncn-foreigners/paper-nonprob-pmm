# coverage 
library(nonprobsvy)
library(sampling)
library(doSNOW)
library(progress)
library(foreach)
library(tidyverse)

seed_for_sim <- str_split(today(), "-") |> unlist() |> as.integer() |> sum()
set.seed(seed_for_sim)

N <- 10000
n_A <- 500
sims <- 500
p <- 50
KK <- 1
alpha_vec1 <- c(-2, 1, 1, 1, 1, rep(0, p - 5))
alpha_vec2 <- c(0, 0, 0, 3, 3, 3, 3, rep(0, p - 7))
beta_vec <- c(1, 0, 0, 1, 1, 1, 1, rep(0, p - 7))

## generate X
X <- cbind(
  "(Intercept)" = 1, 
  matrix(
    rnorm(N * (p - 1)), 
    nrow = N, byrow = TRUE, 
    dimnames = list(NULL, paste0("X", 1:(p - 1)))
  )
)
X_formula  <- as.formula(paste("~", paste0("X", 1:(p - 1), collapse = " + ")))
#X_totals <- colSums(X)
#X_means <- colMeans(X[,-1])

## generate Y
Y_11 <- as.numeric(X %*% beta_vec) +  rnorm(N) ## OM I: linear model
Y_12 <- 1 + exp(3*sin(as.numeric(X %*% beta_vec))) + X[, "X5"] + X[, "X6"] + rnorm(N) ## OM II: nonlinear model
pi_Y_21 <- plogis(as.numeric(X %*% beta_vec)) ## OM III: linear model for binary Y
pi_Y_22 <- plogis(2 - log((X %*% beta_vec)^2) - 2*X[,"X5"] + 2*X[, "X6"]) ## OM IV: nonlinear model for binary Y
Y_21 <- rbinom(N, 1, prob = pi_Y_21)
Y_22 <- rbinom(N, 1, prob = pi_Y_22)

## generate probs
pi_A <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_11), n_A) ## inclusion based on Y_11 only 
pi_B1 <- plogis(as.numeric(X %*% alpha_vec1)) ## PSM I: linear probability

## generate data
pop_data <- data.frame(pi_A, Y_11, Y_12, Y_21, Y_22, X[, 2:p])
head(pop_data)

cores <- 5
cl <- makeCluster(cores)
clusterExport(cl, c("N", "n"))

registerDoSNOW(cl)

pb <- progress_bar$new(total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy", "sampling"),
               .options.snow = opts) %dopar% {
  ## generate sample
  flag_B1 <- rbinom(N, 1, prob = pi_B1)
  flag_A <- UPpoisson(pik = pi_A)
  
  sample_A_svy <- svydesign(ids = ~ 1, probs = ~ pi_A, pps = "brewer", data = pop_data[flag_A == 1, ])
  
  sample_B1 <- pop_data[flag_B1 == 1, ]
  
  ## Mass imputation with pmm, yhat-yhat, no var select
  
  est_mi_pmm_y11_yhat_yhat <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_y12_yhat_yhat <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_y21_yhat_yhat <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_y22_yhat_yhat <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  ## Mass imputation with pmm, yhat-yhat, var select
  
  est_mi_pmm_y11_yhat_yhat_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmm_y12_yhat_yhat_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmm_y21_yhat_yhat_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmm_y22_yhat_yhat_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  
  ## Mass imputation with pmm, yhat-y, no var select
  
  est_mi_pmm_y11_yhat_y <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_y12_yhat_y <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_y21_yhat_y <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  est_mi_pmm_y22_yhat_y <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  ## Mass imputation with pmm, yhat-y, var select
  
  est_mi_pmm_y11_yhat_y_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmm_y12_yhat_y_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmm_y21_yhat_y_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmm_y22_yhat_y_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  ## Mass imputation with glm, no var select
  
  est_mi_glm_y11 <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    family_outcome = "gaussian"
  )
  
  est_mi_glm_y12 <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    family_outcome = "gaussian"
  )
  
  est_mi_glm_y21 <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    family_outcome = "binomial"
  )
  
  est_mi_glm_y22 <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    family_outcome = "binomial"
  )
  
  ## Mass imputation with glm, var select
  
  est_mi_glm_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    family_outcome = "gaussian",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_glm_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    family_outcome = "gaussian",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_glm_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    family_outcome = "binomial",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_glm_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy,
    family_outcome = "binomial",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  cbind(
    # Y11
    est_mi_pmm_y11_yhat_yhat$output,
    est_mi_pmm_y11_yhat_yhat$confidence_interval,
    est_mi_pmm_y11_yhat_yhat_sel$output,
    est_mi_pmm_y11_yhat_yhat_sel$confidence_interval,
    est_mi_pmm_y11_yhat_y$output,
    est_mi_pmm_y11_yhat_y$confidence_interval,
    est_mi_pmm_y11_yhat_y_sel$output,
    est_mi_pmm_y11_yhat_y_sel$confidence_interval,
    est_mi_glm_y11$output,
    est_mi_glm_y11$confidence_interval,
    est_mi_glm_y11_sel$output,
    est_mi_glm_y11_sel$confidence_interval,
    mean(pop_data$Y_11),
    # Y12
    est_mi_pmm_y12_yhat_yhat$output,
    est_mi_pmm_y12_yhat_yhat$confidence_interval,
    est_mi_pmm_y12_yhat_yhat_sel$output,
    est_mi_pmm_y12_yhat_yhat_sel$confidence_interval,
    est_mi_pmm_y12_yhat_y$output,
    est_mi_pmm_y12_yhat_y$confidence_interval,
    est_mi_pmm_y12_yhat_y_sel$output,
    est_mi_pmm_y12_yhat_y_sel$confidence_interval,
    est_mi_glm_y12$output,
    est_mi_glm_y12$confidence_interval,
    est_mi_glm_y12_sel$output,
    est_mi_glm_y12_sel$confidence_interval,
    mean(pop_data$Y_12),
    # Y21
    est_mi_pmm_y21_yhat_yhat$output,
    est_mi_pmm_y21_yhat_yhat$confidence_interval,
    est_mi_pmm_y21_yhat_yhat_sel$output,
    est_mi_pmm_y21_yhat_yhat_sel$confidence_interval,
    est_mi_pmm_y21_yhat_y$output,
    est_mi_pmm_y21_yhat_y$confidence_interval,
    est_mi_pmm_y21_yhat_y_sel$output,
    est_mi_pmm_y21_yhat_y_sel$confidence_interval,
    est_mi_glm_y21$output,
    est_mi_glm_y21$confidence_interval,
    est_mi_glm_y21_sel$output,
    est_mi_glm_y21_sel$confidence_interval,
    mean(pop_data$Y_21),
    # Y22
    est_mi_pmm_y22_yhat_yhat$output,
    est_mi_pmm_y22_yhat_yhat$confidence_interval,
    est_mi_pmm_y22_yhat_yhat_sel$output,
    est_mi_pmm_y22_yhat_yhat_sel$confidence_interval,
    est_mi_pmm_y22_yhat_y$output,
    est_mi_pmm_y22_yhat_y$confidence_interval,
    est_mi_pmm_y22_yhat_y_sel$output,
    est_mi_pmm_y22_yhat_y_sel$confidence_interval,
    est_mi_glm_y22$output,
    est_mi_glm_y22$confidence_interval,
    est_mi_glm_y22_sel$output,
    est_mi_glm_y22_sel$confidence_interval,
    mean(pop_data$Y_22)
  )
}

stopCluster(cl)
res <- as.data.frame(res)

colnames(res) <- c(
  "Y_11-yhat-yhat-match-no-var-sel-est",
  "Y_11-yhat-yhat-match-no-var-sel-se",
  "Y_11-yhat-yhat-match-no-var-sel-lower",
  "Y_11-yhat-yhat-match-no-var-sel-upper",
  "Y_11-yhat-yhat-match-var-sel-est",
  "Y_11-yhat-yhat-match-var-sel-se",
  "Y_11-yhat-yhat-match-var-sel-lower",
  "Y_11-yhat-yhat-match-var-sel-upper",
  "Y_11-yhat-y-match-no-var-sel-est",
  "Y_11-yhat-y-match-no-var-sel-se",
  "Y_11-yhat-y-match-no-var-sel-lower",
  "Y_11-yhat-y-match-no-var-sel-upper",
  "Y_11-yhat-y-match-var-sel-est",
  "Y_11-yhat-y-match-var-sel-se",
  "Y_11-yhat-y-match-var-sel-lower",
  "Y_11-yhat-y-match-var-sel-upper",
  "Y_11-mi-no-var-sel-est",
  "Y_11-mi-no-var-sel-se",
  "Y_11-mi-no-var-sel-lower",
  "Y_11-mi-no-var-sel-upper",
  "Y_11-mi-var-sel-est",
  "Y_11-mi-var-sel-se",
  "Y_11-mi-var-sel-lower",
  "Y_11-mi-var-sel-upper",
  "Y_11_true_mean",
  "Y_12-yhat-yhat-match-no-var-sel-est",
  "Y_12-yhat-yhat-match-no-var-sel-se",
  "Y_12-yhat-yhat-match-no-var-sel-lower",
  "Y_12-yhat-yhat-match-no-var-sel-upper",
  "Y_12-yhat-yhat-match-var-sel-est",
  "Y_12-yhat-yhat-match-var-sel-se",
  "Y_12-yhat-yhat-match-var-sel-lower",
  "Y_12-yhat-yhat-match-var-sel-upper",
  "Y_12-yhat-y-match-no-var-sel-est",
  "Y_12-yhat-y-match-no-var-sel-se",
  "Y_12-yhat-y-match-no-var-sel-lower",
  "Y_12-yhat-y-match-no-var-sel-upper",
  "Y_12-yhat-y-match-var-sel-est",
  "Y_12-yhat-y-match-var-sel-se",
  "Y_12-yhat-y-match-var-sel-lower",
  "Y_12-yhat-y-match-var-sel-upper",
  "Y_12-mi-no-var-sel-est",
  "Y_12-mi-no-var-sel-se",
  "Y_12-mi-no-var-sel-lower",
  "Y_12-mi-no-var-sel-upper",
  "Y_12-mi-var-sel-est",
  "Y_12-mi-var-sel-se",
  "Y_12-mi-var-sel-lower",
  "Y_12-mi-var-sel-upper",
  "Y_12_true_mean",
  "Y_21-yhat-yhat-match-no-var-sel-est",
  "Y_21-yhat-yhat-match-no-var-sel-se",
  "Y_21-yhat-yhat-match-no-var-sel-lower",
  "Y_21-yhat-yhat-match-no-var-sel-upper",
  "Y_21-yhat-yhat-match-var-sel-est",
  "Y_21-yhat-yhat-match-var-sel-se",
  "Y_21-yhat-yhat-match-var-sel-lower",
  "Y_21-yhat-yhat-match-var-sel-upper",
  "Y_21-yhat-y-match-no-var-sel-est",
  "Y_21-yhat-y-match-no-var-sel-se",
  "Y_21-yhat-y-match-no-var-sel-lower",
  "Y_21-yhat-y-match-no-var-sel-upper",
  "Y_21-yhat-y-match-var-sel-est",
  "Y_21-yhat-y-match-var-sel-se",
  "Y_21-yhat-y-match-var-sel-lower",
  "Y_21-yhat-y-match-var-sel-upper",
  "Y_21-mi-no-var-sel-est",
  "Y_21-mi-no-var-sel-se",
  "Y_21-mi-no-var-sel-lower",
  "Y_21-mi-no-var-sel-upper",
  "Y_21-mi-var-sel-est",
  "Y_21-mi-var-sel-se",
  "Y_21-mi-var-sel-lower",
  "Y_21-mi-var-sel-upper",
  "Y_21_true_mean",
  "Y_22-yhat-yhat-match-no-var-sel-est",
  "Y_22-yhat-yhat-match-no-var-sel-se",
  "Y_22-yhat-yhat-match-no-var-sel-lower",
  "Y_22-yhat-yhat-match-no-var-sel-upper",
  "Y_22-yhat-yhat-match-var-sel-est",
  "Y_22-yhat-yhat-match-var-sel-se",
  "Y_22-yhat-yhat-match-var-sel-lower",
  "Y_22-yhat-yhat-match-var-sel-upper",
  "Y_22-yhat-y-match-no-var-sel-est",
  "Y_22-yhat-y-match-no-var-sel-se",
  "Y_22-yhat-y-match-no-var-sel-lower",
  "Y_22-yhat-y-match-no-var-sel-upper",
  "Y_22-yhat-y-match-var-sel-est",
  "Y_22-yhat-y-match-var-sel-se",
  "Y_22-yhat-y-match-var-sel-lower",
  "Y_22-yhat-y-match-var-sel-upper",
  "Y_22-mi-no-var-sel-est",
  "Y_22-mi-no-var-sel-se",
  "Y_22-mi-no-var-sel-lower",
  "Y_22-mi-no-var-sel-upper",
  "Y_22-mi-var-sel-est",
  "Y_22-mi-var-sel-se",
  "Y_22-mi-var-sel-lower",
  "Y_22-mi-var-sel-upper",
  "Y_22_true_mean"
)

true_values <- res[,c("Y_11_true_mean", "Y_12_true_mean", "Y_21_true_mean", "Y_22_true_mean")]

df <- res |> 
  pivot_longer(cols = everything(), names_to = "column", values_to = "value") |>
  mutate(
    y_name   = str_sub(column, start = 1, end = 4),
    est_name = ifelse(
      !is.na(str_match(column, "yhat-yhat")),
      ifelse(
        !is.na(str_match(column, "no-var-sel")),
        "yhat-yhat-no-var-sel",
        "yhat-yhat-var-sel"
      ),
      ifelse(
        is.na(str_match(column, "mi")),
        ifelse(
          !is.na(str_match(column, "no-var-sel")),
          "yhat-y-no-var-sel",
          "yhat-y-var-sel"
        ),
        ifelse(
          !is.na(str_match(column, "no-var-sel")),
          "glm-no-var-sel",
          "glm-var-sel"
        )
      )
    ) |> as.character(),
    column = str_extract(column, "(est|upper|lower|true)")
  ) |> 
  mutate(column = ifelse(is.na(column), "se", column)) |>
  pivot_wider(names_from = column,
              values_from = value,
              values_fn = list) |>
  select(-true) |>
  unnest(cols = c(est, se, lower, upper))

true_values <- true_values[1,] |> unlist()

df$true <- rep(true_values, each = NROW(df) / 4)

pp <- ggplot(data = df, aes(x = est_name, y = est)) + 
  geom_violin(alpha = 0.8, draw_quantiles = 1:9 / 10, scale = "width") +
  geom_hline(aes(yintercept = true), color = "red", linetype = "dashed") +
  stat_summary(fun = function(x) mean(x, na.rm = TRUE), geom = "point") +
  facet_wrap(~ y_name, ncol = 4, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("Estimator name") +
  ylab("Predicted mean")

pp2 <- df |> 
  mutate(covr = lower < true & true < upper) |>
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
  #mutate(est_name = paste0(est_name, " - ", y_name)) |> 
  ggplot(aes(y = est_name, x = mean)) +
  geom_point(col = "blue", size = 5) +
  geom_errorbar(aes(xmin = lower, xmax = upper)) +
  geom_vline(aes(xintercept = .95), color = "red", linetype = "dashed") +
  facet_wrap(~y_name, scale = "free_x") +
  theme_bw() +
  xlab("Coverage") +
  ylab("Estimator")

saveRDS(res, file = "results/yang2020-pmm-500-sims.rds")
ggsave("results/yang2020-pmm-500-sims-plot-errors.png", pp)
ggsave("results/yang2020-pmm-500-sims-plot-coverage.png", pp2)
