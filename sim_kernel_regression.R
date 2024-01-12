#install.packages("gplm")
library(nonprobsvy)
library(gplm)
library(dbscan)
library(doSNOW)
library(progress)

set.seed(321)

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

p1 <- exp(x2)/(1+exp(x2))
p2 <- exp(x1)/(1+exp(x1))
population <- data.frame(
  x1,
  x2,
  y1 = 1 + x1 * .2 + x2 * .1 + epsilon, # linear
  y2 = -1.2 + (x1 - 0.5) ^ 2 + atan(x2) ^ (3 + sin(x1 + x2)) + sin(x1) * cos(x2) + epsilon, # weird
  y3 = x1 * x2 * epsilon, # additive
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
  flag_bd1 <- pmin(
    rbinom(n = 1:N, size = 1, prob = p1),
    epsilon > quantile(epsilon, .8) |
      quantile(epsilon, .2) > epsilon,
    rbinom(n = 1:N, size = 1, prob = p2)
  )
  
  sample_prob <- svydesign(ids= ~1, weights = ~ base_w_srs,
                           data = subset(population, flag_srs == 1))
  
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
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  glm3 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "glm",
    pop_size = N,
    family_outcome = "gaussian"
  )
  
  pmm3 <- nonprob(
    outcome = y3 ~ x1 + x2,
    data = population[flag_bd1 == 1, , drop = FALSE],
    svydesign = sample_prob,
    method_outcome = "pmm",
    pop_size = N,
    family_outcome = "gaussian",
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE)
  )
  
  ## kernel regression
  
  ## y1
  
  yhat_prob <- kreg(
    model.matrix(y1 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    y = population[flag_bd1 == 1, "y1", drop = FALSE], 
    grid = sample_prob$variables[,c("x1", "x2")], 
    kernel = "epanechnikov"
  )
  
  yhat_nons <- kreg(
    model.matrix(y1 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    y = population[flag_bd1 == 1, "y1", drop = FALSE], 
    grid = model.matrix(y1 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    kernel = "epanechnikov"
  )
  
  yhat_prob <- yhat_prob$y[yhat_prob$rearrange]
  yhat_nons <- yhat_nons$y[yhat_nons$rearrange]
  
  nn1 <- kNN(x = cbind(population[flag_bd1 == 1, "y1", drop = FALSE]), k = KK,
             query = cbind(yhat_prob))
  
  nn2 <- kNN(x = cbind(yhat_nons), k = KK,
             query = cbind(yhat_prob))
  
  y_nu1 <- apply(nn1$id, MARGIN = 1,
                 FUN = function(x) mean(population[flag_bd1 == 1, "y1"][x]))
  
  y_nu2 <- apply(nn2$id, MARGIN = 1,
                 FUN = function(x) mean(population[flag_bd1 == 1, "y1"][x]))
  
  (mu_hat_y1   <- weighted.mean(yhat_prob, 1 / sample_prob$allprob[,1]))
  mu_hat_1_y1 <- weighted.mean(y_nu1, 1 / sample_prob$allprob[,1])
  mu_hat_2_y1 <- weighted.mean(y_nu2, 1 / sample_prob$allprob[,1])
  
  ## y2
  
  yhat_prob <- kreg(
    model.matrix(y2 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    y = population[flag_bd1 == 1, "y2", drop = FALSE], 
    grid = sample_prob$variables[,c("x1", "x2")], 
    kernel = "epanechnikov"
  )
  
  yhat_nons <- kreg(
    model.matrix(y2 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    y = population[flag_bd1 == 1, "y2", drop = FALSE], 
    grid = model.matrix(y2 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    kernel = "epanechnikov"
  )
  
  yhat_prob <- yhat_prob$y[yhat_prob$rearrange]
  yhat_nons <- yhat_nons$y[yhat_nons$rearrange]
  
  nn1 <- kNN(x = cbind(population[flag_bd1 == 1, "y2", drop = FALSE]), k = KK,
             query = cbind(yhat_prob))
  
  nn2 <- kNN(x = cbind(yhat_nons), k = KK,
             query = cbind(yhat_prob))
  
  y_nu1 <- apply(nn1$id, MARGIN = 1,
                 FUN = function(x) mean(population[flag_bd1 == 1, "y2"][x]))
  
  y_nu2 <- apply(nn2$id, MARGIN = 1,
                 FUN = function(x) mean(population[flag_bd1 == 1, "y2"][x]))
  
  (mu_hat_y2   <- weighted.mean(yhat_prob, 1 / sample_prob$allprob[,1]))
  mu_hat_1_y2 <- weighted.mean(y_nu1, 1 / sample_prob$allprob[,1])
  mu_hat_2_y2 <- weighted.mean(y_nu2, 1 / sample_prob$allprob[,1])
  
  ## y3
  
  yhat_prob <- kreg(
    model.matrix(y3 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    y = population[flag_bd1 == 1, "y3", drop = FALSE], 
    grid = sample_prob$variables[,c("x1", "x2")], 
    kernel = "epanechnikov"
  )
  
  yhat_nons <- kreg(
    model.matrix(y3 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    y = population[flag_bd1 == 1, "y3", drop = FALSE], 
    grid = model.matrix(y3 ~ x1 + x2 - 1, data = population[flag_bd1 == 1, , drop = FALSE]), 
    kernel = "epanechnikov"
  )
  
  yhat_prob <- yhat_prob$y[yhat_prob$rearrange]
  yhat_nons <- yhat_nons$y[yhat_nons$rearrange]
  
  nn1 <- kNN(x = cbind(population[flag_bd1 == 1, "y3", drop = FALSE]), k = KK,
             query = cbind(yhat_prob))
  
  nn2 <- kNN(x = cbind(yhat_nons), k = KK,
             query = cbind(yhat_prob))
  
  y_nu1 <- apply(nn1$id, MARGIN = 1,
                 FUN = function(x) mean(population[flag_bd1 == 1, "y3"][x]))
  
  y_nu2 <- apply(nn2$id, MARGIN = 1,
                 FUN = function(x) mean(population[flag_bd1 == 1, "y3"][x]))
  
  (mu_hat_y3   <- weighted.mean(yhat_prob, 1 / sample_prob$allprob[,1]))
  mu_hat_1_y3 <- weighted.mean(y_nu1, 1 / sample_prob$allprob[,1])
  mu_hat_2_y3 <- weighted.mean(y_nu2, 1 / sample_prob$allprob[,1])
  
  cbind(
    pmm1$output$mean, glm1$output$mean, mu_hat_y1,
    mu_hat_1_y1, mu_hat_2_y1, mean(population$y1),
    pmm2$output$mean, glm2$output$mean, mu_hat_y2,
    mu_hat_1_y2, mu_hat_2_y2, mean(population$y2),
    pmm3$output$mean, glm3$output$mean, mu_hat_y3,
    mu_hat_1_y3, mu_hat_2_y3, mean(population$y3)
  )
}

stopCluster(cl)

df <- data.frame(
  bias = c(
    mean(res[,  1] - res[,  6], na.rm = TRUE),
    mean(res[,  2] - res[,  6], na.rm = TRUE),
    mean(res[,  3] - res[,  6], na.rm = TRUE),
    mean(res[,  4] - res[,  6], na.rm = TRUE),
    mean(res[,  5] - res[,  6], na.rm = TRUE),
    mean(res[,  7] - res[, 12], na.rm = TRUE),
    mean(res[,  8] - res[, 12], na.rm = TRUE),
    mean(res[,  9] - res[, 12], na.rm = TRUE),
    mean(res[, 10] - res[, 12], na.rm = TRUE),
    mean(res[, 11] - res[, 12], na.rm = TRUE),
    mean(res[, 13] - res[, 18], na.rm = TRUE),
    mean(res[, 14] - res[, 18], na.rm = TRUE),
    mean(res[, 15] - res[, 18], na.rm = TRUE),
    mean(res[, 16] - res[, 18], na.rm = TRUE),
    mean(res[, 17] - res[, 18], na.rm = TRUE)
  ),
  mse = c(
    mean((res[,  1] - res[,  6]) ^ 2, na.rm = TRUE),
    mean((res[,  2] - res[,  6]) ^ 2, na.rm = TRUE),
    mean((res[,  3] - res[,  6]) ^ 2, na.rm = TRUE),
    mean((res[,  4] - res[,  6]) ^ 2, na.rm = TRUE),
    mean((res[,  5] - res[,  6]) ^ 2, na.rm = TRUE),
    mean((res[,  7] - res[, 12]) ^ 2, na.rm = TRUE),
    mean((res[,  8] - res[, 12]) ^ 2, na.rm = TRUE),
    mean((res[,  9] - res[, 12]) ^ 2, na.rm = TRUE),
    mean((res[, 10] - res[, 12]) ^ 2, na.rm = TRUE),
    mean((res[, 11] - res[, 12]) ^ 2, na.rm = TRUE),
    mean((res[, 13] - res[, 18]) ^ 2, na.rm = TRUE),
    mean((res[, 14] - res[, 18]) ^ 2, na.rm = TRUE),
    mean((res[, 15] - res[, 18]) ^ 2, na.rm = TRUE),
    mean((res[, 16] - res[, 18]) ^ 2, na.rm = TRUE),
    mean((res[, 17] - res[, 18]) ^ 2, na.rm = TRUE)
  ),
  mae = c(
    mean(abs(res[,  1] - res[,  6]), na.rm = TRUE),
    mean(abs(res[,  2] - res[,  6]), na.rm = TRUE),
    mean(abs(res[,  3] - res[,  6]), na.rm = TRUE),
    mean(abs(res[,  4] - res[,  6]), na.rm = TRUE),
    mean(abs(res[,  5] - res[,  6]), na.rm = TRUE),
    mean(abs(res[,  7] - res[, 12]), na.rm = TRUE),
    mean(abs(res[,  8] - res[, 12]), na.rm = TRUE),
    mean(abs(res[,  9] - res[, 12]), na.rm = TRUE),
    mean(abs(res[, 10] - res[, 12]), na.rm = TRUE),
    mean(abs(res[, 11] - res[, 12]), na.rm = TRUE),
    mean(abs(res[, 13] - res[, 18]), na.rm = TRUE),
    mean(abs(res[, 14] - res[, 18]), na.rm = TRUE),
    mean(abs(res[, 15] - res[, 18]), na.rm = TRUE),
    mean(abs(res[, 16] - res[, 18]), na.rm = TRUE),
    mean(abs(res[, 17] - res[, 18]), na.rm = TRUE)
  ),
  row.names = c(
    "linear glm-pmm",
    "linear glm",
    "linear kernel imputation",
    "linear kernel-pmm yhat - y",
    "linear kernel-pmm yhat - yhat",
    "non-linear glm-pmm",
    "non-linear glm",
    "non-linear kernel imputation",
    "non-linear kernel-pmm yhat - y",
    "non-linear kernel-pmm yhat - yhat",
    "additive glm-pmm",
    "additive glm",
    "additive kernel imputation",
    "additive kernel-pmm yhat - y",
    "additive kernel-pmm yhat - yhat"
  )
)

saveRDS(df, file = "results/custom-pmm-with-kernel-regression-500-sims.rds")
