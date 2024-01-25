
library(RANN)
library(data.table)

KK <- 5
N <- 100000
n_a <- 500
sims <- 1000
n_b <- 1000
n_b1 <- 0.7 * n_b
n_b2 <- 0.3 * n_b
x <- rnorm(N, 2, 1)
e <- rnorm(N)

y1 <- 1 + 2*x + e
y2 <- 3 + x + 2*e
y3 <- 2.5 + 0.5*x^2 +e

strata <- x <= 2
pop <- data.frame(x, y1, y2, y3, strata)

results_glm <- results <- list()

for (i in 1:sims) {
  if (i %% 100 == 0) print(i)
  set.seed(i)
  sample_a <- pop[sample(1:N, n_a),]
  sample_a$w_a <- N/n_a
  svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_a)
  pop1 <- subset(pop, strata == TRUE)
  pop2 <- subset(pop, strata == FALSE)
  sample_b <- rbind(pop1[sample(1:nrow(pop1), n_b1), ],
                    pop2[sample(1:nrow(pop2), n_b2), ])
  
  lm2 <- lm(y2 ~ x, sample_b) ## change here
  y_nonprob <- fitted(lm2)
  y_prob <- predict(lm2, svy_a$variables)
  res <- nn2(data=y_nonprob, y_prob, k = KK)$nn.idx
  y_imp <- apply(res, 1, function(x) mean(sample_b$y2[x])) # change here
  svy_a <- update(svy_a, y_imp = y_imp)
  #svy_a_jk <- as.svrepdesign(svy_a, type = "JK1")
  y_hat <- svymean(~y_imp, svy_a)
  #v2 <- 1/N^2*sum(weights(svy_a)*(svy_a$variables$y_imp-y_hat)^2) ## centrowanie do weighted.mean
  v1 <- as.numeric(attr(y_hat, "var"))

  v <- v1

  est_mi_glm_y3 <- nonprob(
    outcome = as.formula(y2 ~ x),
    data = sample_b,
    svydesign = svy_a,
    pop_size = N,
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_outcome = controlOut(predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = T)
  )
  
  B <- 1
  boots_var <- numeric(B)
  for (b in 1:B) {
    sample_aa <- sample_a[sample(1:nrow(sample_a), nrow(sample_a)),]
    svy_a <- svydesign(ids= ~1, weights = ~ w_a, data = sample_aa)
    sample_bb <- sample_b[sample(1:nrow(sample_b), nrow(sample_b), T),]
    lm2 <- lm(y2 ~ x, sample_bb) ## change here
    y_nonprob <- fitted(lm2)
    y_prob <- predict(lm2, svy_a$variables)
    res <- nn2(data=y_nonprob, y_prob, k = KK)$nn.idx
    y_imp <- apply(res, 1, function(x) mean(sample_bb$y2[x])) # change here
    svy_a <- update(svy_a, y_imp = y_imp)
    svy_a <- as.svrepdesign(svy_a, type = "JK1")
    boots_var[b] <- svymean(~y_imp, svy_a)
  }
  
  results[[i]] <- data.frame(yhat=y_hat[1], v=v, 
                             lower=y_hat[1]-1.96*sqrt(v), upper=y_hat[1]+1.96*sqrt(v),
                             lower_b=y_hat[1]-1.96*sd(boots_var), upper_b=y_hat[1]+1.96*sd(boots_var))
  
  results_glm[[i]] <- data.frame(yhat=est_mi_glm_y3$output$mean, v=est_mi_glm_y3$output$SE^2, 
                                 lower=est_mi_glm_y3$confidence_interval$lower_bound,
                                 upper=est_mi_glm_y3$confidence_interval$upper_bound)
}

results_df <- rbindlist(results)
results_glm_df <- rbindlist(results_glm)

mean(results_df$lower < mean(y1) & results_df$upper > mean(y1))
mean(results_df$lower_b < mean(y1) & results_df$upper_b > mean(y1))
mean(results_glm_df$lower < mean(y1) & results_glm_df$upper > mean(y1))

summary(results_df)

boxplot(cbind(results_df$yhat, results_glm_df$yhat) - mean(y3))

sqrt((mean(results_df$yhat) - mean(y3))^2 + var(results_df$yhat))
sqrt((mean(results_glm_df$yhat) - mean(y3))^2 + var(results_glm_df$yhat))



