library(nonprobsvy)
library(doSNOW)
library(progress)
library(foreach)
library(xtable)
library(data.table)
library(Rcpp)


# functions ---------------------------------------------------------------

sourceCpp("codes/syspss.cpp") ## randomized systematic PPS
source("codes/functions.R")

# generate data -----------------------------------------------------------

seed_for_sim <- 2024
set.seed(seed_for_sim)

## pop and sample sizes
N <- 100000
n_a <- c(500,1000)
n_b <- 500
n_a1 <- 0.7 * n_a
n_a2 <- 0.3 * n_a
## generate data
x <- rnorm(N, 2, 1)
z <- x + max(x)
e <- rnorm(N)
y1 <- 1 + 2*x + e
y2 <- 3 + x + 2*e
y3 <- 2.5 + 0.5*x^2 +e
strata <- x <= 2
pop <- data.frame(x, z, pi_z = n_b*z/sum(z), y1, y2, y3, strata)

# calculate p_jk ----------------------------------------------------------

# main simulation ---------------------------------------------------------

## setup for parallel computation
sims <- 5000
cores <- 8
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(total = sims)
opts <- list(progress = \(n) pb$tick())

results_simulation1 <- foreach(k=1:sims, .combine = rbind,
                               .packages = c("survey", "nonprobsvy"),
                               .options.snow = opts) %dopar% {
                 
  ## nonprob sample
  pop1 <- subset(pop, strata == TRUE)
  pop2 <- subset(pop, strata == FALSE)
  sample_a_500 <- rbind(pop1[sample(1:nrow(pop1), n_a1[1]), ],
                        pop2[sample(1:nrow(pop2), n_a2[1]), ])
  sample_a_1000 <- rbind(pop1[sample(1:nrow(pop1), n_a1[2]),],
                         pop2[sample(1:nrow(pop2), n_a2[2]),])
  
  ## sample prob 
  sample_b <- pop[sample(1:N, n_b),]
  sample_b$w_b <- N/n_b
  svy_b <- svydesign(ids= ~1, weights = ~ w_b, data = sample_b)
  
  sample_b_pps <- pop[syspss_r(pop$z, n_b),]
  svy_b_pps <- svydesign(ids= ~1, prob = ~ pi_z, data = sample_b_pps, pps = HR())
  
  ## esitmators
  ## true
  trues <- colMeans(pop[, c("y1", "y2", "y3")])
  ## naive
  naive_500 <- colMeans(sample_a_500[, c("y1", "y2", "y3")])
  naive_1000 <- colMeans(sample_a_1000[, c("y1", "y2", "y3")])
  

  # srs ---------------------------------------------------------------------

  ## glm
  mi_glm_500 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b,
                        se = F, family_outcome = "gaussian", method_outcome = "glm")
  
  mi_glm_1000 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b,
                         se = F, family_outcome = "gaussian", method_outcome = "glm")
  
  ## nn with k=5
  mi_nn_500 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b,
                        se = F, family_outcome = "gaussian", method_outcome = "nn", control_outcome = controlOut(k = 5))
  
  mi_nn_1000 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b,
                       se = F, family_outcome = "gaussian", method_outcome = "nn", control_outcome = controlOut(k = 5))
  
  ## hat-hat
  mi_pmm1a_500 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 1, predictive_match = 2))
  ## hat-y
  mi_pmm1b_500 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 1, predictive_match = 1))
  ## hat-hat
  mi_pmm1a_1000 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 1, predictive_match = 2))
  ## hat-y
  mi_pmm1b_1000 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 1, predictive_match = 1))
  
  ## hat-hat
  mi_pmm5a_500 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 5, predictive_match = 2))
  ## hat-y
  mi_pmm5b_500 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 5, predictive_match = 1))
  ## hat-hat
  mi_pmm5a_1000 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b,
                           se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                           control_outcome = controlOut(k = 5, predictive_match = 2))
  ## hat-y
  mi_pmm5b_1000 <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b,
                           se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                           control_outcome = controlOut(k = 5, predictive_match = 1))
  

  # pps ---------------------------------------------------------------------

  ## glm
  mi_glm_500_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b_pps,
                        se = F, family_outcome = "gaussian", method_outcome = "glm")
  
  mi_glm_1000_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b_pps,
                         se = F, family_outcome = "gaussian", method_outcome = "glm")
  
  ## nn with k=5
  mi_nn_500_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b_pps,
                       se = F, family_outcome = "gaussian", method_outcome = "nn", control_outcome = controlOut(k = 5))
  
  mi_nn_1000_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b_pps,
                        se = F, family_outcome = "gaussian", method_outcome = "nn", control_outcome = controlOut(k = 5))
  
  ## hat-hat
  mi_pmm1a_500_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b_pps,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 1, predictive_match = 2))
  ## hat-y
  mi_pmm1b_500_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b_pps,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 1, predictive_match = 1))
  ## hat-hat
  mi_pmm1a_1000_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b_pps,
                           se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                           control_outcome = controlOut(k = 1, predictive_match = 2))
  ## hat-y
  mi_pmm1b_1000_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b_pps,
                           se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                           control_outcome = controlOut(k = 1, predictive_match = 1))
  
  ## hat-hat
  mi_pmm5a_500_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b_pps,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 5, predictive_match = 2))
  ## hat-y
  mi_pmm5b_500_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_500, svydesign = svy_b_pps,
                          se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                          control_outcome = controlOut(k = 5, predictive_match = 1))
  ## hat-hat
  mi_pmm5a_1000_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b_pps,
                           se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                           control_outcome = controlOut(k = 5, predictive_match = 2))
  ## hat-y
  mi_pmm5b_1000_pps <- nonprob(outcome = y1 + y2 + y3 ~ x, data = sample_a_1000, svydesign = svy_b_pps,
                           se = F, family_outcome = "gaussian", method_outcome = "pmm", 
                           control_outcome = controlOut(k = 5, predictive_match = 1))
  
  data.frame(k=k,
             y=c("y1","y2", "y3"), 
             trues=trues,
             naive_5=naive_500, naive_10=naive_1000, 
             ## srs
             glm_5=mi_glm_500$output$mean, glm_10=mi_glm_1000$output$mean,
             nn_5=mi_nn_500$output$mean, nn_10=mi_nn_1000$output$mean,
             pmm1a_5=mi_pmm1a_500$output$mean, pmm1b_5=mi_pmm1b_500$output$mean,
             pmm1a_10=mi_pmm1a_1000$output$mean,pmm1b_10=mi_pmm1b_1000$output$mean,
             pmm5a_5 =mi_pmm5a_500$output$mean, pmm5b_5=mi_pmm5b_500$output$mean,
             pmm5a_10=mi_pmm5a_1000$output$mean,pmm5b_10=mi_pmm5b_1000$output$mean,
             
             ## pps
             glm_5_pps=mi_glm_500_pps$output$mean, glm_10_pps=mi_glm_1000_pps$output$mean,
             nn_5_pps=mi_nn_500_pps$output$mean, nn_10_pps=mi_nn_1000_pps$output$mean,
             pmm1a_5_pps=mi_pmm1a_500_pps$output$mean, pmm1b_5_pps=mi_pmm1b_500_pps$output$mean,
             pmm1a_10_pps=mi_pmm1a_1000_pps$output$mean,pmm1b_10_pps=mi_pmm1b_1000_pps$output$mean,
             pmm5a_5_pps =mi_pmm5a_500_pps$output$mean, pmm5b_5_pps=mi_pmm5b_500_pps$output$mean,
             pmm5a_10_pps=mi_pmm5a_1000_pps$output$mean,pmm5b_10_pps=mi_pmm5b_1000_pps$output$mean)

}

stopCluster(cl)


# processing results ------------------------------------------------------

## processing results
setDT(results_simulation1)

results_simulation1_process <- results_simulation1 |> melt(id.vars = 1:3)
results_simulation1_process[, c("est", "sample", "type"):=tstrsplit(variable, "_")]
results_simulation1_process[is.na(type), type:="srs"]

tab1 <- results_simulation1_process[, .(bias=mean(value)-mean(trues), se = sd(value), 
                                rmse = sqrt((mean(value)-mean(trues))^2 + var(value)) ), .(type, est, sample, y)] |>
  melt(id.vars = c(1, 4,2,3)) |>
  transform(y=paste(y, variable, sep = "_")) |>
  transform(variable=NULL,
            value = value*100) |>
  dcast(... ~ y, value.var = "value")


tab1[, est:=factor(est, c("naive",  "glm", "nn","pmm1a", "pmm1b", "pmm5a", "pmm5b"), 
                   c("Naive", "GLM", "NN", "PMM1A", "PMM1B", "PMM5A", "PMM5B")
                   , ordered = T)]
tab1[, sample:=factor(sample, c(5,10), ordered = T)]
setcolorder(tab1, c("type", "sample", "est", "y1_bias", "y1_se", "y1_rmse",
                    "y2_bias", "y2_se", "y2_rmse", "y3_bias", "y3_se", "y3_rmse"))

## report table1
tab1[order(sample,-type, est),][, ":="(sample=NULL,type=NULL)] |>
  xtable() |>
  print.xtable(include.rownames = F)


## ci coverage
saveRDS(results_simulation1, file = "results/sim1-paper-results.RDS")

