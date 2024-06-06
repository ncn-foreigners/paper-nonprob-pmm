# generate data -----------------------------------------------------------

set.seed(seed_for_sim)

## pop and sample sizes
N <- 100000
n_a <- c(500,1000)
n_b <- 500
n_a1 <- 0.7 * n_a
n_a2 <- 0.3 * n_a
## generate data
x1 <- rnorm(N, 2, 1)
x2 <- rnorm(N, 2, 1)
x3 <- rnorm(N, 2, 1)
e <- rnorm(N)
y1 <- 1 + 2*x1 + e
y2 <- -1 + x1 + x2 + x3 + e
y3 <- -10 + x1^2 + x2^2 + x3^2 +e
strata <- x1 <= 2
pop <- data.frame(x1, x2, x3, y1, y2, y3, strata)

y_formula_1 <- y1 ~ x1
y_formula_2 <- y2 ~ x1 + x2 + x3
y_formula_3_all <- y3 ~ I(x1^2) + I(x2^2) + I(x3^2)
y_formula_3_nn <- y_formula_3_all
y_formula_3_mis <- y3 ~ x1 + x2

# main simulation ---------------------------------------------------------

## setup for parallel computation

cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())


results_simulation1 <- foreach(
  k=1:sims, .combine = rbind, .packages = c("survey", "nonprobsvy"), .options.snow = opts) %dopar% {
    
    ## nonprob sample   
    pop1 <- subset(pop, strata == TRUE)
    pop2 <- subset(pop, strata == FALSE)
    sample_a_500 <- rbind(pop1[sample(1:nrow(pop1), n_a1[1]), ], pop2[sample(1:nrow(pop2), n_a2[1]), ])
    sample_a_1000 <- rbind(pop1[sample(1:nrow(pop1), n_a1[2]), ], pop2[sample(1:nrow(pop2), n_a2[2]), ])
    
    ## sample prob
    sample_b <- pop[sample(1:N, n_b), ]
    sample_b$w_b <- N / n_b
    svy_b <- svydesign(ids = ~ 1, weights = ~ w_b, data = sample_b)
    
    ## estimators
    ## true
    trues <- colMeans(pop[, c("y1", "y2", "y3")])
    ## naive
    naive_500 <- colMeans(sample_a_500[, c("y1", "y2", "y3")])
    naive_1000 <- colMeans(sample_a_1000[, c("y1", "y2", "y3")])
    
    ## glm
    mi_glm_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b)
    
    mi_glm_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b)
    
    ## nn1 (k=1)
    
    mi_nn1_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                             control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn1_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                             control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn1_500_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn1_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_inference = controlInf(pmm_exact_se = TRUE))

    mi_nn1_500b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                             control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn1_500b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                             control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn1_500b_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn1_500b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    mi_nn1_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                              control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn1_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                              control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn1_1000_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn1_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_inference = controlInf(pmm_exact_se = TRUE))

    mi_nn1_1000b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                              control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn1_1000b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                              control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn1_1000b_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn1_1000b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    ## nn5 (k=5)
    
    mi_nn5_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                             control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn5_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                             control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn5_500_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn5_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE))
   
    mi_nn5_500b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                             control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn5_500b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                             control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn5_500b_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn5_500b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    mi_nn5_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                              control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn5_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                              control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn5_1000_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_nn5_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE))
    
    mi_nn5_1000b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                              control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn5_1000b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                              control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn5_1000b_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_nn5_1000b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_outcome = controlOut(k=5), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    ## pnn1a (k=1)
    mi_pmm1a_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1a_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1a_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1a_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))

    mi_pmm1a_500b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1a_500b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1a_500b_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1a_500b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    mi_pmm1a_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1a_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1a_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1a_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    
    mi_pmm1a_1000b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1a_1000b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1a_1000b_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1a_1000b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    ## pmm5a (k=5)
    
    mi_pmm5a_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5a_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5a_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5a_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    
    mi_pmm5a_500b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5a_500b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5a_500b_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5a_500b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    mi_pmm5a_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5a_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5a_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5a_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE))
    
    mi_pmm5a_1000b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5a_1000b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5a_1000b_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5a_1000b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 2), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    ## pnn1b (k=1)
    
    mi_pmm1b_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1b_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1b_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1b_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    
    mi_pmm1b_500b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1b_500b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1b_500b_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1b_500b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    mi_pmm1b_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1b_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1b_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm1b_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    
    mi_pmm1b_1000b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1b_1000b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1b_1000b_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm1b_1000b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    ## pmm5b (k=5)
    
    mi_pmm5b_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5b_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5b_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5b_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    
    mi_pmm5b_500b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5b_500b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5b_500b_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5b_500b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
   
    mi_pmm5b_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                 control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5b_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                 control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5b_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                     control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    mi_pmm5b_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                     control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE))
    
    mi_pmm5b_1000b_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5b_1000b_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5b_1000b_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    mi_pmm5b_1000b_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 1), control_inference = controlInf(pmm_exact_se = TRUE, var_method = "bootstrap"))
    
    
    data.frame(
      k        = k,
      y        = c("y1", "y2", "y3_all", "y3_mis"),
      trues    = trues[c(1,2,3,3)],
      naive_500  = naive_500[c(1,2,3,3)],
      naive_1000 = naive_1000[c(1,2,3,3)],
      ## srs
      glm_500    = c(mi_glm_500_y1$output$mean, mi_glm_500_y2$output$mean, mi_glm_500_y3_all$output$mean, mi_glm_500_y3_mis$output$mean),
      glm_1000   = c(mi_glm_1000_y1$output$mean, mi_glm_1000_y2$output$mean, mi_glm_1000_y3_all$output$mean, mi_glm_1000_y3_mis$output$mean),
      nn1_500    = c(mi_nn1_500_y1$output$mean, mi_nn1_500_y2$output$mean, mi_nn1_500_y3_all$output$mean, mi_nn1_500_y3_mis$output$mean),
      nn1_1000   = c(mi_nn1_1000_y1$output$mean, mi_nn1_1000_y2$output$mean, mi_nn1_1000_y3_all$output$mean, mi_nn1_1000_y3_mis$output$mean),
      nn5_500    = c(mi_nn5_500_y1$output$mean, mi_nn5_500_y2$output$mean, mi_nn5_500_y3_all$output$mean, mi_nn5_500_y3_mis$output$mean), 
      nn5_1000   = c(mi_nn5_1000_y1$output$mean, mi_nn5_1000_y2$output$mean, mi_nn5_1000_y3_all$output$mean, mi_nn5_1000_y3_mis$output$mean),
      pmm1a_500  = c(mi_pmm1a_500_y1$output$mean, mi_pmm1a_500_y2$output$mean, mi_pmm1a_500_y3_all$output$mean, mi_pmm1a_500_y3_mis$output$mean),
      pmm1b_500  = c(mi_pmm1b_500_y1$output$mean, mi_pmm1b_500_y2$output$mean, mi_pmm1b_500_y3_all$output$mean, mi_pmm1b_500_y3_mis$output$mean),
      pmm1a_1000 = c(mi_pmm1a_1000_y1$output$mean, mi_pmm1a_1000_y2$output$mean, mi_pmm1a_1000_y3_all$output$mean, mi_pmm1a_1000_y3_mis$output$mean),
      pmm1b_1000 = c(mi_pmm1b_1000_y1$output$mean, mi_pmm1b_1000_y2$output$mean, mi_pmm1b_1000_y3_all$output$mean, mi_pmm1b_1000_y3_mis$output$mean),
      pmm5a_500  = c(mi_pmm5a_500_y1$output$mean, mi_pmm5a_500_y2$output$mean, mi_pmm5a_500_y3_all$output$mean, mi_pmm5a_500_y3_mis$output$mean),
      pmm5b_500  = c(mi_pmm5b_500_y1$output$mean, mi_pmm5b_500_y2$output$mean, mi_pmm5b_500_y3_all$output$mean, mi_pmm5b_500_y3_mis$output$mean),
      pmm5a_1000 = c(mi_pmm5a_1000_y1$output$mean, mi_pmm5a_1000_y2$output$mean, mi_pmm5a_1000_y3_all$output$mean, mi_pmm5a_1000_y3_mis$output$mean),
      pmm5b_1000 = c(mi_pmm5b_1000_y1$output$mean, mi_pmm5b_1000_y2$output$mean, mi_pmm5b_1000_y3_all$output$mean, mi_pmm5b_1000_y3_mis$output$mean),
      
      glm_500_ci    = c(mi_glm_500_y1$confidence_interval[1] < trues[1] & mi_glm_500_y1$confidence_interval[2] > trues[1], 
                        mi_glm_500_y2$confidence_interval[1] < trues[2] & mi_glm_500_y2$confidence_interval[2] > trues[2],
                        mi_glm_500_y3_all$confidence_interval[1] < trues[3] & mi_glm_500_y3_all$confidence_interval[2] > trues[3],
                        mi_glm_500_y3_mis$confidence_interval[1] < trues[3] & mi_glm_500_y3_mis$confidence_interval[2] > trues[3]),
      glm_1000_ci   = c(mi_glm_1000_y1$confidence_interval[1] < trues[1] & mi_glm_1000_y1$confidence_interval[2] > trues[1], 
                        mi_glm_1000_y2$confidence_interval[1] < trues[2] & mi_glm_1000_y2$confidence_interval[2] > trues[2],
                        mi_glm_1000_y3_all$confidence_interval[1] < trues[3] & mi_glm_1000_y3_all$confidence_interval[2] > trues[3],
                        mi_glm_1000_y3_mis$confidence_interval[1] < trues[3] & mi_glm_1000_y3_mis$confidence_interval[2] > trues[3]),
      nn1_500_ci    = c(mi_nn1_500_y1$confidence_interval[1] < trues[1] & mi_nn1_500_y1$confidence_interval[2] > trues[1], 
                        mi_nn1_500_y2$confidence_interval[1] < trues[2] & mi_nn1_500_y2$confidence_interval[2] > trues[2],
                        mi_nn1_500_y3_all$confidence_interval[1] < trues[3] & mi_nn1_500_y3_all$confidence_interval[2] > trues[3],
                        mi_nn1_500_y3_mis$confidence_interval[1] < trues[3] & mi_nn1_500_y3_mis$confidence_interval[2] > trues[3]),
      nn1_1000_ci   = c(mi_nn1_1000_y1$confidence_interval[1] < trues[1] & mi_nn1_1000_y1$confidence_interval[2] > trues[1], 
                        mi_nn1_1000_y2$confidence_interval[1] < trues[2] & mi_nn1_1000_y2$confidence_interval[2] > trues[2],
                        mi_nn1_1000_y3_all$confidence_interval[1] < trues[3] & mi_nn1_1000_y3_all$confidence_interval[2] > trues[3],
                        mi_nn1_1000_y3_mis$confidence_interval[1] < trues[3] & mi_nn1_1000_y3_mis$confidence_interval[2] > trues[3]),
      nn5_500_ci    = c(mi_nn5_500_y1$confidence_interval[1] < trues[1] & mi_nn5_500_y1$confidence_interval[2] > trues[1], 
                        mi_nn5_500_y2$confidence_interval[1] < trues[2] & mi_nn5_500_y2$confidence_interval[2] > trues[2],
                        mi_nn5_500_y3_all$confidence_interval[1] < trues[3] & mi_nn5_500_y3_all$confidence_interval[2] > trues[3],
                        mi_nn5_500_y3_mis$confidence_interval[1] < trues[3] & mi_nn5_500_y3_mis$confidence_interval[2] > trues[3]),
      nn5_1000_ci   = c(mi_nn5_1000_y1$confidence_interval[1] < trues[1] & mi_nn5_1000_y1$confidence_interval[2] > trues[1], 
                        mi_nn5_1000_y2$confidence_interval[1] < trues[2] & mi_nn5_1000_y2$confidence_interval[2] > trues[2],
                        mi_nn5_1000_y3_all$confidence_interval[1] < trues[3] & mi_nn5_1000_y3_all$confidence_interval[2] > trues[3],
                        mi_nn5_1000_y3_mis$confidence_interval[1] < trues[3] & mi_nn5_1000_y3_mis$confidence_interval[2] > trues[3]),
      pmm1a_500_ci  = c(mi_pmm1a_500_y1$confidence_interval[1] < trues[1] & mi_pmm1a_500_y1$confidence_interval[2] > trues[1], 
                        mi_pmm1a_500_y2$confidence_interval[1] < trues[2] & mi_pmm1a_500_y2$confidence_interval[2] > trues[2],
                        mi_pmm1a_500_y3_all$confidence_interval[1] < trues[3] & mi_pmm1a_500_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm1a_500_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1a_500_y3_mis$confidence_interval[2] > trues[3]),
      pmm1b_500_ci  = c(mi_pmm1b_500_y1$confidence_interval[1] < trues[1] & mi_pmm1b_500_y1$confidence_interval[2] > trues[1], 
                        mi_pmm1b_500_y2$confidence_interval[1] < trues[2] & mi_pmm1b_500_y2$confidence_interval[2] > trues[2],
                        mi_pmm1b_500_y3_all$confidence_interval[1] < trues[3] & mi_pmm1b_500_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm1b_500_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1b_500_y3_mis$confidence_interval[2] > trues[3]),
      pmm1a_1000_ci = c(mi_pmm1a_1000_y1$confidence_interval[1] < trues[1] & mi_pmm1a_1000_y1$confidence_interval[2] > trues[1], 
                        mi_pmm1a_1000_y2$confidence_interval[1] < trues[2] & mi_pmm1a_1000_y2$confidence_interval[2] > trues[2],
                        mi_pmm1a_1000_y3_all$confidence_interval[1] < trues[3] & mi_pmm1a_1000_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm1a_1000_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1a_1000_y3_mis$confidence_interval[2] > trues[3]),
      pmm1b_1000_ci = c(mi_pmm1b_1000_y1$confidence_interval[1] < trues[1] & mi_pmm1b_1000_y1$confidence_interval[2] > trues[1], 
                        mi_pmm1b_1000_y2$confidence_interval[1] < trues[2] & mi_pmm1b_1000_y2$confidence_interval[2] > trues[2],
                        mi_pmm1b_1000_y3_all$confidence_interval[1] < trues[3] & mi_pmm1b_1000_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm1b_1000_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1b_1000_y3_mis$confidence_interval[2] > trues[3]),
      pmm5a_500_ci  = c(mi_pmm5a_500_y1$confidence_interval[1] < trues[1] & mi_pmm5a_500_y1$confidence_interval[2] > trues[1], 
                        mi_pmm5a_500_y2$confidence_interval[1] < trues[2] & mi_pmm5a_500_y2$confidence_interval[2] > trues[2],
                        mi_pmm5a_500_y3_all$confidence_interval[1] < trues[3] & mi_pmm5a_500_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm5a_500_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5a_500_y3_mis$confidence_interval[2] > trues[3]),
      pmm5b_500_ci  = c(mi_pmm5b_500_y1$confidence_interval[1] < trues[1] & mi_pmm5b_500_y1$confidence_interval[2] > trues[1], 
                        mi_pmm5b_500_y2$confidence_interval[1] < trues[2] & mi_pmm5b_500_y2$confidence_interval[2] > trues[2],
                        mi_pmm5b_500_y3_all$confidence_interval[1] < trues[3] & mi_pmm5b_500_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm5b_500_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5b_500_y3_mis$confidence_interval[2] > trues[3]),
      pmm5a_1000_ci = c(mi_pmm5a_1000_y1$confidence_interval[1] < trues[1] & mi_pmm5a_1000_y1$confidence_interval[2] > trues[1], 
                        mi_pmm5a_1000_y2$confidence_interval[1] < trues[2] & mi_pmm5a_1000_y2$confidence_interval[2] > trues[2],
                        mi_pmm5a_1000_y3_all$confidence_interval[1] < trues[3] & mi_pmm5a_1000_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm5a_1000_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5a_1000_y3_mis$confidence_interval[2] > trues[3]),
      pmm5b_1000_ci = c(mi_pmm5b_1000_y1$confidence_interval[1] < trues[1] & mi_pmm5b_1000_y1$confidence_interval[2] > trues[1], 
                        mi_pmm5b_1000_y2$confidence_interval[1] < trues[2] & mi_pmm5b_1000_y2$confidence_interval[2] > trues[2],
                        mi_pmm5b_1000_y3_all$confidence_interval[1] < trues[3] & mi_pmm5b_1000_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm5b_1000_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5b_1000_y3_mis$confidence_interval[2] > trues[3]),
      
      nn1_500b_ci    = c(mi_nn1_500b_y1$confidence_interval[1] < trues[1] & mi_nn1_500b_y1$confidence_interval[2] > trues[1], 
                        mi_nn1_500b_y2$confidence_interval[1] < trues[2] & mi_nn1_500b_y2$confidence_interval[2] > trues[2],
                        mi_nn1_500b_y3_all$confidence_interval[1] < trues[3] & mi_nn1_500b_y3_all$confidence_interval[2] > trues[3],
                        mi_nn1_500b_y3_mis$confidence_interval[1] < trues[3] & mi_nn1_500b_y3_mis$confidence_interval[2] > trues[3]),
      nn1_1000b_ci   = c(mi_nn1_1000b_y1$confidence_interval[1] < trues[1] & mi_nn1_1000b_y1$confidence_interval[2] > trues[1], 
                        mi_nn1_1000b_y2$confidence_interval[1] < trues[2] & mi_nn1_1000b_y2$confidence_interval[2] > trues[2],
                        mi_nn1_1000b_y3_all$confidence_interval[1] < trues[3] & mi_nn1_1000b_y3_all$confidence_interval[2] > trues[3],
                        mi_nn1_1000b_y3_mis$confidence_interval[1] < trues[3] & mi_nn1_1000b_y3_mis$confidence_interval[2] > trues[3]),
      nn5_500b_ci    = c(mi_nn5_500b_y1$confidence_interval[1] < trues[1] & mi_nn5_500b_y1$confidence_interval[2] > trues[1], 
                        mi_nn5_500b_y2$confidence_interval[1] < trues[2] & mi_nn5_500b_y2$confidence_interval[2] > trues[2],
                        mi_nn5_500b_y3_all$confidence_interval[1] < trues[3] & mi_nn5_500b_y3_all$confidence_interval[2] > trues[3],
                        mi_nn5_500b_y3_mis$confidence_interval[1] < trues[3] & mi_nn5_500b_y3_mis$confidence_interval[2] > trues[3]),
      nn5_1000b_ci   = c(mi_nn5_1000b_y1$confidence_interval[1] < trues[1] & mi_nn5_1000b_y1$confidence_interval[2] > trues[1], 
                        mi_nn5_1000b_y2$confidence_interval[1] < trues[2] & mi_nn5_1000b_y2$confidence_interval[2] > trues[2],
                        mi_nn5_1000b_y3_all$confidence_interval[1] < trues[3] & mi_nn5_1000b_y3_all$confidence_interval[2] > trues[3],
                        mi_nn5_1000b_y3_mis$confidence_interval[1] < trues[3] & mi_nn5_1000b_y3_mis$confidence_interval[2] > trues[3]),
      pmm1a_500b_ci  = c(mi_pmm1a_500b_y1$confidence_interval[1] < trues[1] & mi_pmm1a_500b_y1$confidence_interval[2] > trues[1], 
                        mi_pmm1a_500b_y2$confidence_interval[1] < trues[2] & mi_pmm1a_500b_y2$confidence_interval[2] > trues[2],
                        mi_pmm1a_500b_y3_all$confidence_interval[1] < trues[3] & mi_pmm1a_500b_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm1a_500b_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1a_500b_y3_mis$confidence_interval[2] > trues[3]),
      pmm1b_500b_ci  = c(mi_pmm1b_500b_y1$confidence_interval[1] < trues[1] & mi_pmm1b_500b_y1$confidence_interval[2] > trues[1], 
                        mi_pmm1b_500b_y2$confidence_interval[1] < trues[2] & mi_pmm1b_500b_y2$confidence_interval[2] > trues[2],
                        mi_pmm1b_500b_y3_all$confidence_interval[1] < trues[3] & mi_pmm1b_500b_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm1b_500b_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1b_500b_y3_mis$confidence_interval[2] > trues[3]),
      pmm1a_1000b_ci = c(mi_pmm1a_1000b_y1$confidence_interval[1] < trues[1] & mi_pmm1a_1000b_y1$confidence_interval[2] > trues[1], 
                        mi_pmm1a_1000b_y2$confidence_interval[1] < trues[2] & mi_pmm1a_1000b_y2$confidence_interval[2] > trues[2],
                        mi_pmm1a_1000b_y3_all$confidence_interval[1] < trues[3] & mi_pmm1a_1000b_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm1a_1000b_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1a_1000b_y3_mis$confidence_interval[2] > trues[3]),
      pmm1b_1000b_ci = c(mi_pmm1b_1000b_y1$confidence_interval[1] < trues[1] & mi_pmm1b_1000b_y1$confidence_interval[2] > trues[1], 
                        mi_pmm1b_1000b_y2$confidence_interval[1] < trues[2] & mi_pmm1b_1000b_y2$confidence_interval[2] > trues[2],
                        mi_pmm1b_1000b_y3_all$confidence_interval[1] < trues[3] & mi_pmm1b_1000b_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm1b_1000b_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1b_1000b_y3_mis$confidence_interval[2] > trues[3]),
      pmm5a_500b_ci  = c(mi_pmm5a_500b_y1$confidence_interval[1] < trues[1] & mi_pmm5a_500b_y1$confidence_interval[2] > trues[1], 
                        mi_pmm5a_500b_y2$confidence_interval[1] < trues[2] & mi_pmm5a_500b_y2$confidence_interval[2] > trues[2],
                        mi_pmm5a_500b_y3_all$confidence_interval[1] < trues[3] & mi_pmm5a_500b_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm5a_500b_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5a_500b_y3_mis$confidence_interval[2] > trues[3]),
      pmm5b_500b_ci  = c(mi_pmm5b_500b_y1$confidence_interval[1] < trues[1] & mi_pmm5b_500b_y1$confidence_interval[2] > trues[1], 
                        mi_pmm5b_500b_y2$confidence_interval[1] < trues[2] & mi_pmm5b_500b_y2$confidence_interval[2] > trues[2],
                        mi_pmm5b_500b_y3_all$confidence_interval[1] < trues[3] & mi_pmm5b_500b_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm5b_500b_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5b_500b_y3_mis$confidence_interval[2] > trues[3]),
      pmm5a_1000b_ci = c(mi_pmm5a_1000b_y1$confidence_interval[1] < trues[1] & mi_pmm5a_1000b_y1$confidence_interval[2] > trues[1], 
                        mi_pmm5a_1000b_y2$confidence_interval[1] < trues[2] & mi_pmm5a_1000b_y2$confidence_interval[2] > trues[2],
                        mi_pmm5a_1000b_y3_all$confidence_interval[1] < trues[3] & mi_pmm5a_1000b_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm5a_1000b_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5a_1000b_y3_mis$confidence_interval[2] > trues[3]),
      pmm5b_1000b_ci = c(mi_pmm5b_1000b_y1$confidence_interval[1] < trues[1] & mi_pmm5b_1000b_y1$confidence_interval[2] > trues[1], 
                        mi_pmm5b_1000b_y2$confidence_interval[1] < trues[2] & mi_pmm5b_1000b_y2$confidence_interval[2] > trues[2],
                        mi_pmm5b_1000b_y3_all$confidence_interval[1] < trues[3] & mi_pmm5b_1000b_y3_all$confidence_interval[2] > trues[3],
                        mi_pmm5b_1000b_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5b_1000b_y3_mis$confidence_interval[2] > trues[3])
      
    )

}

results_simulation1_no_v2 <- foreach(
  k=1:sims, .combine = rbind, .packages = c("survey", "nonprobsvy"), .options.snow = opts) %dopar% {
 
    ## nonprob sample   
    pop1 <- subset(pop, strata == TRUE)
    pop2 <- subset(pop, strata == FALSE)
    sample_a_500 <- rbind(pop1[sample(1:nrow(pop1), n_a1[1]), ], pop2[sample(1:nrow(pop2), n_a2[1]), ])
    sample_a_1000 <- rbind(pop1[sample(1:nrow(pop1), n_a1[2]), ], pop2[sample(1:nrow(pop2), n_a2[2]), ])
    
    ## sample prob
    sample_b <- pop[sample(1:N, n_b), ]
    sample_b$w_b <- N / n_b
    svy_b <- svydesign(ids = ~ 1, weights = ~ w_b, data = sample_b)
    
    ## estimators
    ## true
    trues <- colMeans(pop[, c("y1", "y2", "y3")])
    ## naive
    naive_500 <- colMeans(sample_a_500[, c("y1", "y2", "y3")])
    naive_1000 <- colMeans(sample_a_1000[, c("y1", "y2", "y3")])
    
    ## glm
    mi_glm_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b)
    mi_glm_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b)
    
    mi_glm_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b)
    mi_glm_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b)
    
    ## nn1 (k=1)
    
    mi_nn1_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "nn")
    mi_nn1_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "nn")
    mi_nn1_500_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_500, svydesign = svy_b, method_outcome = "nn")
    mi_nn1_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "nn")
    
    mi_nn1_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn")
    mi_nn1_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn")
    mi_nn1_1000_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn")
    mi_nn1_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn")
    
    ## nn5 (k=5)
    
    mi_nn5_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", control_outcome = controlOut(k=5))
    mi_nn5_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", control_outcome = controlOut(k=5))
    mi_nn5_500_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_outcome = controlOut(k=5))
    mi_nn5_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "nn", 
                                 control_outcome = controlOut(k=5))
    
    mi_nn5_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", control_outcome = controlOut(k=5))
    mi_nn5_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", control_outcome = controlOut(k=5))
    mi_nn5_1000_y3_all <- nonprob(outcome = y_formula_3_nn, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_outcome = controlOut(k=5))
    mi_nn5_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "nn", 
                                  control_outcome = controlOut(k=5))
    
    ## pnn1a (k=1)
    mi_pmm1a_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 2))
    mi_pmm1a_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 2))
    mi_pmm1a_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 2))
    mi_pmm1a_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 2))
    
    mi_pmm1a_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 2))
    mi_pmm1a_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 2))
    mi_pmm1a_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 2))
    mi_pmm1a_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 2))
    
    ## pmm5a (k=5)
    
    mi_pmm5a_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 2))
    mi_pmm5a_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 2))
    mi_pmm5a_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 2))
    mi_pmm5a_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 2))
    
    mi_pmm5a_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 2))
    mi_pmm5a_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 2))
    mi_pmm5a_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 2))
    mi_pmm5a_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 2))
    
    ## pnn1b (k=1)
    
    mi_pmm1b_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 1))
    mi_pmm1b_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(predictive_match = 1))
    mi_pmm1b_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 1))
    mi_pmm1b_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(predictive_match = 1))
    
    mi_pmm1b_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 1))
    mi_pmm1b_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(predictive_match = 1))
    mi_pmm1b_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 1))
    mi_pmm1b_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(predictive_match = 1))
    
    ## pmm5b (k=5)
    
    mi_pmm5b_500_y1 <- nonprob(outcome = y_formula_1, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 1))
    mi_pmm5b_500_y2 <- nonprob(outcome = y_formula_2, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                               control_outcome = controlOut(k=5, predictive_match = 1))
    mi_pmm5b_500_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 1))
    mi_pmm5b_500_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_500, svydesign = svy_b, method_outcome = "pmm", 
                                   control_outcome = controlOut(k=5, predictive_match = 1))
    
    mi_pmm5b_1000_y1 <- nonprob(outcome = y_formula_1, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 1))
    mi_pmm5b_1000_y2 <- nonprob(outcome = y_formula_2, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                control_outcome = controlOut(k=5, predictive_match = 1))
    mi_pmm5b_1000_y3_all <- nonprob(outcome = y_formula_3_all, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 1))
    mi_pmm5b_1000_y3_mis <- nonprob(outcome = y_formula_3_mis, data = sample_a_1000, svydesign = svy_b, method_outcome = "pmm", 
                                    control_outcome = controlOut(k=5, predictive_match = 1))
    
    
  data.frame(
    k        = k,
    y        = c("y1", "y2", "y3_all", "y3_mis"),
    trues    = trues[c(1,2,3,3)],
    naive_500  = naive_500[c(1,2,3,3)],
    naive_1000 = naive_1000[c(1,2,3,3)],
    ## srs
    glm_500    = c(mi_glm_500_y1$output$mean, mi_glm_500_y2$output$mean, mi_glm_500_y3_all$output$mean, mi_glm_500_y3_mis$output$mean),
    glm_1000   = c(mi_glm_1000_y1$output$mean, mi_glm_1000_y2$output$mean, mi_glm_1000_y3_all$output$mean, mi_glm_1000_y3_mis$output$mean),
    nn1_500    = c(mi_nn1_500_y1$output$mean, mi_nn1_500_y2$output$mean, mi_nn1_500_y3_all$output$mean, mi_nn1_500_y3_mis$output$mean),
    nn1_1000   = c(mi_nn1_1000_y1$output$mean, mi_nn1_1000_y2$output$mean, mi_nn1_1000_y3_all$output$mean, mi_nn1_1000_y3_mis$output$mean),
    nn5_500    = c(mi_nn5_500_y1$output$mean, mi_nn5_500_y2$output$mean, mi_nn5_500_y3_all$output$mean, mi_nn5_500_y3_mis$output$mean), 
    nn5_1000   = c(mi_nn5_1000_y1$output$mean, mi_nn5_1000_y2$output$mean, mi_nn5_1000_y3_all$output$mean, mi_nn5_1000_y3_mis$output$mean),
    pmm1a_500  = c(mi_pmm1a_500_y1$output$mean, mi_pmm1a_500_y2$output$mean, mi_pmm1a_500_y3_all$output$mean, mi_pmm1a_500_y3_mis$output$mean),
    pmm1b_500  = c(mi_pmm1b_500_y1$output$mean, mi_pmm1b_500_y2$output$mean, mi_pmm1b_500_y3_all$output$mean, mi_pmm1b_500_y3_mis$output$mean),
    pmm1a_1000 = c(mi_pmm1a_1000_y1$output$mean, mi_pmm1a_1000_y2$output$mean, mi_pmm1a_1000_y3_all$output$mean, mi_pmm1a_1000_y3_mis$output$mean),
    pmm1b_1000 = c(mi_pmm1b_1000_y1$output$mean, mi_pmm1b_1000_y2$output$mean, mi_pmm1b_1000_y3_all$output$mean, mi_pmm1b_1000_y3_mis$output$mean),
    pmm5a_500  = c(mi_pmm5a_500_y1$output$mean, mi_pmm5a_500_y2$output$mean, mi_pmm5a_500_y3_all$output$mean, mi_pmm5a_500_y3_mis$output$mean),
    pmm5b_500  = c(mi_pmm5b_500_y1$output$mean, mi_pmm5b_500_y2$output$mean, mi_pmm5b_500_y3_all$output$mean, mi_pmm5b_500_y3_mis$output$mean),
    pmm5a_1000 = c(mi_pmm5a_1000_y1$output$mean, mi_pmm5a_1000_y2$output$mean, mi_pmm5a_1000_y3_all$output$mean, mi_pmm5a_1000_y3_mis$output$mean),
    pmm5b_1000 = c(mi_pmm5b_1000_y1$output$mean, mi_pmm5b_1000_y2$output$mean, mi_pmm5b_1000_y3_all$output$mean, mi_pmm5b_1000_y3_mis$output$mean),
    # srs ci
    glm_500_ci    = c(mi_glm_500_y1$confidence_interval[1] < trues[1] & mi_glm_500_y1$confidence_interval[2] > trues[1], 
                   mi_glm_500_y2$confidence_interval[1] < trues[2] & mi_glm_500_y2$confidence_interval[2] > trues[2],
                   mi_glm_500_y3_all$confidence_interval[1] < trues[3] & mi_glm_500_y3_all$confidence_interval[2] > trues[3],
                   mi_glm_500_y3_mis$confidence_interval[1] < trues[3] & mi_glm_500_y3_mis$confidence_interval[2] > trues[3]),
      glm_1000_ci   = c(mi_glm_1000_y1$confidence_interval[1] < trues[1] & mi_glm_1000_y1$confidence_interval[2] > trues[1], 
                     mi_glm_1000_y2$confidence_interval[1] < trues[2] & mi_glm_1000_y2$confidence_interval[2] > trues[2],
                     mi_glm_1000_y3_all$confidence_interval[1] < trues[3] & mi_glm_1000_y3_all$confidence_interval[2] > trues[3],
                     mi_glm_1000_y3_mis$confidence_interval[1] < trues[3] & mi_glm_1000_y3_mis$confidence_interval[2] > trues[3]),
      nn1_500_ci    = c(mi_nn1_500_y1$confidence_interval[1] < trues[1] & mi_nn1_500_y1$confidence_interval[2] > trues[1], 
                     mi_nn1_500_y2$confidence_interval[1] < trues[2] & mi_nn1_500_y2$confidence_interval[2] > trues[2],
                     mi_nn1_500_y3_all$confidence_interval[1] < trues[3] & mi_nn1_500_y3_all$confidence_interval[2] > trues[3],
                     mi_nn1_500_y3_mis$confidence_interval[1] < trues[3] & mi_nn1_500_y3_mis$confidence_interval[2] > trues[3]),
      nn1_1000_ci   = c(mi_nn1_1000_y1$confidence_interval[1] < trues[1] & mi_nn1_1000_y1$confidence_interval[2] > trues[1], 
                     mi_nn1_1000_y2$confidence_interval[1] < trues[2] & mi_nn1_1000_y2$confidence_interval[2] > trues[2],
                     mi_nn1_1000_y3_all$confidence_interval[1] < trues[3] & mi_nn1_1000_y3_all$confidence_interval[2] > trues[3],
                     mi_nn1_1000_y3_mis$confidence_interval[1] < trues[3] & mi_nn1_1000_y3_mis$confidence_interval[2] > trues[3]),
      nn5_500_ci    = c(mi_nn5_500_y1$confidence_interval[1] < trues[1] & mi_nn5_500_y1$confidence_interval[2] > trues[1], 
                     mi_nn5_500_y2$confidence_interval[1] < trues[2] & mi_nn5_500_y2$confidence_interval[2] > trues[2],
                     mi_nn5_500_y3_all$confidence_interval[1] < trues[3] & mi_nn5_500_y3_all$confidence_interval[2] > trues[3],
                     mi_nn5_500_y3_mis$confidence_interval[1] < trues[3] & mi_nn5_500_y3_mis$confidence_interval[2] > trues[3]),
      nn5_1000_ci   = c(mi_nn5_1000_y1$confidence_interval[1] < trues[1] & mi_nn5_1000_y1$confidence_interval[2] > trues[1], 
                     mi_nn5_1000_y2$confidence_interval[1] < trues[2] & mi_nn5_1000_y2$confidence_interval[2] > trues[2],
                     mi_nn5_1000_y3_all$confidence_interval[1] < trues[3] & mi_nn5_1000_y3_all$confidence_interval[2] > trues[3],
                     mi_nn5_1000_y3_mis$confidence_interval[1] < trues[3] & mi_nn5_1000_y3_mis$confidence_interval[2] > trues[3]),
      pmm1a_500_ci  = c(mi_pmm1a_500_y1$confidence_interval[1] < trues[1] & mi_pmm1a_500_y1$confidence_interval[2] > trues[1], 
                     mi_pmm1a_500_y2$confidence_interval[1] < trues[2] & mi_pmm1a_500_y2$confidence_interval[2] > trues[2],
                     mi_pmm1a_500_y3_all$confidence_interval[1] < trues[3] & mi_pmm1a_500_y3_all$confidence_interval[2] > trues[3],
                     mi_pmm1a_500_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1a_500_y3_mis$confidence_interval[2] > trues[3]),
      pmm1b_500_ci  = c(mi_pmm1b_500_y1$confidence_interval[1] < trues[1] & mi_pmm1b_500_y1$confidence_interval[2] > trues[1], 
                     mi_pmm1b_500_y3_all$confidence_interval[1] < trues[3] & mi_pmm1b_500_y3_all$confidence_interval[2] > trues[3],
                     mi_pmm1b_500_y2$confidence_interval[1] < trues[2] & mi_pmm1b_500_y2$confidence_interval[2] > trues[2],
                     mi_pmm1b_500_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1b_500_y3_mis$confidence_interval[2] > trues[3]),
      pmm1a_1000_ci = c(mi_pmm1a_1000_y1$confidence_interval[1] < trues[1] & mi_pmm1a_1000_y1$confidence_interval[2] > trues[1], 
                      mi_pmm1a_1000_y2$confidence_interval[1] < trues[2] & mi_pmm1a_1000_y2$confidence_interval[2] > trues[2],
                      mi_pmm1a_1000_y3_all$confidence_interval[1] < trues[3] & mi_pmm1a_1000_y3_all$confidence_interval[2] > trues[3],
                      mi_pmm1a_1000_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1a_1000_y3_mis$confidence_interval[2] > trues[3]),
      pmm1b_1000_ci = c(mi_pmm1b_1000_y1$confidence_interval[1] < trues[1] & mi_pmm1b_1000_y1$confidence_interval[2] > trues[1], 
                      mi_pmm1b_1000_y2$confidence_interval[1] < trues[2] & mi_pmm1b_1000_y2$confidence_interval[2] > trues[2],
                      mi_pmm1b_1000_y3_all$confidence_interval[1] < trues[3] & mi_pmm1b_1000_y3_all$confidence_interval[2] > trues[3],
                      mi_pmm1b_1000_y3_mis$confidence_interval[1] < trues[3] & mi_pmm1b_1000_y3_mis$confidence_interval[2] > trues[3]),
      pmm5a_500_ci  = c(mi_pmm5a_500_y1$confidence_interval[1] < trues[1] & mi_pmm5a_500_y1$confidence_interval[2] > trues[1], 
                     mi_pmm5a_500_y2$confidence_interval[1] < trues[2] & mi_pmm5a_500_y2$confidence_interval[2] > trues[2],
                     mi_pmm5a_500_y3_all$confidence_interval[1] < trues[3] & mi_pmm5a_500_y3_all$confidence_interval[2] > trues[3],
                     mi_pmm5a_500_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5a_500_y3_mis$confidence_interval[2] > trues[3]),
      pmm5b_500_ci  = c(mi_pmm5b_500_y1$confidence_interval[1] < trues[1] & mi_pmm5b_500_y1$confidence_interval[2] > trues[1], 
                     mi_pmm5b_500_y2$confidence_interval[1] < trues[2] & mi_pmm5b_500_y2$confidence_interval[2] > trues[2],
                     mi_pmm5b_500_y3_all$confidence_interval[1] < trues[3] & mi_pmm5b_500_y3_all$confidence_interval[2] > trues[3],
                     mi_pmm5b_500_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5b_500_y3_mis$confidence_interval[2] > trues[3]),
      pmm5a_1000_ci = c(mi_pmm5a_1000_y1$confidence_interval[1] < trues[1] & mi_pmm5a_1000_y1$confidence_interval[2] > trues[1], 
                      mi_pmm5a_1000_y2$confidence_interval[1] < trues[2] & mi_pmm5a_1000_y2$confidence_interval[2] > trues[2],
                      mi_pmm5a_1000_y3_all$confidence_interval[1] < trues[3] & mi_pmm5a_1000_y3_all$confidence_interval[2] > trues[3],
                      mi_pmm5a_1000_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5a_1000_y3_mis$confidence_interval[2] > trues[3]),
      pmm5b_1000_ci = c(mi_pmm5b_1000_y1$confidence_interval[1] < trues[1] & mi_pmm5b_1000_y1$confidence_interval[2] > trues[1], 
                      mi_pmm5b_1000_y2$confidence_interval[1] < trues[2] & mi_pmm5b_1000_y2$confidence_interval[2] > trues[2],
                      mi_pmm5b_1000_y3_all$confidence_interval[1] < trues[3] & mi_pmm5b_1000_y3_all$confidence_interval[2] > trues[3],
                      mi_pmm5b_1000_y3_mis$confidence_interval[1] < trues[3] & mi_pmm5b_1000_y3_mis$confidence_interval[2] > trues[3])
    
  )
}

stopCluster(cl)


# processing results ------------------------------------------------------

## processing results
setDT(results_simulation1)

results_simulation1_process <- results_simulation1 |> melt(id.vars = 1:3)
results_simulation1_process[, c("est", "sample", "ci"):=tstrsplit(variable, "_")]
results_simulation1_process[sample %in% c("500b","1000b"), ci_boot := TRUE]
results_simulation1_process[sample == "500b", sample := "500"]
results_simulation1_process[sample == "1000b", sample := "1000"]

saveRDS(results_simulation1_process, file = "results/sim1-paper-results.RDS")

# processing results without V2 -------------------------------------------

setDT(results_simulation1_no_v2)

results_simulation1_no_v2_process <- results_simulation1_no_v2 |> melt(id.vars = 1:3)
results_simulation1_no_v2_process[, c("est", "sample", "type"):=tstrsplit(variable, "_")]
results_simulation1_no_v2_process[type == "ci", ci :="ci"]

saveRDS(results_simulation1_no_v2_process, file = "results/sim1-paper-results-no-v2.RDS")

