## generate data

set.seed(seed_for_sim)

N <- 10000
n_A <- 500
p <- 50
KK <- 5
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

## generate Y
Y_11 <- as.numeric(X %*% beta_vec) +  rnorm(N) ## OM I: linear model
Y_12 <- 1 + exp(3*sin(as.numeric(X %*% beta_vec))) + X[, "X5"] + X[, "X6"] + rnorm(N) ## OM II: nonlinear model
pi_Y_21 <- plogis(as.numeric(X %*% beta_vec)) ## OM III: linear model for binary Y
pi_Y_22 <- plogis(2 - log((X %*% beta_vec)^2) - 2*X[,"X5"] + 2*X[, "X6"]) ## OM IV: nonlinear model for binary Y
Y_21 <- rbinom(N, 1, prob = pi_Y_21)
Y_22 <- rbinom(N, 1, prob = pi_Y_22)

## generate probs
pi_A_Y11 <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_11), n_A)
pi_A_Y12 <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_12), n_A)
pi_A_Y21 <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_21), n_A)
pi_A_Y22 <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y_22), n_A)

pi_B1 <- plogis(as.numeric(X %*% alpha_vec1)) ## PSM I: linear probability
pi_B2 <- plogis(3.5 + as.numeric(log(X^2) %*% alpha_vec2) - sin(X[, "X3"] + X[, "X4"]) - X[,"X5"] + X[, "X6"]) ## PSM II: nonlinear 


## generate data
population <- data.frame(pi_A_Y11, pi_A_Y12, pi_A_Y21, pi_A_Y22, Y_11, Y_12, Y_21, Y_22, X[, 2:p])

cl <- makeCluster(cores)
clusterExport(cl, c("N"))

registerDoSNOW(cl)

pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)

opts <- list(progress = \(n) pb$tick())

res <- foreach(k=1:sims, .combine = rbind,
               .packages = c("survey", "nonprobsvy", "sampling"),
               .errorhandling = "stop",
               .options.snow = opts) %dopar% {
  ## generate samples
  flag_B1 <- rbinom(N, 1, prob = pi_B1)
  flag_B2 <- rbinom(N, 1, prob = pi_B2)
  
  flag_A_Y11 <- UPpoisson(pik = pi_A_Y11)
  flag_A_Y12 <- UPpoisson(pik = pi_A_Y12)
  flag_A_Y21 <- UPpoisson(pik = pi_A_Y21)
  flag_A_Y22 <- UPpoisson(pik = pi_A_Y22)
  
  sample_A_svy_Y11 <- svydesign(ids = ~ 1, probs = ~ pi_A_Y11, pps = "brewer", data = population[flag_A_Y11 == 1, ])
  sample_A_svy_Y12 <- svydesign(ids = ~ 1, probs = ~ pi_A_Y12, pps = "brewer", data = population[flag_A_Y12 == 1, ])
  sample_A_svy_Y21 <- svydesign(ids = ~ 1, probs = ~ pi_A_Y21, pps = "brewer", data = population[flag_A_Y21 == 1, ])
  sample_A_svy_Y22 <- svydesign(ids = ~ 1, probs = ~ pi_A_Y22, pps = "brewer", data = population[flag_A_Y22 == 1, ])
  
  sample_B1 <- population[flag_B1 == 1, ]
  sample_B2 <- population[flag_B2 == 1, ]
  
  
  ## sample B1 ---------------------------------------------------------------
  ### Y11 -----------------------------------------------------------------
  est_mi_glm_b1_y11_no <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y11,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )

  est_mi_glm_b1_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y11,
    method_outcome = "glm",
    family_outcome = "gaussian",
    control_inference = controlInf(vars_selection = TRUE)
  )


  est_mi_pmmA_b1_y11_no <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )

  est_mi_pmmA_b1_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_pmmB_b1_y11_no <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmmB_b1_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_nn_b1_y11_no <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "nn"
  )
  
  est_mi_nn_b1_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK),
    family_outcome = "gaussian",
    method_outcome = "nn",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  ### Y12 -----------------------------------------------------------------
  
  est_mi_glm_b1_y12_no <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y12,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )

  est_mi_glm_b1_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y12,
    method_outcome = "glm",
    family_outcome = "gaussian",
    control_inference = controlInf(vars_selection = TRUE)
  )

  est_mi_pmmA_b1_y12_no <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmmA_b1_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_pmmB_b1_y12_no <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )

  est_mi_pmmB_b1_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )

  est_mi_nn_b1_y12_no <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "nn"
  )
  
  est_mi_nn_b1_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK),
    family_outcome = "gaussian",
    method_outcome = "nn",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  ### Y21 -----------------------------------------------------------------
  est_mi_glm_b1_y21_no <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y21,
    method_outcome = "glm",
    family_outcome = "binomial"
  )

  est_mi_glm_b1_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y21,
    method_outcome = "glm",
    family_outcome = "binomial",
    control_inference = controlInf(vars_selection = TRUE)
  )

  est_mi_pmmA_b1_y21_no <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )

  est_mi_pmmA_b1_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )

  est_mi_pmmB_b1_y21_no <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  est_mi_pmmB_b1_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_nn_b1_y21_no <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "nn"
  )
  
  est_mi_nn_b1_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK),
    family_outcome = "binomial",
    method_outcome = "nn",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  ### Y22 -----------------------------------------------------------------
  est_mi_glm_b1_y22_no <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y22,
    method_outcome = "glm",
    family_outcome = "binomial"
  )

  est_mi_glm_b1_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y22,
    method_outcome = "glm",
    family_outcome = "binomial",
    control_inference = controlInf(vars_selection = TRUE)
  )

  est_mi_pmmA_b1_y22_no <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  est_mi_pmmA_b1_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_pmmB_b1_y22_no <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )

  est_mi_pmmB_b1_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE)
  )


  est_mi_nn_b1_y22_no <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "nn"
  )
  
  est_mi_nn_b1_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B1,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK),
    family_outcome = "binomial",
    method_outcome = "nn",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  ## sample B2 ---------------------------------------------------------------
  ### Y11 -----------------------------------------------------------------
  est_mi_glm_b2_y11_no <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y11,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
  
  est_mi_glm_b2_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y11,
    method_outcome = "glm",
    family_outcome = "gaussian",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  
  est_mi_pmmA_b2_y11_no <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmmA_b2_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_pmmB_b2_y11_no <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmmB_b2_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_nn_b2_y11_no <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "nn"
  )
  
  est_mi_nn_b2_y11_sel <- nonprob(
    outcome = as.formula(paste("Y_11 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y11,
    control_outcome = controlOut(k = KK),
    family_outcome = "gaussian",
    method_outcome = "nn",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  ### Y12 -----------------------------------------------------------------
  
  est_mi_glm_b2_y12_no <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y12,
    method_outcome = "glm",
    family_outcome = "gaussian"
  )
  
  est_mi_glm_b2_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y12,
    method_outcome = "glm",
    family_outcome = "gaussian",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmmA_b2_y12_no <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmmA_b2_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_pmmB_b2_y12_no <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "pmm"
  )
  
  est_mi_pmmB_b2_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "gaussian",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_nn_b2_y12_no <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "gaussian",
    method_outcome = "nn"
  )
  
  est_mi_nn_b2_y12_sel <- nonprob(
    outcome = as.formula(paste("Y_12 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y12,
    control_outcome = controlOut(k = KK),
    family_outcome = "gaussian",
    method_outcome = "nn",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  ### Y21 -----------------------------------------------------------------
  est_mi_glm_b2_y21_no <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y21,
    method_outcome = "glm",
    family_outcome = "binomial"
  )
  
  est_mi_glm_b2_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y21,
    method_outcome = "glm",
    family_outcome = "binomial",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmmA_b2_y21_no <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  est_mi_pmmA_b2_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_pmmB_b2_y21_no <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  est_mi_pmmB_b2_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_nn_b2_y21_no <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "nn"
  )
  
  est_mi_nn_b2_y21_sel <- nonprob(
    outcome = as.formula(paste("Y_21 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y21,
    control_outcome = controlOut(k = KK),
    family_outcome = "binomial",
    method_outcome = "nn",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  ### Y22 -----------------------------------------------------------------
  est_mi_glm_b2_y22_no <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y22,
    method_outcome = "glm",
    family_outcome = "binomial"
  )
  
  est_mi_glm_b2_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y22,
    method_outcome = "glm",
    family_outcome = "binomial",
    control_inference = controlInf(vars_selection = TRUE)
  )
  
  est_mi_pmmA_b2_y22_no <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  est_mi_pmmA_b2_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK, predictive_match = 2),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  est_mi_pmmB_b2_y22_no <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "pmm"
  )
  
  est_mi_pmmB_b2_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK, predictive_match = 1),
    family_outcome = "binomial",
    method_outcome = "pmm",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  
  est_mi_nn_b2_y22_no <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK),
    control_inference = controlInf(pmm_exact_se = TRUE),
    family_outcome = "binomial",
    method_outcome = "nn"
  )
  
  est_mi_nn_b2_y22_sel <- nonprob(
    outcome = as.formula(paste("Y_22 ~ ", as.character(X_formula)[2])),
    data = sample_B2,
    svydesign = sample_A_svy_Y22,
    control_outcome = controlOut(k = KK),
    family_outcome = "binomial",
    method_outcome = "nn",
    control_inference = controlInf(vars_selection = TRUE, pmm_exact_se = TRUE)
  )
  
  data.frame(
    k = k,
    sample = rep(c("b1", "b2"), each = 4),
    y = rep(c("y11", "y12", "y21", "y22"), times = 2),
    trues = rep(
      c(mean(population$Y_11), mean(population$Y_12), 
        mean(population$Y_21), mean(population$Y_22)), 2
    ),
    glm_no = c(
      est_mi_glm_b1_y11_no$output$mean, est_mi_glm_b1_y12_no$output$mean,
      est_mi_glm_b1_y21_no$output$mean, est_mi_glm_b1_y22_no$output$mean,
      est_mi_glm_b2_y11_no$output$mean, est_mi_glm_b2_y12_no$output$mean,
      est_mi_glm_b2_y21_no$output$mean, est_mi_glm_b2_y22_no$output$mean
    ),
    pmmA_no = c(
      est_mi_pmmA_b1_y11_no$output$mean, est_mi_pmmA_b1_y12_no$output$mean,
      est_mi_pmmA_b1_y21_no$output$mean, est_mi_pmmA_b1_y22_no$output$mean,
      est_mi_pmmA_b2_y11_no$output$mean, est_mi_pmmA_b2_y12_no$output$mean,
      est_mi_pmmA_b2_y21_no$output$mean, est_mi_pmmA_b2_y22_no$output$mean
    ),
    pmmB_no = c(
      est_mi_pmmB_b1_y11_no$output$mean, est_mi_pmmB_b1_y12_no$output$mean,
      est_mi_pmmB_b1_y21_no$output$mean, est_mi_pmmB_b1_y22_no$output$mean,
      est_mi_pmmB_b2_y11_no$output$mean, est_mi_pmmB_b2_y12_no$output$mean,
      est_mi_pmmB_b2_y21_no$output$mean, est_mi_pmmB_b2_y22_no$output$mean
    ),
    nn_no = c(
      est_mi_nn_b1_y11_no$output$mean, est_mi_nn_b1_y12_no$output$mean,
      est_mi_nn_b1_y21_no$output$mean, est_mi_nn_b1_y22_no$output$mean,
      est_mi_nn_b2_y11_no$output$mean, est_mi_nn_b2_y12_no$output$mean,
      est_mi_nn_b2_y21_no$output$mean, est_mi_nn_b2_y22_no$output$mean
    ),
    glm_sel = c(
      est_mi_glm_b1_y11_sel$output$mean, est_mi_glm_b1_y12_sel$output$mean,
      est_mi_glm_b1_y21_sel$output$mean, est_mi_glm_b1_y22_sel$output$mean,
      est_mi_glm_b2_y11_sel$output$mean, est_mi_glm_b2_y12_sel$output$mean,
      est_mi_glm_b2_y21_sel$output$mean, est_mi_glm_b2_y22_sel$output$mean
    ),
    pmmA_sel = c(
      est_mi_pmmA_b1_y11_sel$output$mean, est_mi_pmmA_b1_y12_sel$output$mean,
      est_mi_pmmA_b1_y21_sel$output$mean, est_mi_pmmA_b1_y22_sel$output$mean,
      est_mi_pmmA_b2_y11_sel$output$mean, est_mi_pmmA_b2_y12_sel$output$mean,
      est_mi_pmmA_b2_y21_sel$output$mean, est_mi_pmmA_b2_y22_sel$output$mean
    ),
    pmmB_sel = c(
      est_mi_pmmB_b1_y11_sel$output$mean, est_mi_pmmB_b1_y12_sel$output$mean,
      est_mi_pmmB_b1_y21_sel$output$mean, est_mi_pmmB_b1_y22_sel$output$mean,
      est_mi_pmmB_b2_y11_sel$output$mean, est_mi_pmmB_b2_y12_sel$output$mean,
      est_mi_pmmB_b2_y21_sel$output$mean, est_mi_pmmB_b2_y22_sel$output$mean
    ),
    nn_sel = c(
      est_mi_nn_b1_y11_sel$output$mean, est_mi_nn_b1_y12_sel$output$mean,
      est_mi_nn_b1_y21_sel$output$mean, est_mi_nn_b1_y22_sel$output$mean,
      est_mi_nn_b2_y11_sel$output$mean, est_mi_nn_b2_y12_sel$output$mean,
      est_mi_nn_b2_y21_sel$output$mean, est_mi_nn_b2_y22_sel$output$mean
    ),
    glm_no_ci = c(
      est_mi_glm_b1_y11_no$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_glm_b1_y11_no$confidence_interval[, 2],
      est_mi_glm_b1_y12_no$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_glm_b1_y12_no$confidence_interval[, 2],
      est_mi_glm_b1_y21_no$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_glm_b1_y21_no$confidence_interval[, 2],
      est_mi_glm_b1_y22_no$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_glm_b1_y22_no$confidence_interval[, 2],
      est_mi_glm_b2_y11_no$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_glm_b2_y11_no$confidence_interval[, 2],
      est_mi_glm_b2_y12_no$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_glm_b2_y12_no$confidence_interval[, 2],
      est_mi_glm_b2_y21_no$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_glm_b2_y21_no$confidence_interval[, 2],
      est_mi_glm_b2_y22_no$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_glm_b2_y22_no$confidence_interval[, 2]
    ),
    pmmA_no_ci = c(
      est_mi_pmmA_b1_y11_no$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_pmmA_b1_y11_no$confidence_interval[, 2],
      est_mi_pmmA_b1_y12_no$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_pmmA_b1_y12_no$confidence_interval[, 2],
      est_mi_pmmA_b1_y21_no$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_pmmA_b1_y21_no$confidence_interval[, 2],
      est_mi_pmmA_b1_y22_no$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_pmmA_b1_y22_no$confidence_interval[, 2],
      est_mi_pmmA_b2_y11_no$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_pmmA_b2_y11_no$confidence_interval[, 2],
      est_mi_pmmA_b2_y12_no$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_pmmA_b2_y12_no$confidence_interval[, 2],
      est_mi_pmmA_b2_y21_no$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_pmmA_b2_y21_no$confidence_interval[, 2],
      est_mi_pmmA_b2_y22_no$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_pmmA_b2_y22_no$confidence_interval[, 2]
    ),
    pmmB_no_ci = c(
      est_mi_pmmB_b1_y11_no$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_pmmB_b1_y11_no$confidence_interval[, 2],
      est_mi_pmmB_b1_y12_no$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_pmmB_b1_y12_no$confidence_interval[, 2],
      est_mi_pmmB_b1_y21_no$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_pmmB_b1_y21_no$confidence_interval[, 2],
      est_mi_pmmB_b1_y22_no$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_pmmB_b1_y22_no$confidence_interval[, 2],
      est_mi_pmmB_b2_y11_no$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_pmmB_b2_y11_no$confidence_interval[, 2],
      est_mi_pmmB_b2_y12_no$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_pmmB_b2_y12_no$confidence_interval[, 2],
      est_mi_pmmB_b2_y21_no$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_pmmB_b2_y21_no$confidence_interval[, 2],
      est_mi_pmmB_b2_y22_no$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_pmmB_b2_y22_no$confidence_interval[, 2]
    ),
    nn_no_ci = c(
      est_mi_nn_b1_y11_no$confidence_interval[, 1] < mean(population$Y_11) & 
        mean(population$Y_11) < est_mi_nn_b1_y11_no$confidence_interval[, 2],
      est_mi_nn_b1_y12_no$confidence_interval[, 1] < mean(population$Y_12) & 
        mean(population$Y_12) < est_mi_nn_b1_y12_no$confidence_interval[, 2],
      est_mi_nn_b1_y21_no$confidence_interval[, 1] < mean(population$Y_21) & 
        mean(population$Y_21) < est_mi_nn_b1_y21_no$confidence_interval[, 2],
      est_mi_nn_b1_y22_no$confidence_interval[, 1] < mean(population$Y_22) & 
        mean(population$Y_22) < est_mi_nn_b1_y22_no$confidence_interval[, 2],
      est_mi_nn_b2_y11_no$confidence_interval[, 1] < mean(population$Y_11) & 
        mean(population$Y_11) < est_mi_nn_b2_y11_no$confidence_interval[, 2],
      est_mi_nn_b2_y12_no$confidence_interval[, 1] < mean(population$Y_12) & 
        mean(population$Y_12) < est_mi_nn_b2_y12_no$confidence_interval[, 2],
      est_mi_nn_b2_y21_no$confidence_interval[, 1] < mean(population$Y_21) & 
        mean(population$Y_21) < est_mi_nn_b2_y21_no$confidence_interval[, 2],
      est_mi_nn_b2_y22_no$confidence_interval[, 1] < mean(population$Y_22) & 
        mean(population$Y_22) < est_mi_nn_b2_y22_no$confidence_interval[, 2]
    ),
    glm_sel_ci = c(
      est_mi_glm_b1_y11_sel$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_glm_b1_y11_sel$confidence_interval[, 2],
      est_mi_glm_b1_y12_sel$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_glm_b1_y12_sel$confidence_interval[, 2],
      est_mi_glm_b1_y21_sel$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_glm_b1_y21_sel$confidence_interval[, 2],
      est_mi_glm_b1_y22_sel$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_glm_b1_y22_sel$confidence_interval[, 2],
      est_mi_glm_b2_y11_sel$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_glm_b2_y11_sel$confidence_interval[, 2],
      est_mi_glm_b2_y12_sel$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_glm_b2_y12_sel$confidence_interval[, 2],
      est_mi_glm_b2_y21_sel$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_glm_b2_y21_sel$confidence_interval[, 2],
      est_mi_glm_b2_y22_sel$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_glm_b2_y22_sel$confidence_interval[, 2]
    ),
    pmmA_sel_ci = c(
      est_mi_pmmA_b1_y11_sel$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_pmmA_b1_y11_sel$confidence_interval[, 2],
      est_mi_pmmA_b1_y12_sel$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_pmmA_b1_y12_sel$confidence_interval[, 2],
      est_mi_pmmA_b1_y21_sel$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_pmmA_b1_y21_sel$confidence_interval[, 2],
      est_mi_pmmA_b1_y22_sel$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_pmmA_b1_y22_sel$confidence_interval[, 2],
      est_mi_pmmA_b2_y11_sel$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_pmmA_b2_y11_sel$confidence_interval[, 2],
      est_mi_pmmA_b2_y12_sel$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_pmmA_b2_y12_sel$confidence_interval[, 2],
      est_mi_pmmA_b2_y21_sel$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_pmmA_b2_y21_sel$confidence_interval[, 2],
      est_mi_pmmA_b2_y22_sel$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_pmmA_b2_y22_sel$confidence_interval[, 2]
    ),
    
    pmmB_sel_ci = c(
      est_mi_pmmB_b1_y11_sel$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_pmmB_b1_y11_sel$confidence_interval[, 2],
      est_mi_pmmB_b1_y12_sel$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_pmmB_b1_y12_sel$confidence_interval[, 2],
      est_mi_pmmB_b1_y21_sel$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_pmmB_b1_y21_sel$confidence_interval[, 2],
      est_mi_pmmB_b1_y22_sel$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_pmmB_b1_y22_sel$confidence_interval[, 2],
      est_mi_pmmB_b2_y11_sel$confidence_interval[, 1] < mean(population$Y_11) &
      mean(population$Y_11) < est_mi_pmmB_b2_y11_sel$confidence_interval[, 2],
      est_mi_pmmB_b2_y12_sel$confidence_interval[, 1] < mean(population$Y_12) &
      mean(population$Y_12) < est_mi_pmmB_b2_y12_sel$confidence_interval[, 2],
      est_mi_pmmB_b2_y21_sel$confidence_interval[, 1] < mean(population$Y_21) &
      mean(population$Y_21) < est_mi_pmmB_b2_y21_sel$confidence_interval[, 2],
      est_mi_pmmB_b2_y22_sel$confidence_interval[, 1] < mean(population$Y_22) &
      mean(population$Y_22) < est_mi_pmmB_b2_y22_sel$confidence_interval[, 2]
    ),
    nn_sel_ci = c(
      est_mi_nn_b1_y11_sel$confidence_interval[, 1] < mean(population$Y_11) & 
        mean(population$Y_11) < est_mi_nn_b1_y11_sel$confidence_interval[, 2],
      est_mi_nn_b1_y12_sel$confidence_interval[, 1] < mean(population$Y_12) & 
        mean(population$Y_12) < est_mi_nn_b1_y12_sel$confidence_interval[, 2],
      est_mi_nn_b1_y21_sel$confidence_interval[, 1] < mean(population$Y_21) & 
        mean(population$Y_21) < est_mi_nn_b1_y21_sel$confidence_interval[, 2],
      est_mi_nn_b1_y22_sel$confidence_interval[, 1] < mean(population$Y_22) & 
        mean(population$Y_22) < est_mi_nn_b1_y22_sel$confidence_interval[, 2],
      est_mi_nn_b2_y11_sel$confidence_interval[, 1] < mean(population$Y_11) & 
        mean(population$Y_11) < est_mi_nn_b2_y11_sel$confidence_interval[, 2],
      est_mi_nn_b2_y12_sel$confidence_interval[, 1] < mean(population$Y_12) & 
        mean(population$Y_12) < est_mi_nn_b2_y12_sel$confidence_interval[, 2],
      est_mi_nn_b2_y21_sel$confidence_interval[, 1] < mean(population$Y_21) & 
        mean(population$Y_21) < est_mi_nn_b2_y21_sel$confidence_interval[, 2],
      est_mi_nn_b2_y22_sel$confidence_interval[, 1] < mean(population$Y_22) & 
        mean(population$Y_22) < est_mi_nn_b2_y22_sel$confidence_interval[, 2]
    )
  )
}

stopCluster(cl)


## processing results
setDT(res)

results_simulation1_process <- res |> melt(id.vars = 1:4)
results_simulation1_process[, c("est", "var", "ci"):=tstrsplit(variable, "_")]

saveRDS(results_simulation1_process, file = "results/sim-appen2-varsel-results.RDS")
