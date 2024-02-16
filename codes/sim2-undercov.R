#remotes::install_github("ncn-foreigners/nonprobsvy@gee_momentfit")

library(nonprobsvy)
library(doSNOW)
library(progress)
library(foreach)
library(xtable)
library(data.table)
library(sampling)
library(jointCalib)

syspps <- function(x,n){ 
  N <- length(x)
  U <- sample(N,N)
  xx <- x[U]
  z <- n*cumsum(xx)/sum(xx)
  r <- runif(1)
  s <- numeric(n)
  k <- 1
  for(i in 1:N){
    if(z[i]>=r){
      s[k] <- U[i]
      r <- r+1
      k <- k+1
    }
  }
  return(s)
}


## sim parameteres
seed_for_sim <- 123
set.seed(seed_for_sim)

eta <- 1
sigma_sim <- 0.50 ## R2 is 0.25
tau <- 0.4

## generate data
N<-20000
n_a <- 1000 ## nonprob svy
n_b <- 500 ## prob svy
x1 <- rnorm(N)
x2 <- rexp(N)
x3 <- rbinom(N,1,0.5)
epsilon <- rnorm(N)
sigma05 <- stats::uniroot(f = function(ss) 
  cor(3+x1+x2+x3-eta*x1^2, 3+x1+x2+x3-eta*x1^2+ss*epsilon) - sigma_sim, c(0,20))$root

y <- 3+x1+x2+x3-eta*x1^2+sigma05*epsilon

## coverage function
p_coverage <- plogis(1-0.6*x1 + 0.5*x2+0.8*x3)
tau_2 <- quantile(p_coverage, tau) ## coverage error

## inclusion into probability sample
c_r <- stats::uniroot(f = function(s) max(s+x2)/min(s+x2) - 50, c(0,2), tol = 1e-10)$root

## dataset
pop_df <- data.frame(id=1:N, x1,x2,x3, y, flag=p_coverage >= tau_2, z=x2 + c_r, pi_b=n_b*(x2 +c_r)/sum(x2 +c_r))
pop_df_cov <- pop_df[pop_df$flag == TRUE,]

theta_r <- with(pop_df_cov, 
                stats::uniroot(f = function(ss) sum(plogis(ss + 0.3*x1-0.3*x2 +0.5*x3))-n_a, c(-10,0), tol = 1e-10)$root)

pop_df_cov$pi_a <- with(pop_df_cov, plogis(theta_r + 0.3*x1-0.3*x2 + 0.5*x3))

N_cov <- nrow(pop_df_cov)
## quantiles
q_probs1 <- c(0.25,0.5,0.75)
q_probs2 <- seq(0.1,0.9,0.1)
## sampling
sims <- 1000
cores <- 8
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())

results_simulation2 <- foreach(k=1:sims, .combine = rbind,
                               .packages = c("survey", "nonprobsvy", "jointCalib"),
                               .options.snow = opts, .errorhandling = "remove") %dopar% {
                                 
  sample_a <- pop_df_cov[which(sampling::UPpoisson(pop_df_cov$pi_a)==1),]  ## nonprob sample: 
  #sample_a <- pop_df_cov[sampling::UPpoisson(pop_df_cov$pi_b) == 1,]  ## nonprob sample - error in the in the code?
  sample_b <- pop_df[syspps(pop_df$z, n_b), ] ## prob sample
  svy_b <- svydesign(ids = ~1, data = sample_b, probs = ~pi_b, pps = HR())
  
  ## quantiles
  quants_est1 <- svyquantile(~x1 + x2, design = svy_b, quantiles = q_probs1)
  quants_est2 <- svyquantile(~x1 + x2, design = svy_b, quantiles = q_probs2)
  ### prob sample
  mat1 <- joint_calib(formula_quantiles = ~x1 + x2,
                      data = svy_b$variables, dweights = weights(svy_b), N = sum(weights(svy_b)),
                      pop_quantiles = quants_est1, method = "linear")
  mat2 <- joint_calib(formula_quantiles = ~x1 + x2,
                      data = svy_b$variables, dweights = weights(svy_b), N = sum(weights(svy_b)),
                      pop_quantiles = quants_est2, method = "linear")
  colnames(mat1$Xs) <- c("i", paste0("x1_", q_probs1), paste0("x2_", q_probs1))
  colnames(mat2$Xs) <- c("i", paste0("x1_", q_probs2), paste0("x2_", q_probs2))
  mat12 <- cbind(mat1$Xs[, -c(1, 3, 6)], mat2$Xs[, -1])
  svy_b$variables <- cbind(svy_b$variables,mat12)
  
  ### non-prob sample
  mat1a <- joint_calib(formula_quantiles = ~x1 + x2,
                      data = sample_a, N = sum(weights(svy_b)),
                      pop_quantiles = quants_est1, method = "linear")
  mat2a <- joint_calib(formula_quantiles = ~x1 + x2,
                      data = sample_a, N = sum(weights(svy_b)),
                      pop_quantiles = quants_est2, method = "linear")
  colnames(mat1a$Xs) <- c("i", paste0("x1_", q_probs1), paste0("x2_", q_probs1))
  colnames(mat2a$Xs) <- c("i", paste0("x1_", q_probs2), paste0("x2_", q_probs2))
  mat12a <- cbind(mat1a$Xs[, -c(1, 3, 6)], mat2a$Xs[, -1])
  sample_a <- cbind(sample_a,mat12a)
  
  ### 
  
  naive <- mean(sample_a$y)
 
  ## standard 
  miglm <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
                   se = F, family_outcome = "gaussian", method_outcome = "glm")
  ipw <- nonprob(selection =  ~ x1 + x2 + x3, target=~y,
                 data = sample_a, svydesign = svy_b, se = F)
  dr <- nonprob(selection =  ~ x1 + x2 + x3, 
                outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b, 
                family_outcome = "gaussian", method_outcome = "glm", se = F)
  
  ## with quartiles
  miglm_q <- nonprob(outcome = y ~ x3 + x1_0.5 + x2_0.5 ,
                     data = sample_a, svydesign = svy_b,
                   se = F, family_outcome = "gaussian", method_outcome = "glm")
  ipw_q <- nonprob(selection =  ~ x3 + x1_0.5 + x2_0.5, 
                   target=~y,
                 data = sample_a, svydesign = svy_b, se = F)
  dr_q <- nonprob(selection =  ~ x3 + x1_0.5 + x2_0.5, 
                outcome = y ~ x3 + x1_0.5 + x2_0.5, 
                data = sample_a, svydesign = svy_b, 
                family_outcome = "gaussian", method_outcome = "glm", se = F)
  ## with deciles
  miglm_d <- nonprob(outcome = y ~ #x1 + x2 + 
                       x3 + x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 + x1_0.7 + x1_0.8 + x1_0.9 + 
                       x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 + x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9, 
                     data = sample_a, svydesign = svy_b,
                     se = F, family_outcome = "gaussian", method_outcome = "glm")
  ipw_d <- nonprob(selection =  ~ #x1 + x2 + 
                     x3 + x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 + x1_0.7 + x1_0.8 + x1_0.9 + 
                     x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 + x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9, 
                   target=~y,
                   data = sample_a, svydesign = svy_b, se = F)
  dr_d <- nonprob(selection =  ~ #x1 + x2 + 
                    x3 + x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 + x1_0.7 + x1_0.8 + x1_0.9 + 
                    x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 + x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9, 
                  outcome = y ~ #x1 + x2 + 
                    x3 + x1_0.1 + x1_0.2 + x1_0.3 + x1_0.4 + x1_0.5 + x1_0.6 + x1_0.7 + x1_0.8 + x1_0.9 + 
                    x2_0.1 + x2_0.2 + x2_0.3 + x2_0.4 + x2_0.5 + x2_0.6 + x2_0.7 + x2_0.8 + x2_0.9, 
                  data = sample_a, svydesign = svy_b, 
                  family_outcome = "gaussian", method_outcome = "glm", se = F)
  
  data.frame(k=k, 
             prob=svymean(~y,svy_b)[1],
             naive=naive,  
             
             miglm=miglm$output$mean, 
             ipw = ipw$output$mean, 
             dr=dr$output$mean,
             
             miglm_q=miglm_q$output$mean, 
             ipw_q = ipw_q$output$mean, 
             dr_q=dr_q$output$mean,
             
             miglm_d=miglm_d$output$mean, 
             ipw_d = ipw_d$output$mean, 
             dr_d=dr_d$output$mean
             #hyb=y_hyb
             #pmm1a=pmm1a$output$mean, pmm1b=pmm1b$output$mean,
             #pmm5a=pmm5a$output$mean, pmm5b=pmm5b$output$mean
             )

}

stopCluster(cl)

#results_simulation2 <- subset(results_simulation2, dr_d > 0 & dr_d < 6)
sort(colMeans((results_simulation2[, -1]-mean(y))/mean(y))*100)
sort(apply(results_simulation2[, -1], 2, FUN = function(x) mean( (x - mean(y))^2))*100)
boxplot(results_simulation2[, -1]-mean(y))
abline(h=0,col='red')

