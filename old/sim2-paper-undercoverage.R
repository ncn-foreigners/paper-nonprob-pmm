#remotes::install_github("ncn-foreigners/nonprobsvy@gee_momentfit")
#remotes::install_github("ncn-foreigners/nonprobsvy@pmm")
library(nonprobsvy)
library(doSNOW)
library(progress)
library(foreach)
library(xtable)
library(data.table)
library(sampling)

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
  cor(3+x1+x2+x3-eta*x1^2, 3+x1+x2+x3-eta*x1^2+ss*epsilon) - sigma_sim, c(-20,0))$root

y <- 3+x1+x2+x3-eta*x1^2+sigma05*epsilon

## coverage function
p_coverage <- plogis(1-0.6*x1 + 0.5*x2+0.8*x3)
tau_2 <- quantile(p_coverage, tau) ## coverage error

## inclusion into probability sample
c_r <- stats::uniroot(f = function(s) max(s+x2)/min(s+x2) - 50, c(0,2), tol = 1e-10)$root

## dataset
pop_df <- data.frame(id=1:N, x1,x2,x3, y, flag=p_coverage >= tau_2, z=x2 + c_r, pi_b=n_b*(x2 +c_r)/sum(x2 +c_r))
aggregate(y~1, pop_df, FUN=mean)
aggregate(y~flag, pop_df, FUN=mean)

pop_df_cov <- pop_df[pop_df$flag == TRUE,]

theta_r <- with(pop_df_cov, 
                stats::uniroot(f = function(ss) sum(plogis(ss + 0.3*x1-0.3*x2 +0.5*x3))-n_a, c(-10,0), tol = 1e-10)$root)

pop_df_cov$pi_a <- with(pop_df_cov, plogis(theta_r + 0.3*x1-0.3*x2 + 0.5*x3))

#pop_df$pi_b <- inclusionprobabilities(pop_df$z, n)
N_cov <- nrow(pop_df_cov)
## sampling
sims <- 5000
cores <- 8
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed: :elapsedfull || Remaining: :eta]",
                       total = sims)
opts <- list(progress = \(n) pb$tick())

results_simulation2 <- foreach(k=1:sims, .combine = rbind,
                               .packages = c("survey", "nonprobsvy"),
                               .options.snow = opts, .errorhandling = 'remove') %dopar% {
                                 
  
  #sample_a <- pop_df_cov[which(sampling::UPpoisson(pop_df_cov$pi_a)==1),]  ## nonprob sample: 
  sample_a <- pop_df_cov[sampling::UPpoisson(pop_df_cov$pi_b) == 1,]  ## nonprob sample - error in the in the code?
  sample_b <- pop_df[syspps(pop_df$z, n_b), ] ## prob sample
  
  ## this is done without the convex hull method -- S_B0
  sample_b0 <- sample_b[sample_b$id %in% pop_df_cov$id, ] ## simple approach
  
  ## this is done without the convex hull method -- S_B1 subsample
  sample_b1 <- sample_b[!sample_b$id %in% pop_df_cov$id, ]
  sample_b1_ss <- sample_b1[sample(1:nrow(sample_b1), round(nrow(sample_b1)*0.2)), ]
  sample_b1_ss$pi_b2 <- nrow(sample_b1_ss)/nrow(sample_b1)
  
  svy_b <- svydesign(ids = ~1, data = sample_b, probs = ~pi_b, pps = HR())
  svy_b0 <- svydesign(ids = ~1, data = sample_b0, probs = ~pi_b, pps = HR())

  
  naive <- mean(sample_a$y)
  
  miglm <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
                   se = F, family_outcome = "gaussian", method_outcome = "glm")
  
  ipw <- nonprob(selection =  ~ x1 + x2 + x3, target=~y,
                 data = sample_a, svydesign = svy_b,
                 se = F)
  
  ipw_b0 <- nonprob(selection =  ~ x1 + x2 + x3, target=~y,
                    data = sample_a, svydesign = svy_b0,
                    se = F)
  
  dr <- nonprob(selection =  ~ x1 + x2 + x3, 
                outcome = y ~ x1 + x2 + x3,
                data = sample_a, svydesign = svy_b,
                family_outcome = "gaussian", method_outcome = "glm",
                se = F)
  
  N_b <- sum(weights(svy_b))
  y_hyb <- 1/N_b*sum(dr$outcome$y$residuals*ipw_b0$weights) + miglm$output$mean
  
  ## subsample
  df_combined <- rbind(sample_a[, -9], sample_b1_ss[, -9])
  m1 <- lm(y ~ x1 + x2 + x3, df_combined)
  m_tilde_a <- predict(m1, sample_a)
  m_tilde_b1 <- predict(m1, sample_b1_ss)
  svy_b <- update(svy_b, m_tilde = predict(m1, svy_b$variables))
  y_ss <- 1/N_b*( sum((sample_a$y - m_tilde_a)/ipw_b0$weights) + 
            sum(1/sample_b1_ss$pi_b*1/sample_b1_ss$pi_b2*(sample_b1_ss$y - m_tilde_b1)) + 
            svytotal(~m_tilde, svy_b)[1] )
  #
  
  # pmm1a <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
  #                 se = F, family_outcome = "gaussian", method_outcome = "pmm",
  #                 control_outcome = controlOut(k = 1, predictive_match = 2))
  # pmm1b <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
  #                 se = F, family_outcome = "gaussian", method_outcome = "pmm",
  #                 control_outcome = controlOut(k = 1, predictive_match = 1))
  # pmm5a <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
  #                  se = F, family_outcome = "gaussian", method_outcome = "pmm",
  #                  control_outcome = controlOut(k = 5, predictive_match = 2))
  # pmm5b <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
  #                  se = F, family_outcome = "gaussian", method_outcome = "pmm",
  #                  control_outcome = controlOut(k = 5, predictive_match = 1))
  
  
  data.frame(k=k, 
             prob=svymean(~y,svy_b)[1],
             naive=naive,  
             miglm=miglm$output$mean, 
             ipw = ipw$output$mean, dr=dr$output$mean,
             hyb=y_hyb, y_ss=y_ss
             #pmm1a=pmm1a$output$mean, pmm1b=pmm1b$output$mean,
             #pmm5a=pmm5a$output$mean, pmm5b=pmm5b$output$mean
             )

}

stopCluster(cl)

#results_simulation2 <- subset(results_simulation2, hyb > 0 & hyb < 10)
colMeans((results_simulation2[, -1]-mean(y))/mean(y))*100
apply(results_simulation2[, -1], 2, FUN = function(x) mean( (x - mean(y))^2))*100
boxplot(results_simulation2[, -1]-mean(y), ylim = c(-2,2))
abline(h=0,col='red')

