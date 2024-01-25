library(nonprobsvy)
library(doSNOW)
library(progress)
library(foreach)
library(xtable)
library(data.table)
library(sampling)


## what about probabilit
syspps <- function(x,n){ 
  N <- length(x)
  U <- sample(N,N)
  xx <- x[U]
  pik <- n*xx/sum(xx)
  z <- n*cumsum(xx)/sum(xx)
  r <- runif(1)
  s <- numeric()
  for(i in 1:N){
    if(z[i]>=r){
      s=c(s,U[i])
      r=r+1
    }
  }
  return(list(s=s[order(s)], pik = pik[s]))
}

## sim parameteres
seed_for_sim <- 2024
set.seed(seed_for_sim)

eta <- 1
sigma_sim <- 0.50 ## R2 is 0.25
tau <- 0.4

## generate data
N<-200000
n_a <- 1000 ## nonprob svy
n_b <- 500 ## prob svy
x1 <- rnorm(N)
x2 <- rexp(N)
x3 <- rbinom(N,1,0.5)
epsilon <- rnorm(N)
sigma05 <- stats::uniroot(f = function(ss) 
  cor(3+x1+x2+x3, 3+x1+x2+x3-eta*x1^2+ss*epsilon) - sigma_sim, c(0,20))$root

y <- 3+x1+x2+x3-eta*x1^2+sigma05*epsilon

## coverage function
p_coverage <- plogis(1-0.6*x1 + 0.5*x2+0.8*x3)
tau_2 <- quantile(p_coverage, tau) ## coverage error

## inclusion into probability sample
c_r <- stats::uniroot(f = function(s) max(s+x2)/min(s+x2) - 50, c(0,2), tol = etol)$root
pi_b <- inclusionprobabilities(x2 +c_r, n_b)

## dataset
pop_df <- data.frame(x1,x2,x3, y, flag=p_coverage >= tau_2, z=x2 +c_r, pi_b=pi_b)
pop_df_cov <- pop_df[pop_df$flag == TRUE,]

c(mean(y), mean(pop_df_cov$y), mean(pop_df$y[pop_df$flag == FALSE]))

theta_r <- with(pop_df_cov, 
                stats::uniroot(f = function(ss) sum(plogis(ss + 0.3*x1-0.3*x2 +0.5*x3))-n_a, c(-10,0))$root)

pop_df_cov$pi_a <- with(pop_df_cov, plogis(theta_r + 0.3*x1-0.3*x2 +0.5*x3))

#pop_df$pi_b <- inclusionprobabilities(pop_df$z, n)
N_cov <- nrow(pop_df_cov)
## sampling
sims <- 1000
cores <- 8
cl <- makeCluster(cores)
registerDoSNOW(cl)
pb <- progress_bar$new(total = sims)
opts <- list(progress = \(n) pb$tick())

results_simulation2 <- foreach(k=1:sims, .combine = rbind,
                               .packages = c("survey", "nonprobsvy"),
                               .options.snow = opts) %dopar% {
                                 
  
  sample_a <- pop_df_cov[which(sampling::UPpoisson(pop_df_cov$pi_a)==1),]  ## nonprob sample
  sample_b <- pop_df[which(sampling::UPrandomsystematic(pop_df$pi_b)==1), ] ## prob sample
  svy_b <- svydesign(ids = ~1, data = sample_b, probs = ~pi_b)
  
  naive <- mean(sample_a$y)
  
  miglm <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
                   se = F, family_outcome = "gaussian", method_outcome = "glm")
  pmm1a <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
                  se = F, family_outcome = "gaussian", method_outcome = "pmm",
                  control_outcome = controlOut(k = 1, predictive_match = 2))
  pmm1b <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
                  se = F, family_outcome = "gaussian", method_outcome = "pmm",
                  control_outcome = controlOut(k = 1, predictive_match = 1))
  pmm5a <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
                   se = F, family_outcome = "gaussian", method_outcome = "pmm",
                   control_outcome = controlOut(k = 5, predictive_match = 2))
  pmm5b <- nonprob(outcome = y ~ x1 + x2 + x3, data = sample_a, svydesign = svy_b,
                   se = F, family_outcome = "gaussian", method_outcome = "pmm",
                   control_outcome = controlOut(k = 5, predictive_match = 1))
  
  data.frame(k=k, prob=svymean(~y,svy_b)[1],naive=naive,  
             miglm=miglm$output$mean, 
             pmm1a=pmm1a$output$mean, pmm1b=pmm1b$output$mean,
             pmm5a=pmm5a$output$mean, pmm5b=pmm5b$output$mean)

}

stopCluster(cl)

colMeans((results_simulation2[, -1]-mean(y))/mean(y))*100
boxplot(results_simulation2[, -1]-mean(y))
abline(h=0,col='red')
