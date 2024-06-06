## packages
library(nonprobsvy)
library(sampling)
library(doSNOW)
library(progress)
library(foreach)
library(data.table)
library(xtable)

## seed
seed_for_sim <- 2024
cores <- 8

start_time <- Sys.time()

## run all
sims <- 500
source("codes/sim-main-paper.R", echo=TRUE)
source("codes/sim-appen-1-choose-k.R", echo=TRUE)
source("codes/sim-appen-2-varsel.R", echo=TRUE)
source("codes/sim-appen-4-positivity.R", echo=TRUE)
source("codes/sim-appen-5-multirobust.R", echo=TRUE)

sims <- 700 # about 30 percent have missing data
source("codes/sim-appen-3-nonparam.R", echo=TRUE)

end_time <- Sys.time() - start_time
