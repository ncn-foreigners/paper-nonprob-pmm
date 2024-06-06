## reportin


# simulation main paper ---------------------------------------------------

results_simulation1_process <- readRDS(file = "results/sim1-paper-results.RDS")
results_simulation1_no_v2_process <- readRDS(file = "results/sim1-paper-results-no-v2.RDS")

tab1_ci <- results_simulation1_process[!is.na(ci) & is.na(ci_boot), .(ci=mean(value)), 
                                       .(type, est, sample, y)] |>
  melt(id.vars = c(1, 4,2,3)) |>
  transform(y=paste(y, variable, sep = "_")) |>
  transform(variable=NULL,
            value = value*100,
            sample=as.character(as.numeric(sample)/100)) |>
  dcast(... ~ y, value.var = "value") 


## stats
tab1 <- results_simulation1_process[is.na(ci), .(bias=mean(value)-mean(trues), se = sd(value), 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))), 
                                    .(type, est, sample, y)] |>
  melt(id.vars = c(1, 4,2,3)) |>
  transform(y=paste(y, variable, sep = "_")) |>
  transform(variable=NULL,
            value = value*100) |>
  dcast(... ~ y, value.var = "value") |>
  transform(type=ifelse(is.na(type), "srs", type))

tab1[tab1_ci, on = c("type","est","sample"), ":="(y1_ci=i.y1_ci,y2_ci=i.y2_ci,y3_ci=i.y3_ci)]
tab1[, y1_ci:=round(y1_ci,1)]
tab1[, y2_ci:=round(y2_ci,1)]
tab1[, y3_ci:=round(y3_ci,1)]

tab1[, est:=factor(est, c("naive",  "glm", "nn1", "nn5", "pmm1a", "pmm1b", "pmm5a", "pmm5b"), 
                   c("Naive", "GLM", "NN1", "NN5", "PMM1A", "PMM1B", "PMM5A", "PMM5B")
                   , ordered = T)]

tab1[, sample:=factor(sample, c(5,10), ordered = T)]

setcolorder(tab1, c("type", "sample", "est", 
                    "y1_bias", "y1_se", "y1_rmse", "y1_ci",
                    "y2_bias", "y2_se", "y2_rmse", "y2_ci",
                    "y3_bias", "y3_se", "y3_rmse", "y3_ci"))


tab1[order(sample,-type, est),][, ":="(sample=NULL,type=NULL)] |>
  xtable() |>
  print.xtable(include.rownames = FALSE)


## table 2 data ------------------------------------------------------------

tab1_ci_no_v2 <- results_simulation1_no_v2_process[!is.na(ci), .(ci=mean(value)), 
                                                   .(type, est, sample, y)] |>
  melt(id.vars = c(1, 4,2,3)) |>
  transform(y=paste(y, variable, sep = "_no_v2")) |>
  transform(variable=NULL,
            value = value*100,
            sample=as.character(as.numeric(sample)/100)) |>
  dcast(... ~ y, value.var = "value")  |>
  subset(est != "glm") |>
  merge(tab1_ci, by = c("type", "est", "sample"))


setcolorder(tab1_ci_no_v2, c("type", "est", "sample",  "y1_no_v2ci", "y1_ci", "y2_no_v2ci", "y2_ci", "y3_no_v2ci", "y3_ci"))


tab1_ci_no_v2[order(-sample,-type, est),][, ":="(sample=NULL,type=NULL)] |>
  xtable() |>
  print.xtable(include.rownames = FALSE)


## table 3 data ------------------------------------------------------------


tab1_ci_b <- results_simulation1_process[!is.na(ci) & !is.na(ci_boot), .(ci=mean(value)), 
                                         .(type, est, sample, y)] |>
  melt(id.vars = c(1, 4,2,3)) |>
  transform(y=paste(y, variable, "b", sep = "_")) |>
  transform(variable=NULL,
            value = value*100,
            sample=as.character(as.numeric(sample)/100)) |>
  dcast(... ~ y, value.var = "value") 

tab3 <- tab1_ci_no_v2[, .(type, est, sample, y1_ci, y2_ci, y3_ci)]
tab3 <- tab3[tab1_ci_b, on = c("type", "est", "sample"), ":="(y1_ci_b=y1_ci_b,y2_ci_b=y2_ci_b,y3_ci_b=y3_ci_b)]

setcolorder(tab3, c("type", "est", "sample",  "y1_ci", "y1_ci_b", "y2_ci", "y2_ci_b", "y3_ci", "y3_ci_b"))



tab3[order(-sample,-type, est),][, ":="(sample=NULL,type=NULL)] |>
  xtable() |>
  print.xtable(include.rownames = FALSE)




# simulation appendix 1 ---------------------------------------------------

results_simulation1_process <- readRDS(file = "results/sim-appen1-choose-k-results.RDS")

tab1 <- results_simulation1_process[is.na(ci) & is.na(k_sel), .(bias=(mean(value)-mean(trues))*100, 
                                                                se = sd(value)*100, 
                                                                rmse = sqrt((mean(value)-mean(trues))^2 + var(value))*100), 
                                    keyby=.(y, est, sample)] 

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(y, est,sample)] 

tab3 <- results_simulation1_process[!is.na(k_sel), .(k = mean(value)), 
                                    keyby=.(y, est,sample)] 


merge(x = tab1[tab2, on = c("y", "est", "sample")],
      y = tab3, 
      by = c("y", "est", "sample"),
      all.x=T) |>
  xtable(digits = 4) |>
  print.xtable(include.rownames = FALSE)



# simulation appendix 2 ---------------------------------------------------

results_simulation1_process <- readRDS(file = "results/sim-appen2-varsel-results.RDS")


tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues))*100, 
                                                 se = sd(value)*100, 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))*100), 
                                    keyby=.(sample, y, est, var)] |>
  melt(id.vars = c(1, 4,2,3)) |>
  transform(sample=paste(sample, variable, sep = "_")) |>
  transform(variable=NULL) |>
  dcast(... ~ sample, value.var = "value") |>
  {\(x) x[order(y, est, var)]}() 

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(sample, y, est, var)]  |>
  dcast(... ~ sample, value.var = "ci")

tab1[tab2, on = c("y", "est", "var")] |>
  setcolorder(c("y", "est", "var", "b1_bias", "b1_se", "b1_rmse", "b1", "b2_bias", "b2_se", "b2_rmse", "b2")) |>
  {\(x) x[,y:=NULL]}() |>
  xtable(digits = 2) |>
  print.xtable(include.rownames = FALSE)


# simulation appendix 3 ---------------------------------------------------

results_simulation1_process <- readRDS(file = "results/sim-appen3-nonparam.RDS")

tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues))*100, 
                                                 se = sd(value)*100, 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))*100), 
                                    keyby=.(y, est, var)]

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(y, est, var)]


tab1[tab2, on = c("y", "est", "var")] |>
  xtable(digits = 2) |>
  print.xtable(include.rownames = F)



# simulation appendix 4 ---------------------------------------------------

results_simulation1_process <- readRDS(file = "results/sim-appen4-positivity.RDS")

tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues))*100, 
                                                 se = sd(value)*100, 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))*100), 
                                    keyby=.(y, type, est)]

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(y, type, est)]


tab1[tab2, on = c("y", "type", "est")] |>
  xtable(digits = 2) |>
  print.xtable(include.rownames = F)


# simulation appendix 5 ---------------------------------------------------

results_simulation1_process <- readRDS(file = "results/sim-appen5-multirobust.RDS")

tab1 <- results_simulation1_process[is.na(ci), .(bias=(mean(value)-mean(trues)), 
                                                 se = sd(value), 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))), 
                                    keyby=.(y, est, var)]

tab2 <- results_simulation1_process[!is.na(ci), .(ci = mean(value)*100), 
                                    keyby=.(y, est, var)]


tab1[tab2, on = c("y", "est", "var")] |>
  xtable(digits = 2) |>
  print.xtable(include.rownames = F)
