## reportin
library(data.table)
library(xtable)

# simulation main paper ---------------------------------------------------

results_simulation1_process <- readRDS(file = "results/sim1-paper-results.RDS")
results_simulation1_no_v2_process <- readRDS(file = "results/sim1-paper-results-no-v2.RDS")

tab1_ci <- results_simulation1_process[!is.na(ci) & is.na(ci_boot), .(ci=mean(value)), 
                                       .(est, sample, y)] |>
  transform(est = ifelse(grepl("all", y), est, est)) |>
  transform(est = ifelse(grepl("mis", y), paste0(est, "_mis"), est),
            y = gsub("_(mis|all)", "", y)) |>
  transform(ci = ci*100) |>
  dcast(... ~ y, value.var = "ci") 


## stats
tab1 <- results_simulation1_process[is.na(ci), .(bias=mean(value)-mean(trues), se = sd(value), 
                                                 rmse = sqrt((mean(value)-mean(trues))^2 + var(value))), 
                                    .(est, sample, y)] |>
  transform(est = ifelse(grepl("all", y), est, est)) |>
  transform(est = ifelse(grepl("mis", y), paste0(est, "_mis"), est),
            y = gsub("_(mis|all)", "", y)) |>
  melt(id.vars = c(1,2,3)) |>
  transform(y=paste(y, variable, sep = "_")) |>
  transform(variable=NULL,
            value = value*100) |>
  dcast(... ~ y, value.var = "value")

tab1[tab1_ci, on = c("est","sample"), ":="(y1_ci=i.y1,y2_ci=i.y2,y3_ci=i.y3)]

setcolorder(tab1, c("sample", "est", 
                    "y1_bias", "y1_se", "y1_rmse", "y1_ci",
                    "y2_bias", "y2_se", "y2_rmse", "y2_ci",
                    "y3_bias", "y3_se", "y3_rmse", "y3_ci"))


tab1[order(sample, est),][, ":="(sample=NULL)] |>
  xtable() |>
  print.xtable(include.rownames = FALSE)


## table 2 data ------------------------------------------------------------

tab1_ci_no_v2 <- results_simulation1_no_v2_process[!is.na(ci), .(ci=mean(value)), 
                                                   .(est, sample, y)] |>
  transform(est = ifelse(grepl("all", y), est, est)) |>
  transform(est = ifelse(grepl("mis", y), paste0(est, "_mis"), est),
            y = gsub("_(mis|all)", "", y)) |>
  transform(ci = ci*100) |>
  dcast(... ~ y, value.var = "ci") 

setnames(tab1_ci_no_v2, names(tab1_ci_no_v2), c("est", "sample", "y1_no_v2ci", "y2_no_v2ci", "y3_no_v2ci"))

tab1_ci[tab1_ci_no_v2][!grepl("glm|mis", est)] |>
  setcolorder(x= _, c("est", "sample", "y1_no_v2ci",  "y1",  "y2_no_v2ci", "y2",  "y3_no_v2ci", "y3")) |> 
  {\(x) x[order(-sample, est)][, ":="(sample=NULL)] }()  |>
  xtable() |>
  print.xtable(include.rownames = FALSE)


## table 3 data ------------------------------------------------------------

tab1_ci <- results_simulation1_process[!is.na(ci) & is.na(ci_boot), .(ci=mean(value)), 
                                       .(est, sample, y)] |>
  transform(est = ifelse(grepl("all", y), est, est)) |>
  transform(est = ifelse(grepl("mis", y), paste0(est, "_mis"), est),
            y = gsub("_(mis|all)", "", y)) |>
  transform(ci = ci*100) |>
  dcast(... ~ y, value.var = "ci") 

tab1_ci_b <- results_simulation1_process[!is.na(ci) & !is.na(ci_boot), .(ci=mean(value)), 
                                         .(est, sample, y)] |>
  transform(est = ifelse(grepl("all", y), est, est)) |>
  transform(est = ifelse(grepl("mis", y), paste0(est, "_mis"), est),
            y = gsub("_(mis|all)", "", y)) |>
  transform(ci = ci*100,
            y=paste0(y,"_b")) |>
  dcast(... ~ y, value.var = "ci") 

tab1_ci[tab1_ci_b, on = c("est", "sample"), ":="(y1_ci_b=y1_b,y2_ci_b=y2_b,y3_ci_b=y3_b)][!grepl("mis|glm", est)] |>
  setcolorder(c("est", "sample", "y1", "y1_ci_b", "y2", "y2_ci_b", "y3", "y3_ci_b")) |>
  {\(x) x[order(-sample, est)][, ":="(sample=NULL)] }()  |>
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
