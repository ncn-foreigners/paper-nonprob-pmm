# Repository with codes for the paper "Data integration of non-probability and probability samples with predictive mean matching"

## Basic info

Repository for the paper:

```         
Chlebicki P., Chrostowski, Ł. & Beręsewicz, M., "Data integration of non-probability and probability samples with predictive mean matching"
```

## Structure of the repo

+ `codes/`
   + `sim-main-paper.R`
   + `sim-appen-1-choose-k.R`
   + `sim-appen-2-varsel.R`
   + `sim-appen-3-nonparam.R`
   + `sim-appen-4-positivity.R`
   + `sim-appen-5-multirobust.R`
+ `results/`
   + `sim1-paper-results-no-v2.RDS`
   + `sim1-paper-results.RDS`
   + `sim-appen1-choose-k-results.RDS`
   + `sim-appen2-varsel-results.RDS`
   + `sim-appen3-nonparam.RDS`
   + `sim-appen4-positivity.RDS`
   + `sim-appen5-multirobust.RDS`

## Basic codes

For this paper we developed two predictive mean matching (PMM) estimators and implemented them into the [`nonprobsvy`](https://github.com/ncn-foreigners/nonprobsvy) package. Below you can find functions to run the

-   PMM estimator $\hat{y}-y$ matching

```r
PMM_A <- nonprob(outcome = y1 ~ x1 + x2,
                 data = sample_nonprob,
                 svydesign = sample_prob,
                 method_outcome = "pmm",
                 family_outcome = "gaussian",
                 control_outcome = controlOut(k = 5, predictive_match = 1),
                 control_inference = controlInf(pmm_exact_se = TRUE))
```

-   PMM estimator $\hat{y}-\hat{y}$ matching

```r
PMM_A <- nonprob(outcome = y1 ~ x1 + x2,
                 data = sample_nonprob,
                 svydesign = sample_prob,
                 method_outcome = "pmm",
                 family_outcome = "gaussian",
                 control_outcome = controlOut(k = 5, predictive_match = 2),
                 control_inference = controlInf(pmm_exact_se = TRUE))
```

## Funding

Work on this package is supported by the National Science Centre, OPUS 22 grant no. 2020/39/B/HS4/00941.
