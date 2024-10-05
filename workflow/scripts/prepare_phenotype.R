library(data.table)
library(magrittr)

###########
## SETUP ##
###########
## read PAGE harmonized phenotype
dt <- fread("source_data/page-harmonizedpheno.topmedimputed.txt")
dt$t2d_status[dt$t2d_status > 0] <- 1 ## recode T2D

## read AA and HCHS/SOL samples
aa_samples <- readLines("source_data/AA_samples.txt")
hchs_sol_samples <- readLines("source_data/HCHS_SOL_samples.txt")

## subset PAGE phenotypes/covariates/samples
phenotypes <- c(
  "height", "bmi", "crp", "diastolic_bp", "pr_interval", "qrs_interval",
  "egfrckdepi_clean", "HbA1c", "hdl", "ldl", "mean_corp_hgb_conc",
  "platelet_cnt", "systolic_bp", "total_cholesterol", "triglycerides",
  "total_wbc_cnt", "chr_kidney_dz", "hypertension", "t2d_status"
)

covars_aa <- c("age", "sex", "study", paste0("ev", 1:19, "_tm"))
covars_hchs_sol <- c("age", "sex", paste0("ev", 1:19, "_tm"))
dt$age2 <- dt$age^2
covars_aa <- c(covars_aa, "age2")
covars_hchs_sol <- c(covars_hchs_sol, "age2")
dt_aa <- dt[page_subject_id %in% aa_samples, ]
dt_hchs_sol <- dt[page_subject_id %in% hchs_sol_samples, ]

#######################
## Adjusting Phenotypes ##
#######################
INT <- function(r) {
  qnorm((rank(r, na.last = "keep") - 0.5) / sum(!is.na(r)))
}

### aa
height_aa <- dt_aa[, .SD, .SDcols = c("height", covars_aa)] %>%
  lm(height ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

bmi_aa <- dt_aa[, .SD, .SDcols = c("bmi", covars_aa)] %>%
  lm(bmi ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

diastolic_bp_aa <- dt_aa[, .SD, .SDcols = c("diastolic_bp", covars_aa)] %>%
  lm(diastolic_bp ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

pr_interval_aa <- dt_aa[, .SD, .SDcols = c("pr_interval", covars_aa)] %>%
  lm(pr_interval ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

qrs_interval_aa <- dt_aa[, .SD, .SDcols = c("qrs_interval", covars_aa)] %>%
  lm(qrs_interval ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

egfrckdepi_clean_aa <- dt_aa[, .SD, .SDcols = c("egfrckdepi_clean", covars_aa)][, c("sex", "study") := NULL] %>%
  lm(egfrckdepi_clean ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

HbA1c_aa <- dt_aa[, .SD, .SDcols = c("HbA1c", covars_aa)] %>%
  lm(HbA1c ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

hdl_aa <- dt_aa[, .SD, .SDcols = c("hdl", covars_aa)] %>%
  lm(hdl ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

ldl_aa <- dt_aa[, .SD, .SDcols = c("ldl", covars_aa)] %>%
  lm(ldl ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

mean_corp_hgb_conc_aa <- dt_aa[, .SD, .SDcols = c("mean_corp_hgb_conc", covars_aa)] %>%
  lm(mean_corp_hgb_conc ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

platelet_cnt_aa <- dt_aa[, .SD, .SDcols = c("platelet_cnt", covars_aa)] %>%
  lm(platelet_cnt ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

systolic_bp_aa <- dt_aa[, .SD, .SDcols = c("systolic_bp", covars_aa)] %>%
  lm(systolic_bp ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

total_cholesterol_aa <- dt_aa[, .SD, .SDcols = c("total_cholesterol", covars_aa)] %>%
  lm(total_cholesterol ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

triglycerides_aa <- dt_aa[, .SD, .SDcols = c("triglycerides", covars_aa)] %>%
  lm(triglycerides ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

total_wbc_cnt_aa <- dt_aa[, .SD, .SDcols = c("total_wbc_cnt", covars_aa)] %>%
  lm(total_wbc_cnt ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

crp_aa <- dt_aa[, .SD, .SDcols = c("crp", covars_aa)] %>%
  lm(crp ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

chr_kidney_dz_aa <- dt_aa[, .SD, .SDcols = c("chr_kidney_dz", covars_aa)] %>%
  glm(chr_kidney_dz ~ ., data = ., na.action = na.exclude, family = binomial) %>%
  resid() %>%
  INT()

hypertension_aa <- dt_aa[, .SD, .SDcols = c("hypertension", covars_aa)] %>%
  glm(hypertension ~ ., data = ., na.action = na.exclude, family = binomial) %>%
  resid() %>%
  INT()

t2d_status_aa <- dt_aa[, .SD, .SDcols = c("t2d_status", covars_aa)] %>%
  glm(t2d_status ~ ., data = ., na.action = na.exclude, family = binomial) %>%
  resid() %>%
  INT()


### HCHS/SOL
height_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("height", covars_hchs_sol)] %>%
  lm(height ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

bmi_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("bmi", covars_hchs_sol)] %>%
  lm(bmi ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

diastolic_bp_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("diastolic_bp", covars_hchs_sol)] %>%
  lm(diastolic_bp ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

pr_interval_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("pr_interval", covars_hchs_sol)] %>%
  lm(pr_interval ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

qrs_interval_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("qrs_interval", covars_hchs_sol)] %>%
  lm(qrs_interval ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

egfrckdepi_clean_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("egfrckdepi_clean", covars_hchs_sol)] %>%
  lm(egfrckdepi_clean ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

HbA1c_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("HbA1c", covars_hchs_sol)] %>%
  lm(HbA1c ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

hdl_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("hdl", covars_hchs_sol)] %>%
  lm(hdl ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

ldl_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("ldl", covars_hchs_sol)] %>%
  lm(ldl ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

mean_corp_hgb_conc_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("mean_corp_hgb_conc", covars_hchs_sol)] %>%
  lm(mean_corp_hgb_conc ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

platelet_cnt_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("platelet_cnt", covars_hchs_sol)] %>%
  lm(platelet_cnt ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

systolic_bp_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("systolic_bp", covars_hchs_sol)] %>%
  lm(systolic_bp ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

total_cholesterol_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("total_cholesterol", covars_hchs_sol)] %>%
  lm(total_cholesterol ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

triglycerides_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("triglycerides", covars_hchs_sol)] %>%
  lm(triglycerides ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

total_wbc_cnt_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("total_wbc_cnt", covars_hchs_sol)] %>%
  lm(total_wbc_cnt ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

crp_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("crp", covars_hchs_sol)] %>%
  lm(crp ~ ., data = ., na.action = na.exclude) %>%
  resid() %>%
  INT()

chr_kidney_dz_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("chr_kidney_dz", covars_hchs_sol)] %>%
  glm(chr_kidney_dz ~ ., data = ., na.action = na.exclude, family = binomial) %>%
  resid() %>%
  INT()

hypertension_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("hypertension", covars_hchs_sol)] %>%
  glm(hypertension ~ ., data = ., na.action = na.exclude, family = binomial) %>%
  resid() %>%
  INT()

t2d_status_hchs_sol <- dt_hchs_sol[, .SD, .SDcols = c("t2d_status", covars_hchs_sol)] %>%
  glm(t2d_status ~ ., data = ., na.action = na.exclude, family = binomial) %>%
  resid() %>%
  INT()


#####################
## Writing Results ##
#####################
pheno_adj_aa <- data.table(
  FID = dt_aa$page_subject_id,
  IID = dt_aa$page_subject_id,
  Height = height_aa,
  WBC = total_wbc_cnt_aa,
  CRP = crp_aa,
  BMI = bmi_aa,
  DBP = diastolic_bp_aa,
  PR = pr_interval_aa,
  QRS = qrs_interval_aa,
  eGFR = egfrckdepi_clean_aa,
  HbA1c = HbA1c_aa,
  HDL = hdl_aa,
  LDL = ldl_aa,
  MCHC = mean_corp_hgb_conc_aa,
  PLT = platelet_cnt_aa,
  SBP = systolic_bp_aa,
  TC = total_cholesterol_aa,
  TG = triglycerides_aa,
  CKD = chr_kidney_dz_aa,
  Hypertension = hypertension_aa,
  T2D = t2d_status_aa
)

fwrite(pheno_adj_aa,
  file = "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_AA.txt",
  sep = " ", na = "NA", quote = FALSE
)

pheno_adj_hchs_sol <- data.table(
  FID = dt_hchs_sol$page_subject_id,
  IID = dt_hchs_sol$page_subject_id,
  Height = height_hchs_sol,
  WBC = total_wbc_cnt_hchs_sol,
  CRP = crp_hchs_sol,
  BMI = bmi_hchs_sol,
  DBP = diastolic_bp_hchs_sol,
  PR = pr_interval_hchs_sol,
  QRS = qrs_interval_hchs_sol,
  eGFR = egfrckdepi_clean_hchs_sol,
  HbA1c = HbA1c_hchs_sol,
  HDL = hdl_hchs_sol,
  LDL = ldl_hchs_sol,
  MCHC = mean_corp_hgb_conc_hchs_sol,
  PLT = platelet_cnt_hchs_sol,
  SBP = systolic_bp_hchs_sol,
  TC = total_cholesterol_hchs_sol,
  TG = triglycerides_hchs_sol,
  CKD = chr_kidney_dz_hchs_sol,
  Hypertension = hypertension_hchs_sol,
  T2D = t2d_status_hchs_sol
)

fwrite(pheno_adj_hchs_sol,
  file = "/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_HCHS_SOL.txt",
  sep = " ", na = "NA", quote = FALSE
)
