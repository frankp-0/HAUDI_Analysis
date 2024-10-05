library(optparse)
library(data.table)
library(HAUDI)
library(genlasso)

option_list <- list(
  make_option("--rds", type = "character"),
  make_option("--info", type = "character"),
  make_option("--pheno_file", type = "character"),
  make_option("--train_file", type = "character"),
  make_option("--validation_file", type = "character"),
  make_option("--test_file", type = "character"),
  make_option("--clump_file", type = "character"),
  make_option("--samples_file", type = "character"),
  make_option("--n_snps", type = "integer"),
  make_option("--gamma", type = "numeric"),
  make_option("--pheno", type = "character"),
  make_option("--model", type = "character"),
  make_option("--results", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

fbm_obj <- readRDS(opt$rds)
dt_info <- fread(opt$info)
dt_pheno <- fread(opt$pheno_file)
train_ids <- readLines(opt$train_file)
test_ids <- readLines(opt$test_file)
validation_ids <- readLines(opt$validation_file)
samples <- readLines(opt$samples_file)

dt_pheno <- dt_pheno[match(samples, IID), ]
y <- dt_pheno[[opt$pheno]]

ind_train <- which(dt_pheno$IID %in% train_ids) |> sort()
ind_test <- which(dt_pheno$IID %in% test_ids) |> sort()
ind_validation <- which(dt_pheno$IID %in% validation_ids) |> sort()

dt_clump <- fread(opt$clump_file)
dt_clump <- setorder(dt_clump, P)
dt_clump <- dt_clump[SNP %in% dt_info$rsid, ]
snps <- dt_clump[seq_len(min(nrow(dt_clump), opt$n_snps)), ]$SNP

time_run <- system.time({
  model <- HAUDI::gaudi(
    fbm_obj = fbm_obj,
    fbm_info = dt_info,
    y = y,
    gamma = opt$gamma,
    ind_train = ind_train,
    snps = snps,
    k = 5
  )
})


X <- HAUDI::construct_gaudi(fbm_obj, dt_info, snps)

y_pred <- predict(model$fit, Xnew = X, lambda = model$best_lambda)$fit[, 1]

r2_validation <- cor(y_pred[ind_validation], y[ind_validation])^2
r2_test <- cor(y_pred[ind_test], y[ind_test])^2

df_result <- data.frame(
  method = "GAUDI",
  cohort = ifelse(grepl("AA", opt$pheno_file), "AA", "HCHS_SOL"),
  phenotype = opt$pheno,
  n_snps = opt$n_snps,
  gamma = opt$gamma,
  r2_validation = r2_validation,
  r2_test = r2_test,
  runtime = time_run[3]
)

write.table(
  x = df_result, file = opt$results, sep = "\t",
  quote = FALSE, row.names = FALSE, col.names = TRUE
)

saveRDS(
  object = model, file = opt$model
)
