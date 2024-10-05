library(data.table)

dt_AA <- fread("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_AA.txt")
dt_HCHS_SOL <- fread("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/adjusted/pheno_adj_HCHS_SOL.txt")

phenotypes <- names(dt_AA)[3:21]

for (sub in c("AA", "HCHS_SOL")) {
  for (pheno in c(phenotypes)) {
    dt <- get(paste0("dt_", sub))[, .SD, .SDcols = c("IID", pheno)]
    ids <- na.omit(dt)$IID
    ids <- sample(ids)
    test_folds <- as.numeric(cut(seq_along(ids), 5))
    for (fold in 1:5) {
      test_ids <- ids[test_folds == fold]
      validation_ids <- sample(ids[test_folds != fold], size = 0.1 * length(ids))
      train_ids <- ids[!(ids %in% c(test_ids, validation_ids))]
      writeLines(
        test_ids,
        con = paste0("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold", fold, "_", pheno, "_", sub, "_test.txt")
      )
      writeLines(
        validation_ids,
        con = paste0("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold", fold, "_", pheno, "_", sub, "_validation.txt")
      )
      writeLines(
        train_ids,
        con = paste0("/work/users/f/r/frocko/HAUDI/PAGE_analysis/data/PAGE/pheno/folds/fold", fold, "_", pheno, "_", sub, "_train.txt")
      )
    }
  }
}
