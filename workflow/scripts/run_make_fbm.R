library(optparse)
library(data.table)

option_list <- list(
  make_option(c("--vcf_file"), type = "character"),
  make_option(c("--fbm_pref"), type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (grepl("AA", opt$fbm_pref)) {
  anc_names <- c("EUR", "AFR")
} else {
  anc_names <- c("EUR", "AFR", "AMR")
}

result <- HAUDI::make_fbm(
  vcf_file = opt$vcf_file,
  fbm_pref = opt$fbm_pref,
  geno_format = "HDS",
  anc_names = anc_names,
  chunk_size = 10000,
  min_ac = 20
)

saveRDS(result$FBM, file = paste0(opt$fbm_pref, ".rds"))
fwrite(result$info, file = paste0(opt$fbm_pref, "_info.txt"))
