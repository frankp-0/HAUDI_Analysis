library(optparse)
library(data.table)
library(bigsnpr)

option_list <- list(
  make_option(c("--input", type = "character")),
  make_option(c("--output", type = "character")),
  make_option(c("--bim_file", type = "character"))
)

opt <- parse_args(OptionParser(option_list = option_list))


dt <- fread(opt$input)
dt <- dt[neglog10_pval_meta > 1, ]

dt <- bigsnpr::snp_modifyBuild(
  info_snp = dt,
  liftOver = "lib/liftOver",
  from = "hg19",
  to = "hg38"
)

dt$p <- 10^(-dt$neglog10_pval_meta)
dt$id <- paste(paste0("chr", dt$chr), dt$pos, dt$ref, dt$alt, sep = ":")
dt$id_swap <- paste(paste0("chr", dt$chr), dt$pos, dt$alt, dt$ref, sep = ":")

plink <- fread(opt$bim_file)

matches <- dt$id %in% plink$V2
matches_swap <- dt$id_swap %in% plink$V2

dt[matches_swap, id := id_swap]

dt <- dt[, .(chr, id, pos, ref, alt, p)]

names(dt) <- c("CHR", "SNP", "BP", "A1", "A2", "P")

fwrite(dt, file = opt$output, sep = "\t")
