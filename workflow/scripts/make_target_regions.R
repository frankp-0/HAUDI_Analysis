library(optparse)
library(data.table)

option_list <- list(
  make_option(c("--flare_snps_file", type = "character")),
  make_option(c("--clumped_snps_file", type = "character")),
  make_option(c("--output", type = "character"))
)

opt <- parse_args(OptionParser(option_list = option_list))

snps_flare <- readLines(opt$flare_snps_file)
snps_clumped <- readLines(opt$clumped_snps_file)

dt_flare <- data.table(
  chr = sapply(strsplit(snps_flare, ":"), "[", 1),
  pos = sapply(strsplit(snps_flare, ":"), "[", 2)
) |> na.omit()

dt_flare$pos <- as.numeric(dt_flare$pos)

dt_clumped <- data.table(
  chr = sapply(strsplit(snps_clumped, ":"), "[", 1),
  pos = sapply(strsplit(snps_clumped, ":"), "[", 2)
) |> na.omit()

dt_clumped$pos <- as.numeric(dt_clumped$pos)

dt_clumped <- dt_clumped[chr %in% dt_flare$chr, ]
dt_clumped$source <- "clump"
dt_flare <- dt_flare[chr %in% dt_clumped$chr, ]
dt_flare$source <- "flare"

dt <- rbind(dt_clumped, dt_flare)
dt <- dt[!duplicated(dt$pos), ]
dt <- setorder(dt, pos)

prev_and_cur_flare <- dt$source[1:(nrow(dt) - 1)] == "flare" &
  dt$source[2:nrow(dt)] == "flare"
keep <- c(TRUE, !prev_and_cur_flare)
dt <- dt[keep, ]
dt$source <- NULL
fwrite(dt, file = opt$output, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
