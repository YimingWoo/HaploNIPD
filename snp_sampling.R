library(data.table)

args <- commandArgs(TRUE)

# read input snps (must have two columns: chrom, pos)
snp_targets <- fread(args[1])

sampling_snps <- rbindlist(lapply(unique(snp_targets$chrom), function(chrN) {
  chrN_dat <- snp_targets[chrom == chrN]
  chrN_st <- min(chrN_dat$pos)
  chrN_ed <- max(chrN_dat$pos)
  wd <- round((chrN_ed - chrN_st) / 25e3)
  chrN_dat$intv <- cut(chrN_dat$pos, breaks = wd)
  downsampling_pos <- sapply(unique(chrN_dat$intv), function(x) {
    sample(chrN_dat[intv == x]$pos, size = 1)
  })
  data.frame(
    chrom = rep(chrN, length(downsampling_pos)),
    pos = downsampling_pos
  )
}))

fwrite(
  sampling_snps, file = "downsampling_snp_targets.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)