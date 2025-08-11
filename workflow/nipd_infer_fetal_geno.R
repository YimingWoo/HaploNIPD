source("/workspace/evan/scripts/nipd/workflow/rhdo.R")
source("/workspace/evan/scripts/nipd/workflow/hmm.R")

#'@title infer_fetal_pat_genotypes
#'@description Infer fetal genotypes based on paternal haplotypes
infer_fetal_pat_genotypes <- function(fam_genodat, dp_thresh = NA, ado1_thresh = 0.55, dist = 200) {
  ff_median <- fam_genodat$qc_info$ff_median
  error_rate <- max(0.0001, fam_genodat$qc_info$seq_error_rate)
  ado1_rate <- fam_genodat$qc_info$ado_rate1

  # paternal inheritance prediction
  if ('pat_snps' %in% names(fam_genodat)) {
    pat_snps <- fam_genodat$pat_snps
    if (nrow(pat_snps) >= 5) {
      if (ado1_rate > ado1_thresh) {
        pat_snps <- pat_snps[M_pRef < 20]
        #pat_snps <- pat_snps[M_pRef > 0.0 & M_pRef < 100.0]
      }
      if (! is.na(dp_thresh)) {
        pat_snps <- pat_snps[cfdna_total_dp >= dp_thresh]
      }
      pat_snps <- pos_dist_filter(pat_snps[snp_type == "type3"], dist)

      # chromosome-wide deduction
      pat_chroms <- unique(pat_snps$chrom)
      pat_hapblks <- vector("list", length(pat_chroms))
      fam_genodat$pat_snps_results_chrom <- data.table::rbindlist(
        lapply(seq_along(pat_chroms), function(i) {
          pat_snps_chrN <- pat_snps[chrom == pat_chroms[i]]
          if (nrow(pat_snps_chrN) < 5) {
            pat_snps_chrN <- NULL
          }
          else {
            pat_snps_chrN <- fetal_paternal_genotyping(
              pat_snps_chrN, ff_median, error_rate
            )
            pat_snps_chrN[, blk_pass := FALSE]
            pat_chrN_blkqc <- paternal_hap_postqc(pat_snps_chrN)
            for (j in seq_len(nrow(pat_chrN_blkqc))) {
              pat_snps_chrN[pat_chrN_blkqc$blk_st[j]:pat_chrN_blkqc$blk_ed[j],
                blk_pass := pat_chrN_blkqc$blk_pass[j]]
            }
            pat_hapblks[[i]] <<- pat_chrN_blkqc
          }
          pat_snps_chrN
        }
      ))
      fam_genodat$pat_hapblks <- rbindlist(pat_hapblks)
    }
  }
  return(fam_genodat)
}


#'@title infer_fetal_pat_genotypes
#'@description Infer fetal genotypes based on maternal haplotypes
infer_fetal_mat_genotypes <- function(fam_genodat, dist = 200, N_thresh = 1000,
  dp_thresh = NA, use_type5 = FALSE, type5_padding_len = 1e6, biased_ref_ratio = NA) {
  ff_median <- fam_genodat$qc_info$ff_median

  # paternal inheritance prediction
  if ('mat_snps' %in% names(fam_genodat)) {
    mat_snps <- fam_genodat$mat_snps
    if (use_type5 == TRUE && "pat_snps_results_chrom" %in% names(fam_genodat))
    {
      hetXhet_snps <- fam_genodat$pat_snps[snp_type == "type5"]
      hetXhet_phased_snps <- data.table::rbindlist(lapply(
        unique(fam_genodat$pat_hapblks$chrom), function(chrN) {
          chrN_phased_info <- fam_genodat$pat_hapblks[chrom == chrN & blk_pass == TRUE & snpnum > 25]
          chrN_hetXhet <- hetXhet_snps[chrom == chrN]
          chrN_hetXhet[, hmm_fetal_hap := -1]

          for (i in seq_len(nrow(chrN_phased_info))) {
            pat_phased_hap <- fam_genodat$pat_snps_results_chrom[
              chrom == chrN & pos > chrN_phased_info$stpos[i] &
                pos < chrN_phased_info$edpos[i]]$hmm_fetal_hap[5] # since snpnum > 25, 5 is a safe choice
            chrN_hetXhet[
              (pos >= chrN_phased_info$stpos[i] + type5_padding_len) &
                (pos <= chrN_phased_info$edpos[i] - type5_padding_len),
              hmm_fetal_hap := pat_phased_hap]
          }
          chrN_hetXhet[hmm_fetal_hap != -1]
      }))
      hetXhet_phased_snps[, F0 := fifelse(hmm_fetal_hap == 0, F0, F1)][, F1 := F0]
      hetXhet_phased_snps[, hmm_fetal_hap := NULL]
      hetXhet_phased_snps_reorder <- cbind(
        hetXhet_phased_snps[, which(! colnames(hetXhet_phased_snps) %in% c("M0", "M1", "F0", "F1", "snp_type")), with = FALSE],
        hetXhet_phased_snps[, c("F0", "F1", "M0", "M1", "snp_type")]
      )
      mat_snps <- rbind(mat_snps[, snp_type := "type4"], hetXhet_phased_snps_reorder)
      mat_snps <- mat_snps[order(chrom, pos)]
    }

    if (! is.na(dp_thresh)) {
      mat_snps <- mat_snps[cfdna_total_dp >= dp_thresh]
    }
    mat_snps <- pos_dist_filter(mat_snps, dist)

    if (nrow(mat_snps) >= 25) {
      # chromosome-wide deduction
      fam_genodat$mat_snps_results_chrom <- data.table::rbindlist(lapply(unique(mat_snps$chrom),
        function(chrN) {
          mat_snps_chrN <- mat_snps[chrom == chrN]
          if (nrow(mat_snps_chrN) < 25) {
            mat_snps_chrN <- NULL
          }
          else if (chrN == "chrX" & fam_genodat$fetal_gender == "male") {
            mat_snps_chrN <- male_fetal_maternal_genotyping_chrX(mat_snps_chrN, ff_median, N = N_thresh)
          }
          else {
            if (! is.na(biased_ref_ratio)) {
              mat_snps_chrN[, ':=' (
                cfdna_ref_dp = round(cfdna_ref_dp / biased_ref_ratio * 50),
                cfdna_alt_dp = round(cfdna_alt_dp / (100 - biased_ref_ratio) * 50)
              )][, cfdna_total_dp := cfdna_ref_dp + cfdna_alt_dp]
            }
            mat_snps_chrN <- fetal_maternal_genotyping(mat_snps_chrN, ff_median, N = N_thresh)
          }
          return(mat_snps_chrN)
        }
      ))
    }
  }
  return(fam_genodat)
}


#'@title pos_dist_filter
#'@description Filter variants by certain distance default
pos_dist_filter <- function(data, dist = 200) {
  data[, distance := c(dist + 100, diff(pos)), by = "chrom"]
  new_data <- data.table::rbindlist(lapply(unique(data$chrom), function(chrN) {
    chrN_dat <- data[chrom == chrN]
    if (nrow(chrN_dat) < 5) {
      chrN_dat
    } else {
      kpt_ind <- rep(0, nrow(chrN_dat))
      kpt_ind[1] <- 1
      chrN_cumsum_dist <- 0
      for (i in seq(2, nrow(chrN_dat), 1)) {
        chrN_cumsum_dist <- chrN_cumsum_dist + chrN_dat$distance[i]
        if (chrN_cumsum_dist >= dist) {
          kpt_ind[i] <- i
          chrN_cumsum_dist <- 0
        }
      }
      kpt_ind <- kpt_ind[which(kpt_ind > 0)]
      chrN_dat[kpt_ind, ]
    }
  }))
  new_data <- new_data[, distance := NULL]
  #new_data <- data[abs(distance) >= dist][, distance := NULL]
  return(new_data)
}


#'@title fetal_paternal_genotyping
#'@description Infer fetal genotypes based on paternal haplotypes
fetal_paternal_genotyping <-
  function(vcfdat, fetal_fraction, error_rate) {
    # infer fetal genotype based on HMM method
    mapinfo <- genetic_mapinfo[chrom == unique(vcfdat$chrom) &
                                 pos >= min(vcfdat$pos) & pos <= max(vcfdat$pos)]
    vcfdat <- HMM_paternal_genotyping(
      vcfdat, fetal_fraction, error_rate, mapinfo
    )

    # infer fetal genotype based on reads-depth method
    vcfdat <- paternal_reads_pass(vcfdat)

    return(vcfdat)
  }


#'@title paternal_reads_pass
#'@description Infer fetal genotypes based on reads-pass criteria ("Naive")
paternal_reads_pass <- function(vcfdat, min_dp_threshold = 2) {
  vcfdat[, naive_fetal_hap := data.table::fifelse(M0 == 0,
            data.table::fifelse(
            F0 == 0,
            data.table::fifelse(cfdna_alt_dp > min_dp_threshold, 1, 0),
            data.table::fifelse(cfdna_alt_dp > min_dp_threshold, 0, 1)
            ),
            data.table::fifelse(
            F0 == 0,
            data.table::fifelse(cfdna_ref_dp > min_dp_threshold, 0, 1),
            data.table::fifelse(cfdna_ref_dp > min_dp_threshold, 1, 0)
            ))]

  # CBS segmentation
  CNA_obj <- DNAcopy::CNA(
    vcfdat$naive_fetal_hap,
    chrom = vcfdat$chrom,
    maploc = vcfdat$pos,
    data.type = "binary"
  )
  seg_CNA_obj <- DNAcopy::segment(CNA_obj, verbose = 0)
  seg_CNA_dat <- seg_CNA_obj$output
  vcfdat <- rbindlist(lapply(seq_len(nrow(seg_CNA_dat)), function(n) {
      seg_vcfdat <- vcfdat[chrom == seg_CNA_dat$chrom[n] &
                             pos >= seg_CNA_dat$loc.start[n] &
                             pos <= seg_CNA_dat$loc.end[n]]
      seg_vcfdat[, naive_fetal_hap := data.table::fifelse(seg_CNA_dat$seg.mean[n] >= 0.5, 1, 0)]
      #return(seg_vcfdat[, naive_fetal_hapName := paste0("F", naive_fetal_hap)])
      return(seg_vcfdat)
    }))

  vcfdat[, naive_fetal_genotype := data.table::fifelse(
            naive_fetal_hap == 0, paste(M0, F0, sep = "/"), paste(M0, F1, sep = "/"))][,
            naive_fetal_genotype := data.table::fifelse(
              naive_fetal_genotype == "1/0", "0/1", naive_fetal_genotype)]

  return(vcfdat)
}


#'@title HMM_paternal_genotyping
#'@description Infer fetal genotypes based on Hidden Markov Model
HMM_paternal_genotyping <-
  function(vcfdat, fetal_fraction, error_rate, mapinfo) {
    uplimit <- 0.95
    lowlimit <- 0.05

    # compute observation probability
    vcfdat[, ':=' (
      F0_prob = data.table::fcase(
        M0 == 0 & F0 == 0, stats::dbinom(cfdna_alt_dp, cfdna_total_dp, error_rate),
        M0 == 0 & F0 == 1, stats::dbinom(cfdna_alt_dp, cfdna_total_dp, 0.5 * fetal_fraction),
        M0 == 1 & F0 == 0, stats::dbinom(cfdna_ref_dp, cfdna_total_dp, 0.5 * fetal_fraction),
        M0 == 1 & F0 == 1, stats::dbinom(cfdna_ref_dp, cfdna_total_dp, error_rate)
      ),
      F1_prob = data.table::fcase(
        M0 == 0 & F0 == 0, stats::dbinom(cfdna_alt_dp, cfdna_total_dp, 0.5 * fetal_fraction),
        M0 == 0 & F0 == 1, stats::dbinom(cfdna_alt_dp, cfdna_total_dp, error_rate),
        M0 == 1 & F0 == 0, stats::dbinom(cfdna_ref_dp, cfdna_total_dp, error_rate),
        M0 == 1 & F0 == 1, stats::dbinom(cfdna_ref_dp, cfdna_total_dp, 0.5 * fetal_fraction)
      )
    )]

    vcfdat[, P0 := F0_prob / (F0_prob + F1_prob)][,
             P0 := pmin(pmax(P0, lowlimit), uplimit)][, P1 := 1 - P0]

    vcfdat <- hmm_genotyping(vcfdat, mapinfo)

    vcfdat[, hmm_fetal_genotype := data.table::fifelse(
      hmm_fetal_hap == 0, paste(F0, M0, sep = "/"), paste(F1, M0, sep = "/"))][,
        hmm_fetal_genotype := data.table::fifelse(
          hmm_fetal_genotype == "1/0", "0/1", hmm_fetal_genotype)]

    #return(vcfdat[, hmm_fetal_hapName := paste0("F", hmm_fetal_hap)])
    return(vcfdat)
  }


#'@title paternal_hap_postqc
#'@description perform postqc on phasing paternal haplotype blocks
paternal_hap_postqc <- function(chrN_dat, min_blksnps = 10) {
  chrN_hmm_hap_switch <- which(diff(chrN_dat$hmm_fetal_hap) != 0)
  chrN_hmm_hap_blks <- c(chrN_hmm_hap_switch, nrow(chrN_dat))
  chrN <- chrN_dat$chrom[1]
  chrN_blkdat <- data.table::rbindlist(lapply(seq_along(chrN_hmm_hap_blks), function(j) {
    hmm_blk_st <- ifelse(j == 1, 1, chrN_hmm_hap_blks[j - 1] + 1)
    hmm_blk_ed <- chrN_hmm_hap_blks[j]
    blk_snpnum <- hmm_blk_ed - hmm_blk_st + 1
    return(list(
      chrom = chrN, blk_st = hmm_blk_st, blk_ed = hmm_blk_ed,
      stpos = chrN_dat$pos[hmm_blk_st], edpos = chrN_dat$pos[hmm_blk_ed],
      snpnum = blk_snpnum,
      blk_pass = ifelse(blk_snpnum < min_blksnps, FALSE, TRUE)
    ))
  }))
  return(chrN_blkdat)
}


#'@title fetal_maternal_genotyping
#'@description Infer fetal genotypes based on maternal haplotypes
fetal_maternal_genotyping <- function(vcfdat, fetal_fraction, N = 1000,
  do_hmm = TRUE, do_rhdo = TRUE) {
  # infer fetal genotype based on HMM method
  mapinfo <- genetic_mapinfo[chrom == unique(vcfdat$chrom) &
                               pos >= min(vcfdat$pos) & pos <= max(vcfdat$pos)]
  if (do_hmm) {
    vcfdat <- HMM_maternal_genotyping(vcfdat, fetal_fraction, mapinfo)
  }

  # infer fetal genotype based on RHDO method (back and forth)
  if (do_rhdo) {
    vcfdat_rev <- RHDO_maternal_genotyping(vcfdat[order(-pos), ], fetal_fraction, N)
    vcfdat <- RHDO_maternal_genotyping(vcfdat, fetal_fraction, N)
    vcfdat[, ':=' (
      rhdo_rev_fetal_hap = vcfdat_rev$rhdo_fetal_hap,
      rhdo_rev_blk = vcfdat_rev$rhdo_blk,
      rhdo_rev_flag = vcfdat_rev$rhdo_flag,
      rhdo_rev_blk_type = vcfdat_rev$rhdo_blk_type,
      rhdo_rev_fetal_genotype = vcfdat_rev$rhdo_fetal_genotype
    )]

    # blk_rev = RHDO_maternal_genotyping(vcfdat[order(-pos), ], fetal_fraction, N)
    # blk_fw <- RHDO_maternal_genotyping(vcfdat, fetal_fraction, N)
    # blks <- rbind(blk_rev[, direction := "rev"], blk_fw[, direction := "fw"])
    # return(blks)
  }
  return(vcfdat)
}

male_fetal_maternal_genotyping_chrX <- function(vcfdat, fetal_fraction, N = 1000,
  do_hmm = TRUE, do_rhdo = TRUE) {
  # infer fetal genotype based on HMM method
  mapinfo <- genetic_mapinfo[chrom == unique(vcfdat$chrom) &
                               pos >= min(vcfdat$pos) & pos <= max(vcfdat$pos)]
  if (do_hmm) {
    vcfdat <- HMM_maternal_genotyping_male_chrX(vcfdat, fetal_fraction, mapinfo)
  }

  # infer fetal genotype based on RHDO method (back and forth)
  if (do_rhdo) {
    vcfdat_rev <- vcfdat[order(-pos), ]
    vcfdat_rev <- RHDO_maternal_genotyping_male_chrX(vcfdat_rev, fetal_fraction, N)
    vcfdat <- RHDO_maternal_genotyping_male_chrX(vcfdat, fetal_fraction, N)
    vcfdat[, ':=' (
      rhdo_rev_fetal_hap = vcfdat_rev$rhdo_fetal_hap,
      rhdo_rev_blk = vcfdat_rev$rhdo_blk,
      rhdo_rev_flag = vcfdat_rev$rhdo_flag,
      rhdo_rev_blk_type = vcfdat_rev$rhdo_blk_type,
      rhdo_rev_fetal_genotype = vcfdat_rev$rhdo_fetal_genotype
    )]
  }
  return(vcfdat)
}


#'@title HMM_maternal_genotyping
#'@description Infer fetal genotypes based on Hidden Markov Model
#'applied on variants data on one chromosome.
#'
#'@param vcfdat vcf data table of maternal specific variants
#'@param ff fetal fraction estimation
#'@param genetic_map genetic maps information
#'
#'@import data.table
#'
#'@return deduced fetal genotypes
#'
#'@export
HMM_maternal_genotyping <- function(vcfdat, fetal_fraction, mapinfo) {
  uplimit <- 0.95
  lowlimit <- 0.05

  vcfdat[, ':=' (
    M0_prob = data.table::fcase(
      F0 == 0 & M0 == 0, dbinom(cfdna_alt_dp, cfdna_total_dp, (1 - fetal_fraction) / 2),
      F0 == 0 & M0 == 1, dbinom(cfdna_alt_dp, cfdna_total_dp, 0.5),
      F0 == 1 & M0 == 0, dbinom(cfdna_ref_dp, cfdna_total_dp, 0.5),
      F0 == 1 & M0 == 1, dbinom(cfdna_ref_dp, cfdna_total_dp, (1 - fetal_fraction) / 2)
    ),
    M1_prob = data.table::fcase(
      F0 == 0 & M0 == 0, dbinom(cfdna_alt_dp, cfdna_total_dp, 0.5),
      F0 == 0 & M0 == 1, dbinom(cfdna_alt_dp, cfdna_total_dp, (1 - fetal_fraction) / 2),
      F0 == 1 & M0 == 0, dbinom(cfdna_ref_dp, cfdna_total_dp, (1 - fetal_fraction) / 2),
      F0 == 1 & M0 == 1, dbinom(cfdna_ref_dp, cfdna_total_dp, 0.5)
    )
  )]

  vcfdat[, P0 := M0_prob / (M0_prob + M1_prob)][,
         P0 := pmin(pmax(P0, lowlimit), uplimit)][, P1 := 1 - P0]

  vcfdat <- hmm_genotyping(vcfdat, mapinfo)

  vcfdat[, hmm_fetal_genotype := data.table::fifelse(hmm_fetal_hap == 0,
            paste(F0, M0, sep = "/"),
            paste(F0, M1, sep = "/"))][,
               hmm_fetal_genotype := data.table::fifelse(hmm_fetal_genotype == "1/0", 
               "0/1", hmm_fetal_genotype)]

  #return(vcfdat[, hmm_fetal_hapName := paste0("M", hmm_fetal_hap)])
  return(vcfdat)
}


HMM_maternal_genotyping_male_chrX <- function(vcfdat, fetal_fraction, mapinfo) {
  uplimit <- 0.95
  lowlimit <- 0.05

  vcfdat[M0 == 0, ':=' (
    M0_prob = dbinom(cfdna_ref_dp, cfdna_total_dp, 1 / (2 - fetal_fraction)),
    M1_prob = dbinom(cfdna_alt_dp, cfdna_total_dp, 1 / (2 - fetal_fraction))
  )]

  vcfdat[M0 == 1, ':=' (
    M0_prob = dbinom(cfdna_alt_dp, cfdna_total_dp, 1 / (2 - fetal_fraction)),
    M1_prob = dbinom(cfdna_ref_dp, cfdna_total_dp, 1 / (2 - fetal_fraction))
  )]

  vcfdat[, P0 := M0_prob / (M0_prob + M1_prob)][,
         P0 := pmin(pmax(P0, lowlimit), uplimit)][, P1 := 1 - P0]

  vcfdat <- hmm_genotyping(vcfdat, mapinfo)

  vcfdat[, hmm_fetal_genotype := data.table::fifelse(hmm_fetal_hap == 0,
            paste(M0, M0, sep = "/"),
            paste(M1, M1, sep = "/"))]

  #return(vcfdat[, hmm_fetal_hapName := paste0("M", hmm_fetal_hap)])
  return(vcfdat)
}


check_rhdo_bks <- function(chrN_dat, left_bk, right_bk) {
  fw_rhdo_haps <- chrN_dat[left_bk:right_bk]$rhdo_fetal_hap
  rev_rhdo_haps <- chrN_dat[left_bk:right_bk]$rhdo_rev_fetal_hap
  fw_rhdo_haptab <- table(fw_rhdo_haps[fw_rhdo_haps > -1])
  rev_rhdo_haptab <- table(rev_rhdo_haps[rev_rhdo_haps > -1])
  ifelse(
    any(c(length(fw_rhdo_haptab), length(rev_rhdo_haptab)) > 1),
    TRUE, FALSE
  )
}

cmp_hmm_rhdo_bks <- function(vcfdat, padding_snps = 25) {
  ret <- data.table::rbindlist(lapply(unique(vcfdat$chrom), function(chrN) {
    chrN_dat <- vcfdat[chrom == chrN]
    chrN_hmm_hap_blks <- c(which(diff(chrN_dat$hmm_fetal_hap) != 0), nrow(chrN_dat))
    chrN_blks <- rbindlist(lapply(seq_along(chrN_hmm_hap_blks), function(j) {
      hap_st <- ifelse(j == 1, 1, chrN_hmm_hap_blks[j - 1] + 1)
      hap_ed <- chrN_hmm_hap_blks[j]
      blk_snpnum <- hap_ed - hap_st + 1
      blk_amn_diff <- chrN_dat[hap_st:hap_ed]$amn_diff
      if (j == 1) {
        left_rhdo_valid <- NA
        bk_ind2 <- hap_ed
        bk_pad_st2 <- ifelse(bk_ind2 - padding_snps < 0, 1, bk_ind2 - padding_snps)
        bk_pad_ed2 <- ifelse(bk_ind2 + padding_snps > nrow(chrN_dat), nrow(chrN_dat), bk_ind2 + padding_snps)
        right_rhdo_valid <- check_rhdo_bks(chrN_dat, bk_pad_st2, bk_pad_ed2)
      } else if (j == length(chrN_hmm_hap_blks)) {
        bk_ind1 <- hap_st
        bk_pad_st1 <- ifelse(bk_ind1 - padding_snps < 0, 1, bk_ind1 - padding_snps)
        bk_pad_ed1 <- ifelse(bk_ind1 + padding_snps > nrow(chrN_dat), nrow(chrN_dat), bk_ind1 + padding_snps)
        left_rhdo_valid <- check_rhdo_bks(chrN_dat, bk_pad_st1, bk_pad_ed1)
        right_rhdo_valid <- NA
      } else {
        bk_ind1 <- hap_st
        bk_pad_st1 <- ifelse(bk_ind1 - padding_snps < 0, 1, bk_ind1 - padding_snps)
        bk_pad_ed1 <- ifelse(bk_ind1 + padding_snps > nrow(chrN_dat), nrow(chrN_dat), bk_ind1 + padding_snps)
        left_rhdo_valid <- check_rhdo_bks(chrN_dat, bk_pad_st1, bk_pad_ed1)

        bk_ind2 <- hap_ed
        bk_pad_st2 <- ifelse(bk_ind2 - padding_snps < 0, 1, bk_ind2 - padding_snps)
        bk_pad_ed2 <- ifelse(bk_ind2 + padding_snps > nrow(chrN_dat), nrow(chrN_dat), bk_ind2 + padding_snps)
        right_rhdo_valid <- check_rhdo_bks(chrN_dat, bk_pad_st2, bk_pad_ed2)
      }


      fw_cmpdat <- chrN_dat[hap_st:hap_ed][rhdo_fetal_hap > -1 & rhdo_flag != "fail"]
      rev_cmpdat <- chrN_dat[hap_st:hap_ed][rhdo_rev_fetal_hap > -1 & rhdo_rev_flag != "fail"]
      blk_rhdo_fw_diff <- length(which(fw_cmpdat$hmm_fetal_hap != fw_cmpdat$rhdo_fetal_hap))
      blk_rhdo_rev_diff <- length(which(rev_cmpdat$hmm_fetal_hap != rev_cmpdat$rhdo_rev_fetal_hap))
      blk_hmm_hap <- chrN_dat[hap_st:hap_ed]$hmm_fetal_hap[2]

      return(list(
        chrom = chrN, stpos = chrN_dat$pos[hap_st], edpos = chrN_dat$pos[hap_ed],
        snpnum = blk_snpnum, len_mb = round((chrN_dat$pos[hap_ed] - chrN_dat$pos[hap_st]) / 1e6, 3),
        concord = round(length(which(blk_amn_diff == FALSE)) / blk_snpnum, 4),
        rhdo_left_valid = left_rhdo_valid, rhdo_right_valid = right_rhdo_valid,
        rhdo_fw_diff = blk_rhdo_fw_diff / nrow(fw_cmpdat), rhdo_rev_diff = blk_rhdo_rev_diff / nrow(rev_cmpdat)
      ))
    }))
    chrN_blks
  }))
  return(ret)
}

# mat_cpt_pred <- function(vcfdat, Q_thresh = 10) {
#   cpt_pred_blkinfo <- data.table::rbindlist(lapply(unique(vcfdat$chrom), function(chrN) {
#     chrN_dat <- vcfdat[chrom == chrN]
#     hap_rel_logratio <- log(chrN_dat$P0 / chrN_dat$P1)
#     cpt_res <- changepoint::cpt.meanvar(
#       hap_rel_logratio, method = "BinSeg", Q = Q_thresh, pen.value = 0.05, penalty="AIC"
#     )
#     pred_blks <- c(changepoint::cpts(cpt_res), nrow(chrN_dat))
#     chrN_blks <- data.table::rbindlist(lapply(seq_along(pred_blks), function(j) {
#       hap_st <- ifelse(j == 1, 1, pred_blks[j-1] + 1)
#       hap_ed <- pred_blks[j]
#       blk_haptab <- table(chrN_dat[hap_st:hap_ed]$hmm_fetal_hap)
#       blk_hap <- as.integer(names(blk_haptab[which.max(blk_haptab)]))
#       return(list(
#         chrom = chrN, stpos = chrN_dat$pos[hap_st], edpos = chrN_dat$pos[hap_ed],
#         hap = blk_hap
#       ))
#     }))
#     chrN_blks
#   }))
#   return(cpt_pred_blkinfo)
# }

# mat_rhdo_blks_merge <- function(vcfdat) {
#   ret <- data.table::rbindlist(lapply(unique(vcfdat$chrom), function(chrN) {
#     chrN_dat <- vcfdat[chrom == chrN]
#     fw_cpt_res <- changepoint::cpt.meanvar(
#       chrN_dat$rhdo_fetal_hap, method = "BinSeg",pen.value = 0.01, Q = 16
#     )
#     type_A <- chrN_dat[rhdo_blk_type == "alpha"]
#     type_B <- chrN_dat[rhdo_blk_type == "beta"]
#     type_A[, rhdo_hap := fcase(
#       rhdo_flag != "fail" & rhdo_rev_flag == "fail", rhdo_fetal_hap,
#       rhdo_flag == "fail" & rhdo_rev_flag != "fail", rhdo_rev_fetal_hap,
#       rhdo_flag != "fail" & rhdo_rev_flag != "fail", (rhdo_fetal_hap + rhdo_rev_fetal_hap) / 2,
#       default = NA
#     )]
#     type_A_bks <- which(diff(type_A[! is.na(rhdo_hap)]$rhdo_hap) != 0)
#     type_B[, rhdo_hap := fcase(
#       rhdo_flag != "fail" & rhdo_rev_flag == "fail", rhdo_fetal_hap,
#       rhdo_flag == "fail" & rhdo_rev_flag != "fail", rhdo_rev_fetal_hap,
#       rhdo_flag != "fail" & rhdo_rev_flag != "fail", (rhdo_fetal_hap + rhdo_rev_fetal_hap) / 2,
#       default = NA
#     )]
#     merged_dat <- rbind(type_A, type_B)
#     merged_dat <- merged_dat[order(pos)]
#     x <- cpt.meanvar(merged_dat$rhdo_hap)
#     pred_blks <- c(changepoint::cpts(cpt_res), nrow(chrN_dat))
#     chrN_blks <- data.table::rbindlist(lapply(seq_along(pred_blks), function(j) {
#       hap_st <- ifelse(j == 1, 1, pred_blks[j-1] + 1)
#       hap_ed <- pred_blks[j]
#       blk_haptab <- table(chrN_dat[hap_st:hap_ed]$hmm_fetal_hap)
#       blk_hap <- as.integer(names(blk_haptab[which.max(blk_haptab)]))
#       return(list(
#         chrom = chrN, stpos = chrN_dat$pos[hap_st], edpos = chrN_dat$pos[hap_ed],
#         hap = blk_hap
#       ))
#     }))
#     chrN_blks
#   }))
#   return(ret)
# }