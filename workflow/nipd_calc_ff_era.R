#'@title preprocess_nipd_info
#'@description Perform preprocess of family-based VCF, estimate
preprocess_nipd_info <- function(fam_genodat, coverage_dir = NA) {

  # basic filtering
  fam_genodat$vcfdat[, tot_rel := (cfdna_ref_dp + cfdna_alt_dp) / cfdna_total_dp]
  vcfdat <- fam_genodat$vcfdat[tot_rel > 0.9][, tot_rel := NULL]

  # calculate sequencing errors and fetal fraction
  vcfdat[, F_pRef := round(100 * father_ref_dp / father_total_dp, 3)]
  if (fam_genodat$has_mom_geno == TRUE) {
    vcfdat[, M_pRef := round(100 * mother_ref_dp / mother_total_dp, 3)]
    fam_genodat$qc_info <- estimate_error_ratio_and_ff(
      vcfdat, rel_mode = fam_genodat$sampmode
    )
  } else {
    vcfdat[, M_pRef := round(100 * cfdna_ref_dp / cfdna_total_dp, 3)]
    fam_genodat$qc_info <- estimate_error_ratio_and_ff(
      vcfdat, rel_mode = fam_genodat$sampmode,
      hom_lower_limit = 15, hom_upper_limit = 85
    )
  }

    # determine fetal gender
  if (! "fetal_gender" %in% names(fam_genodat)) {
    if (is.na(coverage_dir)) {
      fam_genodat$fetal_gender <- NA
    } else {
      dpfile <- sprintf("%s/%s.mosdepth.summary.txt", coverage_dir, fam_genodat$sampname)
      if (! file.exists(dpfile)) {
        fam_genodat$fetal_gender <- NA
      } else {
        # estimate fetus gender
        covdat <- data.table::fread(dpfile)
        cov_regions <- paste0("chr", c(1:18, 20:22, 'Y'), "_region")
        covdat <- covdat[chrom %in% cov_regions]
        cfdna_y_ratio <-
          covdat$bases[nrow(covdat)] / sum(covdat[-nrow(covdat), ]$bases)
        fam_genodat$fetal_gender <-
          estimate_fetal_gender(fam_genodat$qc_info$ff_median, cfdna_y_ratio)
      }
    }
  }
  if (is.na(fam_genodat$fetal_gender)) {
    vcfdat <- vcfdat[chrom != "chrX"]
  }
  
  # trim VCF data to chromosomes of ROI regions
  if ("roi_info" %in% names(fam_genodat)) {
    vcfdat <- vcfdat[chrom %in% fam_genodat$roi_info$chrom]
  }

  fam_genodat$vcfdat <- vcfdat
  return(fam_genodat)
}


#'@title estimate_error_ratio_and_ff
#'@description Estimate global Sequencing Error Ratio from VCF data
estimate_error_ratio_and_ff <-
  function(vcfdat, rel_mode,
           hom_lower_limit = 10, hom_upper_limit = 90,
           het_lower_limit = 35, het_upper_limit = 65, 
           ff_lower_limit = 0.02, ff_upper_limit = 0.4,
           seqerror_threshold = 0.02, cf_dp_thresh = 100) {
    query_dat_autosomal <- vcfdat[chrom != "chrX" & cfdna_total_dp >= cf_dp_thresh &
      father_genotype %in% c("0/0", "1/1")]
    if (! "chrX" %in% vcfdat$chrom) {
      query_dat_chrX <- NULL
    } else {
      query_dat_chrX <- vcfdat[chrom == "chrX" & cfdna_total_dp >= cf_dp_thresh &
        father_genotype %in% c("0/0", "1/1")]
    }

    # snp data for error ratio calculation
    data2a <- query_dat_autosomal[M_pRef >= hom_upper_limit & F_pRef >= hom_upper_limit][,
      obs_error_ratio := cfdna_alt_dp / cfdna_total_dp]
    data2b <- query_dat_autosomal[M_pRef <= hom_lower_limit & F_pRef <= hom_lower_limit][,
      obs_error_ratio := cfdna_ref_dp / cfdna_total_dp]
    # filter abnormal error sites (maybe maternal low fraction variants)
    data2a_flt <- data2a[obs_error_ratio < seqerror_threshold]
    data2b_flt <- data2b[obs_error_ratio < seqerror_threshold]

    data2_seqerror_cov <- sum(data2a_flt$cfdna_alt_dp) + sum(data2b_flt$cfdna_ref_dp)
    data2_total_cov <- sum(data2a_flt$cfdna_total_dp) + sum(data2b_flt$cfdna_total_dp)
    seq_error_rate <- round(data2_seqerror_cov / data2_total_cov, 5)

    # snp data for fetal fraction calculation
    data1a <- query_dat_autosomal[M_pRef >= hom_upper_limit & F_pRef <= hom_lower_limit][,
      ff := 2 * cfdna_alt_dp / cfdna_total_dp]
    data1b <- query_dat_autosomal[M_pRef <= hom_lower_limit & F_pRef >= hom_upper_limit][,
      ff := 2 * cfdna_ref_dp / cfdna_total_dp]
    ff_type1 <- c(data1a$ff, data1b$ff)
    # filter abnormal fetal fraction sites
    ff_type1_flt <- ff_type1[ff_type1 <= ff_upper_limit & ff_type1 >= ff_lower_limit]

    ff_bottom_1perc <- unname(stats::quantile(ff_type1_flt, 0.01))
    ff_top_1perc <- unname(stats::quantile(ff_type1_flt, 0.99))
    ff_type1_flt <- ff_type1_flt[ff_type1_flt > ff_bottom_1perc & ff_type1_flt < ff_top_1perc]
    ff_median <- stats::median(ff_type1_flt)

    ret_info <- list(
      dp_ge100 = nrow(vcfdat[cfdna_total_dp >= 100]) / nrow(vcfdat[cfdna_total_dp >= 80]),
      dp_ge150 = nrow(vcfdat[cfdna_total_dp >= 150]) / nrow(vcfdat[cfdna_total_dp >= 80]),
      seq_error_rate = seq_error_rate, calc_err_snpnum = nrow(data2a_flt) + nrow(data2b_flt),
      ff_mean = mean(ff_type1_flt), ff_median = ff_median,
      #ff_sd = stats::sd(ff_type1_flt), ff_mad = mad(ff_type1_flt),
      ff_snpnum = length(ff_type1_flt),
      ado_rate1 = nrow(data1a[ff < ff_lower_limit & cfdna_alt_dp < 2]) / nrow(data1a),
      ado_rate2 = nrow(data1b[ff < ff_lower_limit & cfdna_ref_dp < 2]) / nrow(data1b)
    )

    # check mendel errors (paternity)
    if (rel_mode == "PAHP") {
      mendel_errors_count <-
        nrow(data2a[father_genotype == "0/0" & proband_genotype %in% c("0/1", "1/1")]) + 
        nrow(data2b[father_genotype == "1/1" & proband_genotype %in% c("0/1", "0/0")]) +
        nrow(data1a[father_genotype == "1/1" & proband_genotype %in% c("0/0", "1/1")]) +
        nrow(data1b[father_genotype == "0/0" & proband_genotype %in% c("0/0", "1/1")])
      mendel_errors_denom <-
        nrow(data2a[father_genotype == "0/0" & proband_genotype != "./."]) +
        nrow(data2b[father_genotype == "1/1" & proband_genotype != "./."]) +
        nrow(data1a[father_genotype == "1/1" & proband_genotype != "./."]) +
        nrow(data1b[father_genotype == "0/0" & proband_genotype != "./."])
      ret_info$proband_mendel_error_ratio <- mendel_errors_count / mendel_errors_denom
      ret_info$proband_pat_chrX_mendel_error_ratio <- NA
      ret_info$proband_mat_chrX_mendel_error_ratio <- NA

      if (! is.null(query_dat_chrX)) {
        query_dat_chrX <- query_dat_chrX[proband_genotype != "./."]
        proband_chrX_genotab <- table(query_dat_chrX$proband_genotype)
        proband_chrX_hets <- proband_chrX_genotab["0/1"] / sum(proband_chrX_genotab)
        if (proband_chrX_hets < 0.2) {
          #proband_gender <- "male"
          ret_info$proband_pat_chrX_mendel_error_ratio <- NA
        } else {
          #proband_gender <- "female"
          pat_chrX_mendel_errors_count <-
            nrow(query_dat_chrX[father_genotype == "0/0" & proband_genotype == "1/1"]) +
            nrow(query_dat_chrX[father_genotype == "1/1" & proband_genotype == "0/0"])
          pat_chrX_mendel_errors_denom <-
            nrow(query_dat_chrX[father_genotype == "0/0"]) + nrow(query_dat_chrX[father_genotype == "1/1"])
          ret_info$proband_pat_chrX_mendel_error_ratio <- pat_chrX_mendel_errors_count / pat_chrX_mendel_errors_denom
        }
        mat_chrX_mendel_errors_count <-
          nrow(query_dat_chrX[M_pRef > hom_upper_limit & proband_genotype == "1/1"]) +
          nrow(query_dat_chrX[M_pRef < hom_lower_limit & proband_genotype == "0/0"])
        mat_chrX_mendel_errors_denom <- 
          nrow(query_dat_chrX[M_pRef > hom_upper_limit]) + nrow(query_dat_chrX[M_pRef < hom_lower_limit])
        ret_info$proband_mat_chrX_mendel_error_ratio <- mat_chrX_mendel_errors_count / mat_chrX_mendel_errors_denom
      }
    } else if (rel_mode == "GAHP")  {
      ret_info$pat_grandpa_mendel_error_ratio <- NA
      ret_info$pat_grandma_mendel_error_ratio <- NA
      ret_info$mat_grandpa_mendel_error_ratio <- NA
      ret_info$mat_grandma_mendel_error_ratio <- NA
      ret_info$mat_grandpa_chrX_mendel_error_ratio <- NA
      ret_info$mat_grandma_chrX_mendel_error_ratio <- NA

      if ("grandpa_pat_genotype" %in% colnames(query_dat_autosomal)) {
        pat_grandpa_mendel_errors_count <-
          nrow(query_dat_autosomal[father_genotype == "0/0" & grandpa_pat_genotype == "1/1"]) +
          nrow(query_dat_autosomal[father_genotype == "1/1"& grandpa_pat_genotype == "0/0"])
        pat_grandpa_mendel_errors_denom <-
          nrow(query_dat_autosomal[father_genotype == "0/0" & grandpa_pat_genotype != "./."]) +
          nrow(query_dat_autosomal[father_genotype == "1/1" & grandpa_pat_genotype != "./."])
        ret_info$pat_grandpa_mendel_error_ratio <- pat_grandpa_mendel_errors_count / pat_grandpa_mendel_errors_denom
      }
      if ("grandma_pat_genotype" %in% colnames(query_dat_autosomal)) {
        pat_grandma_mendel_errors_count <-
          nrow(query_dat_autosomal[father_genotype == "0/0" & grandma_pat_genotype == "1/1"]) +
          nrow(query_dat_autosomal[father_genotype == "1/1" & grandma_pat_genotype == "0/0"])
        pat_grandma_mendel_errors_denom <-
          nrow(query_dat_autosomal[father_genotype == "0/0" & grandma_pat_genotype != "./."]) +
          nrow(query_dat_autosomal[father_genotype == "1/1" & grandma_pat_genotype != "./."])
        ret_info$pat_grandma_mendel_error_ratio <- pat_grandma_mendel_errors_count / pat_grandma_mendel_errors_denom
      }

      if ("grandpa_mat_genotype" %in% colnames(query_dat_autosomal)) {
        mat_grandpa_mendel_errors_count <-
          nrow(query_dat_autosomal[M_pRef > hom_upper_limit & grandpa_mat_genotype == "1/1"]) +
          nrow(query_dat_autosomal[M_pRef < hom_lower_limit & grandpa_mat_genotype == "0/0"])
        mat_grandpa_mendel_errors_denom <-
          nrow(query_dat_autosomal[M_pRef > hom_upper_limit & grandpa_mat_genotype != "./."]) +
          nrow(query_dat_autosomal[M_pRef < hom_lower_limit & grandpa_mat_genotype != "./."])
        ret_info$mat_grandpa_mendel_error_ratio <- mat_grandpa_mendel_errors_count / mat_grandpa_mendel_errors_denom

        if (! is.null(query_dat_chrX)) {
          chrX_mat_grandpa_mendel_errors_count <-
            nrow(query_dat_chrX[M_pRef > hom_upper_limit & grandpa_mat_genotype == "1/1"]) +
            nrow(query_dat_chrX[M_pRef < hom_lower_limit & grandpa_mat_genotype == "0/0"])
          chrX_mat_grandpa_mendel_errors_denom <-
            nrow(query_dat_chrX[M_pRef > hom_upper_limit & grandpa_mat_genotype != "./."]) +
            nrow(query_dat_chrX[M_pRef < hom_lower_limit & grandpa_mat_genotype != "./."])
          ret_info$mat_grandpa_chrX_mendel_error_ratio <- chrX_mat_grandpa_mendel_errors_count / chrX_mat_grandpa_mendel_errors_denom
        }
      }
      if ("grandma_mat_genotype" %in% colnames(query_dat_autosomal)) {
        mat_grandma_mendel_errors_count <-
          nrow(query_dat_autosomal[M_pRef > hom_upper_limit & grandma_mat_genotype == "1/1"]) +
          nrow(query_dat_autosomal[M_pRef < hom_lower_limit & grandma_mat_genotype == "0/0"])
        mat_grandma_mendel_errors_denom <-
           nrow(query_dat_autosomal[M_pRef > hom_upper_limit & grandma_mat_genotype != "./."]) +
          nrow(query_dat_autosomal[M_pRef < hom_lower_limit & grandma_mat_genotype != "./."])
        ret_info$mat_grandma_mendel_error_ratio <- mat_grandma_mendel_errors_count / mat_grandma_mendel_errors_denom

        if (! is.null(query_dat_chrX)) {
          chrX_mat_grandma_mendel_errors_count <-
            nrow(query_dat_chrX[M_pRef > hom_upper_limit & grandma_mat_genotype == "1/1"]) +
            nrow(query_dat_chrX[M_pRef < hom_lower_limit & grandma_mat_genotype == "0/0"])
          chrX_mat_grandma_mendel_errors_denom <-
            nrow(query_dat_chrX[M_pRef > hom_upper_limit & grandma_mat_genotype != "./."]) +
            nrow(query_dat_chrX[M_pRef < hom_lower_limit & grandma_mat_genotype != "./."])
          ret_info$mat_grandma_chrX_mendel_error_ratio <- chrX_mat_grandma_mendel_errors_count / chrX_mat_grandma_mendel_errors_denom
        }
      }
    }

    # generate a depth filter for proceeding robust calculation
    # reads_cutoff <-
    #   round(unname(quantile(
    #     c(data2a$cfdna_total_dp, data2b$cfdna_total_dp), 0.05
    #   )), 0)

    return(ret_info)
  }

#'@title estimate_fetal_gender
#'@description Estimate global Sequencing Error Ratio from VCF data
estimate_fetal_gender <- function(fetal_fraction, cfdna_y_ratio) {
  # predefined gender ratio for NIPDv2 panel
  known_male_y_mean <- 0.000545
  known_male_y_sd <- 2.53e-5
  known_female_y_mean <- 1.41e-6
  known_female_y_sd <- 3.38e-7

  fetal_y_ratio <- cfdna_y_ratio / fetal_fraction
  Zscore_male <-
    (fetal_y_ratio - known_male_y_mean) / known_male_y_sd
  Zscore_female <-
    (fetal_y_ratio - known_female_y_mean) / known_female_y_sd

  return(ifelse(abs(Zscore_male) < abs(Zscore_female), "male", "female"))
}
