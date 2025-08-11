#'@title construct_parental_haplotypes
#'@description Construct parental haplotypes from family-based genotype data
construct_parental_haplotypes <-
  function(fam_genodat, hom_lower_limit = 10, hom_upper_limit = 90,
    het_lower_limit = 35, het_upper_limit = 65, min_pat_alt_dp = 6
  ) {
    # check uninformative reads
    fam_genodat$vcfdat[, tot_rel := (cfdna_ref_dp + cfdna_alt_dp) / cfdna_total_dp]
    vcfdat <- fam_genodat$vcfdat[tot_rel > 0.9][, tot_rel := NULL]

    pat_snps <- vcfdat[(chrom != "chrX") &
      (F_pRef <= het_upper_limit & F_pRef >= het_lower_limit) & (father_genotype == "0/1") & (father_alt_dp > min_pat_alt_dp)]
    mat_snps <- vcfdat[
      (M_pRef <= het_upper_limit & M_pRef >= het_lower_limit) &
      (F_pRef >= hom_upper_limit | F_pRef <= hom_lower_limit) &
      (father_genotype %in% c("0/0", "1/1"))
    ]
    if (fam_genodat$has_mom_geno == FALSE) {
      mat_snps <- mat_snps[cfdna_genotype != "./."]
    } else {
      mat_snps <- mat_snps[mother_genotype != "./." & cfdna_genotype != "./."]
    }
    # remove PAR region on chrX
    if ("chrX" %in% unique(mat_snps$chrom)) {
      mat_snps <- mat_snps[! (chrom == "chrX" & (pos < 2781479 | pos > 156030895))]
    }
    mat_snps[, F0 := data.table::fifelse(F_pRef >= hom_upper_limit, 0, 1)][, F1 := F0]

    if (fam_genodat$sampmode == "PAHP") {
      pat_snps <- PAHP_paternal_phasing(pat_snps)
      if (is.na(fam_genodat$fetal_gender)) {
        mat_snps <- PAHP_maternal_phasing(mat_snps)
      }
      else {
        offsprings <- pedtools::leaves(fam_genodat$pedigree)
        if (fam_genodat$has_mom_geno == FALSE) {
          mat_snps <- PAHP_maternal_phasing(
            mat_snps, pedtools::getSex(fam_genodat$pedigree, offsprings[1])
          )
        } else {
          proband_smpl <- setdiff(offsprings, fam_genodat$sampname)
          mat_snps <- PAHP_maternal_phasing(
            mat_snps, pedtools::getSex(fam_genodat$pedigree, proband_smpl)
          )
        }
      }
    }
    else if (fam_genodat$sampmode == "GAHP") {
      n_pat_gr <- length(grep("^pat_gr", names(fam_genodat)))
      n_mat_gr <- length(grep("^mat_gr", names(fam_genodat)))
      if (n_pat_gr > 0) {
        pat_snps <- GAHP_paternal_phasing(
          pat_snps, n_pat_gr, fam_genodat$index_grandparents
        )
      } else {
        pat_snps <- NULL
      }
      if (n_mat_gr > 0) {
        mat_snps <- GAHP_maternal_phasing(
          mat_snps, n_mat_gr, fam_genodat$index_grandparents
        )
      } else {
        mat_snps <- NULL
      }
      if (length(fam_genodat$index_grandparents) > 1) {
        pat_snps_hetXhet <- pat_snps[snp_type == "type5"]
        pat_snps_hetXhet <- GAHP_maternal_phasing(
          pat_snps_hetXhet, n_mat_gr, fam_genodat$index_grandparents
        )
        pat_snps <- rbind(pat_snps[snp_type == "type3"], pat_snps_hetXhet)
        pat_snps <- pat_snps[order(chrom, pos)]
      }

    }
    fam_genodat$pat_snps <- pat_snps
    fam_genodat$mat_snps <- mat_snps
    return(fam_genodat)
  }

#'@title PAHP_paternal_phasing
#'@description Construct paternal haplotypes from proband information
PAHP_paternal_phasing <- function(pat_snps, hom_lower_limit = 10, hom_upper_limit = 90,
  het_lower_limit = 40, het_upper_limit = 65) {
  pat_snps[, proband_pRef := round(100 * proband_ref_dp / proband_total_dp, 1)]
  pat_snps <- pat_snps[
    (proband_genotype == "0/0" & proband_pRef > hom_upper_limit) |
    (proband_genotype == "1/1" & proband_pRef < hom_lower_limit) |
    (proband_genotype == "0/1" & proband_pRef > het_lower_limit & proband_pRef < het_upper_limit)
  ]

  pat_snps_hetXhom <- pat_snps[M_pRef >= hom_upper_limit | M_pRef <= hom_lower_limit]
  pat_snps_hetXhom[, M0 := data.table::fifelse(M_pRef >= hom_upper_limit, 0, 1)][, M1 := M0]
  pat_snps_hetXhom[, F0 := data.table::fcase(
    proband_genotype == "0/0", 0,
    proband_genotype == "1/1", 1,
    proband_genotype == "0/1", data.table::fifelse(M0 == 0, 1, 0)
  )][, F1 := 1 - F0]
  # remove mendel error snps
  pat_snps_hetXhom <- pat_snps_hetXhom[! (
    (proband_genotype == "0/0" & M0 == 1) | (proband_genotype == "1/1" & M0 == 0)
  )]

  pat_snps_hetXhet <- pat_snps[
    (M_pRef <= het_upper_limit & M_pRef >= het_lower_limit) & (proband_genotype %in% c("0/0", "1/1"))
  ]
  pat_snps_hetXhet[, F0 := data.table::fifelse(proband_genotype == "0/0", 0, 1)][, F1 := 1 - F0]
  pat_snps_hetXhet[, ':=' (M0 = F0, M1 = F1)]

  pat_snps <- rbind(
    pat_snps_hetXhom[, ':=' (proband_pRef = NULL, snp_type = "type3")],
    pat_snps_hetXhet[, ':=' (proband_pRef = NULL, snp_type = "type5")]
  )
  return(pat_snps[order(chrom, pos)])
}

#'@title PAHP_maternal_phasing
#'@description Construct maternal haplotypes from proband information
PAHP_maternal_phasing <- function(mat_snps, proband_gender = 0, hom_lower_limit = 10, hom_upper_limit = 90,
  het_lower_limit = 40, het_upper_limit = 65) {
  mat_snps[, proband_pRef := round(100 * proband_ref_dp / proband_total_dp, 1)]
  mat_snps <- mat_snps[
    (proband_genotype == "0/0" & proband_pRef > hom_upper_limit) |
    (proband_genotype == "1/1" & proband_pRef < hom_lower_limit) |
    (proband_genotype == "0/1" & proband_pRef > het_lower_limit & proband_pRef < het_upper_limit)
  ]
  mat_snps[, M0 := data.table::fcase(
    proband_genotype == "0/0", 0,
    proband_genotype == "1/1", 1,
    proband_genotype == "0/1", data.table::fifelse(F0 == 0, 1, 0)
  )][, M1 := 1 - M0]
  mat_snps <- mat_snps[! (
    ((chrom != "chrX") & (proband_genotype == "0/0" & F0 == 1) | (proband_genotype == "1/1" & F0 == 0))
  )]

  if ("chrX" %in% unique(mat_snps$chrom)) {
    if (proband_gender == 1) {
      mat_snps <- mat_snps[! (chrom == "chrX" & proband_genotype == "0/1")]
    }
    else if (proband_gender == 2) {
      mat_snps <- mat_snps[! (chrom == "chrX" &  (
        (proband_genotype == "0/0" & F0 == 1) | (proband_genotype == "1/1" & F0 == 0)
      ))]
    }
  }
  return(mat_snps[, proband_pRef := NULL])
}

#'@title GAHP_paternal_phasing
#'@description Construct paternal haplotypes from grandparent(s)' information
GAHP_paternal_phasing <- function(pat_snps, n_gr, idx_gr, hom_lower_limit = 10, hom_upper_limit = 90,
  het_lower_limit = 40, het_upper_limit = 65) {

  if (n_gr == 2) {
    # this family includes both paternal parents
    pat_snps <- pat_snps[! (grandpa_pat_genotype == "0/1" & grandma_pat_genotype == "0/1")]
    pat_snps[, grandpa_pRef := round(100 * grandpa_pat_ref_dp / grandpa_pat_total_dp, 1)][,
      grandma_pRef := round(100 * grandma_pat_ref_dp / grandma_pat_total_dp, 1)]
    if ("pat_grandpa" %in% idx_gr) {
      pat_snps <- pat_snps[
        (grandpa_pat_genotype == "0/0" & grandpa_pRef > hom_upper_limit) |
        (grandpa_pat_genotype == "1/1" & grandpa_pRef < hom_lower_limit) |
        (grandpa_pat_genotype == "0/1" & grandpa_pRef > het_lower_limit & grandpa_pRef < het_upper_limit)
      ]
      pat_snps[, F0 := data.table::fcase(
        grandpa_pat_genotype == "0/0", 0,
        grandpa_pat_genotype == "1/1", 1,
        grandpa_pat_genotype == "0/1",
        data.table::fifelse(grandma_pat_genotype == "0/0", 1, 0)
      )][, F1 := 1 - F0]
    }
    else {
      pat_snps <- pat_snps[
        (grandma_pat_genotype == "0/0" & grandma_pRef > hom_upper_limit) |
        (grandma_pat_genotype == "1/1" & grandma_pRef < hom_lower_limit) |
        (grandma_pat_genotype == "0/1" & grandma_pRef > het_lower_limit & grandma_pRef < het_upper_limit)
      ]
      pat_snps[, F0 := data.table::fcase(
        grandma_pat_genotype == "0/0", 0,
        grandma_pat_genotype == "1/1", 1,
        grandma_pat_genotype == "0/1",
        data.table::fifelse(grandpa_pat_genotype == "0/0", 1, 0)
      )][, F1 := 1 - F0]
    }
    pat_snps[, ':=' (grandpa_pRef = NULL, grandma_pRef = NULL)]
  }
  else if (n_gr == 1) {
    # this family includes only one of the paternal grandparents
    if ("grandpa_pat_genotype" %in% colnames(pat_snps)) {
      pat_snps[, grandpa_pRef := round(100 * grandpa_pat_ref_dp / grandpa_pat_total_dp, 1)]
      pat_snps <- pat_snps[
        (grandpa_pat_genotype == "0/0" & grandpa_pRef > hom_upper_limit) |
        (grandpa_pat_genotype == "1/1" & grandpa_pRef < hom_lower_limit)
      ]
      pat_snps[, F0 := data.table::fifelse(grandpa_pat_genotype == "0/0", 0, 1)][, F1 := 1 - F0]
      pat_snps[, grandpa_pRef := NULL]
    } else {
      pat_snps[, grandma_pRef := round(100 * grandma_pat_ref_dp / grandma_pat_total_dp, 1)]
      pat_snps <- pat_snps[
        (grandma_pat_genotype == "0/0" & grandma_pRef > hom_upper_limit) |
        (grandma_pat_genotype == "1/1" & grandma_pRef < hom_lower_limit)
      ]
      pat_snps[, F0 := data.table::fifelse(grandma_pat_genotype == "0/0", 0, 1)][, F1 := 1 - F0]
      pat_snps[, grandma_pRef := NULL]
    }
  }

  pat_snps_hetXhom <- pat_snps[M_pRef >= hom_upper_limit | M_pRef <= hom_lower_limit]
  pat_snps_hetXhom[, M0 := data.table::fifelse(M_pRef >= hom_upper_limit, 0, 1)][, M1 := M0]

  pat_snps_hetXhet <- pat_snps[(M_pRef <= het_upper_limit & M_pRef >= het_lower_limit)]
  pat_snps_hetXhet[, ':=' (M0 = -1, M1 = -1)]

  pat_snps <- rbind(
    pat_snps_hetXhom[, snp_type := "type3"], pat_snps_hetXhet[, snp_type := "type5"]
  )
  return(pat_snps[order(chrom, pos)])
}

#'@title GAHP_maternal_phasing
#'@description Construct maternal haplotypes from grandparent(s)' information
GAHP_maternal_phasing <- function(mat_snps, n_gr, idx_gr, hom_lower_limit = 10, hom_upper_limit = 90,
  het_lower_limit = 40, het_upper_limit = 65) {
  if (n_gr == 2) {
    # this family includes both maternal parents
    mat_snps <- mat_snps[! (grandpa_mat_genotype == "0/1" & grandma_mat_genotype == "0/1")]
    mat_snps[, grandpa_pRef := round(100 * grandpa_mat_ref_dp / grandpa_mat_total_dp, 1)][,
      grandma_pRef := round(100 * grandma_mat_ref_dp / grandma_mat_total_dp, 1)]

    if ("mat_grandpa" %in% idx_gr) {
      mat_snps <- mat_snps[
        (grandpa_mat_genotype == "0/0" & grandpa_pRef > hom_upper_limit) |
        (grandpa_mat_genotype == "1/1" & grandpa_pRef < hom_lower_limit) |
        (grandpa_mat_genotype == "0/1" & grandpa_pRef > het_lower_limit & grandpa_pRef < het_upper_limit)
      ]
      mat_snps <- mat_snps[! (chrom == "chrX" & grandpa_mat_genotype == "0/1")]
      mat_snps[, M0 := data.table::fcase(
        grandpa_mat_genotype == "0/0", 0,
        grandpa_mat_genotype == "1/1", 1,
        grandpa_mat_genotype == "0/1", data.table::fifelse(grandma_mat_genotype == "0/0", 1, 0)
      )][, M1 := 1 - M0]
    } else {
      mat_snps <- mat_snps[
        (grandma_mat_genotype == "0/0" & grandma_pRef > hom_upper_limit) |
        (grandma_mat_genotype == "1/1" & grandma_pRef < hom_lower_limit) |
        (grandma_mat_genotype == "0/1" & grandma_pRef > het_lower_limit & grandma_pRef < het_upper_limit)
      ]
      mat_snps[, M0 := data.table::fcase(
        grandma_mat_genotype == "0/0", 0,
        grandma_mat_genotype == "1/1", 1,
        grandma_mat_genotype == "0/1", data.table::fifelse(grandpa_mat_genotype == "0/0", 1, 0)
      )][, M1 := 1 - M0]
    }
    
    return(mat_snps[, ':=' (grandpa_pRef = NULL, grandma_pRef = NULL)])
  }
  else if (n_gr == 1) {
    # this family includes only one of the maternal grandparents
    if ("grandpa_mat_genotype" %in% colnames(mat_snps)) {
      mat_snps[, grandpa_pRef := round(100 * grandpa_mat_ref_dp / grandpa_mat_total_dp, 1)]
      mat_snps <- mat_snps[
        (grandpa_mat_genotype == "0/0" & grandpa_pRef > hom_upper_limit) |
        (grandpa_mat_genotype == "1/1" & grandpa_pRef < hom_lower_limit)
      ]
      mat_snps[, M0 := data.table::fifelse(grandpa_mat_genotype == "0/0", 0, 1)][, M1 := 1 - M0]
      mat_snps[, grandpa_pRef := NULL]
    } else {
      mat_snps[, grandma_pRef := round(100 * grandma_mat_ref_dp / grandma_mat_total_dp, 1)]
      mat_snps <- mat_snps[
        (grandma_mat_genotype == "0/0" & grandma_pRef > hom_upper_limit) |
        (grandma_mat_genotype == "1/1" & grandma_pRef < hom_lower_limit)
      ]
      mat_snps[, M0 := data.table::fifelse(grandma_mat_genotype == "0/0", 0, 1)][, M1 := 1 - M0]
      mat_snps[, grandma_pRef := NULL]
    }
    return(mat_snps)
  }
}