#'@title summarise_fetal_results
#'@description Summarise cfDNA characteristics and inferred fetal genotypes
summarise_fetal_results <-
  function(fam_genodat, outdir, txt_format = FALSE, unsimplified = FALSE) {
    famName <- fam_genodat$sampname
    if ("pat_snps_results_chrom" %in% names(fam_genodat)) {
      pat_snps <- fam_genodat$pat_snps_results_chrom
    } else {
      pat_snps <- NULL
    }
    if ("mat_snps_results_chrom" %in% names(fam_genodat)) {
      mat_snps <- fam_genodat$mat_snps_results_chrom
    } else {
      mat_snps <- NULL
    }
    if (! "roi_info" %in% names(fam_genodat)) {
      summary_info <- data.frame(
        sample = famName,
        ff_median = round(fam_genodat$qc_info$ff_median, 4),
        seq_error_rate = round(fam_genodat$qc_info$seq_error_rate, 5),
        fetal_gender = ifelse(
          is.na(fam_genodat$fetal_gender),
          "NA", fam_genodat$fetal_gender
        ),
        pat_snps = ifelse(is.null(pat_snps), 0, nrow(pat_snps)),
        mat_snps = ifelse(is.null(mat_snps), 0, nrow(mat_snps))
      )
    } else {
      summary_info <- rbindlist(lapply(seq_len(nrow(fam_genodat$roi_info)), function(i) {
        chrN <- fam_genodat$roi_info$chrom[i]
        roi_st <- fam_genodat$roi_info$start[i]
        roi_ed <- fam_genodat$roi_info$end[i]
        roi_up4mb <- fam_genodat$roi_info$startpos[i]
        roi_down4mb <- fam_genodat$roi_info$endpos[i]
        list(
          sample = famName,
          ff_median = round(fam_genodat$qc_info$ff_median, 4),
          seq_error_rate = round(fam_genodat$qc_info$seq_error_rate, 5),
          fetal_gender = ifelse(
            is.na(fam_genodat$fetal_gender), "NA", fam_genodat$fetal_gender
          ),
          pat_snps_roi = ifelse(
            is.null(pat_snps), 0, nrow(pat_snps[chrom == chrN & pos >= roi_st & pos <= roi_ed])
          ),
          pat_snps_upstream = ifelse(
            is.null(pat_snps), 0, nrow(pat_snps[chrom == chrN & pos >= roi_up4mb & pos < roi_st])
          ),
          pat_snps_downstream = ifelse(
            is.null(pat_snps), 0, nrow(pat_snps[chrom == chrN & pos > roi_ed & pos < roi_down4mb])
          ),
          mat_snps_roi = ifelse(
            is.null(mat_snps), 0, nrow(mat_snps[chrom == chrN & pos >= roi_st & pos <= roi_ed])
          ),
          mat_snps_upstream = ifelse(
            is.null(mat_snps), 0, nrow(mat_snps[chrom == chrN & pos >= roi_up4mb & pos < roi_st])
          ),
          mat_snps_downstream = ifelse(
            is.null(mat_snps), 0, nrow(mat_snps[chrom == chrN & pos > roi_ed & pos < roi_down4mb])
          )
        )
      }))
    }

    ## output genotyped fetal data
    if (txt_format == TRUE) {
      summary_outfile <-
        paste0(outdir, "/", famName, "_nipd_basic_info.tsv")
      utils::write.table(
        summary_info,
        file = summary_outfile,
        sep = "\t",
        row.names = F,
        quote = F
      )

      pat_outfile1 <-
        paste0(outdir, "/", famName, "_nipd_pat_chrom_results.tsv")
      mat_outfile1 <-
        paste0(outdir, "/", famName, "_nipd_mat_chrom_results.tsv")
      if ("pat_snps_results_chrom" %in% names(fam_genodat)) {
        utils::write.table(
          fam_genodat$pat_snps_results_chrom,
          file = pat_outfile1,
          sep = "\t",
          row.names = F,
          quote = F
        )
      }
      if ("mat_snps_results_chrom" %in% names(fam_genodat)) {
        utils::write.table(
          fam_genodat$mat_snps_results_chrom,
          file = mat_outfile1,
          sep = "\t",
          row.names = F,
          quote = F
        )
      }
    } else {
      out_xlsx <- paste0(outdir, "/", famName, "_nipd_results.xlsx")
      wb <- openxlsx::createWorkbook(creator = "yimingwu")
      openxlsx::addWorksheet(wb, sheetName = "Basic Info")
      openxlsx::writeData(wb, sheet = "Basic Info", x = summary_info)

      if ("pat_snps_results_chrom" %in% names(fam_genodat)) {
        if ("roi_info" %in% names(fam_genodat)) {
          for (i in seq_len(nrow(fam_genodat$roi_info))) {
            roi_st <- fam_genodat$roi_info$start[i]
            roi_ed <- fam_genodat$roi_info$end[i]
            roi_up4mb <- fam_genodat$roi_info$startpos[i]
            roi_down4mb <- fam_genodat$roi_info$endpos[i]
            pat_snps[chrom == fam_genodat$roi_info$chrom[i] & pos >= roi_up4mb & pos <= roi_down4mb,
              is_roi := fcase(
                pos >= roi_st & pos <= roi_ed, "genic",
                pos >= roi_up4mb & pos < roi_st, "upstream",
                pos <= roi_down4mb & pos > roi_ed, "downstream"
              )]
          }
        }

        if (unsimplified == FALSE) {
          pat_snps[, ':=' (
            varid = NULL, cfdna_genotype = NULL,
            F_pRef = NULL, M_pRef = NULL, snp_type = NULL,
            F0_prob = NULL, F1_prob = NULL,
            P0 = NULL, P1 = NULL, P0re = NULL, P1re = NULL,
            naive_fetal_hap = NULL, naive_fetal_genotype = NULL
          )]
          geno_cols <- grep(colnames(pat_snps), pattern = "genotype")
          for (i in seq_len(nrow(pat_snps))) {
            for (j in geno_cols) {
              pat_snps[i,j] <- ifelse(
                pat_snps[i,j, with = FALSE] == "0/0", paste0(pat_snps$ref[i], "/", pat_snps$ref[i]),
                ifelse(
                  pat_snps[i,j, with = FALSE] == "0/1",
                  paste0(pat_snps$ref[i], "/", pat_snps$alt[i]),
                  paste0(pat_snps$alt[i], "/", pat_snps$alt[i])
                )
              )
            }
          }
          pat_snps[, F0 := fifelse(F0 == 0, ref, alt)][, F1 := fifelse(F1 == 0, ref, alt)]
          pat_snps[, M0 := fifelse(M0 == 0, ref, alt)][, M1 := M0]
          if (fam_genodat$has_mom_geno == FALSE) {
            pat_snps[, mother_genotype := paste(M0, M1, sep = "/")]
          }
          geno_names <- c(
            "chrom", "pos", "ref", "alt", grep("genotype$", names(pat_snps), value = TRUE),
            "F0", "F1", "M0", "M1", "hmm_fetal_hap"
          )
          pat_snps <- cbind(
            pat_snps[, geno_names, with = FALSE], pat_snps[, !..geno_names]
          )
        }
        openxlsx::addWorksheet(wb, sheetName = "Paternal SNPs")
        openxlsx::writeData(wb,
                  sheet = "Paternal SNPs",
                  x = pat_snps)
      }


      if ("mat_snps_results_chrom" %in% names(fam_genodat)) {
        if ("roi_info" %in% names(fam_genodat)) {
          for (i in seq_len(nrow(fam_genodat$roi_info))) {
            roi_st <- fam_genodat$roi_info$start[i]
            roi_ed <- fam_genodat$roi_info$end[i]
            roi_up4mb <- fam_genodat$roi_info$startpos[i]
            roi_down4mb <- fam_genodat$roi_info$endpos[i]
            mat_snps[chrom == fam_genodat$roi_info$chrom[i] & pos >= roi_up4mb & pos <= roi_down4mb,
              is_roi := fcase(
                pos >= roi_st & pos <= roi_ed, "genic",
                pos >= roi_up4mb & pos < roi_st, "upstream",
                pos <= roi_down4mb & pos > roi_ed, "downstream"
              )]
          }
        }
        if (unsimplified == FALSE) {
          mat_snps[, ':=' (
            varid = NULL, cfdna_genotype = NULL,
            F_pRef = NULL, M_pRef = NULL, snp_type = NULL,
            M0_prob = NULL, M1_prob = NULL,
            P0 = NULL, P1 = NULL, P0re = NULL, P1re = NULL,
            rhdo_fetal_genotype = NULL, rhdo_blk = NULL, rhdo_blk_type = NULL,
            rhdo_rev_blk = NULL, rhdo_rev_blk_type = NULL, rhdo_rev_fetal_genotype = NULL
          )]

          geno_cols <- grep(colnames(mat_snps), pattern = "genotype")
          for (i in seq_len(nrow(mat_snps))) {
            for (j in geno_cols) {
              mat_snps[i,j] <- ifelse(
                mat_snps[i,j, with = FALSE] == "0/0", paste0(mat_snps$ref[i], "/", mat_snps$ref[i]),
                ifelse(
                  mat_snps[i,j, with = FALSE] == "0/1",
                  paste0(mat_snps$ref[i], "/", mat_snps$alt[i]),
                  paste0(mat_snps$alt[i], "/", mat_snps$alt[i])
                )
              )
            }
          }
          mat_snps[, F0 := fifelse(F0 == 0, ref, alt)][, F1 := F0]
          mat_snps[, M0 := fifelse(M0 == 0, ref, alt)][, M1 := fifelse(M1 == 0, ref, alt)]
          if (fam_genodat$has_mom_geno == FALSE) {
            mat_snps[, mother_genotype := paste(M0, M1, sep = "/")]
          }
          geno_names <- c(
            "chrom", "pos", "ref", "alt", grep("genotype$", names(mat_snps), value = TRUE),
            "F0", "F1", "M0", "M1", "hmm_fetal_hap"
          )
          mat_snps <- cbind(
            mat_snps[, geno_names, with = FALSE], mat_snps[, !..geno_names]
          )
        }
        openxlsx::addWorksheet(wb, sheetName = "Maternal SNPs")
        openxlsx::writeData(wb,
                  sheet = "Maternal SNPs",
                  x = mat_snps)
      }
      openxlsx::saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
    }
  }



#'@title plot_fetal_haplotypes
#'@description Plot fetal haplotypes inheritance
plot_fetal_haplotypes <- function(fam_genodat, outdir) {
  ref_chroms <- paste0("chr", c(1:22, "X"))
  chrom_lengths <- data.table::rbindlist(lapply(ref_chroms, function(chrN)
    utils::tail(genetic_mapinfo[chrom == chrN, .(chrom, pos)], 1)))

  famName <- fam_genodat$sampname

  plot_dt_chroms <- data.table()
  if ("pat_snps_results_chrom" %in% names(fam_genodat)) {
    if ("logOR" %in% colnames(fam_genodat$pat_snps_results_chrom))
      plot_dt_chroms <- rbind(plot_dt_chroms,
                              fam_genodat$pat_snps_results_chrom[, .(chrom, pos, logOR)][, class := "pat"])
  }
  if ("mat_snps_results_chrom" %in% names(fam_genodat)) {
    if ("logOR" %in% colnames(fam_genodat$mat_snps_results_chrom))
      plot_dt_chroms <- rbind(plot_dt_chroms,
                              fam_genodat$mat_snps_results_chrom[, .(chrom, pos, logOR)][, class := "mat"])
  }

  plot_dt_roi <- data.table()
  if ("roi_info" %in% names(fam_genodat) & "pat_snps_results_chrom" %in% names(fam_genodat)) {
    pat_snps_roi <- rbindlist(lapply(seq_len(nrow(fam_genodat$roi_info)), function(i) {
      roi_st <- fam_genodat$roi_info$start[i]
      roi_ed <- fam_genodat$roi_info$end[i]
      roi_up4mb <- fam_genodat$roi_info$startpos[i]
      roi_down4mb <- fam_genodat$roi_info$endpos[i]
      fam_genodat$pat_snps_results_chrom[chrom == fam_genodat$roi_info$chrom[i] & pos >= roi_up4mb & pos <= roi_down4mb]
    }))
    plot_dt_roi <- rbind(plot_dt_roi,
                          pat_snps_roi[, .(chrom, pos, logOR)][, class := "pat"])
  
  }
  if ("roi_info" %in% names(fam_genodat) & "mat_snps_results_chrom" %in% names(fam_genodat)) {
    mat_snps_roi <- rbindlist(lapply(seq_len(nrow(fam_genodat$roi_info)), function(i) {
      roi_st <- fam_genodat$roi_info$start[i]
      roi_ed <- fam_genodat$roi_info$end[i]
      roi_up4mb <- fam_genodat$roi_info$startpos[i]
      roi_down4mb <- fam_genodat$roi_info$endpos[i]
      fam_genodat$mat_snps_results_chrom[chrom == fam_genodat$roi_info$chrom[i] & pos >= roi_up4mb & pos <= roi_down4mb]
    }))
    plot_dt_roi <- rbind(plot_dt_roi,
                          mat_snps_roi[, .(chrom, pos, logOR)][, class := "mat"])
    
  }

  if (nrow(plot_dt_chroms) > 0) {
    plot_logOR(plot_dt_chroms, chrom_lengths, outdir, famName, suffix = "_chrom")
  }

  if (nrow(plot_dt_roi) > 0) {
    plot_logOR(plot_dt_roi, chrom_lengths, outdir, famName, suffix = "_ROI")
  }
}

#'@title plot_logOR
#'@description Plot fetal haplotypes inheritance in logOR
plot_logOR <-
  function(plot_dt,
           chrom_lengths,
           outdir,
           sampname,
           suffix = "") {
    pat_col <- "#BB0021FF"
    mat_col <- "#0073C2FF"

    outplot <-
      paste0(outdir,
             "/",
             sampname,
             "_nipd_fetal_haps",
             suffix,
             "_plot.pdf")
    grDevices::pdf(outplot, height = 8, width = 12)
    chroms <- intersect(chrom_lengths$chrom, unique(plot_dt$chrom))
    fig_row <- ifelse(length(chroms) <= 2, 2, 4)
    graphics::par(mfrow = c(fig_row, 1))

    for (chrN in chroms) {
      plot_dt_chrN <- plot_dt[chrom == chrN]
      if (nrow(plot_dt_chrN) == 0) next
      chrN_dist <- max(plot_dt_chrN$pos) - min(plot_dt_chrN$pos)
      if (chrN_dist > 1e7) {
        xlim_lower <- 0
        xlim_upper <-
          round(chrom_lengths[chrom == chrN, pos] / 1e6) + 1
        dist_scale <- 1e6
      } else {
        xlim_lower <- round(min(plot_dt_chrN$pos) / 1e3) - 1
        xlim_upper <- round(max(plot_dt_chrN$pos) / 1e3) + 1
        dist_scale <- 1e3
      }
      ylim_lower <- -10
      ylim_upper <- 10

      fetal_inh <- unique(plot_dt_chrN$class)
      graphics::par(mar = c(5, 5, 2, 2))
      xaxis_lab <-
        ifelse(dist_scale == 1e6, "Position /Mb", "Position /Kb")

      if (dist_scale == 1e6) {
        graphics::plot(
          round(plot_dt_chrN$pos / 1e6, 2),
          plot_dt_chrN$logOR,
          main = chrN, bty = "l",
          ylab = "Log Odds Ratio (hap0/hap1)",
          xlab = xaxis_lab,
          ylim = c(ylim_lower, ylim_upper),
          xlim = c(xlim_lower, xlim_upper),
          xaxs = "i",
          xaxt = "n",
          cex.lab = 1.2,
          cex.axis = 1,
          font.lab = 1.8,
          type = "n",
          lwd = 2.5
        )
        graphics::axis(1, at = seq(xlim_lower, xlim_upper, 10))
      } else {
        graphics::plot(
          round(plot_dt_chrN$pos / 1e6, 2),
          plot_dt_chrN$logOR,
          main = chrN, bty = "l",
          ylab = "Log Odds Ratio (hap0/hap1)",
          xlab = xaxis_lab,
          ylim = c(ylim_lower, ylim_upper),
          xlim = c(xlim_lower, xlim_upper),
          xaxs = "i",
          cex.lab = 1.2,
          cex.axis = 1,
          font.lab = 1.8,
          type = "n",
          lwd = 2.5
        )
      }

      if ("pat" %in% fetal_inh) {
        graphics::lines(
          round(plot_dt_chrN[class == "pat", pos] / dist_scale, 2),
          plot_dt_chrN[class == "pat", logOR],
          col = pat_col,
          lwd = 2.5
        )
      }
      if ("mat" %in% fetal_inh) {
        graphics::lines(
          round(plot_dt_chrN[class == "mat", pos] / dist_scale, 2),
          plot_dt_chrN[class == "mat", logOR],
          col = mat_col,
          lwd = 2.5
        )
      }

      graphics::abline(
        h = 0,
        col = "gray50",
        lwd = 1.2,
        lty = 2
      )
      graphics::legend(
        "topright",
        legend = c("Paternal allele", "Maternal allele"),
        col = c(pat_col, mat_col),
        lty = c(1, 1),
        lwd = 2,
        bty = "n",
        pt.cex = 2,
        cex = 1,
        inset = 0.01,
        text.col = "black",
        horiz = F
      )
    }

    grDevices::dev.off()
  }
