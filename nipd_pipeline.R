#!/usr/bin/Rscript
library(data.table)
library(doParallel)
library(optparse)

load("/workspace/evan/scripts/nipd/data/genetic_mapinfo.rda")

source("/workspace/evan/scripts/nipd/workflow/nipd_vcf_parser.R")
source("/workspace/evan/scripts/nipd/workflow/nipd_calc_ff_era.R")
source("/workspace/evan/scripts/nipd/workflow/nipd_hap_phasing.R")

source("/workspace/evan/scripts/nipd/workflow/nipd_infer_fetal_geno.R")

source("/workspace/evan/scripts/nipd/workflow/nipd_report_utils.R")


## ************************ Part 0. Predefined settings ************************
option_list <- list(
    make_option(c("-i", "--inputVcf"), type="character", default=NULL, action="store", help="combined VCF file"),
    make_option("--ped", type="character", default="fam.ped", action="store", help="pedigree information (at least 5 cols)"),
    make_option(c("-r", "--roi"), type="character", default=NULL, action="store", help="analysis region tsv (2 cols)"),
    make_option("--cutoff", type="integer", default=80, action="store", help="depth cutoff for cfdna"),
    make_option(c("-d","--dist"), type="integer", default=200, action="store", help="minimum inter-variant SNP"),
    make_option(c("-t", "--threads"), type="integer", default=4, action="store", help="num of parallelizing threads (max: 8)"),
    make_option(c("-o", "--outdir"), type="character", default="nipd_analysis", action="store", help="output directory"),
    make_option("--depth_dir", type="character", default="depth", action="store", help="depth directory"),
    make_option("--cfvcf_dir", type="character", default="mutect2_VCF", action="store", help="cfdna VCF directory"),
    make_option("--txt", default=FALSE, action="store_true", help="set output format in txt"),
    make_option("--comprehensive", default=FALSE, action="store_true", help="output comprehensive instead of simplified results") 
)
opt <- parse_args(OptionParser(option_list = option_list, 
    usage="This script performs noninvasive prenatal diagnosis of fetal genotypes in cfDNA."))

# variant filter criteria
mincov_cfdna <- opt$cutoff
stopifnot(mincov_cfdna >= 50)

## **************************** Part 1. read input ***************************
# read family pedigree information from an at least 5-column text
# famid id fid mid sex (index)
pedinfo <- fread(opt$ped)
stopifnot(ncol(pedinfo) == 6)

# 2-column text: famid region
if (! is.null(opt$roi)) {
    roi_regions <- fread(opt$roi)
    stopifnot(ncol(roi_regions) == 2)
    colnames(roi_regions) <- c("famid", "region")
}

# other configuration
stopifnot(opt$dist >= 10 & opt$dist < 5e3)

outdir <- opt$outdir
if (! dir.exists(outdir)) {
    retval <- dir.create(outdir)
    stopifnot(retval == TRUE)
}

# read input vcf
vcf_dat <- fread(opt$inputVcf, skip = "#CHROM")
vcfdat_flt <- vcf_dat[QUAL > 100 & (!grepl(ALT, pattern = ",", fixed = T)) &
    nchar(REF) <= 6 & nchar(ALT) <= 6]
ref_chroms <- paste0("chr", c(1:22, 'X'))
vcfdat_flt <- vcfdat_flt[ref_chroms, on = "#CHROM"]

# register clusters for parallelization
ncores <- opt$threads
if (ncores > 8) ncores <- 8
if (ncores < 2) ncores <- 2
cl <- makeCluster(ncores)
registerDoParallel(cl)

# parse vcf data for each family
if (is.null(opt$roi)) {
    family_genoinfo <- foreach(i=unique(pedinfo$famid), 
        .packages=c("pedtools", "data.table")) %dopar% 
        parse_family_vcf(
            vcfdat = vcfdat_flt, pedigree_dat = pedinfo[famid==i],
            mincov_cf = mincov_cfdna, cfvcf_dir = opt$cfvcf_dir
        )
} else {
    family_genoinfo <- foreach(i=unique(pedinfo$famid), 
        .packages=c("pedtools", "data.table")) %dopar% 
        parse_family_vcf(
            vcfdat = vcfdat_flt, pedigree_dat = pedinfo[famid==i],
            roi = roi_regions[famid==i, region],
            mincov_cf = mincov_cfdna, cfvcf_dir = opt$cfvcf_dir
        )
}
# remove NULL items if exists
family_genoinfo <- family_genoinfo[!sapply(family_genoinfo, is.null)]

## ******** Part 2. calculate fetal fraction, error rate and paternity ********
# vcf filtering and prenatal characteristic calculation
depth_dir <- ifelse(dir.exists(opt$depth_dir), opt$depth_dir, NULL)
family_genoinfo_flt <- foreach(i=seq_along(family_genoinfo), 
    .packages=c("data.table")) %dopar%
    preprocess_nipd_info(family_genoinfo[[i]], coverage_dir = depth_dir)

rm(family_genoinfo); rm(vcf_dat)

## ******************* Part 3. parental haplotype construction ****************
# construct parental haplotypes based on Mendel's law
family_genoinfo_phased <- foreach(i=seq_along(family_genoinfo_flt), 
    .packages=c("data.table")) %dopar%
    construct_parental_haplotypes(family_genoinfo_flt[[i]])

rm(family_genoinfo_flt)

## ********************* Part 4. fetal genotype deduction *********************
# deduct fetal haplotypes
family_genoinfo_inferred <- foreach(i=seq_along(family_genoinfo_phased), 
    .packages=c("data.table", "DNAcopy")) %dopar% {
        fam_genodat <- infer_fetal_pat_genotypes(family_genoinfo_phased[[i]], dp_thresh = mincov_cfdna)
        fam_genodat <- infer_fetal_mat_genotypes(
            fam_genodat, dp_thresh = mincov_cfdna, use_type5 = TRUE, biased_ref_ratio = 51.4
        )
        fam_genodat
    }
#    infer_fetal_genotypes(family_genoinfo_phased[[i]], dist = opt$dist)

rm(family_genoinfo_phased)

## ******************** Part 5. write summary & print output ********************
foreach(i=seq_along(family_genoinfo_inferred), .packages=c("data.table", "openxlsx")) %dopar% {
    summarise_fetal_results(family_genoinfo_inferred[[i]], outdir, txt_format = opt$txt, unsimplified = opt$comprehensive)
    #summarise_fetal_results(family_genoinfo_inferred[[i]], outdir, txt_format = TRUE)
    #summarise_fetal_results(family_genoinfo_inferred[[i]], outdir, txt_format = FALSE)
    plot_fetal_haplotypes(family_genoinfo_inferred[[i]], outdir)
}

# stop registered clusters
stopCluster(cl)
