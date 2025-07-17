library(FarmCPUpp)
library(optparse)
library(vcfR)
library(bigmemory)

#  reads CSV/TSV
read_input_file <- function(filename) {
  if (grepl("\\.tsv$", filename, ignore.case = TRUE)) {
    return(read.delim(filename, header = TRUE))
  } else if (grepl("\\.csv$", filename, ignore.case = TRUE)) {
    return(read.csv(filename, header = TRUE))
  } else {
    stop("Unsupported file format: ", filename)
  }
}

#  creates output folder if it doesn't exist. this can be changed!
if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

#  arguments for CLI input.
option_list <- list(
  make_option("--pheno", type="character", help="Phenotype CSV file (required)"),
  make_option("--vcf",   type="character", help="VCF file (optional)"),
  make_option("--geno",  type="character", help="Genotype CSV file (optional, overrides VCF)"),
  make_option("--map",   type="character", help="Map CSV file (optional, overrides VCF)")
)
opt <- parse_args(OptionParser(option_list = option_list))

#  phenotpye block
if (is.null(opt$pheno)) stop("Error: --pheno is required.")
cat("Loading phenotype:", opt$pheno, "\n")
myY <- read_input_file(opt$pheno)

#  genotype block
if (!is.null(opt$geno) && !is.null(opt$map)) {
  cat("Using provided GD and GM CSV files\n")
  myGD_df <- read_input_file(opt$geno)
  myGM <- read_input_file(opt$map)

  # numeric Check
  myGM$Chromosome <- as.numeric(as.character(myGM$Chromosome))
  if (any(is.na(myGM$Chromosome))) {
    stop("Error: Non-numeric values found in Chromosome column of map file.")
  }

} else if (!is.null(opt$vcf)) {
  cat("Parsing VCF:", opt$vcf, "\n")
  vcf <- read.vcfR(opt$vcf, verbose = FALSE)
  gt_raw <- extract.gt(vcf, element = "GT")

  geno_numeric <- apply(gt_raw, 2, function(genos) {
    sapply(genos, function(gt) {
      if (gt %in% c("0/0","0|0")) return(0)
      else if (gt %in% c("0/1","1/0","0|1","1|0")) return(1)
      else if (gt %in% c("1/1","1|1")) return(2)
      else return(NA)
    })
  })
  geno_numeric <- t(geno_numeric)

  # this is for genotypic data and map files
  myGD_df <- data.frame(Taxa = rownames(geno_numeric), geno_numeric, check.names = FALSE)
  myGM <- data.frame(
    SNP = vcf@fix[, "ID"],
    Chromosome = as.numeric(as.character(vcf@fix[, "CHROM"])),
    Position = as.integer(vcf@fix[, "POS"])
  )

  # Checker, just in case!
  if (any(is.na(myGM$Chromosome))) {
    stop("Error: Non-numeric chromosome values in VCF. Clean or rename chromosomes before running.")
  }

  cat("VCF converted:", nrow(myGD_df), "taxa ×", ncol(myGD_df)-1, "markers\n")

} else {
  stop("Error: Provide either --vcf or both --geno and --map.")
}

#  Taxa validation. Ran into an issue with this once. Here just in case!
if (!all(myY$Taxa %in% myGD_df$Taxa)) {
  warning("Some taxa in phenotype not found in genotype. Consider checking for mismatches.")
}

#  FarmCpupp needs bigmatrix, here's a conversion.
cat("Converting GD to big.matrix...\n")
geno_matrix <- as.matrix(myGD_df[, -1])
myGD <- as.big.matrix(geno_matrix)
options(bigmemory.allow.dimnames = TRUE)
rownames(myGD) <- myGD_df$Taxa

# GWAS block
cat("Running FarmCPUpp\n")
result <- FarmCPUpp::farmcpu(Y = myY, GD = myGD, GM = myGM)

library(ggplot2) #These load plotting libraries
library(qqman)

# Extract results from data
trait_name <- names(result)[1]
gwas_df <- result[[trait_name]]$GWAS

if (!is.null(gwas_df) && nrow(gwas_df) > 0) {
  write.csv(gwas_df, paste0("results/FarmCPUpp_", trait_name, ".csv"), row.names = FALSE)
  gwas_clean <- subset(gwas_df, is.finite(p.value) & !is.na(Chromosome) & !is.na(Position))
  if (nrow(gwas_clean) == 0) {
    cat("⚠No valid GWAS results to plot.\n")
  } else {
    # Manhattan plot
    png(paste0("results/FarmCPUpp_", trait_name, "_manhattan.png"), width = 1000, height = 600)
    manhattan(gwas_clean, chr = "Chromosome", bp = "Position", snp = "SNP", p = "p.value", 
              main = paste("Manhattan Plot -", trait_name))
    dev.off()

    # QQ plot
    png(paste0("results/FarmCPUpp_", trait_name, "_qq.png"), width = 600, height = 600)
    qq(gwas_clean$p.value, main = paste("QQ Plot -", trait_name))
    dev.off()

    cat("Plots generated for", trait_name, "\n")
  }

} else {
  cat("No GWAS results were returned. Check input data.\n")
}

#  Check results 
#cat("Inspecting result object...\n")
#str(result)  # diagnostic: show result contents

cat("Process complete. Please check ./results folder. \n")

trait_name <- colnames(myY)[2]
gwas_result <- result[[trait_name]]$GWAS

#This is an error message for a blank file!
if (is.null(gwas_result) || nrow(gwas_result) == 0) {
  cat("No GWAS results were returned.\n")
  quit(status = 0)
}
