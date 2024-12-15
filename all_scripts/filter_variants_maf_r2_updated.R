if (!"tidyverse" %in% rownames(installed.packages())) {
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
}
if (!"data.table" %in% rownames(installed.packages())) {
  install.packages("data.table", repos = "https://cloud.r-project.org")
}

library(tidyverse)
library(data.table)

arg <- commandArgs(trailingOnly = T)

# define function we will use
filter_info <- function (chr, maf, rsq) {
  # Read the info file of specific chromosome
  data <- paste0("chr", chr,".info") %>% fread()
  # Filter it to maf and rsq
  dat <- data[data$MAF >= maf & data$Rsq >= rsq]
  # Generate chromosome, bp, and range from the "SNP" column
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  dat$range <- paste0(dat$chr, ":", dat$bp, "-", dat$bp)
  # writing files
  dat[,c("SNP","ALT_Frq","Rsq")] %>% fwrite(
    paste0("maf", gsub("0\\.", "", maf), "rsq", gsub("0\\.", "", rsq), "minimums_chr", chr, ".info"),
    row.names = F, quote = F, sep = "\t"
  )
  dat[,c("range")] %>% fwrite(
    paste0("maf", gsub("0\\.", "", maf), "rsq", gsub("0\\.", "", rsq), "minimums_chr", chr, ".txt"),
    col.names = F, row.names = F, quote = F, sep = "\t"
  )
  # report number of SNPs that pass the filter
  paste0("Number of SNPs in chromosome ", chr, " after filter (MAF >= " , maf, ", Rsq >= ", rsq, "): ", nrow(dat)) %>% print()
}

# SET MAF and RSQ filters here 

lapply(1:22, filter_info, maf = arg[1], rsq = arg[2])