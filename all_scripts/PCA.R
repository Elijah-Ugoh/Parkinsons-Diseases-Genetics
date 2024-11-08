# PCA Visualiztion
setwd("~/Documents/cosmos-home")

# load tidyverse package
library(tidyverse)

# read in data
pca <- read_table("MPBC_HRC_PCA.eigenvec", col_names = FALSE)
eigenval <- scan("MPBC_HRC_PCA.eigenval")

# sort out the PCA data and remove the first id column
pca_prunned <- pca[,-1]

#set names
names(pca_prunned)[1] <-"Individuals"

# rename the colmnns
names(pca_prunned)[2:ncol(pca_prunned)] <- paste0("PC", 1:(ncol(pca_prunned) - 1))

# save to a text file
write.table(pca_prunned, file = "pca_prunned_updated.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# plot eigenvalues
# first convert to percentage variance explained
perc_var <- data.frame(PC = 1:20, perc_var = eigenval/sum(eigenval)*100)

# make the plot
plotted = ggplot(perc_var, aes(PC, perc_var)) + geom_bar(stat = "Identity") +
  ylab("Percentage Variance Explained") + theme_light() +
  labs(title = "Percentage Variance Explained by Eigenvalues")

# Save the plot
ggsave(plot = plotted, filename = "percentage_variance_explained.png", width = 8, height = 6, dpi = 300)

# Make the PCA plot using PC1 and PC2
# first, import the merged covariates file since the pca only has individual IDs
covariates <- read.table("mg_covariates.txt", header = TRUE)

# Sort out the PD status and sex
status <- rep(NA, length(covariates$PHENO))
status[grep("0", covariates$PHENO)] <- "Control"
status[grep("1", covariates$PHENO)] <- "Cases"

sex <- rep(NA, length(covariates$SEX))
sex[grep("1", covariates$SEX)] <- "Males"
sex[grep("2", covariates$SEX)] <- "Females"

# combine the variables to plot each in different colours
status_sex <- paste0(status, "_", sex)

# remake data.frame using tibble for easy summary
covs <- as_tibble(data.frame(covariates, status, sex, status_sex))

plotted2 = ggplot(pca_prunned, aes(x = PC1, y = PC2, col = status, shape = sex)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("red", "blue")) +
  coord_equal() + 
  theme_light() +
  xlab(paste0("PC1 (", signif(perc_var$perc_var[1], 3), "%, Variation)")) +
  ylab(paste0("PC2 (", signif(perc_var$perc_var[2], 3), "%, Variation)")) +
  labs(title = "Plot of PC1 vs PC2 on Unlinked Genetic Data Among 1864 Individuals")

ggsave(plot = plotted2, filename = "Plot_of_PC1_vs_PC2.png", width = 8, height = 6, dpi = 300)


