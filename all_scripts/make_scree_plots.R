#!/usr/bin/env Rscript

# Making a Scree plot

# This script only plots variance explained. To see % variance explained, use 'eigenval$VarianceExplained' and 'variance_explained' respectively (variables already defined below) 
# load necessary libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)

# Read in the PCA Eigenvalues and Eigenvectors
print("Read in pca.eigenval files from PLINK")
eigenval <- read.delim("NEW_MPBC_PCA.eigenval", sep ="\t", header = F, stringsAsFactors = F)

# Update column names
colnames(eigenval)[1] <- "Eigenvalues"
eigenval$PC <- as.numeric(rownames(eigenval))

#eigenval$VarianceExplained <- eigenval$Eigenvalues/sum(eigenval$Eigenvalues)*100

# Keeping only the first 10 PCs
eigenval2 <- head(eigenval,10)

# Generating the plot showing only variance explained by the Eigenvalues (not %)
scree <- ggplot(data = eigenval2, aes(x = PC, y = Eigenvalues)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(name = "Principal Components", breaks = seq(0,10,1), limits = c(NA,10)) +
  scale_y_continuous(name = "Variance Explained by Eigenvalues", breaks = seq(0,50,5), limits = c(0,50)) +
  ggtitle("Scree Plot for MPBC Samples (929 Cases; 935 Controls)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Saving out the plot as PDF and JPEG
ggsave("screePlot_MPBC_1-10.jpg", scree, width = 5, height = 3.5, units = "in")


# make another plot to see all 20 principal components ( also not in %)
# Load eigenvalues
eigenvalues <- scan("NEW_MPBC_PCA.eigenval")

# Calculate percentage variance explained
#total_variance <- sum(eigenvalues)
#variance_explained <- (eigenvalues / total_variance) * 100

# Create the scree plot
png("scree_plot_for_20_PCs.png", width=800, height=600)  # Width and height in pixels
plot(1:length(eigenvalues), eigenvalues, type="b", pch=19, col="blue",
     xlab="Principal Component", ylab="Variance Explained by Eigenvalues",
     main="Scree Plot for MPBC Samples (20 PCs)")
grid()
dev.off()
