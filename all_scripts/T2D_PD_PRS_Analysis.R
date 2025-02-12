# PRS Analysis for T2D + PD 
# Download the necessary packages 
if (!require(tidyverse)) install.packages('tidyr')
if (!require(data.table)) install.packages('data.table')
if (!require(dplyr)) install.packages('dplyr')
if (!require(plyr)) install.packages('plyr')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(caret)) install.packages('caret')
# install plotROC from GitHub
remotes::install_github("sachsmc/plotROC")


# Load the necessary packages 
library(tidyr)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(plotROC)
library(caret)

# set working directory
setwd("/Users/Elijah/Documents/Thesis_Project_2024")

# read the calculated scores
Ge_Score <- read.table("PRS_T2D/Ge_PRS.profile", header = T)
Khera_Score <- read.table("PRS_T2D/Khera_PRS.profile", header = T)
Mars_Score <- read.table("PRS_T2D/Mars_PRS.profile", header = T)
Lin_Score <- read.table("PRS_T2D/Lin_PRS.profile", header = T)

# Read the covariates
covariates <- read.table("Scores_PD-DM/new_covariates.txt", header = T)

# Merge scores with covariates
merged_Ge <- merge(Ge_Score, covariates, by = "FID")
merged_Khera <- merge(Khera_Score, covariates, by = "FID")
merged_Mars <- merge(Mars_Score, covariates, by = "FID")
merged_Lin <- merge(Lin_Score, covariates, by = "FID")

# Combine all merged data
combined_scores <- rbind(merged_Ge, merged_Khera, merged_Mars, merged_Lin)

# Function to plot ROC and calculate AUC
plot_roc_auc <- function(data, title) {
  roc_data <- roc(data$T2D, data$PRS)
  auc_value <- auc(roc_data)
  print(paste(title, "AUC:", auc_value))
  plot.roc(roc_data, col = "#1c61b6", lwd = 2, main = paste("ROC Curve for", title))
  return(auc_value)
}

# Plot ROC curves and calculate AUC for each PRS
Ge_AUC <- plot_roc_auc(merged_Ge, "Ge")
Khera_AUC <- plot_roc_auc(merged_Khera, "Khera")
Mars_AUC <- plot_roc_auc(merged_Mars, "Mars")
Lin_AUC <- plot_roc_auc(merged_Lin, "Lin")

# Summary of AUC values
auc_summary <- data.frame(
  PRS = c("Ge", "Khera", "Mars", "Lin"),
  AUC = c(Ge_AUC, Khera_AUC, Mars_AUC, Lin_AUC)
)

print(auc_summary)

# Create Density Plots for each PRS
ggplot(combined_scores, aes(x = PRS, fill = PRS_Name)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of PRS", x = "Polygenic Risk Score", y = "Density") +
  theme_minimal() +
  theme(legend.position = "top")

# Save the density plot
ggsave("PRS_Density_Plot.png", dpi = 300, width = 8, height = 5)























#################################################################################
data <- merge(temp_data, temp_covs, by = "FID") #FID is a common column for both the covariate and plink-generated profile files
data$CASE <- data$PHENO.x - 1 # A new column is added: 1 is subtracted from the PHENO column making cases "1" and controls "0".

# Normalize Score to Z-Score
# Note the mean = 0 and SD = 1
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls

# Perform logistic regression adjusted by covariates
grsTests <- glm(CASE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AGE, family="binomial", data = data)
summary(grsTests)


# GRS versus age at onset
# Subset ONLY cases perform linear regression adjusted by covariates
cases <- subset(data, CASE == 1)
meanPop <- mean(cases$SCORE)
sdPop <- sd(cases$SCORE)
cases$zSCORE <- (cases$SCORE - meanPop)/sdPop
grsTests <- lm(AGE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = cases)
summary(grsTests)


# Data visualization - Violin plots
data$CASE[data$CASE ==0] <- "Controls"
data$CASE[data$CASE ==1] <- "PD"

p <- ggplot(data, aes(x= reorder(as.factor(CASE), zSCORE), y=zSCORE, fill=as.factor(CASE))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange")) + theme_bw() + ylab("PD GRS (Z-transformed)") +xlab("") + theme(legend.position = "none")
ggsave("PD_GRS.jpeg", dpi = 600, units = "in", height = 6, width = 6)

# Quantile plots
# Make quantiles
data$CASE <- data$PHENO.x -1
data$quantile1 <- ifelse(data$zSCORE <= quantile(data$zSCORE)[2], 1, 0)
data$quantile2 <- ifelse(data$zSCORE > quantile(data$zSCORE)[2] & data$zSCORE <= quantile(data$zSCORE)[3], 1, 0)
data$quantile3 <- ifelse(data$zSCORE > quantile(data$zSCORE)[3] & data$zSCORE <= quantile(data$zSCORE)[4], 1, 0)
data$quantile4 <- ifelse(data$zSCORE > quantile(data$zSCORE)[4], 1, 0)
data$quantiles <- 1
data$quantiles[data$quantile2 == 1] <- 2
data$quantiles[data$quantile3 == 1] <- 3
data$quantiles[data$quantile4 == 1] <- 4
quintileTests <- glm(CASE ~ as.factor(data$quantiles) + AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = data)

# Summarize the regression and export a table
summary(quintileTests)
summary_stats <- data.frame(summary(quintileTests)$coef[2:4,1:2])
names(summary_stats) <- c("BETA","SE")
summary_stats$QUANTILE <- c("2nd","3rd","4th")
summary_stats[4,] <- c(0,0,"1st")
summary_stats_sorted <- summary_stats[order(summary_stats$QUANTILE),]
write.table(summary_stats_sorted, "quantile_table.csv", quote = F, row.names = F, sep = ",")


# Make quantile plot
data$CASE <- data$PHENO.x - 1
to_plot <- read.table("quantile_table.csv", header = T, sep = ",")
to_plot$low <- to_plot$BETA - (1.96*to_plot$SE)
to_plot$high <- to_plot$BETA + (1.96*to_plot$SE)
plotted <- ggplot(to_plot, aes(QUANTILE, BETA)) + geom_pointrange(aes(ymin = low, ymax = high))
ggsave(plot = plotted, filename = "plotQuantile.png", width = 4, height = 4, units = "in", dpi = 300)

# Run regression model
Model <- glm(CASE ~ SCORE, data = data, family = 'binomial')

# Make predictions
data$probDisease <- predict(Model, data, type = c("response"))
data$predicted <- ifelse(data$probDisease > 0.5, "DISEASE", "CONTROL")
data$reported <- ifelse(data$CASE == 1, "DISEASE","CONTROL")

# Data visualization - ROC plots
overlayedRocs <- ggplot(data, aes(d = CASE, m = probDisease)) + geom_roc(labels = FALSE) + geom_rocci() + style_roc(theme = theme_gray) + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "plotRoc.png", width = 8, height = 5, units = "in", dpi = 300)


# Show the confusion matrix (specificity and sensitivity)
confMat <- confusionMatrix(data = as.factor(data$predicted), reference = as.factor(data$reported), positive = "DISEASE")
confMat

# Data visualization - Density plots
densPlot <- ggplot(data, aes(probDisease, fill = reported, color = reported)) + geom_density(alpha = 0.5) + theme_bw()
ggsave(plot = densPlot, filename = "plotDensity.png", width = 8, height = 5, units = "in", dpi = 300)




all_factors <- read.table("07_stratification/all_merged_data.txt", header = T)
#rename the FID column
colnames(all_factors)[1] <- "FID"
colnames(all_factors)[2] <- "PD"

covariates <- read.table("covariates.txt", header = T)

# Add only PCs from covariate file to the data
# Select first 5 PCs (columns 10 to 14) and merge with the data
data <- merge(all_factors, covariates[, c(1, 10:14)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(data$SCORE[data$PD == 0])
sdControls <- sd(data$SCORE[data$PD == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls

# Now, logistic Regression model for all three factors covariates
model_covs <- glm(PD ~ zSCORE + Smoking_Ever + Snus_Ever + Pesticides_Ever + Caffeine_Ever + 
                    Caffeine_Consumption_Level + sex + education + Age_Inclusion + 
                    PC1 + PC2 + PC3 + PC4 + PC5, family = binomial, data = data)
summary(model_covs)
