setwd("~/Documents/cosmos-home/04_PRS/")
list.files()


# Download the necessary packages 
if (!require(tidyverse)) install.packages('tidyr')
if (!require(data.table)) install.packages('data.table')
if (!require(dplyr)) install.packages('dplyr')
if (!require(plyr)) install.packages('plyr')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(caret)) install.packages('caret')
if (!require(plotROC)) install.packages("plotROC")


# Load the necessary packages 
library(tidyr)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(plotROC)
library(caret)

# Read in the plink output and merge the data with the covariate file 
temp_data <- read.table("PRS_PD_MPBC.profile", header = T) 
temp_covs <- read.table("01_covariate_file/covariates.txt", header = T)
data <- merge(temp_data, temp_covs, by = "FID")

# stabdardize the phenotype values on the PHENO.X column (from score file)
data$CASE <- data$PHENO.x - 1 #subtract 1 from both the cases and controls to make it easy to handle, and stores the result in a new column CASE.

# Normalization to Z-score
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls

# logistic Regression adjusted by covariates
grsTests <- glm(CASE ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + AAD, family="binomial", data = data)
summary(grsTests)

# GRS against Age at Diagnosis
cases <- subset(data, CASE == 1)
meanPop <- mean(cases$SCORE)
sdPop <- sd(cases$SCORE)
cases$zSCORE <- (cases$SCORE - meanPop)/sdPop

grsTests <- lm(AAD ~ zSCORE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = cases)
summary(grsTests)


# Data visualization
data$CASE[data$CASE ==0] <- "Controls"
data$CASE[data$CASE ==1] <- "PD"

# violin plot
p <- ggplot(data, aes(x= reorder(as.factor(CASE), zSCORE), y=zSCORE, fill=as.factor(CASE))) + geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange")) + theme_bw() + ylab("PD PRS (Z-transformed)") +xlab("") + theme(legend.position = "none")
ggsave("PD_PRS.jpeg", dpi = 600, units = "in", height = 6, width = 6)

# Quantiles 
data$CASE <- data$PHENO.x -1
data$quantile1 <- ifelse(data$zSCORE <= quantile(data$zSCORE)[2], 1, 0)
data$quantile2 <- ifelse(data$zSCORE > quantile(data$zSCORE)[2] & data$zSCORE <= quantile(data$zSCORE)[3], 1, 0)
data$quantile3 <- ifelse(data$zSCORE > quantile(data$zSCORE)[3] & data$zSCORE <= quantile(data$zSCORE)[4], 1, 0)
data$quantile4 <- ifelse(data$zSCORE > quantile(data$zSCORE)[4], 1, 0)
data$quantiles <- 1
data$quantiles[data$quantile2 == 1] <- 2
data$quantiles[data$quantile3 == 1] <- 3
data$quantiles[data$quantile4 == 1] <- 4
quintileTests <- glm(CASE ~ as.factor(data$quantiles) + AAD + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = data)

# Summarize the regression and make a quantile table
summary(quintileTests)
summary_stats <- data.frame(summary(quintileTests)$coef[2:4,1:2])
names(summary_stats) <- c("BETA","SE")
summary_stats$QUANTILE <- c("2nd","3rd","4th")
summary_stats[4,] <- c(0,0,"1st")
summary_stats_sorted <- summary_stats[order(summary_stats$QUANTILE),]
write.table(summary_stats_sorted, "quantile_table.csv", quote = F, row.names = F, sep = ",")

# Make the quantile plot
data$CASE <- data$PHENO.x - 1
to_plot <- read.table("quantile_table.csv", header = T, sep = ",")
to_plot$low <- to_plot$BETA - (1.96*to_plot$SE)
to_plot$high <- to_plot$BETA + (1.96*to_plot$SE)
plotted <- ggplot(to_plot, aes(QUANTILE, BETA)) + geom_pointrange(aes(ymin = low, ymax = high))
ggsave(plot = plotted, filename = "plotQuantile.png", width = 4, height = 4, units = "in", dpi = 300)

# ROC calculation and ROC plots
Model <- glm(CASE ~ SCORE, data = data, family = 'binomial')

# Make predictions
data$probDisease <- predict(Model, data, type = c("response"))
data$predicted <- ifelse(data$probDisease > 0.5, "DISEASE", "CONTROL")
data$reported <- ifelse(data$CASE == 1, "DISEASE","CONTROL")

# Visualizations
overlayedRocs <- ggplot(data, aes(d = CASE, m = probDisease)) + geom_roc(labels = FALSE) + geom_rocci() + style_roc(theme = theme_gray) + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "plotRoc.png", width = 8, height = 5, units = "in", dpi = 300)

# Show the confusion matrix (specificity and sensitivity)
confMat <- confusionMatrix(data = as.factor(data$predicted), reference = as.factor(data$reported), positive = "DISEASE")
confMat

# create Density plots
densPlot <- ggplot(data, aes(probDisease, fill = reported, color = reported)) + geom_density(alpha = 0.5) + theme_bw()
ggsave(plot = densPlot, filename = "plotDensity.png", width = 8, height = 5, units = "in", dpi = 300)
