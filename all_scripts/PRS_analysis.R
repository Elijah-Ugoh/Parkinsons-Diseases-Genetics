rm(list = ls())
setwd("~/Thesis_Project_2024/")
# PRS Analysis using PLINK (90 and 1805 variants) and PRSice outputs (WLD and WO/LD)
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
temp_covs <- read.table("covariates.txt", header = T)
data <- merge(temp_data, temp_covs, by = "FID")

# standardize the phenotype values on the PHENO.X column (from score file)
data$CASE <- data$PHENO.x - 1 #subtract 1 from both the cases and controls to make it easy to handle, and stores the result in a new column CASE.

# Normalization to Z-score
meanControls <- mean(data$SCORE[data$CASE == 0])
sdControls <- sd(data$SCORE[data$CASE == 0])
data$zSCORE <- (data$SCORE - meanControls)/sdControls

# logistic Regression adjusted by covariates
grsTests <- glm(CASE ~ zSCORE + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = data)
summary(grsTests)

# PRS against Age at Diagnosis
cases <- subset(data, CASE == 1)
meanPop <- mean(cases$SCORE)
sdPop <- sd(cases$SCORE)
cases$zSCORE <- (cases$SCORE - meanPop)/sdPop

AADTest <- lm(AAD ~ zSCORE + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = cases)
summary(AADTest)

# Data visualization
data$CASE[data$CASE ==0] <- "Controls"
data$CASE[data$CASE ==1] <- "PD"

# violin plot
p <- ggplot(data, aes(x= reorder(as.factor(CASE), zSCORE), y=zSCORE, fill=as.factor(CASE))) + geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange")) + theme_bw() + ylab("PD PRS (Z-transformed)") + ggtitle("PRS Profile Across Cases and Controls (90 Variants)") + xlab("") + theme(legend.position = "none")
ggsave("90_Variants/PD_PRS_90.jpeg", dpi = 600, units = "in", height = 6, width = 6)

# Quantiles 
# Make Quantile plots using ODDS RATIO
data$CASE <- data$PHENO.x -1
data$quantile1 <- ifelse(data$zSCORE <= quantile(data$zSCORE)[2], 1, 0)
data$quantile2 <- ifelse(data$zSCORE > quantile(data$zSCORE)[2] & data$zSCORE <= quantile(data$zSCORE)[3], 1, 0)
data$quantile3 <- ifelse(data$zSCORE > quantile(data$zSCORE)[3] & data$zSCORE <= quantile(data$zSCORE)[4], 1, 0)
data$quantile4 <- ifelse(data$zSCORE > quantile(data$zSCORE)[4], 1, 0)
data$quantiles <- 1
data$quantiles[data$quantile2 == 1] <- 2
data$quantiles[data$quantile3 == 1] <- 3
data$quantiles[data$quantile4 == 1] <- 4
quintile_PDstatus <- glm(CASE ~ as.factor(data$quantiles) + AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = data)

# Summarize the regression and export a table (to use OR and CI instead of BETA and SE)
summary(quintile_PDstatus)
quintile_PDstatus_summary <- data.frame(exp(cbind(OR = (coef(quintile_PDstatus)), confint(quintile_PDstatus))))
quintile_PDstatus_summary <- quintile_PDstatus_summary[2:4,1:3] #To get only OR and CI for the quantiles
names(quintile_PDstatus_summary) <- c('OR','2.5%','97.5%')
quintile_PDstatus_summary$QUANTILE <- c('2','3','4')
quintile_PDstatus_summary[4,] <- c(0,0,0,'1')
quintile_PDstatus_summary_sorted <- quintile_PDstatus_summary[order(quintile_PDstatus_summary$QUANTILE),]
write.table(quintile_PDstatus_summary_sorted, '90_Variants/quantile_table_90.csv', quote = F, row.names = F, sep = ',')
# Now, make plot
to_plot <- read.table('90_Variants/quantile_table_90.csv', header = T, sep = ',')
to_plot$low <- to_plot$X2.5.
to_plot$high <- to_plot$X97.5.
plotted <- ggplot(to_plot, aes(x=QUANTILE, y=OR)) + geom_pointrange(aes(ymin = low, ymax = high), size=0.7) +
  scale_y_continuous(name = "Odds ratio (95% CI)", breaks = c(0,1,2,3,4,5,6)) + xlab("Quantiles for PRS (Z-standardized)") +
  geom_errorbar(aes(ymin = low, ymax = high, width = 0.05)) +
  ggtitle("OR Against PRS Quantiles, 90 variants") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), panel.background = element_rect(fill = "white"), axis.line = element_line(size = 0.5, colour = "black"), panel.grid.major = element_line(color = "grey90"))
ggsave(plot = plotted, filename = '90_Variants/plotQuantile_90.png', width = 5, height = 6, units = 'in', dpi = 600)


# Create another Quantile plot with just the zScores instead of OR
data <- data %>%
  mutate(quantiles = ntile(zSCORE, 4))  # Divide into 4 groups

# Summarize data for visualization
quantile_summary <- data %>%
  group_by(quantiles, CASE) %>%
  tally(name ="count") %>%
  mutate(
    proportion = count / sum(count),  # Proportion within each quantile
    lower_ci = proportion - qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    upper_ci = proportion + qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    .groups = "drop"  # Ungroup after summarising
  )

# Visualize with geom_pointrange
ggplot(quantile_summary, aes(x = factor(quantiles), y = proportion, color = factor(CASE))) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(0.5)) +
  labs(x = "PD PRS Quantiles", y = "Proportion", color = "PD Status") +
  scale_color_manual(
    values = c("blue", "red"),
    labels = c("0 (Controls)", "1 (Cases)")
  ) +
  theme_minimal() +
  ggtitle("Proportion of Cases and Controls by PD PRS Quantiles")
ggsave(filename = "90_Variants/quantile_plot_by_proportion_90_SNPs.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# ROC calculation and ROC plots
# Plot the main regression model
data$probDisease <- predict(grsTests, type = "response")
ggplot(data, aes(x = zSCORE, y = probDisease, color = as.factor(CASE))) + 
  geom_point(alpha = 0.8) + # Reduce overplotting
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) + # use logistic regression
  labs(x = "Standardized PRS (zSCORE)",
       y = "Predicted Probabilities of PD",
       title = "Predicted Probabilities by zSCORE (90 variants)",
       color = "Disease Status") + 
  scale_color_manual(values = c("blue", "red"), labels = c("Control", "PD")) +
  theme_minimal()
ggsave("90_Variants/Predicted_Probabilities_by_zSCORE_90.jpeg", width = 8, height = 6, dpi = 300) # save plot

# find the best classification cut-off
# First plot histogram if distribution
roc_obj <- roc(data$CASE, data$probDisease)
jpeg("90_Variants/Predicted_Probabilities_for_PD_90.jpeg", width = 1600, height = 1200, res = 300)
par(mar = c(5, 5, 2, 2))  # Adjust margins
hist(data$probDisease, breaks = 30, main = "Predicted Probabilities for PD", col = "skyblue", border = "black")
dev.off()
best_threshold <- coords(roc_obj, "best", ret = "threshold") #use plotROC to find the probability threshold

# Make predictions
data$predicted <- ifelse(data$probDisease > 0.5160492, "DISEASE", "CONTROL")
data$reported <- ifelse(data$CASE == 1, "DISEASE","CONTROL")

# Visualizations
overlayedRocs <- ggplot(data, aes(d = CASE, m = probDisease)) + geom_roc(color ="darkblue", labels = FALSE) + 
  geom_rocci() + style_roc(theme = theme_gray) + ggtitle(paste("ROC Curve based on 90 Variants PRS (AUC =", round(auc(roc_obj), 2), ")")) + 
  geom_abline(linetype = "dashed") + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "90_Variants/plotRoc_90_SNPs.png", width = 8, height = 5, units = "in", dpi = 300)


# Show the confusion matrix (specificity and sensitivity)
confMat <- confusionMatrix(data = as.factor(data$predicted), reference = as.factor(data$reported), positive = "DISEASE")
confMat


# create Density plots
densPlot <- ggplot(data, aes(probDisease, fill = reported, color = reported)) + 
  ggtitle("Overall Distribution of PD Probability Based on PRS from 90 Variants") +
  geom_density(alpha = 0.5) + theme_bw()
ggsave(plot = densPlot, filename = "90_Variants/plotDensity_90_SNPs.png", width = 8, height = 5, units = "in", dpi = 300)
