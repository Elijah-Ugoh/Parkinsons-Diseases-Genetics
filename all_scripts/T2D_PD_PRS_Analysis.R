# PRS Analysis for T2D + PD 
# Download the necessary packages 
if (!require(tidyverse)) install.packages('tidyr')
if (!require(all_data.table)) install.packages('all_data.table')
if (!require(dplyr)) install.packages('dplyr')
if (!require(plyr)) install.packages('plyr')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(caret)) install.packages('caret')
# install plotROC from GitHub
remotes::install_github("sachsmc/plotROC")


# Load the necessary packages 
library(tidyr)
library(all_data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(pROC)
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
colnames(covariates)[6] <- "PD" # Rename the column to PD

# Merge scores with covariates, but only keep FID and SCORE columns. 
# The PHENO column is not needed because it is already in the covariates all_dataframe (T2D)
Ge_Score <- Ge_Score[, c("FID", "SCORE")] # retain only the FID and SCORE columns
colnames(Ge_Score)[colnames(Ge_Score) == "SCORE"] <- "Ge_Score" # Rename the column for easy identification when all all_data is merged
Khera_Score <- Khera_Score[, c("FID", "SCORE")]
colnames(Khera_Score)[colnames(Khera_Score) == "SCORE"] <- "Khera_Score"
Mars_Score <- Mars_Score[, c("FID", "SCORE")]
colnames(Mars_Score)[colnames(Mars_Score) == "SCORE"] <- "Mars_Score"
Lin_Score <- Lin_Score[, c("FID", "SCORE")]
colnames(Lin_Score)[colnames(Lin_Score) == "SCORE"] <- "Lin_Score"

# merge all all_data into one all_dataframe
all_data <- Reduce(function(x, y) merge(x, y, by = "FID", all = TRUE), 
                   list(covariates, Ge_Score, Khera_Score, Mars_Score, Lin_Score))

# Normalize all Scores to Z-Score (using controls, mean = 0 and SD = 1)
# Identify PRS score columns dynamically (all columns ending in "_Score")
prs_columns <- grep("_Score$", colnames(all_data), value = TRUE)

# Loop over each PRS column and calculate Z-score
for (score in prs_columns) {
  meanControls <- mean(all_data[[score]][all_data$T2D == 0], na.rm = TRUE)
  sdControls <- sd(all_data[[score]][all_data$T2D == 0], na.rm = TRUE)
  
  # Create new column with Z-score
  z_score <- paste0("z_", score)
  all_data[[z_score]] <- (all_data[[score]] - meanControls) / sdControls
}


# Further classification: Modify the dataset to include a new column that categorizes individuals into one of the four groups

# Create a new classification variable
all_data <- all_data %>%
  mutate(
    Disease_Group = case_when(
      PD == 1 & T2D == 0 ~ "PD Only",
      PD == 0 & T2D == 1 ~ "T2D Only",
      PD == 1 & T2D == 1 ~ "PD + T2D",
      PD == 0 & T2D == 0 ~ "Healthy Controls"
    )
  )

# Convert to a factor
all_data$Disease_Group <- factor(all_data$Disease_Group, 
                                 levels = c("Healthy Controls", "PD Only", "T2D Only", "PD + T2D"))
# check the new data column
table(all_data$Disease_Group)
# Save to file
write.table(all_all_data, "Scores_PD-DM/merged_all_scores.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Make plots
# Perform logistic regression adjusted by covariates
Ge_model <- glm(T2D ~ z_Ge_Score + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5, 
                family="binomial", data = all_data)
summary(Ge_model)
Khera_model <- glm(T2D ~ z_Khera_Score + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5, 
                family="binomial", data = all_data)
summary(Khera_model)
Mars_model <- glm(T2D ~ z_Mars_Score + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5, 
                   family="binomial", data = all_data)
summary(Mars_model)
Lin_model <- glm(T2D ~ z_Lin_Score + SEX + AGE + PC1 + PC2 + PC3 + PC4 + PC5, 
                  family="binomial", data = all_data)
summary(Lin_model)


# all_data visualization
####################################################################################################
# Ge, 2022
# Violin plots
# First, create a new column
all_data$CASE_Label <- ifelse(all_data$T2D == 0, "CONTROL", "T2D")

p <- ggplot(all_data, aes(x= reorder(as.factor(Disease_Group), z_Ge_Score), y=z_Ge_Score, fill=as.factor(Disease_Group))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange", "blue", "red")) + 
  theme_bw() + ylab("T2D PRS (Z-transformed)") +
  xlab("Disease Status") + 
  theme(legend.position = "none") +
  ggtitle("Distribution of T2D PRS (Ge Score) by Disease Status")
ggsave("Scores_PD-DM/T2D_Ge_PRS.jpeg", dpi = 600, units = "in", height = 6, width = 6)

# Make Quantile plots using ODDS RATIO
all_data$quantile1 <- ifelse(all_data$z_Ge_Score <= quantile(all_data$z_Ge_Score)[2], 1, 0)
all_data$quantile2 <- ifelse(all_data$z_Ge_Score > quantile(all_data$z_Ge_Score)[2] & all_data$z_Ge_Score <= quantile(all_data$z_Ge_Score)[3], 1, 0)
all_data$quantile3 <- ifelse(all_data$z_Ge_Score > quantile(all_data$z_Ge_Score)[3] & all_data$z_Ge_Score <= quantile(all_data$z_Ge_Score)[4], 1, 0)
all_data$quantile4 <- ifelse(all_data$z_Ge_Score > quantile(all_data$z_Ge_Score)[4], 1, 0)
all_data$quantiles <- 1
all_data$quantiles[all_data$quantile2 == 1] <- 2
all_data$quantiles[all_data$quantile3 == 1] <- 3
all_data$quantiles[all_data$quantile4 == 1] <- 4
quintile_T2Dstatus <- glm(as.factor(CASE_Label) ~ as.factor(all_data$quantiles) + AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5, family='binomial', data = all_data)

# Summarize the regression and export a table (to use OR and CI instead of BETA and SE)
summary(quintile_T2Dstatus)

quintile_T2Dstatus_summary <- data.frame(exp(cbind(OR = (coef(quintile_T2Dstatus)), confint(quintile_T2Dstatus))))
quintile_T2Dstatus_summary <- quintile_T2Dstatus_summary[2:4,1:3] #To get only OR and CI for the quantiles
names(quintile_T2Dstatus_summary) <- c('OR','2.5%','97.5%')
quintile_T2Dstatus_summary$QUANTILE <- c('2','3','4')
quintile_T2Dstatus_summary[4,] <- c(0,0,0,'1')
quintile_T2Dstatus_summary_sorted <- quintile_T2Dstatus_summary[order(quintile_T2Dstatus_summary$QUANTILE),]
write.table(quintile_T2Dstatus_summary_sorted, 'Scores_PD-DM/quantile_table_MPBC_T2D_PRS_OR_20250214_ge.csv', quote = F, row.names = F, sep = ',')
# Now, make plot
to_plot <- read.table('Scores_PD-DM/quantile_table_MPBC_T2D_PRS_OR_20250214_ge.csv', header = T, sep = ',')
to_plot$low <- to_plot$X2.5.
to_plot$high <- to_plot$X97.5.
plotted <- ggplot(to_plot, aes(x=QUANTILE, y=OR)) + geom_pointrange(aes(ymin = low, ymax = high), size=0.7) +
  scale_y_continuous(name = "Odds ratio (95% CI)", breaks = c(0,1,2,3,4,5,6)) + xlab("Quantiles for PRS (Z-standardized)") +
  geom_errorbar(aes(ymin = low, ymax = high, width = 0.05)) +
  ggtitle("Ge, 2022") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), panel.background = element_rect(fill = "white"), axis.line = element_line(size = 0.5, colour = "black"), panel.grid.major = element_line(color = "grey90"))
ggsave(plot = plotted, filename = 'Scores_PD-DM/plotQuantile_PRS_MPBC_OR_T2Dstatus_ge.png', width = 5, height = 6, units = 'in', dpi = 600)

# Density plot T2D Cases and controls only
ggplot(all_data, aes(x = z_Ge_Score, fill = factor(T2D))) +
  geom_density(alpha = 0.5, size = 0.1) +
  labs(x = "Polygenic Risk Score (PRS)", fill = "T2D") +
  scale_fill_manual(values = c("cyan", "purple"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  ggtitle("PRS Distribution for T2D in MPBC Cohort (GWAS - Ge, 2022)")
ggsave(filename = "Scores_PD-DM/plotDensity_ge.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Use the four-group classification created above
# Density plot for PRS across disease groups
ggplot(all_data, aes(x = z_Ge_Score, fill = Disease_Group)) +
  geom_density(alpha = 0.5, size = 0.1) +
  labs(x = "Polygenic Risk Score (PRS)", fill = "Disease Group") +
  scale_fill_manual(values = c("skyblue", "salmon", "purple", "green")) + # Customize colors
  theme_minimal() +
  ggtitle("PRS Distribution Across Disease Groups (Ge, 2022)")

# Save the plot
ggsave(filename = "Scores_PD-DM/plotDensity_groups_ge.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Make predictions to access model performance
# First, run a simple regression model
Model_ge <- glm(T2D ~ Ge_Score + AGE + SEX, data = all_data, family = 'binomial')
all_data$probDisease <- predict(Model_ge, all_data, type = c("response"))
# find the best classification cut-off
roc_obj <- roc(all_data$T2D, all_data$probDisease)
hist(all_data$probDisease, breaks = 30, main = "Predicted Probabilities (GWAS - Ge, 2022)") # use histogram of predicted probabilities
best_threshold <- coords(roc_obj, "best", ret = "threshold") #use pROC to find the probability threshold
all_data$predicted <- ifelse(all_data$probDisease > 0.0740, "DISEASE", "CONTROL")
all_data$reported <- ifelse(all_data$T2D == 1, "DISEASE","CONTROL")

# ROC plots
overlayedRocs <- ggplot(all_data, aes(d = T2D, m = probDisease)) + geom_roc(color ="darkblue", labels = FALSE) +
  geom_rocci() + style_roc(theme = theme_gray) + ggtitle(paste("ROC Curve based on Ge, 2022 GWAS (AUC =", round(auc(roc_obj), 2), ")")) + 
  geom_abline(linetype = "dashed") + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "Scores_PD-DM/plotRoc_ge.png", width = 8, height = 5, units = "in", dpi = 300)

# Show the confusion matrix (specificity and sensitivity)
confMat <- confusionMatrix(data = as.factor(all_data$predicted), reference = as.factor(all_data$reported), positive = "DISEASE")
confMat
# save it to an image
library(gridExtra)
# Convert confusion matrix to a table format
confMat_text <- capture.output(print(confMat))
# Save as an image
png("Scores_PD-DM/confusion_matrix_ge.png", width = 800, height = 600)
grid.table(confMat_text)  # Display text as an image
dev.off()

# Using more visualizations
# Boxplot of PRS scores by T2D status
ggplot(all_data, aes(x = reorder(as.factor(T2D), z_Ge_Score), y = z_Ge_Score, fill = as.factor(T2D))) +
  geom_boxplot() +
  labs(x = "T2D Status", y = "Polygenic Risk Score (PRS)", fill = "T2D Status") +
  scale_fill_manual(values = c("skyblue", "salmon"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  ggtitle("PRS Distribution by T2D Status (GWAS - Ge, 2022)")
ggsave(filename = "Scores_PD-DM/boxplot_ge.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Create another Quantile plot with just the zScores instead of OR
all_data <- all_data %>%
  mutate(quantiles = ntile(z_Ge_Score, 4))  # Divide into 4 groups

# Summarize data for visualization
quantile_summary <- all_data %>%
  group_by(quantiles, T2D) %>%
  tally(name ="count") %>%
  mutate(
    proportion = count / sum(count),  # Proportion within each quantile
    lower_ci = proportion - qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    upper_ci = proportion + qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    .groups = "drop"  # Ungroup after summarising
  )

# Visualize with geom_pointrange
ggplot(quantile_summary, aes(x = factor(quantiles), y = proportion, color = factor(T2D))) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(0.5)) +
  labs(x = "T2D PRS Quantiles", y = "Proportion", color = "T2D Status") +
  scale_color_manual(
    values = c("blue", "red"),
    labels = c("0 (Controls)", "1 (Cases)")
  ) +
  theme_minimal() +
  ggtitle("Proportion of Cases and Controls by T2D PRS Quantiles (GWAS - Ge, 2022)")
ggsave(filename = "Scores_PD-DM/quantile_plot_by_proportion_ge.jpeg", width = 8, height = 5, units = "in", dpi = 300)

####################################################################################################
# Khera, 2018
# Violin plots
# Use the CASE_Label column created above
p <- ggplot(all_data, aes(x= reorder(as.factor(Disease_Group), z_Khera_Score), y=z_Khera_Score, fill=as.factor(Disease_Group))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange", "blue", "red")) + theme_bw() + 
  ylab("T2D PRS (Z-transformed)") +
  xlab("Disease Status") + 
  theme(legend.position = "none") +
  ggtitle("Distribution of T2D PRS (Khera Score) by Disease Status")
ggsave("Scores_PD-DM/T2D_Khera_PRS.jpeg", dpi = 600, units = "in", height = 6, width = 6)

# Make Quantile plots using ODDS RATIO
all_data$quantile1 <- ifelse(all_data$z_Khera_Score <= quantile(all_data$z_Khera_Score)[2], 1, 0)
all_data$quantile2 <- ifelse(all_data$z_Khera_Score > quantile(all_data$z_Khera_Score)[2] & all_data$z_Khera_Score <= quantile(all_data$z_Khera_Score)[3], 1, 0)
all_data$quantile3 <- ifelse(all_data$z_Khera_Score > quantile(all_data$z_Khera_Score)[3] & all_data$z_Khera_Score <= quantile(all_data$z_Khera_Score)[4], 1, 0)
all_data$quantile4 <- ifelse(all_data$z_Khera_Score > quantile(all_data$z_Khera_Score)[4], 1, 0)
all_data$quantiles <- 1
all_data$quantiles[all_data$quantile2 == 1] <- 2
all_data$quantiles[all_data$quantile3 == 1] <- 3
all_data$quantiles[all_data$quantile4 == 1] <- 4
quintile_T2Dstatus <- glm(as.factor(CASE_Label) ~ as.factor(all_data$quantiles) + AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5, family='binomial', data = all_data)

# Summarize the regression and export a table (to use OR and CI instead of BETA and SE)
summary(quintile_T2Dstatus)

quintile_T2Dstatus_summary <- data.frame(exp(cbind(OR = (coef(quintile_T2Dstatus)), confint(quintile_T2Dstatus))))
quintile_T2Dstatus_summary <- quintile_T2Dstatus_summary[2:4,1:3] #To get only OR and CI for the quantiles
names(quintile_T2Dstatus_summary) <- c('OR','2.5%','97.5%')
quintile_T2Dstatus_summary$QUANTILE <- c('2','3','4')
quintile_T2Dstatus_summary[4,] <- c(0,0,0,'1')
quintile_T2Dstatus_summary_sorted <- quintile_T2Dstatus_summary[order(quintile_T2Dstatus_summary$QUANTILE),]
write.table(quintile_T2Dstatus_summary_sorted, 'Scores_PD-DM/quantile_table_MPBC_T2D_PRS_OR_20250214_khera.csv', quote = F, row.names = F, sep = ',')
# Now, make plot
to_plot <- read.table('Scores_PD-DM/quantile_table_MPBC_T2D_PRS_OR_20250214_khera.csv', header = T, sep = ',')
to_plot$low <- to_plot$X2.5.
to_plot$high <- to_plot$X97.5.
plotted <- ggplot(to_plot, aes(x=QUANTILE, y=OR)) + geom_pointrange(aes(ymin = low, ymax = high), size=0.7) +
  scale_y_continuous(name = "Odds ratio (95% CI)", breaks = c(0,1,2,3,4,5,6)) + xlab("Quantiles for PRS (Z-standardized)") +
  geom_errorbar(aes(ymin = low, ymax = high, width = 0.05)) +
  ggtitle("Khera, 2018") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), panel.background = element_rect(fill = "white"), axis.line = element_line(size = 0.5, colour = "black"), panel.grid.major = element_line(color = "grey90"))
ggsave(plot = plotted, filename = 'Scores_PD-DM/plotQuantile_PRS_MPBC_OR_T2Dstatus_khera.png', width = 5, height = 6, units = 'in', dpi = 600)

# Density plot T2D cases and controls only
ggplot(all_data, aes(x = z_Khera_Score, fill = factor(T2D))) +
  geom_density(alpha = 0.5, size = 0.1) +
  labs(x = "Polygenic Risk Score (PRS)", fill = "T2D") +
  scale_fill_manual(values = c("cyan", "purple"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  ggtitle("PRS Distribution for T2D in MPBC Cohort (GWAS - Khera, 2018)")
ggsave(filename = "Scores_PD-DM/plotDensity_khera.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Use the four-group classification created above
# Density plot for PRS across disease groups
ggplot(all_data, aes(x = z_Khera_Score, fill = Disease_Group)) +
  geom_density(alpha = 0.5, size = 0.1) +
  labs(x = "Polygenic Risk Score (PRS)", fill = "Disease Group") +
  scale_fill_manual(values = c("skyblue", "salmon", "purple", "green")) + # Customize colors
  theme_minimal() +
  ggtitle("PRS Distribution Across Disease Groups (Khera, 2018)")
# Save the plot
ggsave(filename = "Scores_PD-DM/plotDensity_groups_khera.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Make predictions to access model performance
# First, run a simple regression model
Model_khera <- glm(T2D ~ Khera_Score + AGE + SEX, data = all_data, family = 'binomial')
all_data$probDisease <- predict(Model_khera, all_data, type = c("response"))
# find the best classification cut-off
roc_obj <- roc(all_data$T2D, all_data$probDisease)
hist(all_data$probDisease, breaks = 30, main = "Predicted Probabilities (GWAS - Khera, 2018)") # use histogram of predicted probabilities
best_threshold <- coords(roc_obj, "best", ret = "threshold") #use pROC to find the probability threshold
all_data$predicted <- ifelse(all_data$probDisease > 0.061, "DISEASE", "CONTROL")
all_data$reported <- ifelse(all_data$T2D == 1, "DISEASE","CONTROL")

# ROC plots
overlayedRocs <- ggplot(all_data, aes(d = T2D, m = probDisease)) + geom_roc(color ="darkblue", labels = FALSE) +
  geom_rocci() + style_roc(theme = theme_gray) + ggtitle(paste("ROC Curve based on Khera, 2018 GWAS (AUC =", round(auc(roc_obj), 2), ")")) + 
  geom_abline(linetype = "dashed") + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "Scores_PD-DM/plotRoc_khera.png", width = 8, height = 5, units = "in", dpi = 300)

# Show the confusion matrix (specificity and sensitivity)
confMat <- confusionMatrix(data = as.factor(all_data$predicted), reference = as.factor(all_data$reported), positive = "DISEASE")
confMat
# save to an image
library(gridExtra)
# Convert confusion matrix to a table format
confMat_text <- capture.output(print(confMat))
# Save as an image
png("Scores_PD-DM/confusion_matrix_khera.png", width = 800, height = 600)
grid.table(confMat_text)  # Display text as an image
dev.off()

# Using more visualizations
# Boxplot of PRS scores by T2D status
ggplot(all_data, aes(x = reorder(as.factor(T2D), z_Khera_Score), y = z_Khera_Score, fill = as.factor(T2D))) +
  geom_boxplot() +
  labs(x = "T2D Status", y = "Polygenic Risk Score (PRS)", fill = "T2D Status") +
  scale_fill_manual(values = c("skyblue", "salmon"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  ggtitle("PRS Distribution by T2D Status (Khera, 2018)")
ggsave(filename = "Scores_PD-DM/boxplot_khera.jpeg", width = 8, height = 5, units = "in", dpi = 300)


# Create another Quantile plot with just the zScores instead of OR
all_data <- all_data %>%
  mutate(quantiles = ntile(z_Khera_Score, 4))  # Divide into 4 groups

# Summarize data for visualization
quantile_summary <- all_data %>%
  group_by(quantiles, T2D) %>%
  tally(name ="count") %>%
  mutate(
    proportion = count / sum(count),  # Proportion within each quantile
    lower_ci = proportion - qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    upper_ci = proportion + qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    .groups = "drop"  # Ungroup after summarising
  )

# Visualize with geom_pointrange
ggplot(quantile_summary, aes(x = factor(quantiles), y = proportion, color = factor(T2D))) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(0.5)) +
  labs(x = "T2D PRS Quantiles", y = "Proportion", color = "T2D Status") +
  scale_color_manual(
    values = c("blue", "red"),
    labels = c("0 (Controls)", "1 (Cases)")
  ) +
  theme_minimal() +
  ggtitle("Proportion of Cases and Controls by T2D PRS Quantiles (GWAS - Khera, 2018)")
ggsave(filename = "Scores_PD-DM/quantile_plot_by_proportion_khera.jpeg", width = 8, height = 5, units = "in", dpi = 300)

#######################################################################################################
# Mars, 2020
# Violin plots
# Use the CASE_Label column created above
p <- ggplot(all_data, aes(x= reorder(as.factor(Disease_Group), z_Mars_Score), y=z_Mars_Score, fill=as.factor(Disease_Group))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange", "blue", "red")) + theme_bw() + 
  ylab("T2D PRS (Z-transformed)") +
  xlab("Disease Status") + 
  theme(legend.position = "none") +
  ggtitle("Distribution of T2D PRS (Mars Score) by Disease Status")
ggsave("Scores_PD-DM/T2D_Mars_PRS.jpeg", dpi = 600, units = "in", height = 6, width = 6)

# Make Quantile plots using ODDS RATIO
all_data$quantile1 <- ifelse(all_data$z_Mars_Score <= quantile(all_data$z_Mars_Score)[2], 1, 0)
all_data$quantile2 <- ifelse(all_data$z_Mars_Score > quantile(all_data$z_Mars_Score)[2] & all_data$z_Mars_Score <= quantile(all_data$z_Mars_Score)[3], 1, 0)
all_data$quantile3 <- ifelse(all_data$z_Mars_Score > quantile(all_data$z_Mars_Score)[3] & all_data$z_Mars_Score <= quantile(all_data$z_Mars_Score)[4], 1, 0)
all_data$quantile4 <- ifelse(all_data$z_Mars_Score > quantile(all_data$z_Mars_Score)[4], 1, 0)
all_data$quantiles <- 1
all_data$quantiles[all_data$quantile2 == 1] <- 2
all_data$quantiles[all_data$quantile3 == 1] <- 3
all_data$quantiles[all_data$quantile4 == 1] <- 4
quintile_T2Dstatus <- glm(as.factor(CASE_Label) ~ as.factor(all_data$quantiles) + AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5, family='binomial', data = all_data)

# Summarize the regression and export a table (to use OR and CI instead of BETA and SE)
summary(quintile_T2Dstatus)

quintile_T2Dstatus_summary <- data.frame(exp(cbind(OR = (coef(quintile_T2Dstatus)), confint(quintile_T2Dstatus))))
quintile_T2Dstatus_summary <- quintile_T2Dstatus_summary[2:4,1:3] #To get only OR and CI for the quantiles
names(quintile_T2Dstatus_summary) <- c('OR','2.5%','97.5%')
quintile_T2Dstatus_summary$QUANTILE <- c('2','3','4')
quintile_T2Dstatus_summary[4,] <- c(0,0,0,'1')
quintile_T2Dstatus_summary_sorted <- quintile_T2Dstatus_summary[order(quintile_T2Dstatus_summary$QUANTILE),]
write.table(quintile_T2Dstatus_summary_sorted, 'Scores_PD-DM/quantile_table_MPBC_T2D_PRS_OR_20250214_Mars.csv', quote = F, row.names = F, sep = ',')
# Now, make plot
to_plot <- read.table('Scores_PD-DM/quantile_table_MPBC_T2D_PRS_OR_20250214_Mars.csv', header = T, sep = ',')
to_plot$low <- to_plot$X2.5.
to_plot$high <- to_plot$X97.5.
plotted <- ggplot(to_plot, aes(x=QUANTILE, y=OR)) + geom_pointrange(aes(ymin = low, ymax = high), size=0.7) +
  scale_y_continuous(name = "Odds ratio (95% CI)", breaks = c(0,1,2,3,4,5,6)) + xlab("Quantiles for PRS (Z-standardized)") +
  geom_errorbar(aes(ymin = low, ymax = high, width = 0.05)) +
  ggtitle("Mars, 2020") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), panel.background = element_rect(fill = "white"), axis.line = element_line(size = 0.5, colour = "black"), panel.grid.major = element_line(color = "grey90"))
ggsave(plot = plotted, filename = 'Scores_PD-DM/plotQuantile_PRS_MPBC_OR_T2Dstatus_Mars.png', width = 5, height = 6, units = 'in', dpi = 600)

# Density plot
ggplot(all_data, aes(x = z_Mars_Score, fill = factor(T2D))) +
  geom_density(alpha = 0.5, size = 0.1) +
  labs(x = "Polygenic Risk Score (PRS)", fill = "T2D") +
  scale_fill_manual(values = c("cyan", "purple"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  ggtitle("PRS Distribution for T2D in MPBC Cohort (GWAS - Mars, 2020)")
ggsave(filename = "Scores_PD-DM/plotDensity_Mars.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Use the four-group classification created above
# Density plot for PRS across disease groups
ggplot(all_data, aes(x = z_Mars_Score, fill = Disease_Group)) +
  geom_density(alpha = 0.5, size = 0.1) +
  labs(x = "Polygenic Risk Score (PRS)", fill = "Disease Group") +
  scale_fill_manual(values = c("skyblue", "salmon", "purple", "green")) + # Customize colors
  theme_minimal() +
  ggtitle("PRS Distribution Across Disease Groups (Mars, 2020)")
# Save the plot
ggsave(filename = "Scores_PD-DM/plotDensity_groups_Mars.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Make predictions to access model performance
# First, run a simple regression model
Model_mars <- glm(T2D ~ Mars_Score + AGE + SEX, data = all_data, family = 'binomial')
all_data$probDisease <- predict(Model_mars, all_data, type = c("response"))
# find the best classification cut-off
roc_obj <- roc(all_data$T2D, all_data$probDisease)
hist(all_data$probDisease, breaks = 30, main = "Predicted Probabilities (GWAS - Mars, 2020)") # use histogram of predicted probabilities
best_threshold <- coords(roc_obj, "best", ret = "threshold") #use pROC to find the probability threshold
all_data$predicted <- ifelse(all_data$probDisease > 0.075, "DISEASE", "CONTROL")
all_data$reported <- ifelse(all_data$T2D == 1, "DISEASE","CONTROL")

# ROC plots
overlayedRocs <- ggplot(all_data, aes(d = T2D, m = probDisease)) + geom_roc(color ="darkblue", labels = FALSE) +
  geom_rocci() + style_roc(theme = theme_gray) + ggtitle(paste("ROC Curve based on Mars, 2020 GWAS (AUC =", round(auc(roc_obj), 2), ")")) + 
  geom_abline(linetype = "dashed") + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "Scores_PD-DM/plotRoc_Mars.png", width = 8, height = 5, units = "in", dpi = 300)

# Show the confusion matrix (specificity and sensitivity)
confMat <- confusionMatrix(data = as.factor(all_data$predicted), reference = as.factor(all_data$reported), positive = "DISEASE")
confMat
# save to an image
# Convert confusion matrix to a table format
confMat_text <- capture.output(print(confMat))
# Save as an image
png("Scores_PD-DM/confusion_matrix_Mars.png", width = 800, height = 600)
grid.table(confMat_text)  # Display text as an image
dev.off()

# Using more visualizations
# Boxplot of PRS scores by T2D status
ggplot(all_data, aes(x = reorder(as.factor(T2D), z_Mars_Score), y = z_Mars_Score, fill = as.factor(T2D))) +
  geom_boxplot() +
  labs(x = "T2D Status", y = "Polygenic Risk Score (PRS)", fill = "T2D Status") +
  scale_fill_manual(values = c("skyblue", "salmon"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  ggtitle("PRS Distribution by T2D Status (Mars, 2020)")
ggsave(filename = "Scores_PD-DM/boxplot_Mars.jpeg", width = 8, height = 5, units = "in", dpi = 300)


# Create another Quantile plot with just the zScores instead of OR
all_data <- all_data %>%
  mutate(quantiles = ntile(z_Mars_Score, 4))  # Divide into 4 groups

# Summarize data for visualization
quantile_summary <- all_data %>%
  group_by(quantiles, T2D) %>%
  tally(name ="count") %>%
  mutate(
    proportion = count / sum(count),  # Proportion within each quantile
    lower_ci = proportion - qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    upper_ci = proportion + qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    .groups = "drop"  # Ungroup after summarising
  )

# Visualize with geom_pointrange
ggplot(quantile_summary, aes(x = factor(quantiles), y = proportion, color = factor(T2D))) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(0.5)) +
  labs(x = "T2D PRS Quantiles", y = "Proportion", color = "T2D Status") +
  scale_color_manual(
    values = c("blue", "red"),
    labels = c("0 (Controls)", "1 (Cases)")
  ) +
  theme_minimal() +
  ggtitle("Proportion of Cases and Controls by T2D PRS Quantiles (GWAS - Mars, 2020)")
ggsave(filename = "Scores_PD-DM/quantile_plot_by_proportion_Mars.jpeg", width = 8, height = 5, units = "in", dpi = 300)

#######################################################################################################
# Lin, 2023
# Violin plots
# Use the CASE_Label column created above
p <- ggplot(all_data, aes(x= reorder(as.factor(Disease_Group), z_Lin_Score), y=z_Lin_Score, fill=as.factor(Disease_Group))) +
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("lightblue", "orange", "blue", "red")) + theme_bw() + 
  ylab("T2D PRS (Z-transformed)") +
  xlab("Disease Status") + 
  theme(legend.position = "none") +
  ggtitle("Distribution of T2D PRS (Lin Score) by Disease Status")
ggsave("Scores_PD-DM/T2D_Lin_PRS.jpeg", dpi = 600, units = "in", height = 6, width = 6)

# Make Quantile plots using ODDS RATIO
all_data$quantile1 <- ifelse(all_data$z_Lin_Score <= quantile(all_data$z_Lin_Score)[2], 1, 0)
all_data$quantile2 <- ifelse(all_data$z_Lin_Score > quantile(all_data$z_Lin_Score)[2] & all_data$z_Lin_Score <= quantile(all_data$z_Lin_Score)[3], 1, 0)
all_data$quantile3 <- ifelse(all_data$z_Lin_Score > quantile(all_data$z_Lin_Score)[3] & all_data$z_Lin_Score <= quantile(all_data$z_Lin_Score)[4], 1, 0)
all_data$quantile4 <- ifelse(all_data$z_Lin_Score > quantile(all_data$z_Lin_Score)[4], 1, 0)
all_data$quantiles <- 1
all_data$quantiles[all_data$quantile2 == 1] <- 2
all_data$quantiles[all_data$quantile3 == 1] <- 3
all_data$quantiles[all_data$quantile4 == 1] <- 4
quintile_T2Dstatus <- glm(as.factor(CASE_Label) ~ as.factor(all_data$quantiles) + AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5, family='binomial', data = all_data)

# Summarize the regression and export a table (to use OR and CI instead of BETA and SE)
summary(quintile_T2Dstatus)

quintile_T2Dstatus_summary <- data.frame(exp(cbind(OR = (coef(quintile_T2Dstatus)), confint(quintile_T2Dstatus))))
quintile_T2Dstatus_summary <- quintile_T2Dstatus_summary[2:4,1:3] #To get only OR and CI for the quantiles
names(quintile_T2Dstatus_summary) <- c('OR','2.5%','97.5%')
quintile_T2Dstatus_summary$QUANTILE <- c('2','3','4')
quintile_T2Dstatus_summary[4,] <- c(0,0,0,'1')
quintile_T2Dstatus_summary_sorted <- quintile_T2Dstatus_summary[order(quintile_T2Dstatus_summary$QUANTILE),]
write.table(quintile_T2Dstatus_summary_sorted, 'Scores_PD-DM/quantile_table_MPBC_T2D_PRS_OR_20250214_Lin.csv', quote = F, row.names = F, sep = ',')
# Now, make plot
to_plot <- read.table('Scores_PD-DM/quantile_table_MPBC_T2D_PRS_OR_20250214_Lin.csv', header = T, sep = ',')
to_plot$low <- to_plot$X2.5.
to_plot$high <- to_plot$X97.5.
plotted <- ggplot(to_plot, aes(x=QUANTILE, y=OR)) + geom_pointrange(aes(ymin = low, ymax = high), size=0.7) +
  scale_y_continuous(name = "Odds ratio (95% CI)", breaks = c(0,1,2,3,4,5,6)) + xlab("Quantiles for PRS (Z-standardized)") +
  geom_errorbar(aes(ymin = low, ymax = high, width = 0.05)) +
  ggtitle("Lin, 2023") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16), panel.background = element_rect(fill = "white"), axis.line = element_line(size = 0.5, colour = "black"), panel.grid.major = element_line(color = "grey90"))
ggsave(plot = plotted, filename = 'Scores_PD-DM/plotQuantile_PRS_MPBC_OR_T2Dstatus_Lin.png', width = 5, height = 6, units = 'in', dpi = 600)

# Density plot
ggplot(all_data, aes(x = z_Lin_Score, fill = factor(T2D))) +
  geom_density(alpha = 0.5, size = 0.1) +
  labs(x = "Polygenic Risk Score (PRS)", fill = "T2D") +
  scale_fill_manual(values = c("cyan", "purple"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  ggtitle("PRS Distribution for T2D in MPBC Cohort (GWAS - Lin, 2023)")
ggsave(filename = "Scores_PD-DM/plotDensity_Lin.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Use the four-group classification created above
# Density plot for PRS across disease groups
ggplot(all_data, aes(x = z_Lin_Score, fill = Disease_Group)) +
  geom_density(alpha = 0.5, size = 0.1) +
  labs(x = "Polygenic Risk Score (PRS)", fill = "Disease Group") +
  scale_fill_manual(values = c("skyblue", "salmon", "purple", "green")) + # Customize colors
  theme_minimal() +
  ggtitle("PRS Distribution Across Disease Groups (Lin, 2023)")
# Save the plot
ggsave(filename = "Scores_PD-DM/plotDensity_groups_Lin.jpeg", width = 8, height = 5, units = "in", dpi = 300)

# Make predictions to access model performance
# First, run a simple regression model
Model_lin <- glm(T2D ~ Lin_Score + AGE + SEX, data = all_data, family = 'binomial')
all_data$probDisease <- predict(Model_lin, all_data, type = c("response"))
# find the best classification cut-off
roc_obj <- roc(all_data$T2D, all_data$probDisease)
hist(all_data$probDisease, breaks = 30, main = "Predicted Probabilities (GWAS - Lin, 2023)") # use histogram of predicted probabilities
best_threshold <- coords(roc_obj, "best", ret = "threshold") #use pROC to find the probability threshold
all_data$predicted <- ifelse(all_data$probDisease > 0.0678, "DISEASE", "CONTROL")
all_data$reported <- ifelse(all_data$T2D == 1, "DISEASE","CONTROL")

# ROC plots
overlayedRocs <- ggplot(all_data, aes(d = T2D, m = probDisease)) + geom_roc(color ="darkblue", labels = FALSE) +
  geom_rocci() + style_roc(theme = theme_gray) + ggtitle(paste("ROC Curve based on Lin, 2023 GWAS (AUC =", round(auc(roc_obj), 2), ")")) + 
  geom_abline(linetype = "dashed") + theme_bw() + scale_fill_brewer(palette="Spectral")
ggsave(plot = overlayedRocs, filename = "Scores_PD-DM/plotRoc_Lin.png", width = 8, height = 5, units = "in", dpi = 300)

# Show the confusion matrix (specificity and sensitivity)
confMat <- confusionMatrix(data = as.factor(all_data$predicted), reference = as.factor(all_data$reported), positive = "DISEASE")
confMat
# save to an image
# Convert confusion matrix to a table format
confMat_text <- capture.output(print(confMat))
# Save as an image
png("Scores_PD-DM/confusion_matrix_Lin.png", width = 800, height = 600)
grid.table(confMat_text)  # Display text as an image
dev.off()

# Using more visualizations
# Boxplot of PRS scores by T2D status
ggplot(all_data, aes(x = reorder(as.factor(T2D), z_Lin_Score), y = z_Lin_Score, fill = as.factor(T2D))) +
  geom_boxplot() +
  labs(x = "T2D Status", y = "Polygenic Risk Score (PRS)", fill = "T2D Status") +
  scale_fill_manual(values = c("skyblue", "salmon"), labels = c("Controls", "Cases")) +
  theme_minimal() +
  ggtitle("PRS Distribution by T2D Status (Lin, 2023)")
ggsave(filename = "Scores_PD-DM/boxplot_Lin.jpeg", width = 8, height = 5, units = "in", dpi = 300)


# Create another Quantile plot with just the zScores instead of OR
all_data <- all_data %>%
  mutate(quantiles = ntile(z_Lin_Score, 4))  # Divide into 4 groups

# Summarize data for visualization
quantile_summary <- all_data %>%
  group_by(quantiles, T2D) %>%
  tally(name ="count") %>%
  mutate(
    proportion = count / sum(count),  # Proportion within each quantile
    lower_ci = proportion - qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    upper_ci = proportion + qnorm(0.975) * sqrt((proportion * (1 - proportion)) / sum(count)),
    .groups = "drop"  # Ungroup after summarising
  )

# Visualize with geom_pointrange
ggplot(quantile_summary, aes(x = factor(quantiles), y = proportion, color = factor(T2D))) +
  geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(0.5)) +
  labs(x = "T2D PRS Quantiles", y = "Proportion", color = "T2D Status") +
  scale_color_manual(
    values = c("blue", "red"),
    labels = c("0 (Controls)", "1 (Cases)")
  ) +
  theme_minimal() +
  ggtitle("Proportion of Cases and Controls by T2D PRS Quantiles (GWAS - Lin, 2023)")
ggsave(filename = "Scores_PD-DM/quantile_plot_by_proportion_Lin.jpeg", width = 8, height = 5, units = "in", dpi = 300)

