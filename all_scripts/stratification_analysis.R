setwd("/Users/Elijah/Documents/Thesis_Project_2024")
list.files()


# Load the necessary packages 
library(tidyr)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(plotROC)
library(caret)

# Logistic regressions for all environmental factors
#########################################################################
# Smoking
# Read in the smoking data with the covariate file 
smoking_data <- read.table("07_stratification/smoking_ever_never.txt", header = T)
#rename the FID column
colnames(smoking_data)[1] <- "FID"
colnames(smoking_data)[2] <- "PD"
covariates <- read.table("covariates.txt", header = T)
# Add only PCs from covariate file to the smoking_data
# Select PCs (columns 10 to 14) and merge with smoking_data
smoking <- merge(smoking_data, covariates[, c(1, 10:14)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(smoking$SCORE[smoking$PD == 0])
sdControls <- sd(smoking$SCORE[smoking$PD == 0])
smoking$zSCORE <- (smoking$SCORE - meanControls)/sdControls

# logistic Regression model for smoking adjusted by covariates
smoking_model_covs <- glm(PD ~ zSCORE + Smoking_Ever + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                          family="binomial", data = smoking)
summary(smoking_model_covs)

# Modeling the impact of interaction between PRS and smoking on PD diagnosis
smoking_model_interaction <- glm(PD ~ zSCORE * Smoking_Ever + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                          family="binomial", data = smoking)
summary(smoking_model_interaction)

# test of model fitness
anova(smoking_model_covs, smoking_model_interaction, test= "Chisq")
# Using McFadden’s Pseudo-R-squared
library(pscl)
pR2(smoking_model_covs)
pR2(smoking_model_interaction)
# check confidence interval for interaction term
confint(smoking_model_interaction)

#########################################################################
# Snus usage
# Read in the snus data
snus_data <- read.table("07_stratification/snus_ever_never.txt", header = T)
#rename the FID column
colnames(snus_data)[1] <- "FID"
colnames(snus_data)[2] <- "PD"

# Add only PCs from covariate file to the snus_data
# Select PCs (columns 10 to 14) and merge with snus_data
snus <- merge(snus_data, covariates[, c(1, 10:14)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(snus$SCORE[snus$PD == 0])
sdControls <- sd(snus$SCORE[snus$PD == 0])
snus$zSCORE <- (snus$SCORE - meanControls)/sdControls

# logistic Regression model for snus usage adjusted by covariates
model_snus_covs <- glm(PD ~ zSCORE + Snus_Ever + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5,
                       family = binomial, data = snus)
summary(model_snus_covs)

# Modeling the impact of interaction between PRS and snus usage on PD diagnosis
model_snus_interaction <- glm(PD ~ zSCORE * Snus_Ever + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                                 family="binomial", data = snus)
summary(model_snus_interaction)

# test of model fitness
anova(model_snus_covs, model_snus_interaction, test= "Chisq")
# Using McFadden’s Pseudo-R-squared
pR2(model_snus_covs)
pR2(model_snus_interaction)
# check confidence interval for interaction term
confint(model_snus_interaction)

#########################################################################
# Pesticide exposure
# Read in the pesticides exposure data 
pesticides_data <- read.table("07_stratification/pesticides_ever_never.txt", header = T)
#rename the FID column
colnames(pesticides_data)[1] <- "FID"
colnames(pesticides_data)[2] <- "PD"

# Add only PCs from covariate file to the pesticides_data
# Select PCs (columns 10 to 14) and merge with pesticides_data
pesticides <- merge(pesticides_data, covariates[, c(1, 10:14)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(pesticides$SCORE[pesticides$PD == 0])
sdControls <- sd(pesticides$SCORE[pesticides$PD == 0])
pesticides$zSCORE <- (pesticides$SCORE - meanControls)/sdControls

# logistic Regression model for pesticides exposure adjusted by covariates
model_pesticides_covs <- glm(PD ~ zSCORE + Pesticides_Ever + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5,
                       family = binomial, data = pesticides)
summary(model_pesticides_covs)

# Modeling the impact of interaction between PRS and pesticides exposure on PD diagnosis
model_pesticides_interaction <- glm(PD ~ zSCORE * Pesticides_Ever + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                              family="binomial", data = pesticides)
summary(model_pesticides_interaction)

# test of model fitness
anova(model_pesticides_covs, model_pesticides_interaction, test= "Chisq")
# Using McFadden’s Pseudo-R-squared
pR2(model_pesticides_covs)
pR2(model_pesticides_interaction)
# check confidence interval for interaction term
confint(model_pesticides_interaction)

#########################################################################
# Caffeine consumption
# Read in the caffeine consumption data 
caffeine_consumption_data <- read.table("07_stratification/caffeine_ever_never.txt", header = T)
#rename the FID column
colnames(caffeine_consumption_data)[1] <- "FID"
colnames(caffeine_consumption_data)[2] <- "PD"

# Add only PCs from covariate file to the caffeine_consumption_data
# Select PCs (columns 10 to 14) and merge with caffeine_consumption_data
caffeine1 <- merge(caffeine_consumption_data, covariates[, c(1, 10:14)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(caffeine1$SCORE[caffeine1$PD == 0])
sdControls <- sd(caffeine1$SCORE[caffeine1$PD == 0])
caffeine1$zSCORE <- (caffeine1$SCORE - meanControls)/sdControls

# logistic Regression model for caffeine consumption adjusted by covariates
model_caffeine1_covs <- glm(PD ~ zSCORE + Caffeine_Ever + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5,
                             family = binomial, data = caffeine1)
summary(model_caffeine1_covs)

# Modeling the impact of interaction between PRS and caffeine consumption on PD diagnosis
model_caffeine1_interaction <- glm(PD ~ zSCORE * Caffeine_Ever + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                                    family="binomial", data = caffeine1)
summary(model_caffeine1_interaction)

# test of model fitness
anova(model_caffeine1_covs, model_caffeine1_interaction, test= "Chisq")
# Using McFadden’s Pseudo-R-squared
pR2(model_caffeine1_covs)
pR2(model_caffeine1_interaction)
# check confidence interval for interaction term
confint(model_caffeine1_interaction)

#########################################################################
# Caffeine consumption Level
# Read in the caffeine consumption data 
caffeine_consumption_level <- read.table("07_stratification/caffeine_low_vs_high.txt", header = T)
#rename the FID column
colnames(caffeine_consumption_level)[1] <- "FID"
colnames(caffeine_consumption_level)[2] <- "PD"

# Add only PCs from covariates file to the caffeine_consumption_data
# Select PCs (columns 10 to 14) and merge with caffeine_consumption_data
caffeine2 <- merge(caffeine_consumption_level, covariates[, c(1, 10:14)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(caffeine2$SCORE[caffeine2$PD == 0])
sdControls <- sd(caffeine2$SCORE[caffeine2$PD == 0])
caffeine2$zSCORE <- (caffeine2$SCORE - meanControls)/sdControls

# Convert Caffeine_Consumption_Level to a factor
caffeine2$Caffeine_Consumption_Level <- factor(caffeine2$Caffeine_Consumption_Level)
# Relevel to set "LOW" as the reference level
caffeine2$Caffeine_Consumption_Level <- relevel(caffeine2$Caffeine_Consumption_Level, ref = "LOW")

# Now, logistic Regression model for caffeine consumption adjusted by covariates
model_caffeine2_covs <- glm(PD ~ zSCORE + Caffeine_Consumption_Level + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5,
                            family = binomial, data = caffeine2)
summary(model_caffeine2_covs)

# Modeling the impact of interaction between PRS and caffeine consumption on PD diagnosis
model_caffeine2_interaction <- glm(PD ~ zSCORE * Caffeine_Consumption_Level + sex + education + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                                   family="binomial", data = caffeine2)
summary(model_caffeine2_interaction)

# test of model fitness
anova(model_caffeine2_covs, model_caffeine2_interaction, test= "Chisq")
# Using McFadden’s Pseudo-R-squared
pR2(model_caffeine2_covs)
pR2(model_caffeine2_interaction)
# check confidence interval for interaction term
confint(model_caffeine2_interaction)

#########################################################################
# Combining all three factors
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

# Filter the whole dataset to ensure all individuals (patients and controls) have complete data across all
# relevant factors
filtered_data <- data[complete.cases(data[, c("PD", 
                                              "Smoking_Ever", 
                                              "Snus_Ever", 
                                              "Pesticides_Ever", 
                                              "Caffeine_Ever", 
                                              "Caffeine_Consumption_Level")]), ]


# check the number of individuals removed
nrow(data) - nrow(filtered_data)
# Convert Caffeine_Consumption_Level to a factor
data$Caffeine_Consumption_Level <- factor(data$Caffeine_Consumption_Level)
# Relevel to set "LOW" as the reference level
data$Caffeine_Consumption_Level <- relevel(data$Caffeine_Consumption_Level, ref = "LOW")


# Now, logistic Regression model for all three factors covariates
model_covs <- glm(PD ~ zSCORE + Smoking_Ever + Snus_Ever + Pesticides_Ever + Caffeine_Ever + 
                    Caffeine_Consumption_Level + sex + education + Age_Inclusion + 
                    PC1 + PC2 + PC3 + PC4 + PC5, family = binomial, data = data)
summary(model_covs)

# Three way interact model for all factors
three_way_interact_model <- glm(PD ~ zSCORE * Smoking_Ever * Pesticides_Ever * Caffeine_Consumption_Level 
                                + sex + education + Age_Inclusion + PC3, family = binomial, data = data)
summary(three_way_interact_model)

# Three way interact model for all protective factors
protective_interact_model <- glm(PD ~ zSCORE * Smoking_Ever * Caffeine_Consumption_Level 
                                + sex + education + Age_Inclusion + PC3, family = binomial, data = data)
summary(protective_interact_model)

# Checking for multicollinearity
library(car)

# test of model fitness using McFadden’s Pseudo-R-squared
pR2(model_covs)
pR2(three_way_interact_model)
pR2(protective_interact_model)
# check confidence interval for interaction term
confint(model_covs)

# Visualizations
#########################################################################
# Visualizations
#Smoking
# Step 1: Create a new column combining snus usage and PD status
smoking$Smoking_PD_Status <- with(smoking, paste(Smoking_Ever, ifelse(PD == 1, "PD", "Control"), sep = "_"))

# Step 2: Plot the violin plot with the new variable
p <- ggplot(smoking, aes(x = reorder(as.factor(Smoking_PD_Status), zSCORE), y = zSCORE, 
                         fill = as.factor(Smoking_PD_Status))) + geom_violin(trim = FALSE)
p2 <- p + geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("PRS Distribution by Smoking and PD Status") +
  scale_fill_manual(values = c("blue", "red", "blue", "red")) + 
  theme_bw() + 
  ylab("PD PRS (zSCORE)") + 
  xlab("Smoking and PD Status") + 
  theme(legend.position = "none")
# Display the plot
p2
ggsave("Smoking_and_PD_Status.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot

# visualize interaction between PRS and smoking in PD
ggplot(smoking, aes(x = zSCORE, y = PD, color = Smoking_Ever)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  labs(title = "Interaction Between PRS and Smoking", x = "PRS Score", y = "Probability of PD") +
  theme_minimal()
ggsave("Smoking_interaction.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot

# Alternatively, use the interaction.plot() function
interaction.plot(x.factor = as.factor(smoking$PD), 
                 trace.factor = as.factor(smoking$Smoking_Ever), 
                 response = smoking$zSCORE, 
                 fun = mean, xlab="PD Status", 
                 ylab="PRS (zSCORE)", trace.label="Ever Smoked?", 
                 col=c("green","red"), 
                 lty=1, lwd=2 )

#Snus
# Step 1: Create a new column combining snus usage and PD status
snus$Snus_PD_Status <- with(snus, paste(Snus_Ever, ifelse(PD == 1, "PD", "Control"), sep = "_"))

# Step 2: Plot the violin plot with the new variable
p <- ggplot(snus, aes(x = reorder(as.factor(Snus_PD_Status), zSCORE), y = zSCORE, 
                         fill = as.factor(Snus_PD_Status))) + geom_violin(trim = FALSE)
p2 <- p + geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("PRS Distribution by Snus and PD Status") +
  scale_fill_manual(values = c("blue", "red", "blue", "red")) + 
  theme_bw() + 
  ylab("PD PRS (zSCORE)") + 
  xlab("Snus Usage and PD Status") + 
  theme(legend.position = "none")
# Display the plot
p2
ggsave("Snus_and_PD_Status.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot

# visualize interaction between PRS and snus usage in PD
ggplot(snus, aes(x = zSCORE, y = PD, color = Snus_Ever)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  labs(title = "Interaction Between PRS and Snus Usage", x = "PRS Score", y = "Probability of PD") +
  theme_minimal()
ggsave("Snus_interaction.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot

# Alternatively, use the interaction.plot() function
interaction.plot(x.factor = as.factor(snus$PD), 
                 trace.factor = as.factor(snus$Snus_Ever), 
                 response = snus$zSCORE, 
                 fun = mean, xlab="PD Status", 
                 ylab="PRS (zSCORE)", trace.label="Ever Used Snus?", 
                 col=c("green","red"), 
                 lty=1, lwd=2 )

#Pesticides Exposure
# Step 1: Create a new column combining snus usage and PD status
pesticides$Pestices_PD_Status <- with(pesticides, paste(Pesticides_Ever, ifelse(PD == 1, "PD", "Control"), sep = "_"))

# Step 2: Plot the violin plot with the new variable
p <- ggplot(pesticides, aes(x = reorder(as.factor(Pestices_PD_Status), zSCORE), y = zSCORE, 
                      fill = as.factor(Pestices_PD_Status))) + geom_violin(trim = FALSE)
p2 <- p + geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("PRS Distribution by Pesticide Exposure and PD Status") +
  scale_fill_manual(values = c("blue", "red", "blue", "red")) + 
  theme_bw() + 
  ylab("PD PRS (zSCORE)") + 
  xlab("Pesticides Exposure and PD Status") + 
  theme(legend.position = "none")
# Display the plot
p2
ggsave("Pesticides_Exposure_and_PD_Status.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot

# visualize interaction between PRS and Pesticides Exposure in PD
ggplot(pesticides, aes(x = zSCORE, y = PD, color = Pesticides_Ever)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  labs(title = "Interaction Between PRS and Pesticides Exposure", x = "PRS Score", y = "Probability of PD") +
  theme_minimal()
ggsave("Pesticides_interaction.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot


# Alternatively, use the interaction.plot() function
interaction.plot(x.factor = as.factor(pesticides$PD), 
                 trace.factor = as.factor(pesticides$Pesticides_Ever), 
                 response = pesticides$zSCORE, 
                 fun = mean, xlab="PD Status", 
                 ylab="PRS (zSCORE)", trace.label="Ever Used Pesticides?", 
                 col=c("green","red"), 
                 lty=1, lwd=2 )

#Caffeine Consumption
# Step 1: Create a new column combining Caffeine Consumption and PD status
caffeine1$Caffeine_PD_Status <- with(caffeine1, paste(Caffeine_Ever, ifelse(PD == 1, "PD", "Control"), sep = "_"))

# Step 2: Plot the violin plot with the new variable
p <- ggplot(caffeine1, aes(x = reorder(as.factor(Caffeine_PD_Status), zSCORE), y = zSCORE, 
                            fill = as.factor(Caffeine_PD_Status))) + geom_violin(trim = FALSE)
p2 <- p + geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("PRS Distribution by Caffeine Consumption and PD Status") +
  scale_fill_manual(values = c("blue", "red", "blue", "red")) + 
  theme_bw() + 
  ylab("PD PRS (zSCORE)") + 
  xlab("Caffeine Consumption and PD Status") + 
  theme(legend.position = "none")
# Display the plot
p2
ggsave("Caffeine_Consumption_and_PD_Status.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot

# visualize interaction between PRS and Caffeine Consumption in PD
ggplot(caffeine1, aes(x = zSCORE, y = PD, color = Caffeine_Ever)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  labs(title = "Interaction Between PRS and Caffeine Consumption", x = "PRS Score", y = "Probability of PD") +
  theme_minimal()
ggsave("Caffeine_Consumption_interaction.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot


# Alternatively, use the interaction.plot() function
interaction.plot(x.factor = as.factor(caffeine1$PD), 
                 trace.factor = as.factor(caffeine1$Caffeine_Ever), 
                 response = caffeine1$zSCORE, 
                 fun = mean, xlab="PD Status", 
                 ylab="PRS (zSCORE)", trace.label="Ever Drank Coffee/Tea?", 
                 col=c("green","red"), 
                 lty=1, lwd=2 )

#Caffeine Consumption Level
# Step 1: Create a new column combining Caffeine Consumption Level and PD status
caffeine2$Caffeine_Level_PD_Status <- with(caffeine2, paste(Caffeine_Consumption_Level, ifelse(PD == 1, "PD", "Control"), sep = "_"))

# Step 2: Plot the violin plot with the new variable
p <- ggplot(caffeine2, aes(x = reorder(as.factor(Caffeine_Level_PD_Status), zSCORE), y = zSCORE, 
                           fill = as.factor(Caffeine_Level_PD_Status))) + geom_violin(trim = FALSE)
p2 <- p + geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("PRS Distribution by Caffeine Consumption Level and PD Status") +
  scale_fill_manual(values = c("blue", "red", "blue", "red")) + 
  theme_bw() + 
  ylab("PD PRS (zSCORE)") + 
  xlab("Caffeine Consumption Level and PD Status") + 
  theme(legend.position = "none")
# Display the plot
p2
ggsave("Caffeine_Consumption_Level_and_PD_Status.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot

# visualize interaction between PRS and Caffeine Consumption in PD
ggplot(caffeine2, aes(x = zSCORE, y = PD, color = Caffeine_Consumption_Level)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  labs(title = "Interaction Between PRS and Caffeine Consumption Level", x = "PRS Score", y = "Probability of PD") +
  theme_minimal()
ggsave("Caffeine_Consumption_Level_interaction.jpeg", dpi = 600, units = "in", height = 6, width = 6) # save plot

# Alternatively, use the interaction.plot() function
interaction.plot(x.factor = as.factor(caffeine2$PD), 
                 trace.factor = as.factor(caffeine2$Caffeine_Consumption_Level), 
                 response = caffeine2$zSCORE, 
                 fun = mean, xlab="PD Status", 
                 ylab="PRS (zSCORE)", trace.label="Do You Drink Lots of Coffee/Tea?", 
                 col=c("green","red"), 
                 lty=1, lwd=2 )


# Make tables to know the number of individuals in each group
table(covariates$SEX[covariates$PHENO == 0]) # sex of controls
table(covariates$SEX[covariates$PHENO == 1]) #sex of cases
table(smoking$Smoking_PD_Status)
table(snus$Snus_PD_Status)
table(pesticides$Pestices_PD_Status)
table(caffeine1$Caffeine_PD_Status)
table(caffeine2$Caffeine_Level_PD_Status)

# Mean of ages
# Calculate the mean age for cases (PHENO == 1)
mean(covariates$AGE[covariates$PHENO == 1], na.rm = TRUE)
# Get a summary of AGE for cases
summary(covariates$AGE[covariates$PHENO == 1])
# Calculate the mean age for controls(PHENO == 0)
mean(covariates$AGE[covariates$PHENO == 0], na.rm = TRUE)
# Get a summary of AGE for controls
summary(covariates$AGE[covariates$PHENO == 0])

# Calculate the mean age at diagnoses in the study for cases(PHENO == 1)
mean(covariates$AAD[covariates$PHENO == 1], na.rm = TRUE)
# Get a summary of AGE for cases
summary(covariates$AAD[covariates$PHENO == 1])
# Standard deviation of the age at diagnoses for the cases
sd(covariates$AAD[covariates$PHENO == 1], na.rm = TRUE)


# Standard deviation of ages
# Standard deviation of AGE for cases (PHENO == 1)
sd(covariates$AGE[covariates$PHENO == 1], na.rm = TRUE)
# Standard deviation of AGE for controls (PHENO == 0)
sd(covariates$AGE[covariates$PHENO == 0], na.rm = TRUE)

# to investigate the mean PRS score in groups of individuals with/without PD that smoke, use snus, drink coffee, and are exposed to pesticides
mean(caffeine2$zSCORE[caffeine2$Caffeine_Consumption_Level == "HIGH" & caffeine2$PD == 1])
mean(caffeine2$zSCORE[caffeine2$Caffeine_Consumption_Level == "LOW" & caffeine2$PD == 1])
mean(caffeine2$zSCORE[caffeine2$Caffeine_Consumption_Level == "HIGH" & caffeine2$PD == 0])
mean(caffeine2$zSCORE[caffeine2$Caffeine_Consumption_Level == "LOW" & caffeine2$PD == 0])

