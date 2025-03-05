# Stratification of T2D PRS based on environmental exposures (Smoking, Snus, Caffeine, and Pesticides)
# Since the PRS from the Ge, 2022 GWAS produces a better-performing model across all metrics, we will use it for the stratification analysis
# Install and load dplyr
install.packages("dplyr")
library(dplyr)
library(ggplot2)
# set working directory
setwd("/Users/Elijah/Documents/Thesis_Project_2024")
#########################################################################
# Smoking
# Read in the smoking data and covariate file
smoking_data <- read.table("07_stratification/smoking_ever_never.txt", header = T)
#rename the first and second columns
colnames(smoking_data)[1:2] <- c("FID", "PD")
# read the Ge, 2022 polygenic risk scores
Ge_Score <- read.table("PRS_T2D/Ge_PRS.profile", header = T)
covariates <- read.table("Scores_PD-DM/new_covariates.txt", header = T)
colnames(covariates)[6] <- "PD" # Rename the column to PD
# Select PCs (columns 10 to 14) from covariates file and merge with smoking_data and Ge_Score
smoking <- smoking_data[, c(1:7,9)] %>%
  inner_join(Ge_Score[, c(1, 6)], by = "FID") %>%
  inner_join(covariates[, c(1, 10:14, 20)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(smoking$SCORE[smoking$T2D == 0])
sdControls <- sd(smoking$SCORE[smoking$T2D == 0])
smoking$zSCORE <- (smoking$SCORE - meanControls)/sdControls

# logistic Regression model for smoking adjusted by covariates
smoking_model_covs <- glm(T2D ~ zSCORE + Smoking_Ever + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                          family="binomial", data = smoking)
summary(smoking_model_covs)

# Modeling the impact of interaction between PRS and smoking on T2D diagnosis
smoking_model_interaction <- glm(T2D ~ zSCORE * Smoking_Ever + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                                 family="binomial", data = smoking)
summary(smoking_model_interaction)

#######################################################################
# Snus
# Read in the snus data and covariate file
snus_data <- read.table("07_stratification/snus_ever_never.txt", header = T)
#rename the first and second columns
colnames(snus_data)[1:2] <- c("FID", "PD")
# Select PCs (columns 10 to 14) from covariates file and merge with snus_data and Ge_Score
snus <- snus_data[, c(1:7,9)] %>%
  inner_join(Ge_Score[, c(1, 6)], by = "FID") %>%
  inner_join(covariates[, c(1, 10:14, 20)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(snus$SCORE[snus$T2D == 0])
sdControls <- sd(snus$SCORE[snus$T2D == 0])
snus$zSCORE <- (snus$SCORE - meanControls)/sdControls

# logistic Regression model for smoking adjusted by covariates
snus_model_covs <- glm(T2D ~ zSCORE + Snus_Ever + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                          family="binomial", data = snus)
summary(snus_model_covs)

# Modeling the impact of interaction between PRS and smoking on T2D diagnosis
snus_model_interaction <- glm(T2D ~ zSCORE * Snus_Ever + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                                 family="binomial", data = snus)
summary(snus_model_interaction)

#######################################################################
# Pesticides exposure
# Read in the pesticides data and covariate file
pesticides_data <- read.table("07_stratification/pesticides_ever_never.txt", header = T)
#rename the first and second columns
colnames(pesticides_data)[1:2] <- c("FID", "PD")
# Select PCs (columns 10 to 14) from covariates file and merge with snus_data and Ge_Score
pesticides <- pesticides_data[, c(1:7,9)] %>%
  inner_join(Ge_Score[, c(1, 6)], by = "FID") %>%
  inner_join(covariates[, c(1, 10:14, 20)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(pesticides$SCORE[pesticides$T2D == 0])
sdControls <- sd(pesticides$SCORE[pesticides$T2D == 0])
pesticides$zSCORE <- (pesticides$SCORE - meanControls)/sdControls

# logistic Regression model for pesticides exposure adjusted by covariates
model_pesticides_covs <- glm(T2D ~ zSCORE + Pesticides_Ever + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5,
                             family = binomial, data = pesticides)
summary(model_pesticides_covs)

# Modeling the impact of interaction between PRS and pesticides exposure on T2D diagnosis
model_pesticides_interaction <- glm(T2D ~ zSCORE * Pesticides_Ever + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                                    family="binomial", data = pesticides)
summary(model_pesticides_interaction)

#######################################################################
# Caffeiene consumption
# Read in the Caffeiene consumption data
caffeine1_data <- read.table("07_stratification/caffeine_ever_never.txt", header = T)
#rename the first and second columns
colnames(caffeine1_data)[1:2] <- c("FID", "PD")
# Select PCs (columns 10 to 14) from covariates file and merge with snus_data and Ge_Score
caffeine1<- caffeine1_data[, c(1:4, 14:15,17)] %>%
  inner_join(Ge_Score[, c(1, 6)], by = "FID") %>%
  inner_join(covariates[, c(1, 10:14, 20)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(caffeine1$SCORE[caffeine1$T2D == 0])
sdControls <- sd(caffeine1_data$SCORE[caffeine1$T2D == 0])
caffeine1$zSCORE <- (caffeine1$SCORE - meanControls)/sdControls

# logistic Regression model for caffeine consumption adjusted by covariates
model_caffeine1_covs <- glm(T2D ~ zSCORE + Caffeine_Ever + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5,
                            family = binomial, data = caffeine1)
summary(model_caffeine1_covs)

# Modeling the impact of interaction between PRS and caffeine consumption on PD diagnosis
model_caffeine1_interaction <- glm(T2D ~ zSCORE * Caffeine_Ever + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                                   family="binomial", data = caffeine1)
summary(model_caffeine1_interaction)

#######################################################################
# Caffeiene consumption level 
# Read in the Caffeiene consumption level data 
caffeine2_data <- read.table("07_stratification/caffeine_low_vs_high.txt", header = T)
#rename the first and second columns
colnames(caffeine2_data)[1:2] <- c("FID", "PD")
# Select PCs (columns 10 to 14) from covariates file and merge with snus_data and Ge_Score
caffeine2<- caffeine2_data[, c(1:4, 14:15,17)] %>%
  inner_join(Ge_Score[, c(1, 6)], by = "FID") %>%
  inner_join(covariates[, c(1, 10:14, 20)], by = "FID")

# Normalization of PRS to Z-score
meanControls <- mean(caffeine2$SCORE[caffeine2$T2D == 0])
sdControls <- sd(caffeine2$SCORE[caffeine2$T2D == 0])
caffeine2$zSCORE <- (caffeine2$SCORE - meanControls)/sdControls

# Convert Caffeine_Consumption_Level to a factor
caffeine2$Caffeine_Consumption_Level <- factor(caffeine2$Caffeine_Consumption_Level)
# Relevel to set "LOW" as the reference level
caffeine2$Caffeine_Consumption_Level <- relevel(caffeine2$Caffeine_Consumption_Level, ref = "LOW")

# logistic Regression model for caffeine consumption adjusted by covariates
model_caffeine2_covs <- glm(T2D ~ zSCORE + Caffeine_Consumption_Level + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5,
                            family = binomial, data = caffeine2)
summary(model_caffeine2_covs)

# Modeling the impact of interaction between PRS and caffeine consumption on PD diagnosis
model_caffeine2_interaction <- glm(T2D ~ zSCORE * Caffeine_Consumption_Level + sex + Age_Inclusion + PC1 + PC2 + PC3 + PC4 + PC5, 
                                   family="binomial", data = caffeine2)
summary(model_caffeine2_interaction)

######################################################################
# Get statistics sex
table(covariates$SEX[covariates$PD == 1 & covariates$T2D == 0]) # sex for PD only
table(covariates$SEX[covariates$PD == 0 & covariates$T2D == 1]) # sex for T2D only
table(covariates$SEX[covariates$PD == 1 & covariates$T2D == 1]) # sex for PD + T2D
table(covariates$SEX[covariates$PD == 0 & covariates$T2D == 0]) # sex for healthy controls

# for smoking
table(smoking$Smoking_Ever[smoking$PD == 1 & smoking$T2D == 0])
table(smoking$Smoking_Ever[smoking$PD == 0 & smoking$T2D == 1])
table(smoking$Smoking_Ever[smoking$PD == 1 & smoking$T2D == 1])
table(smoking$Smoking_Ever[smoking$PD == 0 & smoking$T2D == 0])

# for snusing
table(snus$Snus_Ever[snus$PD == 1 & snus$T2D == 0])
table(snus$Snus_Ever[snus$PD == 0 & snus$T2D == 1])
table(snus$Snus_Ever[snus$PD == 1 & snus$T2D == 1])
table(snus$Snus_Ever[snus$PD == 0 & snus$T2D == 0])

# for pesticides
table(pesticides$Pesticides_Ever[pesticides$PD == 1 & pesticides$T2D == 0])
table(pesticides$Pesticides_Ever[pesticides$PD == 0 & pesticides$T2D == 1])
table(pesticides$Pesticides_Ever[pesticides$PD == 1 & pesticides$T2D == 1])
table(pesticides$Pesticides_Ever[pesticides$PD == 0 & pesticides$T2D == 0])

# for caffeine consumption
table(caffeine1$Caffeine_Ever[caffeine1$PD == 1 & caffeine1$T2D == 0])
table(caffeine1$Caffeine_Ever[caffeine1$PD == 0 & caffeine1$T2D == 1])
table(caffeine1$Caffeine_Ever[caffeine1$PD == 1 & caffeine1$T2D == 1])
table(caffeine1$Caffeine_Ever[caffeine1$PD == 0 & caffeine1$T2D == 0])

# for caffeine consumption level
table(caffeine2$Caffeine_Consumption_Level[caffeine2$PD == 1 & caffeine2$T2D == 0])
table(caffeine2$Caffeine_Consumption_Level[caffeine2$PD == 0 & caffeine2$T2D == 1])
table(caffeine2$Caffeine_Consumption_Level[caffeine2$PD == 1 & caffeine2$T2D == 1])
table(caffeine2$Caffeine_Consumption_Level[caffeine2$PD == 0 & caffeine2$T2D == 0])

# get stats for AAD (PD only)
mean(covariates$AAD[covariates$PD == 1 & covariates$T2D == 0], na.rm = TRUE) 
sd(covariates$AAD[covariates$PD == 1 & covariates$T2D == 0], na.rm = TRUE)
summary(covariates$AAD[covariates$PD == 1 & covariates$T2D == 0], na.rm = TRUE)

# AAD for (T2D only)
mean(covariates$AAD[covariates$PD == 0 & covariates$T2D == 1], na.rm = TRUE) 
sd(covariates$AAD[covariates$PD == 0 & covariates$T2D == 1], na.rm = TRUE)
summary(covariates$AAD[covariates$PD == 0 & covariates$T2D == 1], na.rm = TRUE)

# AAD for (T2D + PD)
mean(covariates$AAD[covariates$PD == 1 & covariates$T2D == 1], na.rm = TRUE) 
sd(covariates$AAD[covariates$PD == 1 & covariates$T2D == 1], na.rm = TRUE)
summary(covariates$AAD[covariates$PD == 1 & covariates$T2D == 1], na.rm = TRUE)

# Now, AGE
# PD only
mean(covariates$AGE[covariates$PD == 1 & covariates$T2D == 0], na.rm = TRUE) 
sd(covariates$AGE[covariates$PD == 1 & covariates$T2D == 0], na.rm = TRUE)
summary(covariates$AGE[covariates$PD == 1 & covariates$T2D == 0], na.rm = TRUE)

# AGE for (T2D only)
mean(covariates$AGE[covariates$PD == 0 & covariates$T2D == 1], na.rm = TRUE) 
sd(covariates$AGE[covariates$PD == 0 & covariates$T2D == 1], na.rm = TRUE)
summary(covariates$AGE[covariates$PD == 0 & covariates$T2D == 1], na.rm = TRUE)

# AGE for (T2D + PD)
mean(covariates$AGE[covariates$PD == 1 & covariates$T2D == 1], na.rm = TRUE) 
sd(covariates$AGE[covariates$PD == 1 & covariates$T2D == 1], na.rm = TRUE)
summary(covariates$AGE[covariates$PD == 1 & covariates$T2D == 1], na.rm = TRUE)

# AGE for (Healthy Controls)
mean(covariates$AGE[covariates$PD == 0 & covariates$T2D == 0], na.rm = TRUE) 
sd(covariates$AGE[covariates$PD == 0 & covariates$T2D == 0], na.rm = TRUE)
summary(covariates$AGE[covariates$PD == 0 & covariates$T2D == 0], na.rm = TRUE)


# Visualizations for stratifications
#########################################################################
# Smoking
# Step 1: Create the grouping variable for Smoking, T2D, PD, and Healthy controls
smoking$Group_Status <- with(smoking, ifelse(T2D == 1 & PD == 1, "PD + T2D",
                                             ifelse(T2D == 1 & PD == 0, "T2D Only",
                                                    ifelse(T2D == 0 & PD == 1, "PD Only", "Healthy Controls"))))

# Step 2: Combine Smoking Status with the Group Status
smoking$Smoking_Group_Status <- with(smoking, paste(Smoking_Ever, Group_Status, sep = "_"))

# Step 3: Plot the violin plot with the new variable
p <- ggplot(smoking, aes(x = reorder(as.factor(Smoking_Group_Status), zSCORE), y = zSCORE, 
                         fill = as.factor(Smoking_Group_Status))) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("T2D PRS Distribution by Smoking and Group Status in MPBC Cohort") +
  scale_fill_manual(values = c("blue", "red", "green", "purple", "blue", "red", "green", "purple")) + 
  theme_bw() + 
  ylab("T2D PRS (zSCORE)") + 
  xlab("Smoking and Group Status") + 
  theme(legend.position = "none")

# Display the plot
print(p)

# Save the plot
ggsave("Scores_PD-DM/Smoking_Group_Status_PRS_Distribution.jpeg", plot = p, dpi = 600, units = "in", height = 6, width = 6)

# Snus
# Step 1: Create the grouping variable for Smoking, T2D, PD, and Healthy controls
snus$Group_Status <- with(snus, ifelse(T2D == 1 & PD == 1, "PD + T2D",
                                             ifelse(T2D == 1 & PD == 0, "T2D Only",
                                                    ifelse(T2D == 0 & PD == 1, "PD Only", "Healthy Controls"))))

# Step 2: Combine Smoking Status with the Group Status
snus$Snus_Group_Status <- with(snus, paste(Snus_Ever, Group_Status, sep = "_"))

# Step 3: Plot the violin plot with the new variable
p <- ggplot(snus, aes(x = reorder(as.factor(Snus_Group_Status), zSCORE), y = zSCORE, 
                         fill = as.factor(Snus_Group_Status))) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("T2D PRS Distribution by Snus Usage and Group Status in MPBC Cohort") +
  scale_fill_manual(values = c("blue", "red", "green", "purple", "blue", "red", "green", "purple")) + 
  theme_bw() + 
  ylab("T2D PRS (zSCORE)") + 
  xlab("Snus and Group Status") + 
  theme(legend.position = "none")

# Display the plot
print(p)

# Save the plot
ggsave("Scores_PD-DM/Snus_Group_Status_PRS_Distribution.jpeg", plot = p, dpi = 600, units = "in", height = 6, width = 6)


# Caffeine Consu,ption Level
# Step 1: Create the grouping variable for Smoking, T2D, PD, and Healthy controls
caffeine2$Group_Status <- with(caffeine2, ifelse(T2D == 1 & PD == 1, "PD + T2D",
                                       ifelse(T2D == 1 & PD == 0, "T2D Only",
                                              ifelse(T2D == 0 & PD == 1, "PD Only", "Healthy Controls"))))

# Step 2: Combine Smoking Status with the Group Status
caffeine2$caffeine2_Group_Status <- with(caffeine2, paste(Caffeine_Consumption_Level, Group_Status, sep = "_"))

# Step 3: Plot the violin plot with the new variable
p <- ggplot(caffeine2, aes(x = reorder(as.factor(caffeine2_Group_Status), zSCORE), y = zSCORE, 
                      fill = as.factor(caffeine2_Group_Status))) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("T2D PRS Distribution by Caffeine Consumption Level and Group Status in MPBC Cohort") +
  scale_fill_manual(values = c("blue", "red", "green", "purple", "blue", "red", "green", "purple")) + 
  theme_bw() + 
  ylab("T2D PRS (zSCORE)") + 
  xlab("Caffeine Consumption Level and Group Status") + 
  theme(legend.position = "none")

# Display the plot
print(p)

# Save the plot
ggsave("Scores_PD-DM/caffeine2_Group_Status_PRS_Distribution.jpeg", plot = p, dpi = 600, units = "in", height = 6, width = 6)

# Pesticide Exposure
# Step 1: Create the grouping variable for Smoking, T2D, PD, and Healthy controls
pesticides$Group_Status <- with(pesticides, ifelse(T2D == 1 & PD == 1, "PD + T2D",
                                                 ifelse(T2D == 1 & PD == 0, "T2D Only",
                                                        ifelse(T2D == 0 & PD == 1, "PD Only", "Healthy Controls"))))

# Step 2: Combine Smoking Status with the Group Status
pesticides$pesticides_Group_Status <- with(pesticides, paste(Pesticides_Ever, Group_Status, sep = "_"))

# Step 3: Plot the violin plot with the new variable
p <- ggplot(pesticides, aes(x = reorder(as.factor(pesticides_Group_Status), zSCORE), y = zSCORE, 
                           fill = as.factor(pesticides_Group_Status))) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.4, fill = "white", alpha = 0.7) + 
  theme_minimal() + 
  ggtitle("T2D PRS Distribution by Pesticides Exposure and Group Status in MPBC Cohort") +
  scale_fill_manual(values = c("blue", "red", "green", "purple", "blue", "red", "green", "purple")) + 
  theme_bw() + 
  ylab("T2D PRS (zSCORE)") + 
  xlab("Pesticides Exposure and Group Status") + 
  theme(legend.position = "none")

# Display the plot
print(p)

# Save the plot
ggsave("Scores_PD-DM/pesticides_Group_Status_PRS_Distribution.jpeg", plot = p, dpi = 600, units = "in", height = 6, width = 6)


table()
