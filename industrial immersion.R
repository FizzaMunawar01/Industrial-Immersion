
### II. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# use of the constraint 'set-to-zero' for ANOVAs ## will see later in this script
options(contrasts=c('contr.treatment','contr.poly'))
#can also use 'contr.sum' for a 'sum-to-zero' constraint


#############################################################################################################

## Loading of the R packages needed for the analysis.
library(ggplot2) # Needed for some graphs (e.g. bwplots)
library(tidyverse)
library(agricolae) # For multiple mean comparisons
library(car) # for Levene's test
library(lme4)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(lmerTest)
library(caret)


#Set working directory
maindir <- here::here()
setwd(maindir)


##############################################################################################################



#Get the data
pod_shattering <- read.csv("pod_shattering.csv", header = TRUE, sep = ",")
str(pod_shattering)
summary(pod_shattering)

#mutate
pod_shattering$Location <- as.factor(pod_shattering$Location)
pod_shattering$Genotype <- as.factor(pod_shattering$Genotype)
pod_shattering$Year <- as.factor(pod_shattering$Year)
pod_shattering$ID <- toupper(pod_shattering$ID)
pod_shattering$pod_shattering <- as.numeric(pod_shattering$pod_shattering)

# unite columns 3 and 4 and store the result in a new column named 'united_column'
pod_shattering <- unite(pod_shattering, 
                        col = "Location_Year", 
                        Location, Year, 
                        sep = "_")  
head(pod_shattering)



##############################################################################################################

# Data transformation for normality
# Apply the arcsine square root transformation  
# As percentage values range from 0-100 other type of transformations won't work with zero

pod_shattering$normalized_pod_shattering <- asin(sqrt(pod_shattering$pod_shattering/100))

# Now 'normalized_pod_shattering' column contains the transformed values


##############################################################################################################


# Fit a linear mixed-effects model
model <- lmer(normalized_pod_shattering ~ Location_Year + (1 | Genotype), 
              data = pod_shattering)

# Get BLUPs for the random effect (Genotype)
BLUPs <- ranef(model)$Genotype

# Extract Genotype levels from the fitted model
genotype_levels <- levels(pod_shattering$Genotype)

# Combine BLUPs with Genotype levels in a data frame
BLUPs_with_variety <- data.frame(Genotype = genotype_levels, BLUP = BLUPs)

# View the BLUPs with variety names
print(BLUPs_with_variety)

# Write BLUPs to a CSV file named "BLUPs.csv"
#write.csv(BLUPs_with_variety, file = "BLUPs_pod.csv", row.names = FALSE)


###############################################################################################################

#check assumptions for linear model on transformed data

# Get the residuals from the fitted model
residuals_model <- resid(model)

# Diagnostic plot for normality of residuals (Q-Q plot)
qqnorm(residuals_model)
qqline(residuals_model)

# Diagnostic plot for constant variance (Residuals vs. Fitted values plot)
plot(fitted(model), residuals_model,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted Values")

# Diagnostic plot for linearity (Observed vs. Fitted values plot)
# Get the fitted values from the model
fitted_values <- fitted(model)

# Get the observed 'pod_shattering' values
observed_values <- pod_shattering$normalized_pod_shattering

# Remove any rows with missing values from both fitted values and observed values
complete_cases <- complete.cases(fitted_values, observed_values)
fitted_values <- fitted_values[complete_cases]
observed_values <- observed_values[complete_cases]

# Create the plot for Observed vs. Fitted values
plot(fitted_values, observed_values,
     xlab = "Fitted Values",
     ylab = "Observed Pod Shattering",
     main = "Observed vs. Fitted Values")



#### result #####

# QQplot looks okay with asci transformation but the other two aren't okay.


################




#############################################################################################################
############################################## Data file 2###################################################


#Get the data
veg_period <- read.csv("vegetation_period.csv", 
                       header = TRUE, sep = ",")
str(veg_period)
summary(veg_period)
# unite columns 3 and 4 and store the result in a new column
veg_period <- unite(veg_period, 
                        col = "Location_Year", 
                        Location, Year, 
                        sep = "_")  
#mutate
veg_period$Genotype <- as.factor(veg_period$Genotype)
veg_period$ID <- toupper(veg_period$ID)
veg_period$Veg_period <- as.numeric(veg_period$Veg_period)
head(veg_period)


###########################################################################################################


# Apply the log transformation on 'Veg_period'
veg_period$log_Veg_period <- log(veg_period$Veg_period)

# Now 'log_Veg_period' column contains the log-transformed values


##########################################################################################################

# Fit a linear mixed-effects model
model2 <- lmer(log_Veg_period ~ Location_Year + (1 | Genotype), data = veg_period)

# Get BLUPs for the random effect (Genotype)
BLUPs2 <- ranef(model2)$Genotype

# Extract Genotype levels from the fitted model
genotype_levels <- levels(veg_period$Genotype)

# Combine BLUPs with Genotype levels in a data frame
BLUPs_with_variety2 <- data.frame(Genotype = genotype_levels, BLUP = BLUPs2)

# View the BLUPs with variety names
print(BLUPs_with_variety2)

# Write BLUPs to a CSV file named "BLUPs.csv"
#write.csv(BLUPs_with_variety2, file = "BLUPs_veg_period.csv", row.names = FALSE)


#########################################################################################################


#check assumptions for linear model on transformed data

# Get the residuals from the fitted model
residuals_model <- resid(model2)

# Diagnostic plot for normality of residuals (Q-Q plot)
qqnorm(residuals_model)
qqline(residuals_model)

# Diagnostic plot for constant variance (Residuals vs. Fitted values plot)
plot(fitted(model2), residuals_model,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted Values")

# Diagnostic plot for linearity (Observed vs. Fitted values plot)
plot(fitted(model2), veg_period,
     xlab = "Fitted Values",
     ylab = "Observed Veg period",
     main = "Observed vs. Fitted Values")

# Get the fitted values from the model
fitted_values <- fitted(model2)

# Get the observed 'veg_period' values
observed_values <- veg_period$log_Veg_period

# Remove any rows with missing values from both fitted values and observed values
complete_cases <- na.exclude(cbind(fitted_values, observed_values))
fitted_values <- complete_cases[, 1]
observed_values <- complete_cases[, 2]

# Create the plot for Observed vs. Fitted values
plot(fitted_values, observed_values,
     xlab = "Fitted Values",
     ylab = "Observed Veg period",
     main = "Observed vs. Fitted Values")


#### result #####

# data is normalized for veg period using log transformation.

#######################################################################################################


#descriptive statistics

# Data frame 'pod_shattering' has a response variable 'pod_shattering'
# and predictor variable 'Location_Year'

# Calculate summary statistics for the 'pod_shattering' variable
summary(pod_shattering$pod_shattering)

# Calculate mean and median for 'pod_shattering'
mean_value <- mean(pod_shattering$pod_shattering)
median_value <- median(pod_shattering$pod_shattering)
cat("Mean:", mean_value, "\n")
cat("Median:", median_value, "\n")

# Calculate standard deviation for 'pod_shattering'
sd_value <- sd(pod_shattering$pod_shattering)
cat("Standard Deviation:", sd_value, "\n")

# Calculate minimum and maximum for 'pod_shattering'
min_value <- min(pod_shattering$pod_shattering)
max_value <- max(pod_shattering$pod_shattering)
cat("Minimum:", min_value, "\n")
cat("Maximum:", max_value, "\n")

# Calculate range for 'pod_shattering'
range_value <- range(pod_shattering$pod_shattering)
cat("Range:", range_value[2] - range_value[1], "\n")




