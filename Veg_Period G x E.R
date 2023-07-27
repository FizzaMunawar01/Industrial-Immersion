
rm(list = ls())
gc()

library(tidyverse)
library(lme4)
library(lmerTest)
library(pracma)
library(ggplot2)
library(multcomp)
library(lattice)
library(car)    # Levene test
library(agricolae)   # Newman-Keuls & Tukey test
library(ggplot2)
library(gridExtra)
library(dplyr)
library(MASS)#BoxCox


#Set working directory
maindir <- here::here()
setwd(maindir)

##############################################################################################################

#Get the data
veg_period <- read.csv("vegetation_period.csv", 
                       header = TRUE, sep = ",")
str(veg_period)
summary(veg_period)
# unite columns 3 and 4 and store the result in a new column
#veg_period <- unite(veg_period, 
#col = "Location_Year", 
#Location, Year, 
#sep = "_") 

#mutate
veg_period$Genotype <- as.factor(veg_period$Genotype)
veg_period$ID <- toupper(veg_period$ID)
veg_period$Veg_period <- as.numeric(veg_period$Veg_period)
head(veg_period)


##############################################################################################################

# plots 

x11()
(graf2a <- ggplot(veg_period, aes(x =Location , y = Veg_period, 
                                  fill = Location)) +
    geom_boxplot())


#x11()
#(graf1b <- ggplot(pod_shattering, aes(x = Genotype, y = pod_shattering, fill = Genotype)) +
#geom_boxplot())


## 2. A variant: Violin plot
x11()
(graf2 <- ggplot(veg_period, aes(x = Location, y = Veg_period, fill = Location)) +
    geom_violin() 
)
# facet_wrap( ~ Genotype) useless bcz too many genotypes


## Calculation  of means per site and per genotype, and measure numbers
Summaries <- veg_period%>%
  group_by(Location, Genotype) %>%
  summarise(avgveg_period = mean(Veg_period, na.rm = TRUE),
            sdveg_period = sd(Veg_period, na.rm = TRUE),
            nbData = n()
  ) %>%
  arrange(desc(avgveg_period)) %>%
  print(n = Inf)  # to see all data

# NAs are produced as some of the genotypes have only 1 observation across allyears and thus 1 datapoint can't calculate the sd.

# In summary, the code calculates the mean, standard deviation, and count of data points 
#for the "pod_shattering" trait within each combination of "Location" and "Genotype." 
#It then arranges the summarized data in descending order based on the mean "pod_shattering" values and displays the results. 
#This provides a clear overview of how the "pod_shattering" trait varies across different "Location" and "Genotype" combinations, highlighting the groups with higher mean values.


## Barplot classical but not very informative compared to boxplot
x11()
(graf2 <- ggplot( Summaries, aes( x = Location, y = avgveg_period, fill = Location)) +
    geom_bar(stat = "identity"))

# too many genotypes facet wrap with genotype isn't producing any meaningful result


## classical for the study of the interaction between 2 factors
x11()
interaction.plot(veg_period$Location,veg_period$Genotype,
                 veg_period$Veg_period,
                 type="b",
                 xlab="Lieux",ylab="veg_period",trace.label="Genotypes",
                 pch=c(1:10),cex=2, lty=1,lwd=2)

# too many genotypes no meaningful result


# Apply the log transformation on 'Veg_period'
veg_period$log_Veg_period <- log(veg_period$Veg_period)



#######################################################################################################################

# ANOVA 2 factors with interaction

# Perform two-way ANOVA
veg_period_anova1 <- aov(log_Veg_period ~ Genotype * Location, 
                         data = veg_period)

# Print ANOVA summary
print(summary(veg_period_anova1))


# Perform three-way ANOVA 
veg_period_anova2 <- aov(log_Veg_period ~ Genotype * Location +
                           Year, data = veg_period)

# Print ANOVA summary
print(summary(veg_period_anova2))


# Check ANOVA pre-requisites
x11()
par(mfrow = c(2, 2))
# Diagnostic plots for the model
plot(veg_period_anova2)


# Normality of the residuals
shapiro.test(veg_period_anova2$residuals)  

# shapiro test is okay we can go ahead

# Homogeneity of variances
## attention, in this plan with nested blocks, it is necessary to "forget" a factor, otherwise no test possible
summary(aov(abs(veg_period_anova2$residuals) ~ veg_period$Gen:veg_period$Location:veg_period$Year))

## so if one consider Rep not a pb
summary(aov(abs(veg_period_anova2$residuals) ~ veg_period$Gen:veg_period$Location))   

#homoegenity of variance isn't okay we need to do box cox transformation


############ Data transformation
bc <- boxcox(veg_period_anova2)
x11()
boxcox(veg_period_anova2)
trans <-bc$x[which.max(bc$y)]
trans

###### anova with transformed data

modele2trans <- aov( log_Veg_period^trans ~ Genotype * Location + Year, 
                     data = veg_period)
summary(modele2trans)

# The Shapiro-Wilk test for normality is produced by:
shapiro.test(resid(modele2trans))

## Variances homogeneity for the residuals
leveneTest(log_Veg_period^trans ~ Genotype * Location, data= veg_period)

# Get the residuals from the fitted model
residuals_model <- resid(modele2trans)

# Diagnostic plot for normality of residuals (Q-Q plot)
qqnorm(residuals_model)
qqline(residuals_model)

# leven test is still not okay may be bcz unbalnced data?

# no significant interaction between G X E thus no need to do multiple mean comparison test.

