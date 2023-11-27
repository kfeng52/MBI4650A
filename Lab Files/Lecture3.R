#Lecture 3 Live Coding Exercise

#install.packages("Lock5Data", repos="http://R-Forge.R-project.org")
#install.packages("lme4",
# repos=c("http://lme4.r-forge.r-project.org/repos",
#  getOption("repos")[["CRAN"]]))

#Clear your environment
rm(list = ls())

library(Lock5Data)
library(lattice)
library(ggplot2)
library(lme4)

#TOPIC 1: Regression Models

#Simple Regression
#Load the data
data(SleepStudy)
help(SleepStudy)

#Look at your data
SleepStudy <- as.data.frame(SleepStudy)
dim(SleepStudy)
names(SleepStudy)
str(SleepStudy)

#One thing we notice right away is that some of the variables are coded as integers when they should be factors or numeric

#Rename gender to sex
colnames(SleepStudy)[1] <- "Sex"
names(SleepStudy)

#We can rename the dataframe to make it easier to type
sleep <- SleepStudy

#Convert integer columns to factors
colsF <- c("Sex", "ClassYear", "NumEarlyClass", "EarlyClass", "AllNighter")
colsN <- c("PoorSleepQuality", "DepressionScore", "AnxietyScore", "StressScore", "Happiness","DASScore", "Drinks", "ClassesMissed")
sleep[colsF] <- lapply(sleep[colsF], factor)  ## as.factor() could also be used
sleep[colsN] <- lapply(sleep[colsN], as.numeric)  ## as.numeric() could also be used
str(sleep)

#Simple linear regression (outcome is continuous)
#The term ‘simple’ refers to the fact that we only have one explanatory variable, so the resulting linear model will only have one intercept term and one slope.
#One question we can ask: is the average amount of sleep associated with an individual's GPA?
#We regress GPA (the response variable, dependent variable, y) on AverageSleep (the explanatory variable or predictor, independent variable, x)
##lm(dependent~independent)
summary(lm(sleep$GPA  ~  sleep$AverageSleep))

#You can also set your lm call as an object
model1 <- lm(sleep$GPA  ~  sleep$AverageSleep)

#Some variables in R are of specific types, meaning R knows what they are, and knows how to handle them in certain circumstances. 
#For example, our mod1 variable is of class lm, which comes with some useful commands such as summary, abline, and plot. 
#This makes it easy to print a nice summary for the regression and also add the fitted line in red to a plot of our original data. 
summary(model1)
plot(GPA ~ AverageSleep, data=sleep, pch=19, cex=0.5)
abline(model1, col="red")
abline(lm(sleep$GPA ~ sleep$AverageSleep), col="red")

#Let's look at the p-value of the association, and the effect size or estimate, sometimes also called the beta value
#In this example the p-value is 0.3368 and the effect size is -0.02541
#The units of the effect size are always the same as the units for the dependent variable (in this case GPA)
#The easiest way to interpret this is to say that for every ONE unit increase in the independent variable the dependent variable is expected to change by -0.02541. 
#In this case, this effect is not statistically significant

#Simple logistic regression (outcome is categorical)
summary(glm(sleep$AllNighter~sleep$ClassesMissed, family=binomial))

#Some general insights:
##The ‘Residual standard error’ reported for the model fit is the estimated standard deviation of dependent variable values, at fixed values of the independent variable
##The adjusted R-squared value indicates the strength of the trend (like cor) 
##Interpreting coefficients
#**Intercept** If the value of [explanatory variable] was zero, our model would predict [response variable] to be [intercept value].
#**Slope** For a one-[unit] increase in [explanatory variable], our model would predict a [slope value]-[unit] [increase/decrease] in [response variable].

#Simple Plots
#Before proceeding with fitting any model it is always important to look at a plot of the underlying data, first to ensure a linear model framework is appropriate and second to have an understanding of what the association looks like

#Density plot
ggplot(sleep,aes(x=AverageSleep)) + geom_density(alpha=0.25) + geom_rug() + xlab("AverageSleep")

ggplot(sleep,aes(x=GPA)) + geom_density(alpha=0.25) + geom_rug() + xlab("GPA")

#Scatter plot with a regression line
ggplot(sleep, aes(y=GPA, x=AverageSleep)) + geom_smooth(method="lm",color="black") + geom_point(shape=19, size=1)

#colour by sex
ggplot(sleep, aes(y=GPA, x=AverageSleep)) + geom_smooth(method="lm",color="black") + geom_point(shape=19, size=1, aes(colour=Sex))

#TOPIC 2:Multiple Regression

#Now let’s fit a multiple linear regression so we have more than one explanatory variable in our model. To include additional individual predictors, use the ‘+’ sign.
summary(lm(GPA  ~  AverageSleep, data=sleep))
summary(lm(GPA  ~  AverageSleep + CognitionZscore, data=sleep))
#How did the estimate change?

#The most common reason to use multiple regression is to adjust for covariates or confounding variables in a study, however multiple regression also allows for multiple variables to be in your prediction model
#Let's plot the multiple regression

ggplot(sleep, aes(x=AverageSleep, y=GPA)) + geom_point(shape=19, size=1, aes(colour=CognitionZscore)) +geom_smooth(method=lm, # Add linear regression lines 
                                                                                                                   se=FALSE, # Do not add shaded confidence region    
                                                                                                                   fullrange=TRUE) # Extend regression lines
#Colinearity
#multicollinearity occurs when independent variables in a regression model are correlated. This correlation is a problem because independent variables should be independent. 
#If the degree of correlation between variables is high enough, it can cause problems when you fit the model and interpret the results.

#Colinearity can lead to highly variable and uninterpretable results, as a general rule it is your job to make sure that multiple independent variables in a regression are not highly correlated

#Calculating Residuals
#We can ask for the residuals of our model as lists
M1R <- resid(model1)
head(M1R)

#I like to append residuals to my original dataframe
sleep$resids_model1 <- resid(model1)

M2R <- resid(lm(GPA  ~  AverageSleep + CognitionZscore, data=sleep))
head(M2R)
#This can be very useful for generating variables with other effects 'removed'
#In this case you might actually run this instead
GPA_CogResids <- resid(lm(GPA  ~  CognitionZscore, data=sleep))

#This would allow you to use the variable GPA_CogResids in another model as an independent variable without having to worry about the effect of cognitionZscore on GPA (it will have been regressed out)

#TOPIC 3:Mixed Effects Regression Models
#Use package lme4 and command 'lmer'
#Let's add a random variable, in this case we will randomly assign null data that reflects which researcher performed the survey for each individual
sleep$Researcher <- sample(x=1:5, size=253, replace=TRUE)

summary(lm(StressScore  ~  PoorSleepQuality, data=sleep))

#Let's add in the researcher effect on StressScore
#In this case research effect would be a random variable (batch effect)
#Because we have mixed effects (random and fixed), we need to use the lmer command instead of lm (lmer is from package lme4)
summary(lmer(StressScore  ~  PoorSleepQuality + (1|Researcher), data=sleep))
#If you see a message that indicates 'singular fit' the model is overfitted – that is, the random effects structure is too complex to be supported by the data

#lmer does not return a pvalue, we can use the t-value and convert it to a p-value using the following formula
pchisq(5.493**2,1,lower.tail=F)

#END