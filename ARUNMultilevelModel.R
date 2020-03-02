#Loading required modules
library(nlme)
library(ggplot2)

# set working directory and load data
setwd("/Users/ArunMahadevan/Documents/BBL/studies/alpraz/scripts/ControlCode/Results/06-Feb-2020_avge_FD_thresh_0.5_parcelCoverageThresh_0.5_allNodes_QA")
data <- read.csv("allControlEnergies_emotionrec_contrast1_threatcorrectStd.csv")

# converting categorical variables to factors
data$genderFactor <- as.factor(data$gender)
data$drugFactor <- as.factor(data$drug)
data$groupFactor <- as.factor(data$group)

# run main model
hist(data$persistence_allNodes)
persistenceMainModel <- lme(fixed = persistence_DMN ~ drugFactor + groupFactor + drugFactor*groupFactor + STAI_TRAIT*drugFactor, # + SISTOTAL + genderFactor + age, 
                data = data,
                random =~ 1 | subjectID,
                na.action = na.omit)

summary(persistenceMainModel)

persistenceMainModel2 <- lme(fixed = persistence_allNodes ~ drugFactor + groupFactor + drugFactor*groupFactor + STAI_TRAIT*drugFactor + SISTOTAL + genderFactor + age, 
                             data = data,
                             random =~ 1 + drugFactor | subjectID,
                             na.action = na.omit)

summary(persistenceMainModel2)

anova(persistenceMainModel, persistenceMainModel2) # check whether random slopes model is significantly better than random intercepts model

# accuracy versus persistence
hist(data$pctcorr_threat, breaks=25)
accuracyModel <- lme(fixed = pctcorr_threat ~ persistence_allNodes + drugFactor + groupFactor + genderFactor + age,
                        data = data,
                        random =~ 1 | subjectID,
                        na.action = na.omit)
summary(accuracyModel)

# median reaction time versus persistence
hist(data$rtmdn_threatcorr, breaks=25)
reactionTimeModel <- lme(fixed = rtmdn_threatcorr ~ persistence_allNodes + drugFactor + groupFactor + persistence_allNodes*drugFactor*groupFactor + genderFactor + age,
                         data = data,
                         random =~ 1 | subjectID,
                         na.action = na.omit)
summary(reactionTimeModel)
