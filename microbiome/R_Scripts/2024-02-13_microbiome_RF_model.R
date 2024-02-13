# Microbiome Random Forest Model
# RJH
# 2024-02-13

# ---- library ----

library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(ranger)
library(Boruta)

# ---- set working directory ----
getwd()
# setwd("2024_SIO278_collab-main/microbiome/")

# ---- set seed for reproducible results ----
set.seed(1234)

# ---- read in data ----

asv.meta.df.before <- readRDS("2024-02-13_asv_meta_before.rds")
asv.meta.df.during <- readRDS("2024-02-13_asv_meta_during.rds")

# ---- split into training and testing data ----

# microbiome.df <- asv.meta.df.before
microbiome.df <- asv.meta.df.during
microbiome.df <- microbiome.df[,c(ncol(microbiome.df), 1:(ncol(microbiome.df)-1))]

meta.cols <- c(1:50)
asv.cols <- c(51:ncol(microbiome.df))

index <- sample(x = c(1:nrow(microbiome.df)), size = round(0.8*nrow(microbiome.df)), replace = F)

scfa.na.remove <- which(is.na(microbiome.df$Acetate) | is.na(microbiome.df$Butyrate) | is.na(microbiome.df$Propionate))
microbiome.df <- microbiome.df[-scfa.na.remove,]

microbiome.df[,asv.cols] <- microbiome.df[,asv.cols] %>% mutate_all(~ifelse(is.nan(.), NA, .))
microbiome.df[,asv.cols] <- microbiome.df[,asv.cols] %>% replace(is.na(.), 0)

low.abundance.index <- which(colSums(microbiome.df[,asv.cols]) <= 10) # remove taxa with low abundance
microbiome.df <- microbiome.df[,-low.abundance.index]
asv.cols <- c(51:ncol(microbiome.df))

microbiome.df[,asv.cols] <- microbiome.df[,asv.cols]/colSums(microbiome.df[,asv.cols]) # normalize

train <- microbiome.df[index,]
test <- microbiome.df[-index,]

# ---- (OPTIONAL) Boruta predictor selection algorithm ----

## this removes precitor variables (ASVs) that predict output worse than a randomized version of themselves

boruta.df <- Boruta(x = train[,asv.cols], y = train$Butyrate) # running for butyrate

boruta.results.df <- as.data.frame(boruta.df[["finalDecision"]])
boruta.index <- rownames(boruta.results.df)[which(boruta.results.df$`boruta.df[["finalDecision"]]` != "Rejected")]

# ---- set up random forest model ----

response <- "Butyrate"

## set predictors as only the ASVs that Boruta has determined to be decent predictors 
# predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)]
predictors <- colnames(microbiome.df[,asv.cols])


m1 <- ranger(as.formula(paste(response, '.', sep = '~')),
             data = train[,c(response, predictors)])

## two different ways of calculating RMSE
sqrt(mean((m1$predictions - train.train$aop.corrected)^2))
sqrt(m1$prediction.error)

plot(m1$predictions ~ train.train[,response])



















