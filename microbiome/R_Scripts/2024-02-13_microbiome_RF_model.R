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

asv.meta.df <- readRDS("2024-03-05_asv_meta.rds")
asv.meta.df.before <- readRDS("2024-02-13_asv_meta_before.rds")
asv.meta.df.during <- readRDS("2024-02-13_asv_meta_during.rds")


# ---- take the difference between before and during samples ----

asv.meta.df <- asv.meta.df[,c(which(colnames(asv.meta.df) %in% c("Acetate", "Butyrate", "Propionate", "host_subject_id", "status")), 51:ncol(asv.meta.df))]
asv.meta.df <- na.omit(asv.meta.df)
asv.meta.df <- asv.meta.df[,-ncol(asv.meta.df)]
asv.meta.df$host_subject_id <- factor(asv.meta.df$host_subject_id)
asv.meta.df$status <- factor(asv.meta.df$status)
asv.meta.df <- asv.meta.df %>% group_by(host_subject_id, status) %>% summarize_all(mean)

asv.meta.df[,c(6:ncol(asv.meta.df))] <- asv.meta.df[,c(6:ncol(asv.meta.df))]/rowSums(asv.meta.df[,c(6:ncol(asv.meta.df))])

asv.meta.df.diff <- asv.meta.df %>% group_by(host_subject_id) %>% filter(n() > 1)
index <- which(duplicated(asv.meta.df.diff$host_subject_id) == FALSE)
my.df <- asv.meta.df.diff[1,]
my.df <- my.df[,-2]

for(i in index){
  
  temp <- asv.meta.df.diff[c(i:(i+1)),]
  subject <- as.character(temp$host_subject_id[1])
  temp <- temp[,-c(1:2)]
  temp <- temp[2,] - temp[1,]
  temp$host_subject_id <- subject
  temp <- temp[,c(ncol(temp), 1:ncol(temp)-1)]
  my.df <- rbind(my.df, temp)
  
}
my.df <- my.df[-1,]

saveRDS(my.df, file = "2024-03-05_asv_meta_df_diffs.rds")

microbiome.df <- readRDS("2024-03-05_asv_meta_df_diffs.rds")

microbiome.df <- as.data.frame(microbiome.df)

index <- sample(x = c(1:nrow(microbiome.df)), size = round(0.7*nrow(microbiome.df)), replace = F)
train <- microbiome.df[index,]
test <- microbiome.df[-index,]

# response <- "Butyrate"
response <- "Acetate"
# response <- "Propionate"

asv.cols <- c(5:ncol(microbiome.df))
meta.cols <- c(1:4)

# ---- split into training and testing data (SKIP  IF TAKING DIFFERENCE BETWEEN BEFORE/DURING SAMPLES) ----

# microbiome.df <- asv.meta.df.before 
microbiome.df <- asv.meta.df.during
microbiome.df <- microbiome.df[,c(ncol(microbiome.df), 1:(ncol(microbiome.df)-1))]

meta.cols <- c(1:50)
asv.cols <- c(51:ncol(microbiome.df))

scfa.na.remove <- which(is.na(microbiome.df$Acetate) | is.na(microbiome.df$Butyrate) | is.na(microbiome.df$Propionate))
microbiome.df <- microbiome.df[-scfa.na.remove,]

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
microbiome.df[is.nan(microbiome.df)] <- 0
microbiome.df[is.na(microbiome.df)] <- 0


low.abundance.index <- which(colSums(microbiome.df[,asv.cols]) <= 10) # remove taxa with low abundance
microbiome.df <- microbiome.df[,-low.abundance.index]
asv.cols <- c(51:ncol(microbiome.df))

microbiome.df[,asv.cols] <- microbiome.df[,asv.cols]/rowSums(microbiome.df[,asv.cols]) # normalize

index <- sample(x = c(1:nrow(microbiome.df)), size = round(0.7*nrow(microbiome.df)), replace = F)
train <- microbiome.df[index,]
test <- microbiome.df[-index,]

response <- "Butyrate"
# response <- "Acetate"
# response <- "Propionate"


# ---- (OPTIONAL) Boruta predictor selection algorithm ----

## this removes predictor variables (ASVs) that predict output worse than a randomized version of themselves

# boruta.df <- Boruta(x = train[,asv.cols], y = train$Butyrate) # running for butyrate
boruta.df <- Boruta(x = train[,asv.cols], y = train$Acetate) # running for acetate
# boruta.df <- Boruta(x = train[,asv.cols], y = train$Propionate) # running for propionate

boruta.results.df <- as.data.frame(boruta.df[["finalDecision"]])
boruta.index <- rownames(boruta.results.df)[which(boruta.results.df$`boruta.df[["finalDecision"]]` != "Rejected")]

microbiome.df <- microbiome.df[,c(meta.cols,which(colnames(microbiome.df) %in% boruta.index))]
train <- microbiome.df[index,]
test <- microbiome.df[-index,]
asv.cols <- c(51:ncol(microbiome.df))
# asv.cols <- c(5:ncol(microbiome.df))


# ---- set up random forest model ----

## set predictors as only the ASVs that Boruta has determined to be decent predictors 
predictors <- boruta.index[order(colSums(train[asv.cols])[colnames(train[asv.cols]) %in% boruta.index], decreasing = T)]
# predictors <- colnames(microbiome.df[,asv.cols])
# predictors <- c(colnames(microbiome.df[,asv.cols]), "Acetate", "Propionate")


m1 <- ranger(as.formula(paste(response, '.', sep = '~')),
             data = train[,c(response, predictors)])

## two different ways of calculating RMSE
sqrt(mean((m1$predictions - train$Butyrate)^2))
sqrt(mean((m1$predictions - train$Acetate)^2))
sqrt(mean((m1$predictions - train$Propionate)^2))
sqrt(m1$prediction.error)

# plot training data
plot(m1$predictions ~ train[,response])
abline(a = 0, b = 1, col = "red")

## assess model performance by testing model on witheld data (test data)
scfa.predict <- predict(m1, test) 

# sqrt(mean((test$Butyrate - scfa.predict$predictions)^2))
sqrt(mean((test$Acetate - scfa.predict$predictions)^2))
# sqrt(mean((test$Propionate - scfa.predict$predictions)^2))


plot(scfa.predict$predictions ~ test[,response],
     ylab = 'Observed',
     xlab = 'Predicted')
abline(a = 0, b = 1, col = "red")


m1.lm <- lm(scfa.predict$predictions ~ test[,response])
# abline(0, 1, lty = 2)
abline(m1.lm)
summary(m1.lm)



saveRDS(m1, "2024-02-27_microbiome_RF_m1_diff.rds")

# ---- parameter optimization ----

## define the parameter space
## these are all of the little settings in the model that we will try and adjust to find the best fit for the data

hyper.grid <- expand.grid(
  n.edges = seq(100, 3000, 100), # aka n.trees
  mtry       = seq(10, 30, by = 2), # 
  node_size  = seq(3, 9, by = 2), 
  sample_size = c(.55, .632, .70, .80), # internal
  OOB_RMSE   = 0 # internal
)

for(i in 1:nrow(hyper.grid)){ ## AKA for every combination of parameter settings
  
  # predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)] # cannot use more n.edges than boruta predictors
  
  try({ ## try clause necessary because some parameter combinations are incompatible
    
    model <- ranger(
      formula = as.formula(paste(response, '.', sep = '~')),
      data = train[,c(response, predictors)], 
      num.trees       = 500,
      mtry            = hyper.grid$mtry[i],
      min.node.size   = hyper.grid$node_size[i],
      sample.fraction = hyper.grid$sample_size[i],
      seed            = 123
    )
    
    ## add OOB error to grid
    hyper.grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
    
    ## From the internet: 
    ## OOB (out-of-bag) score is a performance metric for a machine learning model, 
    ## specifically for ensemble models such as random forests. 
    ## It is calculated using the samples that are not used in the training of the model, 
    ## which is called out-of-bag samples.
    ## The OOB_score is computed as the number of correctly predicted rows from the out-of-bag sample. 
    ## OOB Error is the number of wrongly classifying the OOB Sample.
    
  }, silent = F)
  
  print(paste(i, 'out of', nrow(hyper.grid), hyper.grid$OOB_RMSE[i]))
  
}

hyper.grid$OOB_RMSE[hyper.grid$OOB_RMSE == 0] <- NA
hyper.grid <- na.omit(hyper.grid)

hist(hyper.grid$OOB_RMSE, breaks = 100)

## define selected optimal parameters for the model
selected.params <- hyper.grid[which.min(hyper.grid$OOB_RMSE),]

# ---- create final model ----

# predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)]

## create second model using optimal selected paramters
m2 <- ranger(
  formula = as.formula(paste(response, '.', sep = '~')),
  data = train[,c(response, predictors)],
  num.trees       = 500,
  mtry            = selected.params$mtry,
  min.node.size   = selected.params$node_size,
  sample.fraction = selected.params$sample_size,
  seed            = 123,
  importance = 'permutation',
  oob.error = T
)

saveRDS(m2, "2024-02-27_microbiome_RF_m2_diff.rds")

## compare m1 and m2 with and without parameter optimization
scfa.predict.1 <- predict(m1, test)
scfa.predict.2 <- predict(m2, test)

plot(scfa.predict$predictions ~ test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(scfa.predict.1$predictions ~ test[,response]), col = "blue")
abline(lm(scfa.predict.2$predictions ~ test[,response]), col = "red")

summary(lm(scfa.predict.1$predictions ~ test[,response]))
summary(lm(scfa.predict.2$predictions ~ test[,response]))

sqrt(m1$prediction.error)
sqrt(m2$prediction.error)


# ---- final model ----

m3 <- ranger(
  formula = as.formula(paste(response, '.', sep = '~')),
  data = microbiome.df[,c(response, predictors)],
  num.trees       = 500,
  mtry            = selected.params$mtry,
  min.node.size   = selected.params$node_size,
  sample.fraction = selected.params$sample_size,
  seed            = 123,
  importance = 'permutation',
  oob.error = T
)


saveRDS(m3, "2024-02-27_microbiome_RF_m3_final_diff.rds")

## compare m1 and m2 with and without parameter optimization
scfa.predict.1 <- predict(m1, test)
scfa.predict.2 <- predict(m2, test)
scfa.predict.3 <- predict(m3, microbiome.df)

plot(scfa.predict.1$predictions ~ test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(scfa.predict.1$predictions ~ test[,response]), col = "blue")
summary(lm(scfa.predict.1$predictions ~ test[,response]))


plot(scfa.predict.2$predictions ~ test[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(scfa.predict.2$predictions ~ test[,response]), col = "red")
summary(lm(scfa.predict.2$predictions ~ test[,response]))


plot(scfa.predict.3$predictions ~ microbiome.df[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(scfa.predict.3$predictions ~ microbiome.df[,response]), col = "green")
summary(lm(scfa.predict.3$predictions ~ microbiome.df[,response]))

index <- which(is.infinite((scfa.predict.3$predictions - microbiome.df[,response])^2/microbiome.df[,response])) 
mean((scfa.predict.3$predictions[-index] - microbiome.df[-index,response])^2/microbiome.df[-index,response])









