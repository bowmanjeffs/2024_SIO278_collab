# Download Required Libraries----
library(randomForest)
library(dplyr)
library(ggplot2)

#Read in your Data, includes the environmental parameters and ASV.
rf_data <- read.csv("PnB_ASV_RA_Env_Final_QC.csv", header=TRUE, row.name=1)


#Partition data 70:30 ratio, 70% data for training (train) and 30% data for validation (test)
set.seed(123)
ind <- sample(2, nrow(rf_data), replace = TRUE, prob = c(0.7, 0.3)) #why do they start from 2 ??
train <- rf_data[ind==1,]
test <- rf_data[ind==2,]



#Output Data Frame 3 columsn (1) Env Pred, (2) R^2 Training Data- Testing Model Validity, (3) R^2 Testing Data- Predicting Env. Parameter 

r_sq_dataframe<- matrix(NA, nrow=length(6:55), ncol =3) #empty output matrix 

colnames(r_sq_dataframe)<-c("Variable","Model Validation","Predictive Power")


for (i in 43:55){ # 6 through 55 represents the columns in which environmental data is stored in rf_data

  train_env<- subset(train, select= -c(1:(i-1), (i+1):58)) #Modify training data to exclude all dependent variables except the variable of interest 
  train_env<-na.omit(train_env)
  
  #save the formula that will be inputted into the RF regression model 
  env_var<-colnames(train_env[1])
  form_input<-paste0(env_var,"~.")
  formula<-as.formula(form_input)
  list<-list(formula)
 
 #running rf
set.seed(4444)
rf <- randomForest(list[[1]], data=train_env, #here column[i] is the dependent variable and the ASV relative abundances are the independent variables
                   ntree = 300, #number of trees will depend on the error plot
                   importance = TRUE, #evaluates importance of a predictor
                   proximity = TRUE) #calculates the proximity measure among the rows


plot(rf)

#predicting the dependent variable in the validation set using generated RF model

p1 <- predict(rf, train_env)

train_env$pred<-p1

#plotting predicted vs actual 

plot<-ggplot(data=train_env, aes(x=get(env_var), y=pred)) +
  geom_point()+
  geom_smooth(method = "lm")+
  xlab(paste("Actual", env_var, "Concentration"))+
  ylab(paste("Predicted", env_var ,"Concentration"))+
  ggtitle("RF Training Data Validation")

print(plot)
ggplot2::ggsave(filename=paste0("plot_val_", env_var,".png"), plot, height = 3, width = 5)

#Saving the r^2 values into a data frame 

model_train_val<-lm(get(env_var)~pred,data=train_env)

r_sq_dataframe[i-5, 1]<-env_var
r_sq_dataframe[i-5, 2]<-summary(model_train_val)$adj.r.squared



#predicting the dependent variable in the validation set using generated RF model
test_env<- subset(test, select= -c(1:(i-1), (i+1):58)) #Modify testing data to exclude all dependent variables except the variable of interest 
test_env<-na.omit(test_env)

p2 <- predict(rf, test_env)

test_env$pred<-p2

#plotting predicted vs actual for the dependent variable 
plot_test<-ggplot(data=test_env, aes(x=get(env_var), y=pred)) +
  geom_point()+
  geom_smooth(method = "lm")+
  xlab(paste("Actual", env_var, "Concentration"))+
  ylab(paste("Predicted", env_var ,"Concentration"))+
  ggtitle("RF Test (Non-Training) Data Validation")

print(plot_test)
ggplot2::ggsave(filename=paste0("plot_rf_test_", env_var,".png"), plot_test, height = 3, width = 5)

#Saving the r^2 values into a data frame 


model_test_val<-lm(get(env_var)~pred,data=test_env)

r_sq_dataframe[i-5, 3]<-summary(model_test_val)$adj.r.squared

}


#saving outputted R^2 values to determine best model performance 
write.csv(r_sq_dataframe, file= "PnB_Random_Forest_Results.csv")





















































#predicting the dependent variable in the validation set using generated RF model
p1 <- predict(rf, train_env)


#ggplot(x=train_env$PO4, y=p1)
train$PO4pred<-p1
ggplot(data=train, aes(x=PO4, y=PO4pred)) +
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Actual PO4 Concentration")+
  ylab("Predicted PO4 Concentration")+
  ggtitle("Random Forest (Training Data Validation")

model_train_val<-lm(PO4~PO4pred,data=train_env)
summary(model_train_val)
adj_r_sq_train<-summary(model_train_val)$adj.r.squared



#predicting the dependent variable in the validation set using generated RF model

p2 <- predict(rf, test)
test$PO4pred<-p2
ggplot(data=train, aes(x=PO4, y=PO4pred)) +
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Actual PO4 Concentration")+
  ylab("Predicted PO4 Concentration")+
  ggtitle("Random Forest Non-Training Data Validation")

model_test_val<-lm(PO4~PO4pred,data=test)
adj_r_sq_train<-summary(model_test_val)$adj.r.squared


































# Download Required Libraries----
library(randomForest)
library(dplyr)
library(ggplot2)

#Read in your Data, includes the environmental parameters and ASV.
rf_data <- read.csv("PnB_ASV_RA_Env_Final_QC.csv", header=TRUE, row.name=1)


#Partition data 70:30 ratio, 70% data for training (train) and 30% data for validation (test)
set.seed(123)
ind <- sample(2, nrow(rf_data), replace = TRUE, prob = c(0.7, 0.3)) #why do they start from 2 ??
train <- rf_data[ind==1,]
test <- rf_data[ind==2,]



#Output Data Frame 3 columsn (1) Env Pred, (2) R^2 Training Data- Testing Model Validity, (3) R^2 Testing Data- Predicting Env. Parameter 

r_sq_dataframe<- matrix(NA, nrow=length(6:55), ncol =3) #empty output matrix 

colnames(r_sq_dataframe)<-c("Variable","Model Validation","Predictive Power")


for (i in 6:55){ # 6 through 55 represents the columns in which environmental data is stored in rf_data
  
  train_env<- subset(train, select= -c(1:7, 9:58)) #Modify training data to exclude all dependent variables except the variable of interest 
  
  
  #running rf
  set.seed(4444)
  rf <- randomForest(PO4~., data=train_env, #here column[i] is the dependent variable and the ASV relative abundances are the independent variables
                     ntree = 300, #number of trees will depend on the error plot
                     importance = TRUE, #evaluates importance of a predictor
                     proximity = TRUE) #calculates the proximity measure among the rows
  
  plot(rf)
  
}




#predicting the dependent variable in the validation set using generated RF model
p1 <- predict(rf, train_env)

ggplot(x=train_env$PO4, y=p1)
train$PO4pred<-p1
ggplot(data=train, aes(x=PO4, y=PO4pred)) +
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Actual PO4 Concentration")+
  ylab("Predicted PO4 Concentration")+
  ggtitle("Random Forest Training Data Validation")

model_train_val<-lm(PO4~PO4pred,data=train_env)
summary(model_train_val)
adj_r_sq_train<-summary(model_train_val)$adj.r.squared



#predicting the dependent variable in the validation set using generated RF model

p2 <- predict(rf, test)
test$PO4pred<-p2
ggplot(data=train, aes(x=PO4, y=PO4pred)) +
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("Actual PO4 Concentration")+
  ylab("Predicted PO4 Concentration")+
  ggtitle("Random Forest Non-Training Data Validation")

model_test_val<-lm(PO4~PO4pred,data=test)
adj_r_sq_train<-summary(model_test_val)$adj.r.squared





  
  
  
  