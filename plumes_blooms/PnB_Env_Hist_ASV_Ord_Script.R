# Review and Quality Control of Environmental and ASV Data 

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

#Quality Control of data----

#Are there any extreme outliers found in the environmental parameters? 
#Use histograms to visualize the data:


#Quality Control of data----

#read in data
seq_env_comb<-read.csv("PnB_ASV_RA_Env_Final.csv")

#forloop, to examine histograms of environmental data 
for(i in 5:55){
  env_data<-seq_env_comb[i]
  hist(env_data[,], main=paste("Histogram of", colnames(env_data)), xlab="Value")
  box(lty = "solid")
}

#Remove GYRO (only zero values)

seq_env_comb<-subset(seq_env_comb, select= -c(Gyro))

#Downloading the final organized dataframe----

write.csv(seq_env_comb, file = "PnB_ASV_RA_Env_Final_QC.csv")


#Other way to visualize specific Histograms with axis labels

#Salinity 
hist(seq_env_comb$btl_sal, xlab="Salinity", main="Histogram of Salinity")


#Temperature
hist(seq_env_comb$btl_temp, xlab="Temperature (C)", main="Histogram of Temp (C)")


#Phosphate
hist(seq_env_comb$PO4, xlab="[PO4]", main="Histogram of PO4")


#Silica
hist(seq_env_comb$SiO2, xlab="[SiO2]", main="Histogram of Silicon Dioxide")


#Nitrite 
hist(seq_env_comb$NO2, xlab="Nitrite NO2", main="Histogram of Nitrite")


#NO3NO2
hist(seq_env_comb$NO3NO2, xlab="NO3NO2", main="Histogram of NO3NO2")

#POC
hist(seq_env_comb$POC, xlab="POC", main="Histogram of POC")


#PON

hist(seq_env_comb$PON, xlab="PON", main="Histogram of PON")



#Ordination----

# ASVs

pnb_ASV <- seq_env_comb[,c(60:ncol(seq_env_comb)-1)] # keep only the numeric ASV data

#Remove NA values 
pnb_ASV <- na.omit(pnb_ASV)

#pnb_ASV <- pnb_ASV[which(rowSums(pnb_ASV) > 4000),]
#pnb_ASV <- pnb_ASV[,which(colSums(pnb_ASV) > 10)] # remove ASVs with less than 10 reads across all samples

asv.mat <- as.matrix(pnb_ASV) #convert to matrix form


NMS <- metaMDS(asv.mat, distance = "bray", k = 3)
#goodness(NMS)
stressplot(NMS)# Plot stressplot

asv.nmds <- as.data.frame(NMS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
# asv.nmds <- cbind(asv.nmds, na.omit(pnb_ASV))

#Plot ordination 
ggplot(data = asv.nmds) +
  geom_point(aes(x = MDS1, y = MDS2)) +
  theme_bw()
