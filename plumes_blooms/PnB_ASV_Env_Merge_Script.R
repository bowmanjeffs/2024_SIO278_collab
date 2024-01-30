#PNB: downloading and merging PnB 18s and environmental data from 2011-2014

#load packages ----
library(tidyverse)
library(dplyr)

#choose working directory, choose in files and using settings in bottom right box 
#pivoting ASV data 
#observations by row(samples) and asv (columns) and environmental data 

#loading in the data----

P2B_ra <- read.csv('PnB_18sV9_protist_sequence_counts.csv')
env_data <- read.csv('Catlett_PnB_composite_discrete_sample_2020.csv')
tax_data <- read.csv('PnB_18sV9_protist_taxonomy_trophicMode.csv')


# making sample ID that matches env data to P2B_ra----
  
env_data$date <- sub("^", "PB", env_data$date)
env_data$station <- sub("^", "St", env_data$station)
env_data$depth <- paste(env_data$depth, "m", sep="")
env_data$cruise<-paste(env_data$date,env_data$station, env_data$depth, sep="") 
colnames(env_data)[1] <- c("sample_id")
 
#pivoting the P2B data table so sample ID is matched with ASV relative abundance and ASV is the column title---- 
P2B_ra_pivot <- P2B_ra%>%
pivot_wider(names_from='ASV', values_from = 'rel_seq_abundance')

#combine reformatted data into final table----
seq_env_comb <-right_join(env_data, P2B_ra_pivot, by=c('sample_id'))

  seq_env_comb$depth <- gsub("m","",as.character(seq_env_comb$depth))
  seq_env_comb$station <- gsub("St", "",as.character(seq_env_comb$station))
  seq_env_comb$date<- gsub("PB", "", as.character(seq_env_comb$date))

#convert final table into csv----

write.csv(seq_env_comb, file = "PnB_ASV_RA_Env_Final.csv")


  
  
