# merge abundance data with metadata
# RJH
# 2024-01-30

# ---- library ----

library(tidyverse)
library(readxl)

# ---- read in data ----

getwd()
setwd("2024_SIO278_collab-main/microbiome/")

metadata <- read_excel("metadata_SRP128128.xlsx")

asv.abund.df.full <- read.csv("microbiome_ASVs.csv")
asv.abund.df.reduced <- read.csv("microbiome_ASVs_reduced.csv")

# ---- format and merge ----

# asv.abund.df <- asv.abund.df.full
asv.abund.df <- asv.abund.df.reduced # pick which file to work with

asv.abund.df <- t(asv.abund.df)
colnames(asv.abund.df) <- asv.abund.df[1,]
asv.abund.df <- as.data.frame(asv.abund.df[-1,])

asv.abund.df$Run <- rownames(asv.abund.df)

combined.df <- merge(metadata, asv.abund.df, by = "Run")

# ---- save as RDS ----
# saveRDS(combined.df, file = "2024-01-30_combined_abund_metadata.rds")
saveRDS(combined.df, file = "2024-01-30_combined_abund_metadata_reduced.rds")


