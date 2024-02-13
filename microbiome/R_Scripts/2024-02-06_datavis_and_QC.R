# QC and data visualization
# 2024-02-06
# RJH

# ---- library ----

library(tidyverse)
library(lubridate)
library(vegan)
library(plotly)
library(goeveg)
library(patchwork)

# ---- read in data ----
getwd()
# setwd("2024_SIO278_collab-main/microbiome/")

asv.meta.df <- readRDS("2024-01-30_combined_abund_metadata_reduced.rds")


# ---- format data ----

asv.meta.df$Butyrate <- as.numeric(asv.meta.df$Butyrate)
asv.meta.df$Acetate <- as.numeric(asv.meta.df$Acetate)
asv.meta.df$Propionate <- as.numeric(asv.meta.df$Propionate)

asv.meta.df$supplement <- as.factor(asv.meta.df$supplement)

my.asv.cols <- asv.meta.df[,c(50:ncol(asv.meta.df))]
asv.meta.df <- asv.meta.df[,-c(50:ncol(asv.meta.df))]
my.asv.cols <- mutate_all(my.asv.cols, function(x) as.numeric(as.character(x)))
asv.meta.df <- cbind(asv.meta.df, my.asv.cols)

asv.meta.df$subject <- substr(asv.meta.df$sample_title, 1, 4)

# ---- quick look at data ----

asv.meta.df.before <- asv.meta.df[which(asv.meta.df$status == "before"),]
asv.meta.df.during <- asv.meta.df[which(asv.meta.df$status == "during"),]

saveRDS(asv.meta.df.before, file = "2024-02-13_asv_meta_before.rds")
saveRDS(asv.meta.df.during, file = "2024-02-13_asv_meta_during.rds")


# ---- histograms ----

# butyrate
hist(asv.meta.df.before$Butyrate)
hist(asv.meta.df.during$Butyrate)

# acetate
hist(asv.meta.df.before$Acetate)
hist(asv.meta.df.during$Acetate)

# propionate
hist(asv.meta.df.before$Propionate)
hist(asv.meta.df.during$Propionate)

# test ASV
hist(asv.meta.df.before$ATCGGAATTACTGGGCGTAAAGGGTGCGCAGGCGGTTGAGTAAGACAGATGTGAAATCCCCGAGCTTAACTCGGGAATGGCATATGTGACTGCTCGACTAGAGTGTGTCAGAGGGAGGTGGAATTCCACGTGTAGCAGTGAAATGCGTAGATATGTGGAAGAACACCGATGGCGAAGGCAGCCTCCTGGGACATAACTGACGCTCAGGCACGA)


# ---- ordinations ----

my.scfa.names <- c("Acetate", "Butyrate", "Propionate")
scfa.df.before <- asv.meta.df.before[,which(colnames(asv.meta.df.before) %in% c(my.scfa.names, "sample_title", "subject"))]
scfa.df.during <- asv.meta.df.during[,which(colnames(asv.meta.df.during) %in% c(my.scfa.names, "sample_title", "subject"))]


scfa.mat <- as.matrix(scfa.df.during[,c(1:3)])
scfa.mat <- na.omit(scfa.mat)

# dimcheckMDS(scfa.mat)

NMS <- metaMDS(scfa.mat, distance = "bray", k = 3)
# goodness(NMS)
stressplot(NMS)

scfa.nmds <- as.data.frame(NMS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scfa.nmds <- cbind(scfa.nmds, na.omit(scfa.df.during))

ggplot(data = scfa.nmds) +
  geom_point(aes(x = MDS1, y = MDS2)) +
  theme_bw()


# ASVs

asv.df.during <- asv.meta.df.during[,c(50:ncol(asv.meta.df.during)-1)] # only numeric ASV data

# asv.df.during <- mutate_all(asv.df.during, function(x) as.numeric(as.character(x)))
asv.df.during <- na.omit(asv.df.during)

asv.df.during <- asv.df.during[which(rowSums(asv.df.during) > 4000),]
asv.df.during <- asv.df.during[,which(colSums(asv.df.during) > 10)] # remove ASVs with less than 10 reads across all samples

asv.mat <- as.matrix(asv.df.during)

# dimcheckMDS(asv.mat)

NMS <- metaMDS(asv.mat, distance = "bray", k = 3)
# goodness(NMS)
stressplot(NMS)

asv.nmds <- as.data.frame(NMS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
# asv.nmds <- cbind(asv.nmds, na.omit(asv.df.during))


ggplot(data = asv.nmds) +
  geom_point(aes(x = MDS1, y = MDS2)) +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_bw()










