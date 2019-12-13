library(dplyr)
library(tidyverse)
library(magrittr)

GI_df <- read.table('GI_short.txt', header = T)
GI_df[,1] <- as.character(GI_df[,1])
GI_df[,2] <- as.character(GI_df[,2])
PPI_df <- read.table('PPI_short.txt', header = T)
PPI_df[,1] <- as.character(PPI_df[,1])
PPI_df[,2] <- as.character(PPI_df[,2])
CE_df <- read.table('CE_short.txt', header = T)
CE_df[,1] <- as.character(CE_df[,1])
CE_df[,2] <- as.character(CE_df[,2])

# get merged values for both gene A and gene B
GI_df$merged_A <- paste(GI_df$Query_Strain, GI_df$Array_Strain, sep='_')
GI_df$merged_B <- paste(GI_df$Array_Strain, GI_df$Query_Strain, sep='_')

CE_df$merged_A <- paste(CE_df$Gene_A, CE_df$Gene_B, sep='_')
CE_df$merged_B <- paste(CE_df$Gene_B, CE_df$Gene_A, sep='_')

PPI_df$merged_A <- paste(PPI_df$Interactor_A, PPI_df$Interactor_B, sep='_')
PPI_df$merged_B <- paste(PPI_df$Interactor_B, PPI_df$Interactor_A, sep='_')

# get non-duplicated values for GI, CE and PPI

non_dup_GI <- select(GI_df, merged_A, merged_B) %>% gather() %>% select(value) %>%
  arrange(value) %>% 
  filter(!duplicated(.))

non_dup_CE <- select(CE_df, merged_A, merged_B) %>% gather() %>% select(value) %>%
  arrange(value) %>% 
  filter(!duplicated(.))

non_dup_PPI <- select(PPI_df, merged_A, merged_B) %>% gather() %>% select(value) %>%
  arrange(value) %>% 
  filter(!duplicated(.))

# get non-duplicated values for merged GI, CE and PPI
non_dup_all <- rbind(non_dup_CE, non_dup_GI, non_dup_PPI) %>%
  filter(!duplicated(.)) %>% set_colnames("merged")
sum(duplicated(non_dup_all))

# create strange table

# NB: filter(PPI_df, merged_A == 'YCR011C_YFL031W')
# NB: filter(PPI_df, merged_B == 'YCR011C_YFL031W')
# I have "duplicated" merged values in merged_A and merged_B columns, however PPI_score are identical

non_dup_all$Interactor_A <- NA
non_dup_all$Interactor_B <- NA
non_dup_all$GI_score <- NA
non_dup_all$CE_score <- NA
non_dup_all$PPI_score <- NA

non_dup_all$Interactor_A <- gsub('_.*', '', non_dup_all$merged) # use regex for getting  Interactor_A
non_dup_all$Interactor_B <- gsub('.*_', '', non_dup_all$merged)

# because I have identical score values for "duplicated" merged_A and merged_B, this code works correctly

non_dup_all$PPI_score <- PPI_df$PPI_score[match(non_dup_all$merged, PPI_df$merged_A)]
non_dup_all[is.na(non_dup_all$PPI_score) == T,]$PPI_score <- PPI_df$PPI_score[match(non_dup_all[is.na(non_dup_all$PPI_score) == T,]$merged, PPI_df$merged_B)]

#non_dup_all$CE_score <- CE_df$LLS[match(non_dup_all$merged, CE_df$merged_A) || match(non_dup_all$merged, CE_df$merged_B)]
non_dup_all$CE_score <- CE_df$LLS[match(non_dup_all$merged, CE_df$merged_A)]
non_dup_all[is.na(non_dup_all$CE_score) == T,]$CE_score <- CE_df$LLS[match(non_dup_all[is.na(non_dup_all$CE_score) == T,]$merged, CE_df$merged_B)]

non_dup_all$GI_score <- GI_df$GI_score[match(non_dup_all$merged, GI_df$merged_A)]
non_dup_all[is.na(non_dup_all$GI_score) == T,]$GI_score <- GI_df$GI_score[match(non_dup_all[is.na(non_dup_all$GI_score) == T,]$merged, GI_df$merged_B)]
