library(dplyr)

non_dup_all <- read.csv('non_duplicated_values.csv', header = T, stringsAsFactors = F)[,-1]
test <- non_dup_all %>% select(Interactor_A, Interactor_B, GI_score, PPI_score, CE_score) %>% arrange(Interactor_A) %>% 
  mutate(PPI_score = abs(PPI_score), GI_score = abs(GI_score), CE_score = abs(CE_score))
test[is.na(test)] <- 0
sub_test <- test[test$PPI_score > 5 | test$GI_score > 0.2 | test$CE_score > 2, ]

overlap_inter <- sub_test %>%  select(Interactor_A, Interactor_B)
overlap_inter$over_all <- NA
str(overlap_inter)


for (i in seq_along(1:nrow(overlap_inter))){
  overlap_inter[i,3] <- sum(sort(overlap_inter[overlap_inter$Interactor_A == overlap_inter[i,1],]$Interactor_B) %in% 
                              sort(overlap_inter[overlap_inter$Interactor_A == overlap_inter[i,2],]$Interactor_B))
}


str(overlap_inter)
overlap_inter <- overlap_inter %>% select(Interactor_A, Interactor_B, over_all, over_GI, over_PPI, over_CE)
overlap_inter$over_GI <- NA
overlap_inter$over_PPI <- NA
overlap_inter$over_CE <- NA

str(sub_test)
GI_df <- sub_test %>% select(Interactor_A, Interactor_B, GI_score)
PPI_df <- sub_test %>% select(Interactor_A, Interactor_B, PPI_score)
CE_df <- sub_test %>% select(Interactor_A, Interactor_B, CE_score)
GI_df <- na.omit(GI_df)
PPI_df <- na.omit(PPI_df)
CE_df <- na.omit(CE_df)


for (i in seq_along(1:nrow(overlap_inter))){
  overlap_inter[i,4] <- sum(sort(GI_df[GI_df$Interactor_A == overlap_inter[i,1] & GI_df$GI_score > 0.2,]$Interactor_B) %in% 
                              sort(GI_df[GI_df$Interactor_A == overlap_inter[i,2] & GI_df$GI_score > 0.2,]$Interactor_B))
}

for (i in seq_along(1:nrow(overlap_inter))){
  overlap_inter[i,5] <- sum(sort(PPI_df[PPI_df$Interactor_A == overlap_inter[i,1] & PPI_df$PPI_score > 5,]$Interactor_B) %in% 
                              sort(PPI_df[PPI_df$Interactor_A == overlap_inter[i,2] & PPI_df$PPI_score > 5,]$Interactor_B))
}


for (i in seq_along(1:nrow(overlap_inter))){
  overlap_inter[i,6] <- sum(sort(CE_df[CE_df$Interactor_A == overlap_inter[i,1] & CE_df$CE_score > 2,]$Interactor_B) %in% 
                              sort(CE_df[CE_df$Interactor_A == overlap_inter[i,2] & CE_df$CE_score > 2,]$Interactor_B))
}