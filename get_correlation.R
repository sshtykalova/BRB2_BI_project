library(dplyr)

non_dup_all <- read.csv('non_duplicated_values.csv', header = T, stringsAsFactors = F)[,-1]
overlap_table <- read.csv('overlapping_corrected.csv', header = T, stringsAsFactors = F)[,-1]

test <- non_dup_all %>% select(Interactor_A, Interactor_B, GI_score, PPI_score, CE_score) %>% arrange(Interactor_A) %>% 
  mutate(PPI_score = abs(PPI_score), GI_score = abs(GI_score), CE_score = abs(CE_score))
test[is.na(test)] <- 0
sub_test <- test[test$PPI_score > 5 | test$GI_score > 0.2 | test$CE_score > 2, ] #filtered total table to save only relevant values
sub_test[sub_test$PPI_score == 0, ]$PPI_score <- NA
sub_test[sub_test$GI_score == 0, ]$GI_score <- NA
sub_test[sub_test$CE_score == 0, ]$CE_score <- NA


#get correlation values for each interactor

all_cor <- cbind(unique(sub_test$Interactor_A))
all_cor <- as.data.frame(all_cor)
all_cor[,1] <- as.character(all_cor[,1])
colnames(all_cor) <- 'Interactor'
all_cor$PPI_num <- cbind(overlap_table$PPI_num)
all_cor$GI_num <- cbind(overlap_table$GI_num)
all_cor$CE_num <- cbind(overlap_table$CE_num)
all_cor$GI_PPI <- NA
all_cor$GI_PPI_pval <- NA
all_cor$GI_CE <- NA
all_cor$GI_CE_pval <- NA
all_cor$PPI_CE <- NA
all_cor$PPI_CE_pval <- NA



for (i in seq_along(1:nrow(all_cor))) {
  if (nrow(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, GI_score))) > 5) {
    all_cor[i, 5] <- cor.test(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, GI_score))$PPI_score, na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, GI_score))$GI_score,  method = "spearman")$estimate[[1]]
    all_cor[i, 6] <- cor.test(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, GI_score))$PPI_score, na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, GI_score))$GI_score,  method = "spearman")$p.value[[1]]
  }
  if (nrow(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, GI_score, CE_score))) > 5) {
    all_cor[i, 7] <- cor.test(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, GI_score, CE_score))$GI_score, na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, GI_score, CE_score))$CE_score,  method = "spearman")$estimate[[1]]
    all_cor[i, 8] <- cor.test(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, GI_score, CE_score))$GI_score, na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, GI_score, CE_score))$CE_score,  method = "spearman")$p.value[[1]]
  }
  if (nrow(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, CE_score))) > 5) {
    all_cor[i, 9] <- cor.test(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, CE_score))$PPI_score, na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, CE_score))$CE_score,  method = "spearman")$estimate[[1]]
    all_cor[i, 10] <- cor.test(na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, CE_score))$PPI_score, na.omit(sub_test[sub_test$Interactor_A == all_cor[i,1],] %>% select(Interactor_A, PPI_score, CE_score))$CE_score,  method = "spearman")$p.value[[1]]
  }
}
