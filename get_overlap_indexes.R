library(dplyr)

non_dup_all <- read.csv('non_duplicated_values.csv', header = T, stringsAsFactors = F)[,-1]

non_dup_all <- non_dup_all %>% arrange(merged)
test <- non_dup_all %>% select(Interactor_A, Interactor_B, GI_score, PPI_score, CE_score) %>% arrange(Interactor_A) %>% 
  mutate(PPI_score = abs(PPI_score), GI_score = abs(GI_score), CE_score = abs(CE_score))

test[is.na(test)] <- 0
sub_test <- test[test$PPI_score > 5 | test$GI_score > 0.2 | test$CE_score > 2, ]
sub_test[sub_test$PPI_score == 0, ]$PPI_score <- NA
sub_test[sub_test$GI_score == 0, ]$GI_score <- NA
sub_test[sub_test$CE_score == 0, ]$CE_score <- NA

#get number of neighbours
overlap_table <- cbind(unique(sub_test$Interactor_A))
overlap_table <- as.data.frame(overlap_table)
overlap_table$neighbours_num <- NA
overlap_table$PPI_num <- NA
overlap_table$GI_num <- NA
overlap_table$CE_num <- NA
overlap_table[,1] <- as.character(overlap_table[,1])
colnames(overlap_table) <- c('Interactor_A', 'neighbours_num', 'PPI_num', 'GI_num', 'CE_num')

str(sub_test)
str(overlap_table)

for (i in seq_along(1:nrow(overlap_table))){
  overlap_table[i,2] <- nrow(select(filter(sub_test, Interactor_A == overlap_table[i,1]), Interactor_B))
  overlap_table[i,3] <- nrow(na.omit(select(filter(sub_test, Interactor_A == overlap_table[i,1]), PPI_score)))
  overlap_table[i,4] <- nrow(na.omit(select(filter(sub_test, Interactor_A == overlap_table[i,1]), GI_score)))
  overlap_table[i,5] <- nrow(na.omit(select(filter(sub_test, Interactor_A == overlap_table[i,1]), CE_score)))
}

#get number of overlap in pairs of webs (GI_PPI, GI_CE, PPI_CE)

#for PPI_GI overlap

PPI_GI_num_log <- cbind(sub_test$Interactor_A)
PPI_GI_num_log <- as.data.frame(PPI_GI_num_log)
for (i in seq_along(1:nrow(sub_test))){
  PPI_GI_num_log[i,2:3] <- is.na(select(sub_test[i,], PPI_score, GI_score))
}
for (i in seq_along(1:nrow(PPI_GI_num_log))){
  if (sum(PPI_GI_num_log[i,]$V2, PPI_GI_num_log[i,]$V3) == 0){
    PPI_GI_num_log[i,4] <- 1
  }else{
    PPI_GI_num_log[i,4] <- 0
  }
}
colnames(PPI_GI_num_log) <- c('Interactor_A', 'PPI_log', 'GI_log', 'sum')
PPI_GI_num_log[,1] <- as.character(PPI_GI_num_log[,1])

for (i in seq_along(1:nrow(overlap_table))){
  overlap_table[i,6] <- sum(PPI_GI_num_log[PPI_GI_num_log$Interactor_A == overlap_table[i,1],]$sum)
}

#for GI_CE overlap

GI_CE_num_log <- cbind(sub_test$Interactor_A)
GI_CE_num_log <- as.data.frame(GI_CE_num_log)
for (i in seq_along(1:nrow(sub_test))){
  GI_CE_num_log[i,2:3] <- is.na(select(sub_test[i,], GI_score, CE_score))
}
for (i in seq_along(1:nrow(GI_CE_num_log))){
  if (sum(GI_CE_num_log[i,]$V2, GI_CE_num_log[i,]$V3) == 0){
    GI_CE_num_log[i,4] <- 1
  }else{
    GI_CE_num_log[i,4] <- 0
  }
}
colnames(GI_CE_num_log) <- c('Interactor_A', 'GI_log', 'CE_log', 'sum')
GI_CE_num_log[,1] <- as.character(GI_CE_num_log[,1])

for (i in seq_along(1:nrow(overlap_table))){
  overlap_table[i,7] <- sum(GI_CE_num_log[GI_CE_num_log$Interactor_A == overlap_table[i,1],]$sum)
}

#for CE_PPI overlap

CE_PPI_num_log <- cbind(sub_test$Interactor_A)
CE_PPI_num_log <- as.data.frame(CE_PPI_num_log)
for (i in seq_along(1:nrow(sub_test))){
  CE_PPI_num_log[i,2:3] <- is.na(select(sub_test[i,], CE_score, PPI_score))
}
for (i in seq_along(1:nrow(CE_PPI_num_log))){
  if (sum(CE_PPI_num_log[i,]$V2, CE_PPI_num_log[i,]$V3) == 0){
    CE_PPI_num_log[i,4] <- 1
  }else{
    CE_PPI_num_log[i,4] <- 0
  }
}
colnames(CE_PPI_num_log) <- c('Interactor_A', 'CE_log', 'PPI_log', 'sum')
CE_PPI_num_log[,1] <- as.character(CE_PPI_num_log[,1])

for (i in seq_along(1:nrow(overlap_table))){
  overlap_table[i,8] <- sum(CE_PPI_num_log[CE_PPI_num_log$Interactor_A == overlap_table[i,1],]$sum)
}

#three overlaps

all_num_log <- cbind(sub_test$Interactor_A)
all_num_log <- as.data.frame(all_num_log)
for (i in seq_along(1:nrow(sub_test))){
  all_num_log[i,2:4] <- is.na(select(sub_test[i,], CE_score, PPI_score, GI_score))
}
for (i in seq_along(1:nrow(all_num_log))){
  if (sum(all_num_log[i,]$V2, all_num_log[i,]$V3, all_num_log[i,]$V4) == 0){
    all_num_log[i,5] <- 1
  }else{
    all_num_log[i,5] <- 0
  }
}
colnames(all_num_log) <- c('Interactor_A', 'CE_log', 'PPI_log', 'GI_score', 'sum')
all_num_log[,1] <- as.character(all_num_log[,1])
overlap_table$all_overlap_num <- NA

for (i in seq_along(1:nrow(overlap_table))){
  overlap_table[i,9] <- sum(all_num_log[all_num_log$Interactor_A == overlap_table[i,1],]$sum)
}

colnames(overlap_table) <- c('Interactor', 'neighbours_num', 'PPI_num', 'GI_num', 'CE_num', 'PPI_GI_num',
                             'GI_CE_num', 'CE_PPI_num', 'all_overlap_num')

#count the overlap index for all of three molecular webs for each gene
overlap_table$PPI_GI_OI <- NA
overlap_table$GI_CE_OI <- NA
overlap_table$CE_PPI_OI <- NA
overlap_table$THREE_OI <- NA

for (i in seq_along(1:nrow(overlap_table))){
  overlap_table[i,10] <- (overlap_table[i,]$PPI_GI_num)/min(overlap_table[i,]$PPI_num, overlap_table[i,]$GI_num)
  overlap_table[i,11] <- (overlap_table[i,]$GI_CE_num)/min(overlap_table[i,]$GI_num, overlap_table[i,]$CE_num)
  overlap_table[i,12] <- (overlap_table[i,]$CE_PPI_num)/min(overlap_table[i,]$CE_num, overlap_table[i,]$PPI_num)
  overlap_table[i,13] <- (overlap_table[i,]$all_overlap_num)/min(overlap_table[i,]$CE_num, overlap_table[i,]$PPI_num, overlap_table[i,]$GI_num)
  
}
overlap_table[is.nan(overlap_table$PPI_GI_OI),]$PPI_GI_OI <- 0
overlap_table[is.nan(overlap_table$GI_CE_OI),]$GI_CE_OI <- 0
overlap_table[is.nan(overlap_table$CE_PPI_OI),]$CE_PPI_OI <- 0
overlap_table[is.nan(overlap_table$THREE_OI),]$THREE_OI <- 0

