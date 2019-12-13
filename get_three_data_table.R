#=============MERGED===========================

#Gene_interaction_data

GI_df <- read.table('GI_short.txt', header = T)

#PEARSONS_CORRELATION_COEF

GI_df[,1] <- as.character(GI_df[,1])
GI_df[,2] <- as.character(GI_df[,2])
#hist(abs(GI_df$GI_score))
str(GI_df)

#===========================================================================

#Protein-protein_interaction_data

PPI_df <- read.table('PPI_short.txt', header = T)

#Quantitative Score. This will be a positive for negative value recorded by the original 
#publication depicting P-Values, Confidence Score, SGA Score,etc

PPI_df[,1] <- as.character(PPI_df[,1])
PPI_df[,2] <- as.character(PPI_df[,2])
#hist(PPI_df$PPI_score)
str(PPI_df)

#==========================================================================

#Co-expression_data

CE_df <- read.table('CE_short.txt', header = T)
CE_df[,1] <- as.character(CE_df[,1])
CE_df[,2] <- as.character(CE_df[,2])
#hist(CE_df$LLS)
str(CE_df)



#========plot_MERGED===============================
#colnames(CE_df) <- c('Gene_A','Gene_B','Score')
#colnames(PPI_df) <- c('Interactor_A', 'Interactor_B', 'Score')
#colnames(GI_df) <- c('Query_Strain', 'Array_Strain', 'Score')


par(mfrow=c(1,3)) 
hist(GI_df$GI_score, main = "Gene interaction", xlab = "Pearson's coefficient", ylab = "№ of interactions")
hist(PPI_df$PPI_score, main = "Protein-protein interaction", xlab = "Quantitative Score", ylab = "№ of interactions")
hist(CE_df$LLS, main = "Gene coexpression", xlab = "LLS", ylab = "№ of interactions")