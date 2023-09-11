setwd("/Users/matthewgruen/Desktop/Graeber Lab /R exercises/IVYGAP")
library(dplyr)


TPM_matrix <- read.table("IVYGAP_TPM_matrix.txt", sep = "\t", header = T, stringsAsFactors = F, fill =T, row.names =1, check.names = F)
annotations <- read.table("IVYGAP_sample_annotations.txt", sep = "\t", header = T, stringsAsFactors = F, fill =T, row.names =1)

setwd("/Users/matthewgruen/Desktop/Graeber Lab /R exercises/Immune Response Project")

AB_signature <- read.table("groupAB_signature.txt", sep = "\t", header = T, stringsAsFactors = F, fill =T, row.names =1)[,1]


#do sum of z score analysis and save result 

#restrict TPM_matrix to just gene names in immune gene signature 

TPM_matrix_restricted <- TPM_matrix[AB_signature,]

#restrict TPM matrix to new subsetted signatures
TPM_matrix_restricted_immune_sub <- TPM_matrix[immune_sig_subsetted,]
TPM_matrix_restricted_vasculature_sub <- TPM_matrix[vasculature_sig_subsetted,]
TPM_matrix_restricted_ecm_sub <- TPM_matrix[ecm_sig_subsetted,]


#take log2(x+1) of all values in TPM_matrix 
TPM_matrix_restricted_log <- log2(TPM_matrix_restricted+1)

TPM_matrix_restricted_immune_sub_log <- log2(TPM_matrix_restricted_immune_sub+1)
TPM_matrix_restricted_vasculature_sub_log <- log2(TPM_matrix_restricted_vasculature_sub+1)
TPM_matrix_restricted_ecm_sub_log <- log2(TPM_matrix_restricted_ecm_sub+1)

#scale values of matrix to get z scores by centering and dividing by standard deviation of the column
#do transposing because we want to normalize across the different samples for a given gene, not across the genes for a given sample
TPM_matrix_zscored <- t(scale(t(TPM_matrix_restricted_log)))

TPM_matrix_zscored_immune_sub <- t(scale(t(TPM_matrix_restricted_immune_sub_log)))
TPM_matrix_zscored_vasculature_sub <- t(scale(t(TPM_matrix_restricted_vasculature_sub_log)))
TPM_matrix_zscored_ecm_sub <- t(scale(t(TPM_matrix_restricted_ecm_sub_log)))

#create new row at bottom of matrix with sum of z-scores
sums <- colSums(TPM_matrix_zscored, na.rm = T)
TPM_matrix_zscored <- rbind(TPM_matrix_zscored, sums)

sums_immune_sub <- colSums(TPM_matrix_zscored_immune_sub, na.rm = T)
means_immune_sub <- colMeans(x=TPM_matrix_zscored_immune_sub, na.rm = T)
TPM_matrix_zscored_immune_sub <- rbind(TPM_matrix_zscored_immune_sub, means_immune_sub)

sums_vasculature_sub <- colSums(TPM_matrix_zscored_vasculature_sub, na.rm = T)
means_vasculature_sub <- colMeans(x=TPM_matrix_zscored_vasculature_sub, na.rm = T)
TPM_matrix_zscored_vasculature_sub <- rbind(TPM_matrix_zscored_vasculature_sub, means_vasculature_sub)

sums_ecm_subbed <- colSums(TPM_matrix_zscored_ecm_sub, na.rm = T)
means_ecm_subbed <- colMeans(x=TPM_matrix_zscored_ecm_sub, na.rm = T)
TPM_matrix_zscored_ecm_sub <- rbind(TPM_matrix_zscored_ecm_sub, means_ecm_subbed)

#remove all NA rows
TPM_matrix_zscored_immune_sub <- na.omit(TPM_matrix_zscored_immune_sub)
TPM_matrix_zscored_vasculature_sub <- na.omit(TPM_matrix_zscored_vasculature_sub)
TPM_matrix_zscored_ecm_sub <- na.omit(TPM_matrix_zscored_ecm_sub)

write.table(TPM_matrix_zscored, file = "IVYGAP_scored.txt", sep = "\t", quote = F, row.names = T, col.names = T)
setwd("/Users/matthewgruen/Desktop/Graeber Lab /R exercises/IVYGAP")
write.table(TPM_matrix_zscored_immune_sub, file = "IVYGAP_scored_immune_sub.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(TPM_matrix_zscored_vasculature_sub, file = "IVYGAP_scored_vasculature_sub.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(TPM_matrix_zscored_ecm_sub, file = "IVYGAP_scored_ecm_sub.txt", sep = "\t", quote = F, row.names = T, col.names = T)


