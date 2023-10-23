library(survival)
library(maxstat)
library(glmnet)
library(dplyr)


parse.id <- function(id){
	first <- toupper(substr(id, 6, 7))	
	second <- toupper(substr(id, 9, 12))
	final <- paste(first, second, sep = "")
	return(final)
}

setwd("/Users/matthewgruen/Desktop/Graeber Lab /Survival Analysis/")

#Seq data

# load in expression data, reformat TCGA IDs, find sample ovelap with survival data
tcga_exp_mat <- read.table("TOIL_TCGA_RSEM_raw_DESeq2_norm_CPM_20200609.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
ids <- parse.id(colnames(tcga_exp_mat))
shared_ids <- intersect(ids, surv_df$id)

# filter both survival and expression datasets to overlapping samples
filt_exp_mat <- tcga_exp_mat[,na.omit(match(shared_ids, ids))]
filt_surv_df <- surv_df[na.omit(match(shared_ids, surv_df$id)),]

surv_df <- read.table("TCGA GBM survival data CBioPortal.tsv", sep = "\t", header = T, stringsAsFactors = F)

# load in expression data, reformat TCGA IDs, find sample ovelap with survival data
tcga_exp_mat <- read.table("TOIL_TCGA_RSEM_raw_DESeq2_norm_CPM_20200609.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
ids <- parse.id(colnames(tcga_exp_mat))
shared_ids <- intersect(ids, surv_df$id)

# filter both survival and expression datasets to overlapping samples
filt_exp_mat <- tcga_exp_mat[,na.omit(match(shared_ids, ids))]
filt_surv_df <- surv_df[na.omit(match(shared_ids, surv_df$id)),]

#read in column of immune scores for TCGA samples
immune_signature <- read.table("groupAB_signature.txt", sep = "\t", stringsAsFactors = F, header = T, fill = NA)[,1]
#z-score analysis


#restrict TPM_matrix to just gene names in immune gene signature 

filt_exp_mat <- filt_exp_mat[immune_signature,]

#take log2(x+1) of all values
filt_exp_mat_log <- log2(filt_exp_mat+1)

#scale values of matrix to get z scores by centering and dividing by standard deviation of the column
#do transposing because we want to normalize across the different samples for a given gene, not across the genes for a given sample
mat_zscored <- t(scale(t(filt_exp_mat_log)))

#create new row at bottom of matrix with sum of z-scores
sums <- colSums(mat_zscored, na.rm = T)
mat_zscored <- rbind(mat_zscored, sums)
#means <- colMeans(mat_zscored[-1,], na.rm = T)
#mat_zscored <- rbind(mat_zscored, means)
#mat_zscored <- na.omit(mat_zscored) --do I need to cut out any sample that has NA for any of the genes. 

immune_score <- sums

filt_surv_df$immune_score <- rep(0, nrow(filt_surv_df))
filt_surv_df$immune_score <- immune_score


#run log rank test with calculated cutoff for immune score


surv_model <- Surv(time = filt_surv_df$last_contact, event = as.numeric(as.factor(filt_surv_df$vital)))

# find the best way to split our dataset into two groups based on Log Rank Test significance
# this is an example, other ways of splitting the data into groups (by median, top third vs. bottom third)
cutoff <- maxstat.test(surv_model ~ filt_surv_df$immune_score, data = filt_surv_df, smethod = "LogRank", pmethod = "condMC")$estimate
#binary_labels <- ifelse(filt_surv_df$immune_score >= cutoff, 1, 0)
# or example splitting by median expression of EGFR (leave commented out if you prefer optimal cutoff)
binary_labels <- ifelse(filt_surv_df$immune_score >= median(filt_surv_df$immune_score), 1, 0)
# get number of samples in each group
num_samples <- table(binary_labels)

# perform log rank test, extract chi-squared statistic and convert to p-value
log_rank_test_immune <- survdiff(surv_model ~ binary_labels)
chisq <- log_rank_test_immune[[5]]
p_val <- pchisq(chisq, df = 1, lower.tail = FALSE)

#splitting into two groups by median leads to p-value of .160. 
#transition into single cell to get better signature

#greatest difference between 20 and 50 months for survival based on immune score. -add

# even with the best possible cutoff, EGFR groups are not quite significant (p = 0.06)

# create Kaplan-Meier plot 
png(file="tcga_survival.png", units = "in", type = "cairo",
    width= 8, height=6 , res = 1200)
plot(survfit(surv_model ~ binary_labels), col = c("black", "red"), xlab = "Time (Months)", ylab = "Survival Probability", main = "TCGA GBM\nKaplan-Meier Curve")
legend("topright", legend = c(paste0("high (", num_samples[2],")"), paste0("low (", num_samples[1],")")), pch = 3, col = c("red", "black"), title = "Immune Score")
#add log rank and cox p-values to bottom of plot
dev.off()
#add p values to plot

#run cox proportional hazards model
plot(surv_model)
# model with Cox proportional hazards model
cox_model <- coxph(Surv(last_contact, as.numeric(as.factor(vital))) ~ immune_score, data = filt_surv_df)
summary(cox_model)


cox_model <- coxph(Surv(last_contact, as.numeric(as.factor(vital))) ~ gender + immune_score, data = filt_surv_df)
summary(cox_model)
