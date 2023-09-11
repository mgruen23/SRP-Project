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

###################
#### Load Data ####
###################

setwd("/Users/matthewgruen/Desktop/Graeber Lab /Survival Analysis/")

#Seq data

#bio_surv_data <- read.table("gbm_tcga_clinical_data.tsv", sep = "\t", stringsAsFactors = F, header = T)
#surv_df <- data.frame(id = as.character(apply(cbio_surv_data, 1, function(x){
#	parse.id(x[2])
#})), last_contact = cbio_surv_data[,36], birth = cbio_surv_data[,4], vital = cbio_surv_data[,67], gender = cbio_surv_data[,55], stringsAsFactors = F)
#surv_df <- surv_df[!apply(surv_df, 1, function(x) any(is.na(x))),]
#write.table(surv_df, "TCGA GBM survival data CBioPortal.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
surv_df <- read.table("TCGA GBM survival data CBioPortal.tsv", sep = "\t", header = T, stringsAsFactors = F)

# load in expression data, reformat TCGA IDs, find sample ovelap with survival data
tcga_exp_mat <- read.table("TOIL_TCGA_RSEM_raw_DESeq2_norm_CPM_20200609.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
ids <- parse.id(colnames(tcga_exp_mat))
shared_ids <- intersect(ids, surv_df$id)

# filter both survival and expression datasets to overlapping samples
filt_exp_mat <- tcga_exp_mat[,na.omit(match(shared_ids, ids))]
filt_surv_df <- surv_df[na.omit(match(shared_ids, surv_df$id)),]

######################
#### Add Variable ####
######################

# as an example we can see if there is a significant prognostic benefit to the GCIMP methylation subtype
# we will use this in the "Cox Model" section below
# read in data, clean up ids, match to survival data
gbm_class <- read.table("TCGA_2013_GBM_classification.txt", sep = "\t", stringsAsFactors = F, header = T, fill = NA)
gbm_class$id <- parse.id(gbm_class$Case.ID)
m_ind <- match(filt_surv_df$id, gbm_class$id)
filt_surv_df$gcimp <- gbm_class$G.CIMP..methylation[m_ind]
# remove samples without data for our variable of interest
filt_surv_df <- filt_surv_df[!is.na(filt_surv_df$gcimp) & !filt_surv_df$gcimp == "",]

final_exp_mat <- filt_exp_mat[,!is.na(filt_surv_df$gcimp) & !filt_surv_df$gcimp == ""]


#######################
#### Log Rank Test ####
#######################

# for continuous variables, if we would like to plot a Kaplan-Meier curve showing the significance of
# our variable we need to split the samples into two groups based on our variable
# in this example,we will use the continuous expression values of EGFR for TCGA GBM samples
# here I add the log2(x+1) expression value to the `filt_surv_df` data.frame
# since final_exp_mat is already in the same order, we don't need to worry about that
filt_surv_df$EGFR <- log2(as.numeric(final_exp_mat["EGFR",])+1)

# create Survival object for survival analysis
surv_model <- Surv(time = filt_surv_df$last_contact, event = as.numeric(as.factor(filt_surv_df$vital)))

# find the best way to split our dataset into two groups based on Log Rank Test significance
# this is an example, other ways of splitting the data into groups (by median, top third vs. bottom third)
cutoff <- maxstat.test(surv_model ~ filt_surv_df$EGFR, data = filt_surv_df, smethod = "LogRank", pmethod = "condMC")$estimate
binary_labels <- ifelse(filt_surv_df$EGFR >= cutoff, 1, 0)
# or example splitting by median expression of EGFR (leave commented out if you prefer optimal cutoff)
#binary_labels <- ifelse(filt_surv_df$EGFR >= median(filt_surv_df$EGFR), 1, 0)
# get number of samples in each group
num_samples <- table(binary_labels)

# perform log rank test, extract chi-squared statistic and convert to p-value
log_rank_test_EGFR <- survdiff(surv_model ~ binary_labels)
chisq <- log_rank_test_EGFR[[5]]
p_val <- pchisq(chisq, df = 1, lower.tail = FALSE)
# even with the best possible cutoff, EGFR groups are not quite significant (p = 0.06)

# create Kaplan-Meier plot 
plot(survfit(surv_model ~ binary_labels), col = c("red", "black"), xlab = "Time (Months)", ylab = "Survival Probability", main = "TCGA GBM\nKaplan-Meier Curve")
legend("topright", legend = c(paste0("high (", num_samples[2],")"), paste0("low (", num_samples[1],")")), pch = 3, col = c("red", "black"), title = "EGFR Expression")

###################
#### Cox Model ####
###################

### G-CIMP ###


#plots the survival data (with confidence interval by default
plot(surv_model)
# model with Cox proportional hazards model
cox_model <- coxph(Surv(last_contact, as.numeric(as.factor(vital))) ~ gcimp, data = filt_surv_df)
summary(cox_model)
# in the summary, we pay attention to `exp(coef)` this is the hazard ratio
# values greater than 1 indicate increased risk, values below indicate decreased risk
# p-value determines whether this increase or decrease in risk is significant
# in this example being non G-CIMP is associated with a 3.4x increase in risk of death

# the advantage of the Cox model is that we can include covariates (like age and gender) to assess
# the significance of our variable of interest (G-CIMP status) while controlling for potential confounders
# In general you should always account for potential confounders, but it can potentially worsen or improve p-values
# for your variable of interest
cox_model <- coxph(Surv(last_contact, as.numeric(as.factor(vital))) ~ gender + birth + gcimp, data = filt_surv_df)
summary(cox_model)
# Age (coded as `birth`) has the most significant impact on survival time and G-CIMP status is no longer significant 
# when accounting for age (and gender), as well it's hazard ratio has greatly decreased

# GCIMP status is strongly associated with age, so when we control for age it became less significant in our Cox Model
boxplot(birth ~ gcimp, data = filt_surv_df)
t.test(birth ~ gcimp, data = filt_surv_df)

# Kaplan-Meier plot by G-CIMP status
surv_fit <- survfit(surv_model ~ filt_surv_df$gcimp)
plot(surv_fit, col = c("red", "black"), xlab = "Time (Months)", ylab = "Survival Probability", main = "TCGA GBM\nKaplan-Meier Curve")
legend("topright", legend = c("Positive (9)", "Negative (146)"), pch = 3, col = c("red", "black"), title = "G-CIMP Status")

# notice that log rank test is more similar to results from univariate Cox model
# HOWEVER they are not identical!! different stats
log_rank_test_gcimp <- survdiff(surv_model ~ filt_surv_df$gcimp)
chisq <- log_rank_test_gcimp[[5]]
p_val <- pchisq(chisq, df = 1, lower.tail = FALSE)

### EGFR ###

# For continuous variables like EGFR, we can plug them in directly without having to convert our data into groups
cox_model_EGFR_multi <- coxph(Surv(last_contact, as.numeric(as.factor(vital))) ~ gender + birth + EGFR, data = filt_surv_df)
summary(cox_model_EGFR_multi)
# We get more significant results for EGFR expression when accounting for age and gender
cox_model_EGFR_uni <- coxph(Surv(last_contact, as.numeric(as.factor(vital))) ~ EGFR, data = filt_surv_df)
summary(cox_model_EGFR_uni)











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

# # filter both survival and expression datasets to overlapping samples
# 
# ids <- parse.id(colnames(mat_zscored))
# shared_ids <- intersect(ids, surv_df$id)
# filt_surv_df <- surv_df[na.omit(match(shared_ids, surv_df$id)),]

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