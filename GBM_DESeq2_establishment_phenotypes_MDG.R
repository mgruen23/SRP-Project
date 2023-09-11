library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(RColorBrewer)

###################
#### Read Data ####
###################

### Count matrix for DE ###
setwd("/Users/matthewgruen/Desktop/Graeber Lab /Manuscript")

count_matrix <- read.table("Nathanson_RSEM_all+bbseal_xenos_raw_counts_20220830.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1)

tmed_anno <- read.table("Nathanson TMED annotation table from Box.tsv", sep = "\t", header = T)
rownames(tmed_anno) <- tmed_anno$PT.ID

# read in annotation information
annotation <- read.table(file ="Sequencing Metadata batch15.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "", fileEncoding="latin1")
rna_anno <- annotation[!is.na(annotation$"RNA.Batch.."),]
rownames(rna_anno) <- rna_anno$Short.ID




# As I mentioned, I've written code for reading in the data
#and there is old code for running DESeq2 that grabs our comparison of interest (TME.D vs. TME.I).
#If you can update the deseq2 code to use our DESeq_functions.R script and then add in the FGSEA analysis
#using the fgsea_functions.R that would be fantastic!


###################
#### Filter and process count and annotation matrices ####
###################


#Warning message:
# Problem with `mutate()` column `^^--arrange_quosure_1`.
# ℹ `^^--arrange_quosure_1 = as.numeric(Line..)`.
# ℹ NAs introduced by coercion 
filt_anno <- rna_anno %>% 
	arrange(as.numeric(Line..),Line..,Sample.Type, Short.ID) %>%
	filter(!FLAG.RNA %in% c("terminate", "relabel", "remove")) %>%
	filter(Proj.ID == 3 | Cross.Project == "Y")

pair_ind <- match(colnames(count_matrix), rna_anno$Short.ID)
nat_anno <- rna_anno[pair_ind,]

pt_count_matrix <- count_matrix[,nat_anno$Sample.Type %in% c("PT")]

#changed to pt_count_matrix
tmed_anno <- tmed_anno[colnames(pt_count_matrix),]

#use table to check number of entries by parameter
filt_tmed_anno <- tmed_anno[tmed_anno$EP %in% c("TME.I", "TME.D"),]

filt_pt_count_matrix <- pt_count_matrix[rowSums(pt_count_matrix > 50) > 10,rownames(filt_tmed_anno)]

filt_tmed_anno$RNA.Batch <- nat_anno[rownames(filt_tmed_anno),"RNA.Batch.."]

filt_tmed_anno <- data.frame(lapply(filt_tmed_anno, factor))



#load in functions
source("/Users/matthewgruen/Desktop/Graeber Lab /Manuscript/DESeq_functions.R") 
source("/Users/matthewgruen/Desktop/Graeber Lab /Manuscript/FGSEA_functions.R") 

#Generate DESEQ data set
pt_dds <- DESeqDataSetFromMatrix(countData = filt_pt_count_matrix, colData = filt_tmed_anno, design = ~ RNA.Batch + EP)
#Run DESEQ
seq <- DESeq(pt_dds)
#check results names
resultsNames(seq)
## processing steps and then extract results table
#establishment/engraftment phenotype = EP
pt_res <- get.contrast.results(seq, contrast = c("EP", "TME.D", "TME.I"), lfc_shrink = T, shrink_type = "ashr", filename = F)

 #change ord to just pt_res

sig_pt_res <- pt_res[ord_pt_res$padj < 0.05 & abs(ord_pt_res$log2FoldChange) > log2(2),]


#add logp column for GSEA ranking
signed_sig_log_p <- add.signed.log.p(sig_pt_res, col_name = "padj", sign_name = "stat")
signed_res_log_p <- add.signed.log.p(pt_res, col_name = "padj", sign_name = "stat")
  
# Plots -log of padj vs. log2fc. Colors by sig decrease, no change, or sig increase
#should look like a v shape
#called volcano plot
sig <- as.numeric(abs(ord_pt_res$log2FoldChange) > 1 & ord_pt_res$padj < 0.05) * sign(ord_pt_res$log2FoldChange)
plot(ord_pt_res$log2FoldChange, -log10(ord_pt_res$padj), col = as.numeric(factor(sig)), xlim = c(-8,8))
abline(h = -log10(0.05))

#print out DESEQ results table
write.table(signed_res_log_p, "UCLA TME.D vs. TME.I DESeq2 lfcS results batch15.tsv", sep = "\t", quote = F, row.names = T, col.names = NA)

#all genes with increased expression in TME.D are labelled as TME.D and all genes with increased expression in TME I are
#labelled as  labelled as TME.I. 
labels <- ifelse(sig_pt_res$log2FoldChange > 0, "TME.D", "TME.I")
names(labels) <- rownames(sig_pt_res)
#check row, column names
write.table(labels, "UCLA TME.D vs. TME.I DESeq2 lfcS signature batch15.tsv", sep = "\t", quote = F, row.names = T, col.names = T)


#FGSEA
#should this be only sig or all ordered results?
rnk <- get.rnk.vector(signed_res_log_p, column_name = "signed.log.p")
msig_df <- load.MSigDB(species = "Homo sapiens")
genesets <- get.MSigDB.genesets(msig_df, genesets = c("KEGG", "REACTOME", "H$"))
gsea_res <- run.FGSEA(rnk, genesets, nproc = 2, minGenes = 3, maxGenes = 5000, reformat = T, filename = F, minP = 1e-30)

#check row, column names
write.table(gsea_res, "UCLA TME.D vs. TME.I FGSEA results batch15.tsv", sep = "\t", quote = F, row.names = T, col.names = NA)
#file named based on DESEQ comparison