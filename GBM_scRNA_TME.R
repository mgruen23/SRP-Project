setwd("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project")
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)


#####read in expression data
exp_tbl_SMARTSEQ <- data.frame(fread("GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv"), row.names = 1)

# create seurat object
seur_dar_SMARTSEQ <- CreateSeuratObject(exp_tbl_SMARTSEQ, project = "Neftel SMARTSEQ", assay = "RNA")

#read in metadata
meta_data <- read.csv(file = "GSE131928_single_cells_tumor_name_and_adult_or_peidatric.csv", sep = ",", header = T)
subject_names <- meta_data[1:7930,"tumour.name"] #because only the first 7930 are smartseq samples

#set orig.ident for each cell to their sample name
seur_dar_SMARTSEQ@meta.data$orig.ident <- subject_names


Idents(seur_dar_SMARTSEQ) <- "orig.ident"
######QC STEPS


#check percentage of mitochondrial genes (if too high means cell is dying or dead and maybe should be eliminated from analysis)
seur_dar_SMARTSEQ[["percent.mt"]] <- PercentageFeatureSet(seur_dar_SMARTSEQ, pattern = "^MT-")

# QC Violin Plots
VlnPlot(seur_dar_SMARTSEQ, features = "nFeature_RNA")
VlnPlot(seur_dar_SMARTSEQ, features = "nCount_RNA")
VlnPlot(seur_dar_SMARTSEQ, features = "percent.mt")

feature_cutoff <- 10000


seur_dar_SMARTSEQ <- subset(seur_dar_SMARTSEQ, subset = nFeature_RNA < feature_cutoff)


#Normalize data -   normalizes individual gene expression counts for a given gene by overall counts for that cell. This is done
#to account for the fact that some cells simply have more reads overall and this has nothing to do with the biology. 
seur_dar <- NormalizeData(seur_dar, normalization.method = "LogNormalize", scale.factor = 10000)
 

#nfeatures?
seur_dar_SMARTSEQ <- FindVariableFeatures(seur_dar_SMARTSEQ, selection.method = "vst")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seur_dar_SMARTSEQ), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seur_dar_SMARTSEQ, log=T) #without labels
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) #with labels
plot1 + plot2

#scaling - shift and scale gene expresion data so that mean expression for a given gene across cells is 0
#and the variability across cells is 1. This is neccesary to do prior to PCA so that highly overexpressed
#genes don't dominate and overly influence the top principal components. 

all.genes <- rownames(seur_dar_SMARTSEQ)
seur_dar_SMARTSEQ <- ScaleData(seur_dar_SMARTSEQ, features = all.genes)


saveRDS(seur_dar_SMARTSEQ, file = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Clustered Seurat by sample/clustered_Neftel_SMARTSEQ.RDS")


#function to run workflow through a set of samples


seur_dar_SMARTSEQ <- readRDS(file = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Clustered Seurat by sample/clustered_Neftel_SMARTSEQ.RDS")

#add parameters for platform, res

#Platform = "SMARTSEQ" or "10X"
subset_workflow <- function(main_seurat, platform, resolution = "", UMAP_folder, 
                            clustered_object_folder, findmarkers_folder, annotations_folder){
  
  res = resolution
  if(platform == "SMARTSEQ"){
    dimens = 15
    elbow_dims = 30
    if (res == ""){
      res = 0.8
    }
  }
  if(platform == "10X"){
    dimens = 25
    elbow_dims = 50
    if (res == ""){
      res = 1.0
    }
  }
  

  
  sample_ids <-unique(main_seurat$orig.ident)
  
  for(sample in  sample_ids){
    sample_seurat <- subset(main_seurat, subset = orig.ident == sample)
    sample_seurat <- RunPCA(sample_seurat, features = VariableFeatures(object = sample_seurat))
    sample_seurat <- FindNeighbors(sample_seurat, dims = 1:dimens)
    sample_seurat <- FindClusters(sample_seurat, resolution = res, algorithm = 1)
    sample_seurat <- RunUMAP(sample_seurat, dims = 1:dimens)
    
    sample_seurat.markers <- FindAllMarkers(sample_seurat, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.0)
    sample_seurat.markers %>%
      group_by(cluster) %>%
      slice_max(n = 2, order_by = avg_log2FC)

    write.table(sample_seurat.markers, file=paste(findmarkers_folder, "find_all_markers_", platform, "_" , sample, ".txt"),
                row.names = F, col.names = T, quote =F, sep = "\t")

    saveRDS(sample_seurat, file=paste(clustered_object_folder, "clustered_Neftel_", platform, "_", sample, ".RDS"))


    #save umap plot
    png(file=paste(UMAP_folder, "umap_", platform, "_", sample, ".png"), units = "in", type = "cairo",
        width= 8, height=6 , res = 1200)
    d <- DimPlot(sample_seurat, reduction = "umap")
    plot(d)
    dev.off()


    #save elbow plot
    png(file=paste(UMAP_folder, "elbow_plot_", platform, "_", sample, ".png"), units = "in", type = "cairo",
        width= 8, height=6 , res = 1200)
    e <- ElbowPlot(sample_seurat, ndims=elbow_dims)
    plot(e)
    dev.off()

#run 
    old_wd <- getwd()
    setwd(annotations_folder)

    test_template <- read.table("Allen_humanareas_clustermarkers_v3_for_celltype_annotation.txt", sep = "\t", header = T, stringsAsFactors = F)
    test_template2 <- read.table("Nowakowski_2017_Science_cell_type_DE_table_for_celltype_annotation.txt", sep = "\t", header = T, stringsAsFactors = F)
    test_template3 <- read.table("Bhaduri_2020_oRG_Cell_Type_DE_for_celltype_annotation.txt", sep = "\t", header = T, stringsAsFactors = F)
    test_template4 <-read.table("GBM Seurat_v4 gbm_filtered_batch3 final_integrated_STACAS clustering Cell.Type MAST DE results.txt", sep = "\t", header = T, stringsAsFactors = F)[,-1] #remove row names
    test_template5 <-read.table("GBM Seurat_v4 gbm_filtered_batch3 final_integrated_STACAS clustering seurat clusters MAST DE results annotated.txt", sep = "\t", header = T, stringsAsFactors = F)[,-1] #remove row names
    #switching "cluster" col with "Annotation" column for template 5
    test_template5$cluster <- test_template5$Annotation

    # # read in un-annotated cluster DE results
    test_results <- sample_seurat.markers

        # filter  DE results, adjust pct.1 pct.2
    filt_test_results <- prep.cluster.DE(test_results, p_val_thresh = 0.1, logfc_thresh = 0.25)
    filt_template4 <- prep.cluster.DE(test_template4)
    filt_template5 <- prep.cluster.DE(test_template5)


    # get annotation results for each cluster with the given template
    cor_res <- get.cluster.correlations(de_results = filt_test_results, template_results = test_template)
    cor_res2 <- get.cluster.correlations(de_results = filt_test_results, template_results = test_template2)
    cor_res3 <- get.cluster.correlations(de_results = filt_test_results, template_results = test_template3)
    cor_res4 <- get.cluster.correlations(de_results = filt_test_results, template_results = filt_template4)
    cor_res5 <- get.cluster.correlations(de_results = filt_test_results, template_results = filt_template5)

    #save results using all 5 templates
    final_res <- combine.results(list(allen = cor_res, nowakowski = cor_res2, bhaduri = cor_res3, Bayley = cor_res4, Bayley_annotations = cor_res5))

    prefix <- "Unsupervised clustering"

    for(i in 1:length(final_res)){
      temp_name <- names(final_res)[i]
      write.table(final_res[[i]], paste(prefix, platform, sample, temp_name, "cluster correlations.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
    }

    setwd(old_wd)

  }
  
}
  





  #plot UMAP and save to separate folder
  png(file=paste(annotated_UMAP_folder, "annotated_umap_", platform, "_", sample, ".png"), units = "in", type = "cairo",
      width= 8, height=6 , res = 1200)
  d <- DimPlot(clustered_object, reduction = "umap")
  plot(d)
  dev.off()
  
}
#script call

subset_workflow(main_seurat = seur_dar_SMARTSEQ, platform = "SMARTSEQ",
                UMAP_folder = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Test/UMAP_test/",
                clustered_object_folder = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Test/clustered_object_test/", 
                findmarkers_folder =  "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Test/findmarkers_test/", 
                annotations_folder =  "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Test/annotations_test/")

rewrite_annotations(clustered_R_file = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Test/clustered_object_test/ clustered_Neftel_ SMARTSEQ _ MGH101 .RDS",
                    cluster_names = c("immune", "RG", "NPC", "OPC"))
