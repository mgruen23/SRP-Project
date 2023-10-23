setwd("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project")
#load libraries

library("CellChat")
library("patchwork")
options(stringsAsFactors = FALSE)

#cellchat processing to perform on each sample's seurat object
cellchat_processing <- function(clustered_seurat_folder, platform, cellchat_folder){
  old_dir <- getwd()
  setwd(clustered_seurat_folder)
  files <- list.files(clustered_seurat_folder)
  

  #create dataframe to be iteratively added to, with data from each subject
  all_subjects_overall_interactions <-data_frame()
  colnames(all_subjects_overall_interactions) <- c("Sample","Source","Target","Number","Weight")
  all_subjects_pathways_probabilities <-data_frame()
  colnames(all_subjects_pathways_probabilities) <- c("Sample", "Pathway","Source", "Target", "Probability")
  
  for (file in files){
    seur <- readRDS(file)
    #pull out sample name from seurat object - make proper file name
    cellchat <- createCellChat(seur, group.by = "ident")
    
    CellChatDB <- CellChatDB.human
    CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 4) # do parallel
    
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    df.net <- subsetCommunication(cellchat)
    df.netP <- subsetCommunication(cellchat, slot.name = 'netP')
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    #save cellchat object to desired folder with proper naming scheme
    sample <- cellchat@meta$orig.ident[1]
    cellchat_file_name <- paste(cellchat_folder, "cellchat_", platform, "_", sample, ".RDS")
    saveRDS(cellchat, cellchat_file_name)
    
    #add any other results tables that should be saved
    count_matrix <- cellchat@net$count #returns matrix 
    weight_matrix <- cellchat@net$weight
    count_dataframe <- reshape2::melt(count_matrix)   #var1 is sender, var2 is receiver
    weight_dataframe <- reshape2::melt(weight_matrix) 
    
    #create circle plots showing signaling to and from immune cluster
    if(!all(cellchat@idents != "N-Macrophage")){
      png(file=paste("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/cellchat_save/",sample ,".png"), units = "in", type = "cairo",
          width= 8, height=6 , res = 1200)
      groupSize <- as.numeric(table(cellchat_test@idents))
      par(mfrow = c(1,2), xpd=TRUE)
      netVisual_circle(cellchat_test@net$count, targets.use = c("N-Macrophage"), vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste(sample, "Number of interactions"))
      dev.off()
    }
  
  }

}



#create circle plot for exisitng cellchat object signaling to a particular group
circle_plot <- function(cellchat_path, save_folder, target = "N-Macrophage"){
  cellchat <- readRDS(cellchat_path)
  
  sample <- cellchat@meta$orig.ident[1]
  if(!all(cellchat@idents != "N-Macrophage")){
    png(file=paste(save_folder,sample,".png"), units = "in", type = "cairo",
        width= 8, height=6 , res = 1200)
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = sample)
    dev.off()
  }

}

circle_plot("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/SMARTSEQ/cellchat_SMARTSEQ/ cellchat_ SMARTSEQ _ BT771 .RDS","/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/SMARTSEQ/cellchat_SMARTSEQ/")


sample <- cellchat@meta$orig.ident[1]
if(!all(cellchat@idents != "N-Macrophage")){
  png(file=paste(save_folder,sample,".png"), units = "in", type = "cairo",
      width= 8, height=6 , res = 1200)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, targets.use = c("N-Macrophage"), vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = sample)
  dev.off()
}

if(!all(cellchat@idents != "N-Macrophage")){
  print("test")
}



#calling cellchat function

cellchat_processing(clustered_seurat_folder = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/SMARTSEQ/annotated_object_SMARTSEQ/", 
                    platform = "SMARTSEQ", 
                    cellchat_folder = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/SMARTSEQ/cellchat_SMARTSEQ/")


cellchat_processing(clustered_seurat_folder = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/10X/annotated_object_10X/", 
                    platform = "10X", 
                    cellchat_folder = "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/10X/cellchat_10X/")
