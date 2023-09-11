setwd("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project")
#load libraries

library("CellChat")
library("patchwork")
options(stringsAsFactors = FALSE)

#Create CellChat object from Seurat object

seur_dar_SMARTSEQ_102 <- readRDS(file="/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Clustered Seurat by sample/clustered_Neftel_SMARTSEQ_102.RDS")
cellchat <- createCellChat(seur_dar_SMARTSEQ_102, group.by = "ident")
#test <-readRDS("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Clustered Seurat by sample/clustered_Neftel_SMARTSEQ_124.RDS")
#cellchat_124 <- createCellChat(seur_dar_124, group.by = "ident")
#cellchat_101 <- createCellChat(seur_dar_curr, group.by = "ident")
cellchat_102 <- createCellChat(seur_dar_102, group.by = "ident")
cellchat <- cellchat_102



#do cellchat processing on each sample's seurat object
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
      
    
    
    
    
  #combine count and weight into single dataframe 
#     count_dataframe[,4] <- weight_dataframe[,3]
#     overall_interactions <- count_dataframe
#     sample_names <- rep(sample, nrow(overall_interactions))
#     overall_interactions <- cbind(sample_names ,overall_interactions)
#     colnames(overall_interactions) <- c("Sample", "Source", "Target", "Number", "Weight")
# # add to large dataframe for all samples
#     all_subjects_overall_interactions <- rbind(all_subjects_overall_interactions, overall_interactions)
#     
#     sample_probabilities <- data_frame()
#     colnames(sample_probabilities) <- c("Pathway","Source", "Target", "Probability")
#     num_pathways <- length(cellchat@netP$pathways)
#     for (i in 1:num_pathways){
#       prob <- cellchat@netP$prob[,,i]
#       prob_df <- reshape2::melt(prob)
#       #add name of pathway as first column 
#       pathway_name <- rep(cellchat@netP$pathways[i], nrow(prob_df))
#       prob_df <- cbind(pathway_name, prob_df)
#       #add to dataframe for the entire sample
#       sample_probabilities <- rbind(sample_probabilities, prob_df)
#     }
#     sample_names <- rep(sample, nrow(sample_probabilities))
#     sample_probabilities <- cbind(sample_names, sample_probabilities)
#     colnames(sample_probabilities) <- c("Sample", "Pathway","Source", "Target", "Probability")
#   
#     #add probability data for entire sample to large dataframe for all samples
#     all_subjects_pathways_probabilities <- rbind(all_subjects_pathways_probabilities, sample_probabilities)
#     
    
    
  }
  setwd(old_dir)
  #save the two dataframes
  #write.table(all_subjects_overall_interactions, paste(cellchat_folder, "group2group_interactions.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
  #write.table(all_subjects_pathways_probabilities, paste(cellchat_folder, "pathway_probabilities.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)
  
  
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

circle_plot("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/SMARTSEQ/cellchat_SMARTSEQ/ cellchat_ SMARTSEQ _ MGH100 .RDS","/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/SMARTSEQ/cellchat_SMARTSEQ/")
cellchat_path <- "/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/10X/cellchat_10X/ cellchat_ 10X _ MGH102 .RDS"
cellchat <- readRDS(cellchat_path)


circle_plot("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/10X/cellchat_10X/ cellchat_ 10X _ MGH102 .RDS","/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/10X/cellchat_10X/")




circle_plot("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/SMARTSEQ/cellchat_SMARTSEQ/ cellchat_ SMARTSEQ _ MGH101 .RDS","/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/SMARTSEQ/cellchat_SMARTSEQ/")


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







cellchat_test <- readRDS(paste("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/cellchat_save/",  "cellchat_ 10X _ MGH102 .RDS"))
rankNet(cellchat_test, slot.name = "netP")


png(file=paste("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/cellchat_save/", "102_circle plots_num_interactions.png"), units = "in", type = "cairo",
    width= 8, height=6 , res = 1200)

groupSize <- as.numeric(table(cellchat_test@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_test@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()








#Data to be grabbed from the cellchat object

#net number(count) or strength(weight) of interaction between any two cell groups
#can also isolate individual groups as source or target(ex. target is immune cluster)
#matrices with rows =source and cols = target

count_matrix <- cellchat_test@net$count #returns matrix 
weight_matrix <- cellchat_test@net$weight
#
#
count_dataframe <- reshape2::melt(count_matrix)

cellchat <- cellchat_test

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

netVisual_circle(weight_matrix, targets.use = "N-immune", vertex.weight = groupSize, weight.scale = T, title.name = "Tumor to Immune Signaling Weights", vertex.label.cex = 2)


















#setting database of ligand-receptor interactions
CellChatDB <- CellChatDB.human

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


#skip projecting onto Protein to Protein Interaction Network
#USERS can also skip this step and set raw.use = TRUE in the function computeCommunProb().

#computes probabilility of interactions. Uses triMean method by default which has 25% truncated mean. Can change to lower values
#by setting set type = "truncatedMean" and trim = 0.1
#CellChat can also consider the effect of cell proportion in each cell group in the probability calculation.
#USER can set population.size = TRUE.

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

#The function computeAveExpr can help to check the average expression of signaling genes of interest, 
#e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1).

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
#Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# Note: The inferred intercellular communication network of each ligand-receptor pair is stored in slot 'net' 
#and each signaling pathway is stored in the slot 'netP'

df.net <- subsetCommunication(cellchat)
df.netP <- subsetCommunication(cellchat, slot.name = 'netP')

#The following gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
#df.net <- subsetCommunication(cellchat, sources.use = c(1), targets.use = c(1,5)) 

#The following gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 

cellchat <- computeCommunProbPathway(cellchat)

#calculate aggregated cell-cell communication network
#USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use
cellchat <- aggregateNet(cellchat)
#cellchat <- aggregateNet(cellchat, sources.use = "Neoplastic", targets.use = "Immune cell", remove.isolate = F)

#use Circle plots to visualize aggregate network of cell communications
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


rankNet(cellchat, slot.name = "netP")



#visualizing interactions one cell group at a time
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#pathway specific visual

cellchat_test@netP$pathways #gives the most significant pathways - PTN is number 1. FN1 is also high
pathways.show <- c("EGF") 
#EGF, MHC 1,2, CXCL, CCL, CSF (promotes macrophage proliferation) TNF, TGFb (immune signaling pathways) VEGF, WNT, NOTCH, BMP (differentiation), NRXN, NCAM, NEGR (neuron)
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 

#hierarchgy plot same as circle plot
vertex.receiver = seq(1,4) # a numeric vector. 
#levels(cellchat@idents) gives index for each cell type
netVisual_aggregate(cellchat_test, signaling = cellchat_test@netP$pathways, vertex.receiver = seq(1,4))


ht1 <- netAnalysis_signalingRole_heatmap(cellchat_test, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = cellchat@netP$pathways, color.heatmap = "Reds")
#> Do heatmap based on a single object

#contribution of each single receptor-ligand pair to all signaling intections within a given pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
#visualize this single receptor-ligand pair. 
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair , change number to change pair
#do plot for specific cell types
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)



#show network centrality heatmap and circle map
#Try netvisual heatmap
#Network centrality scores




#comparing signaling between tumor samples - compareInteractions and RankNET














#Create CellChat object from Seurat object

#test <-readRDS("/Users/matthewgruen/Desktop/Graeber Lab /SRP Final Project/Clustered Seurat by sample/clustered_Neftel_SMARTSEQ_124.RDS")
cellchat_124 <- createCellChat(seur_dar_124, group.by = "ident")

#setting database of ligand-receptor interactions
CellChatDB <- CellChatDB.human

CellChatDB.use <- CellChatDB # simply use the default CellChatDB

cellchat_124@DB <- CellChatDB.use

cellchat_124 <- subsetData(cellchat_124) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel

cellchat_124 <- identifyOverExpressedGenes(cellchat_124)
cellchat_124 <- identifyOverExpressedInteractions(cellchat_124)

cellchat_124 <- computeCommunProb(cellchat_124, raw.use = TRUE)

#The function computeAveExpr can help to check the average expression of signaling genes of interest, 
#e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1).

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_124 <- filterCommunication(cellchat_124, min.cells = 5)

#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
#Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# Note: The inferred intercellular communication network of each ligand-receptor pair is stored in slot 'net' 
#and each signaling pathway is stored in the slot 'netP'

df.net_124 <- subsetCommunication(cellchat_124)
df.netP_124 <- subsetCommunication(cellchat_124, slot.name = 'netP')

#The following gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
#df.net <- subsetCommunication(cellchat, sources.use = c(1), targets.use = c(1,5)) 

#The following gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 

cellchat_124 <- computeCommunProbPathway(cellchat_124)

#calculate aggregated cell-cell communication network
#USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use
cellchat_124 <- aggregateNet(cellchat_124)
#cellchat <- aggregateNet(cellchat, sources.use = "Neoplastic", targets.use = "Immune cell", remove.isolate = F)

#use Circle plots to visualize aggregate network of cell communications
groupSize <- as.numeric(table(cellchat_124@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_124@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_124@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#visualizing interactions one cell group at a time
mat <- cellchat_124@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#pathway specific visual

cellchat_124@netP$pathways #gives the most significant pathways - PTN is number 1. FN1 is also high
pathways.show <- c("EGF") 
#EGF, MHC 1,2, CXCL, CCL, CSF (promotes macrophage proliferation) TNF, TGFb (immune signaling pathways) VEGF, WNT, NOTCH, BMP (differentiation), NRXN, NCAM, NEGR (neuron)
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 

#hierarchgy plot same as circle plot
vertex.receiver = seq(1,4) # a numeric vector. 
#levels(cellchat@idents) gives index for each cell type
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

#contribution of each single receptor-ligand pair to all signaling intections within a given pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
#visualize this single receptor-ligand pair. 
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair , change number to change pair
#do plot for specific cell types
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)



#show network centrality heatmap and circle map
#Try netvisual heatmap
#Network centrality scores




#comparing signaling between tumor samples - compareInteractions and RankNET