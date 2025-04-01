## Libraries ####
library(R.utils)
library(svglite)
library(tidyverse)
library(data.table)
library(readxl)
library(openxlsx)
library(dplyr)
library(Seurat)
library(openxlsx)
library(patchwork)
library(car)
library(sm)
library(reticulate)
library(anndata)
library(harmony)
library(renv)
library(radiant.data)
library(BUSpaRse) #knee plot
library("loupeR")


#STEP 1: Load Seurat Object of single-nuclei data #####
snDRG <- readRDS("C:/Users/axw131330/Documents/DRG_snRNAseq/all_samples.labeled.rds")

#STEP 2: Get Enriched genes per cluster #####

DefaultAssay(snDRG) <- 'SCT'
snDRG <- PrepSCTFindMarkers(snDRG, assay = "SCT")
cluster.markers <- Seurat::FindAllMarkers(snDRG, only.pos = T, logfc.threshold = 1.0, 
                                          min.pct = 0.5, assay = 'SCT', slot = 'data')
cluster.markers <- cluster.markers[which(cluster.markers$p_val_adj<1e-20),]
enrGenes <- cluster.markers$gene
which(duplicated(enrGenes)) #check
enrGenes <- enrGenes[!duplicated(enrGenes)]
which(duplicated(enrGenes)) #check
write.xlsx(cluster.markers, 'markers.xlsx', rowNames=T )

#STEP 3: Get the count matrix for just the sig matrix genes #####
count_matrix <- snDRG[["RNA"]]$counts
count_matrix <- data.frame(count_matrix[which(rownames(count_matrix) %in% enrGenes),])
count_matrix <- rownames_to_column(count_matrix)
write_tsv(count_matrix, "C:/Users/axw131330/Documents/DRG_snRNAseq/sn_matrix.tsv")
write_tsv(data.frame(enrGenes), "C:/Users/axw131330/Documents/DRG_snRNAseq/sn_matrixGenes.tsv")


#STEP 4: Get the cell type annotations per barcode #####
cellstypelist <- data.frame(Idents(snDRG))
cellstypelist <- rownames_to_column(cellstypelist)
colnames(cellstypelist) <- c('cellname','celltype')
cellstypelist$cellname <- sub("-", ".",cellstypelist$cellname)
write_tsv(cellstypelist, "C:/Users/axw131330/Documents/DRG_snRNAseq/sn_labels.tsv", col_names = T)


#STEP 5: Prepare spatial data input #####
spatial_dir <- "C:/Users/axw131330/Documents/DRG_snRNAseq/spatial/"
samples <- list.files(spatial_dir)

for(sample in samples){
  
  print(paste0("sample is ", sample))
  sampName <- sample
  cur_path <- paste0(spatial_dir, sample)
  
  #matrix w/ sig genes
  h5_file <- list.files(path = cur_path, pattern = "h5")#, full.names = TRUE)
  h5_path <- paste0(cur_path, "/", h5_file)
  mat_values <- Read10X_h5(h5_path)
  mat_values <- mat_values[which(rownames(mat_values) %in% enrGenes[,1]),]
  mat_values2 <- data.frame(mat_values)
  mat_values2 <- rownames_to_column(mat_values2)
  write_tsv(mat_values2, paste0(spatial_dir,sampName,"_vis_matrix.tsv"))  

  #Get Barcode Coords
  csv_file <- list.files(path = cur_path, pattern = "csv")#, full.names = TRUE)
  csv_path <- paste0(cur_path, "/", csv_file)
  vis_coords <- read.table(csv_path, sep = ",")
  vis_coords <- vis_coords[,c(1,6,5)]
  colnames(vis_coords) <- c("barcodes", "xcoord", "ycoord")
  vis_coords$barcodes <- sub("-", ".",vis_coords$barcodes)
  print(paste0("at end of spatial pipeline ", sampName))
  write.table(vis_coords, paste0(spatial_dir,sampName,'_vis_coords.tsv'), sep="\t", row.names = F)
  
}

