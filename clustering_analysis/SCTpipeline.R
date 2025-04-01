#### SCT Pipeline ####

#This script runs Seurat's SCT Pipeline with Harmony Integration
#INPUT: a Seurat Object with only RNA assay and counts slot,
#with each sample in a separate layer in counts slot

#following sctransform tutorial (pub. 10/31/2023): https://satijalab.org/seurat/articles/sctransform_vignette



###### DEFINE PLOTTING FUNCTIONS ######

qc_violinsRNA <- function( seuratobj, name ) {
  
  #with dots
  qcplotname <- paste0("QCViolinRNA_",name,".png")
  violin <- VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.3)
  png(qcplotname, width = 16, height = 8, units = "in", res = 180)
  print(violin)
  dev.off()
  
  #w/out dots
  qcplotname <- paste0("QCViolinRNANODOTS_", name,".svg")
  violin <- VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
  svg(qcplotname, width = 16, height = 8)
  print(violin)
  dev.off()
}

qc_violinsSCT <- function( seuratobj, name ) {
  
  #with dots
  qcplotname <- paste0("QCViolinSCT_",name,".png")
  violin <- VlnPlot(seuratobj, features = c("nFeature_SCT", "nCount_SCT", "percent.mt"), ncol = 3, pt.size = 0.3)
  png(qcplotname, width = 16, height = 8, units = "in", res = 180)
  print(violin)
  dev.off()
  
  #w/out dots
  qcplotname <- paste0("QCViolinSCTNODOTS_", name,".svg")
  violin <- VlnPlot(seuratobj, features = c("nFeature_SCT", "nCount_SCT", "percent.mt"), ncol = 3, pt.size = 0)
  svg(qcplotname, width = 16, height = 8)
  print(violin)
  dev.off()
}

qc_featRNA <- function( seuratobj, name, reduc ) {
  
  prev_assay <- DefaultAssay(seuratobj)
  
  #RNA Assay
  DefaultAssay(seuratobj) <- "RNA"
  featplotname <- paste0("QCFeatRNA_",name,".svg")
  plotheight <- ceiling(3/2)*10
  svg(featplotname, width = 20, height = plotheight)
  print(FeaturePlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    reduction = reduc, pt.size = 0.1, ncol = 2,
                    cols = c("lightyellow","firebrick4"),label = T, raster = F))
  dev.off()
  
  DefaultAssay(seuratobj) <- prev_assay
}

qc_featSCT <- function( seuratobj, name, reduc ) {
  
  prev_assay <- DefaultAssay(seuratobj)
  
  #SCT Assay
  DefaultAssay(seuratobj) <- "SCT"
  featplotname <- paste0("QCFeatSCT_",name,".svg")
  plotheight <- ceiling(3/2)*10
  svg(featplotname, width = 20, height = plotheight)
  print(FeaturePlot(seuratobj, features = c("nFeature_SCT", "nCount_SCT", "percent.mt"),
                    reduction = reduc, pt.size = 0.2, ncol = 2,
                    cols = c("lightyellow","firebrick4"),label = T, raster = F))
  dev.off()
  
  DefaultAssay(seuratobj) <- prev_assay
}

umaps <- function( seuratobj, name, reduc ) {
  
  #UMAP plot w/ active ident
  umapplotname <- paste0("UMAP_", name,".svg")
  umap_plot <- DimPlot(
    seuratobj, label = F, pt.size = 0.1,
    reduction = reduc, raster = F,
    combine = FALSE, label.size = 4)
  svg(umapplotname, width = 10, height = 10)
  print(umap_plot)
  dev.off()
  
  #UMAP plot w/ active ident + LABEL
  umapplotname <- paste0("UMAPLABEL_", name,".svg")
  umap_plot <- DimPlot(
    seuratobj, label = T, pt.size = 0.1,
    reduction = reduc, raster = F,
    combine = FALSE, label.size = 4)
  svg(umapplotname, width = 10, height = 10)
  print(umap_plot)
  dev.off()
  
  #UMAP plot w/ sample
  umapplotname <- paste0("UMAPsamp_", name,".svg")
  umap_plot <- DimPlot(
    seuratobj, label = F, group.by = "sample",
    reduction = reduc, raster = F, pt.size = 0.1,
    combine = FALSE)
  svg(umapplotname, width = 10, height = 10)
  print(umap_plot)
  dev.off()
  
  #UMAP plot w/ active ident but split by sample
  umapplotname <- paste0("UMAPsplit_", name,".svg")
  umap_plot <- DimPlot(
    seuratobj, label = F, pt.size = 0.1,
    reduction = reduc, raster = F,
    split.by = "sample",
    combine = FALSE, label.size = 2)
  svg(umapplotname, width = 65, height = 5)
  print(umap_plot)
  dev.off()
}


umap3d <- function( seuratobj, name, metaSlot ) {
  
  #3D UMAP plot w/ ident specified in metaSlot
  #https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0
  seuratobj$clus3d <- seuratobj[[metaSlot]]
  head(Embeddings(object = seuratobj, reduction = "umap.3d"))
  UMAP_1 <- seuratobj[["umap.3d"]]@cell.embeddings[,1]
  UMAP_2 <- seuratobj[["umap.3d"]]@cell.embeddings[,2]
  UMAP_3 <- seuratobj[["umap.3d"]]@cell.embeddings[,3]
  plot.data <- FetchData(object = seuratobj, vars = c("umap3d_1", "umap3d_2", "umap3d_3", "clus3d"))
  plot.data$label <- paste(rownames(plot.data))
  fig <- plot_ly(data = plot.data, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                 color = ~clus3d, 
                 colors = c("red","orange", "gold", "limegreen", 'darkgreen', 'turquoise',
                            'skyblue', 'blue','purple', 'orchid1', 'magenta'),
                 type = "scatter3d", mode = "markers", marker = list(size = 5, width=2), # controls size of points
                 text=~label, #This is that extra column we made earlier for which we will use for cell ID
                 hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
  axx <- list(nticks = 5, range = c(-10,10)) # Before you plot, set the ranges of the axis you desire. 
  axy <- list(nticks = 5, range = c(-10,10))
  axz <- list(nticks = 5, range = c(-10,10))
  #fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
  save_html(fig_cube, paste0("umap3D_", name,".html"))
  
}


feat_plots <- function( seuratobj, name, reduc, seurassay, seurslot, feats ) {
  
  prev_assay <- DefaultAssay(seuratobj)
  
  #Feature
  DefaultAssay(seuratobj) <- seurassay
  featplotname <- paste0("Feat",seurassay,seurslot,"_",name,".png")
  plotheight <- ceiling(length(feats)/2)*10
  png(featplotname, width = 20, height = plotheight, units = "in", res = 180)
  print(FeaturePlot(seuratobj, features = feats, cols = c("lightyellow","firebrick4"),
                    reduction = reduc, pt.size = 0.1, ncol = 2, order = T,
                    label = T, raster = F, slot = seurslot))
  dev.off()
  
  #Unscaled Dot
  dotplotname <- paste0("Dot",seurassay,"data_",name,".svg")
  plotheight <- 2+length(unique(Idents(seuratobj)))/3
  plotwidth <- 2+length(feats)/3
  svg(dotplotname, width = plotwidth, height = plotheight)
  print(DotPlot(seuratobj, features = feats, assay = seurassay, scale = F) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic")))
  dev.off()
  
  #Scaled Dot
  dotplotname <- paste0("DotScaled",seurassay,"data_",name,".svg")
  plotheight <- 2+length(unique(Idents(seuratobj)))/3
  plotwidth <- 2+length(feats)/3
  svg(dotplotname, width = plotwidth, height = plotheight)
  print(DotPlot(seuratobj, features = feats, assay = seurassay, scale = T) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic")))
  dev.off()
  
  #Violin No Dots
  vlnplotname <- paste0("VlnNoDots",seurassay,seurslot,"_",name,".svg")
  plotheight <- 2+ceiling(length(feats)/4)*4
  plotwidth <- 2+length(unique(Idents(seuratobj)))
  svg(vlnplotname, width = plotwidth, height = plotheight)
  print(VlnPlot(seuratobj, features = feats, assay = seurassay, 
                layer = seurslot, pt.size = 0, ncol=4))
  dev.off()
  
  #Violin w/ Dots
  vlnplotname <- paste0("Vln",seurassay,seurslot,"_",name,".png")
  plotheight <- 2+ceiling(length(feats)/4)*4
  plotwidth <- 2+length(unique(Idents(seuratobj)))
  png(vlnplotname, width = plotwidth, height = plotheight, units = "in", res = 180)
  print(VlnPlot(seuratobj, features = feats, assay = seurassay, 
                layer = seurslot, pt.size = 0.5, ncol=4))
  dev.off()
  
  DefaultAssay(seuratobj) <- prev_assay
}

markers.gen <- c("SNAP25", "PRPH", "PIRT", "RBFOX3", "TMSB10", "TMSB4X", 
                 "PKD1", "PDIA2", "MRGPRE", "NTRK1", "NTRK3", "HTR3A", "CHRNA3",
                 "CCK", "GFRA2", "MRGPRX1", "FABP7", "GJA1", 
                 "SCN7A", "NCAM1","MBP", "POU3F1", "NKX2-2", "PI16", 
                 "ANPEP", "PDGFRA", "ACTA2", "DES", "TEK", "VWF",
                 "FABP4", "ADIPOQ", "PLCB2", "TRAC", "CD34", "IL6", "OSM",
                 "CDKN1A", "CDKN2A", "PECAM1", "CLDN5" )

markers.STMneurons <- c('PVALB','KCNS1','NTRK3','SFRP1','NTRK2',
                        'ANKRD37','SCN10A','TRPM8','FAM49A','CPNE6',
                        'NTRK1', 'NGFR','PENK','GZMB','TRPA1',
                        'GAL','CHRNA3', 'CCK','GFRA2','NPPB','CSF3R',
                        'MRGPRF', 'TRPV1', 'HTR3A')

markers.neural <- c("SNAP25", "PRPH", "PIRT", "RBFOX3", "HEPACAM",
                    "FABP7", "SOX10", "GJA1", "KCNJ10", "NCAM1", 
                    "SCN7A", "MPZ", "MBP", "EGR2", "DHH", "POU3F1", "NKX2-2",
                    "TMSB10", "TMSB4X")





###### DEFINE SCT PIPELINE FXN ###### 

runSCT <- function( inSeur, numPC=30, seurRes=0.8 ) {
  
  #Create dir for output
  outDir <- getwd()
  if(!dir.exists("SCT_Harmony")){dir.create("SCT_Harmony")}
  setwd("SCT_Harmony")
  
  #Processing pipeline
  print( c( 'Starting SCT Pipeline with Harmony Integration', format(Sys.time(), "%H:%M:%S")))
  inSeur <- SCTransform( inSeur ) %>% 
    RunPCA(npcs = numPC) %>%
    IntegrateLayers(method = HarmonyIntegration, normalization.method = 'SCT', 
                    orig.reduction = 'pca', new.reduction = 'harmony', verbose = TRUE)  %>%
    FindNeighbors(reduction = 'harmony', dims = 1:numPC) %>%
    FindClusters(cluster.name = 'harmony_clusters', resolution = seurRes) %>%
    RunUMAP(reduction = 'harmony', dims = 1:numPC, reduction.name = 'umap.3d', n.components = 3L)  %>%
    RunUMAP(reduction = 'harmony', dims = 1:numPC, reduction.name = 'umap.harmony')
  inSeur[["RNA"]]  <- JoinLayers(inSeur[["RNA"]]) 
  inSeur <- PrepSCTFindMarkers(inSeur)
  
  #Plots
  Idents(inSeur) <- "sample"
  qc_violinsRNA(inSeur, "samp_postHarmony")
  qc_violinsSCT(inSeur, "samp_postHarmony")
  Idents(inSeur) <- "seurat_clusters"
  umaps(inSeur, "postHarmony", "umap.harmony" )
  umap3d(inSeur, "postHarmony", "seurat_clusters" )
  feat_plots(inSeur, "postHarmony", "umap.harmony", "SCT", "data", markers.gen )
  qc_violinsRNA(inSeur, "clus_postHarmony")
  qc_violinsSCT(inSeur, "clus_postHarmony")
  qc_featSCT(inSeur, "clus_postHarmony", "umap.harmony")
  qc_featRNA(inSeur, "clus_postHarmony", "umap.harmony")
  write.xlsx(data.table(table(Idents(inSeur))), "ClusterCounts.xlsx")
  
  # #Save a Loupe object
  # copySeur <- inSeur
  # copySeur$origBarcode <- NULL
  # create_loupe_from_seurat( copySeur )
  
  #Find markers for each cluster
  snDRG.markers <- FindAllMarkers(inSeur, only.pos = TRUE)
  snDRG.markers %>%  group_by(cluster)
  write.xlsx(snDRG.markers, "snDRG_FindAllMarkers.xlsx")
  
  #Enrichr on clusters
  require(enrichR) # https://github.com/wjawaid/enrichR
  cellTypes <- c('')
  dbs <- c("PanglaoDB_Augmented_2021", "CellMarker_2024")
  
  for( clus in 0:length(unique(inSeur$seurat_clusters))) {
    
    clusMarkers <- snDRG.markers[ which(snDRG.markers$cluster == clus), ]
    clusMarkers <- clusMarkers[ which(snDRG.markers$p_val_adj < 0.01), ]
    clusMarkers <- clusMarkers[ which(snDRG.markers$avg_log2FC > 1), ]
    clusMarkers <- clusMarkers[ which(snDRG.markers$pct.1 > 0.5), ]
    clusMarkers <- clusMarkers$gene
    clusMarkers <- clusMarkers[!is.na(clusMarkers)]
    
    sink('clusMarkers.txt', append=T)
    cat('\n\nMarkers for cluster', clus)
    cat('\nFindAllMarkers output that has padj<0.01, LFC>1, pct>0.5\n')
    print(clusMarkers)
    sink()
    
    enriched <- enrichr( clusMarkers, dbs )
    pang <- enriched[['PanglaoDB_Augmented_2021']]
    pang <- pang$Term[1]
    cm <- enriched[['CellMarker_2024']]
    cm <- cm$Term[1]
    
    sink('clusMarkers.txt', append=T)
    cat('Predicted cell type with PanglaoDB: ', pang)
    cat('\nPredicted cell type with CellMarker: ', cm)
    sink()
    
    cellTypes <- append(cellTypes, pang)
    
  }
  
  #Print File with Predicted cell types
  cellTypes <- cellTypes[-1]
  cellTypes <- data.frame( cellTypes = cellTypes )
  cellTypes$clusNum <- rownames(cellTypes)
  cellTypes$clusNum <- as.integer(cellTypes$clusNum) - 1
  write.xlsx(cellTypes, "predictedCellTypes.xlsx")
  write.xlsx(cellTypes, "finalCellTypes.xlsx")
  
  print( c( 'Finished SCT Pipeline with Harmony Integration', format(Sys.time(), "%H:%M:%S")))
  setwd(outDir)
  return(inSeur)
  
}