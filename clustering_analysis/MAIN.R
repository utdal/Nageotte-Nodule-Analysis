
#This is the main script for nageotte nodule gene expression analysis
#INSTRUCTIONS:
#1) Save Spaceranger filtered matrix h5 files from each sample in a directory, which should have no other contents
#   - give the path to this dir to the variable 'inputFolder' below
#2) Save annotated barcode csv files from each sample in a directory, which should have no other contents
#   - this file has 2 columns, the first has barcode IDs and the second has a label for each barcodes (ex. Nageotte Nodules)
#   - the csv file is genereated by LoupeBrowser when barcodes are labeled on the platform
#   - give the path to this dir to the variable 'barcodes' below
#3) Rename the files so that the h5 and csv files for the same sample have the same name (before the file extension)
#   - Lines 81-87 must also be adjusted to reflect the sample names
#4) Set 'project_dir' to the folder where you want output to be saved 
#5) The scripts "LoadFiles.R" and "SCTpipeline.R" should also be present in the 'project_dir' folder


## Libraries ####
library(svglite)
library(tidyverse)
library(data.table) #data table methods
library(readxl) #read in excel files
library(openxlsx) #write excel files
library(dplyr)
library(glmGamPoi) #seurat recommends
library(Seurat)
library(SeuratWrappers)
library(car)
library(sm)
library(harmony)
library(loupeR)
library(enrichR)
library(plotly) #3D umap
library(htmltools) #3D umap
library(viridis)


###### SET UP ######

#Set wd
project_dir <- "" #set path to working dir
setwd(project_dir)

#Session Info
sink('sessioninfo.txt', append = T)
sessionInfo()
sink()
loadedPkgs <- sessionInfo()
loadedPkgs <- names(loadedPkgs$otherPkgs)

#Pipeline fxn scripts
source("SCTpipeline.R")
source("LoadFiles.R")

#Set Input Dir
inputFolder <- "../direct_output" #set path to dir with h5 files
barcodes <- "../nageotte_barcodes" #set path to dir with csv files

#Set Working Dir
setwd(project_dir)



###### LOAD FILES - ######

seurObjs <- loadFiles( inputFolder, barcodes, 'nageotte_nodules' )
ngn <- seurObjs[[1]]
neu <- seurObjs[[2]]


###### GET/LABEL HIGH LEVEL CLUSTERS - ######

#Seurat SCT Pipeline
ngn <- runSCT( ngn )
feat_plots(ngn, 'ngn.STMneurons', 'umap.harmony', 'SCT', 'data', markers.STMneurons )

#Backup data
saveRDS(ngn, "ngn.rds")



###### MORE METADATA ######

ngn$barcodes <- colnames(ngn)
ngn$donor <- ifelse((ngn$orig.ident == 'ngn1' | ngn$orig.ident == 'ngn2'), 'donor14', 
                    ifelse((ngn$orig.ident == 'ngn3' | ngn$orig.ident == 'ngn4'), 'donor1',
                           ifelse((ngn$orig.ident == 'ngn5' | ngn$orig.ident == 'ngn6'), 'donor2',
                                  ifelse((ngn$orig.ident == 'ngn7' | ngn$orig.ident == 'ngn8'), 'donor2',
                                         ifelse((ngn$orig.ident == 'ngn9' | ngn$orig.ident == 'ngn10' | ngn$orig.ident == 'ngn11'), 'donor3',
                                                ifelse((ngn$orig.ident == 'ngn12' | ngn$orig.ident == 'ngn13' | ngn$orig.ident == 'ngn14'), 'donor6', 'donor4'))))))
sampleKey <- data.frame(sampleID = c('14.1.1','14.1.2','1.1.1','1.1.2','2.1.1',
                                     '2.1.2','2.2.1','2.2.2','3.1.1','3.1.2',
                                     '3.1.3','6.1.1','6.1.2','6.1.3','4.1.1','4.1.2'))
sampleKey$sample <- paste0('ngn', rownames(sampleKey))
sampleKey$batch <- c(rep(1,times=2), rep(2,times=6), rep(3,times=8))
sampleKey$slide <- c(4,3,2,1,1,2,1,2,5,8,7,5,6,7,6,8)
ngn@meta.data <- left_join(ngn@meta.data, sampleKey, by='sample')
rownames(ngn@meta.data) <- ngn$barcodes
p1 <- DimPlot(ngn, reduction = 'umap.harmony')
p2 <- DimPlot(ngn, reduction = 'umap.harmony', group.by = 'slide')
p3 <- DimPlot(ngn, reduction = 'umap.harmony', group.by = 'donor')
p1 + p2 + p3



###### GENE EXPRESSION TABLES FOR ALL GENES (SEP) - ######

#Get gene expression table (in CPM) per ngn cluster
Idents(ngn) <- 'seurat_clusters'
geneExpression <- AggregateExpression(ngn,return.seurat=T, normalization.method = 'RC', assays = 'RNA')
geneExpression <- data.frame(geneExpression[['RNA']]$data)
geneExpression <- geneExpression*100
colSums(geneExpression)

#Get gene expression table (in CPM) per all ngn
ngn$bcType <- 'ngn'
Idents(ngn) <- 'bcType'
geneExpression2 <- data.frame(AggregateExpression(ngn,return.seurat=F, assays = 'RNA'))
colnames(geneExpression2)[1] <- 'all'
geneExpression2$all <- 1000000 * geneExpression2$all / sum(geneExpression2$all)
sum(geneExpression2$all)
geneExpression$all <- geneExpression2$all
colnames(geneExpression) <- paste0( colnames(geneExpression), '.ngn')

#Get gene expression table (in CPM) per all neu
neu$bcType <- 'neu'
Idents(neu) <- 'bcType'
geneExpression4 <- data.frame(AggregateExpression(neu,return.seurat=F, assays = 'SCT'))
geneExpression4$all <- 1000000 * geneExpression4$all / sum(geneExpression4$all)
sum(geneExpression4$all)
geneExpression3$all <- geneExpression4$all
colnames(geneExpression4) <- paste0( colnames(geneExpression4), '.neu')

#Write
geneExpression$Gene <- rownames(geneExpression)
geneExpression4$Gene <- rownames(geneExpression4)
geneExpression <- left_join( geneExpression, geneExpression4, by = 'Gene')
write.xlsx( geneExpression, 'ngn_cpm.xlsx', rowNames=T)



###### CYTOKINES BAR PLOT ######

#Cytokines list
cytokines <- c('SPP1',     'CCL3',     'CCL4',     'CCL4L2',   'BMP8B',    'IL34',     'TNFSF13',  'CCL14',   
               'CCL2',     'CSF1',     'TNFSF12',  'IK',       'TGFB1',    'BMP7',     'IL18',     'GDF10',   
               'CXCL16',   'GDF11',    'TNFSF10',  'IL16',     'FLT3LG',   'CX3CL1',   'IL1B',     'WNT5B',   
               'IL33',     'WNT6',     'CLCF1',    'KITLG',    'EDA',      'TNFSF13B', 'BMP3',     'IL32',    
               'IL15',     'IL18BP',   'IL17B',    'TNFSF9',   'CCL28',    'WNT5A',    'IL1RN',    'CYTL1',   
               'OSM',      'IL17D',    'WNT10B',   'BMP2',     'CCL5',     'CNTF',     'BMP5',     'CXCL10',  
               'BMP6',     'WNT2B',    'WNT11',    'INHBA',    'CCL13',    'BMP8A',    'WNT7B',    'INHA',    
               'TNFSF14',  'TGFB3',    'GDF9',     'CCL17',    'INHBB',    'CCL27',    'GDF5',     'CXCL13',  
               'CXCL8',    'CCL8',     'CXCL9',    'CCL19',    'IL6',      'BMP4',     'WNT2',     'CCL23',   
               'TNF',      'WNT3',     'IL12A',    'TNFSF15',  'WNT9A',    'IL24',     'PF4',      'LTB',     
               'CXCL11',   'WNT4',     'CCL15',    'XCL1',     'IL23A',    'IL10',     'TNFSF18',  'TNFSF8',  
               'IL1A',     'CCL1',     'CXCL2',    'CCL22',    'CXCL1',    'IL27',     'IL7',      'LEFTY1',  
               'TNFSF4',   'WNT9B',    'CCL21',    'WNT7A',    'IL11',     'GDF3',     'LIF',      'LEP',     
               'CCL16',    'IL5',      'MSTN',     'CCL25',    'GDF15',    'CXCL6',    'IL25',     'LTA',     
               'PF4V1',    'CD70',     'CCL11',    'CCL18',    'CCL20',    'CCL26',    'CCL7',     'CXCL3',   
               'CXCL5',    'EPO',      'GDF6',     'GDF7',     'IFNA17',   'IL37',     'LEFTY2',   'NODAL',   
               'SCGB1A1',  'TGFB2',    'WNT10A',   'WNT16',    'XCL2',     'IL17C',    'IL17F',    'IL26',    
               'CD40LG',   'FASLG',    'IFNL3' )
  
#Get expression of cytokines in ngn and top 20 cytokine genes
#Cytokines identified by PANTHER protein class of gene product (those labeled as 'cytokines' and 'interleukin)
cytoNGN <- ngn[ which(rownames(ngn) %in% cytokines$symbol), ]
cytoNGN$type <- 'nodule'
Idents(cytoNGN) <- 'type'
cytokineExpression <- data.frame(AggregateExpression(cytoNGN,return.seurat=F, assays = 'RNA'))
cytokineExpression <- rownames_to_column(cytokineExpression)
cytokineExpression <- cytokineExpression[ order(cytokineExpression$RNA, decreasing = T), ]
cytokinesTop <- cytokineExpression$rowname[1:20]
write.xlsx( cytokineExpression, 'cytokineExpression.xlsx')

#Set Gene List for plotting
genes <- cytokinesTop

#Prep for Plot
ngn2 <- ngn
DefaultAssay(ngn2) <- 'RNA'
ngn2 <- NormalizeData( ngn2, normalization.method = 'RC', scale.factor = 1e6)
ngnData <- ngn2[['RNA']]$data
head(colSums(ngn2))
ngnSubset <- data.frame(ngnData[ which( rownames(ngnData) %in% genes), ])

#Reorder by total expression of genes
ngnSubset$TOTALS <- rowSums(ngnSubset)
ngnSubset <- ngnSubset[ order( ngnSubset$TOTALS ) , ]
genes <- rownames(ngnSubset) 
ngnSubset <- ngnSubset[ , -c(length(colnames(ngnSubset)))]

#format data for plot
ngnSubset <- t(ngnSubset)
ngnPlot <- data.frame( values = 0 , gene = '')
ngnPlot$barcode <- ''
ngnPlot <- ngnPlot[-c(1),]
for( gene in 1:length(colnames(ngnSubset))){
  ngnGene <- data.frame(values = ngnSubset[ , gene], gene = colnames(ngnSubset)[gene])
  ngnGene$barcode <- rownames(ngnGene)
  ngnPlot <- rbind(ngnPlot, ngnGene)
}
ngnPlot$values <- log10(ngnPlot$values) #log transform
ngnPlot$values <- ifelse(ngnPlot$values == '-Inf', 0, ngnPlot$values) #log transform

#extra steps for cytokine plot
ngnClus <- data.frame( barcode = colnames(ngn), type = ngn$seurat_clusters)
ngnClus$barcode <- str_replace_all( ngnClus$barcode, '-', '.')
ngnPlot <- left_join( ngnPlot, ngnClus, by = 'barcode')

#barplot
numGenes <- 20
colors <- viridis(numGenes, alpha =  0.5, option = 'B')
colors <- data.frame( colors = colors, genes = genes)
colors <- colors[ order(colors$genes), ]
ngnPlot$ptColors <- ifelse( ngnPlot$type == '0', 'blue', 
                    ifelse( ngnPlot$type == '1', 'green',
                            ifelse( ngnPlot$type == '2', 'magenta',
                            ifelse( ngnPlot$type == '3', 'red', 'cyan') ) ) )
ptColors <- ngnPlot$ptColors
boxPlotNgn <- ggplot(ngnPlot, aes(x = factor(gene, level=genes), y = values, fill = gene)) + 
  stat_summary(geom = "col", fun = mean) + ylim(-0.1, 4.1)+ 
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) + coord_flip() +
  theme(legend.position="none") +
  scale_fill_viridis(discrete = TRUE, alpha=0.5, option = 'B') +
  scale_fill_manual(values=colors$colors) + 
  geom_jitter(position=position_jitter(0.25, 0.05), color= 'black', fill = 'black' , size=0.5, alpha=0.3, pch=21) +
  labs(x = "", y = "Log10(CPM)", title = "Top 20 Cytokine Genes in the Nodule Barcodes",
       subtitle = "Mean expression across all Nodule barcodes") +
  theme(axis.text.y = element_text(face = "bold.italic", size = 14), 
        axis.text.x = element_text(size = 14), 
        axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 18, face = 'bold'),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle=element_text(size=16, hjust=0.5, face="bold.italic", color="black"))
svg( 'cytokines.svg', 8, 8)
print(boxPlotNgn)
dev.off()

