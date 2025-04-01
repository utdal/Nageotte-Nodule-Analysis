

###### DEFINE loadFiles FXN ###### 

loadFiles <- function( inDir, barcodeDir, projName ) {

  #Get lıst of fıles that need to be processed with SoupX
  samplesFolder <- inDir
  samples <- list.files(samplesFolder)
  
  #Make empty list to add seurat objects
  ngn <- list()
  neu <- list()
  
  #Run for each sample
  for(sampNum in 1:length(samples)){
    
    #Setup for sample output dir
    sampName <- samples[sampNum]
    sampName <- substr(sampName, 1, (nchar(sampName) - 3) )
    print(paste('Running sample:', sampName))
    
    #read in cellranger output and barcodes
    filtMat <- Read10X_h5(paste0(samplesFolder, "/", sampName, '.h5'))
    barcodes <- read.table(paste0(barcodeDir, "/", sampName, '.csv'), sep = ',')
    ngnBC <- barcodes$V1[ which(barcodes$V2 == 'Nageotte Nodule' | 
                                  barcodes$V2 == 'Nageotte Nodules')]
    neuBC <- barcodes$V1[ which(barcodes$V2 == 'Multiple neurons near nageotte' | 
                                  barcodes$V2 == 'Multiple neurons near Nageotte' |
                                  barcodes$V2 == 'Multiple Neurons near Nageotte' |
                                  barcodes$V2 == 'Single neuron near nageotte' |
                                  barcodes$V2 == 'Single neuron near Nageotte' |
                                  barcodes$V2 == 'Single neuron near nageottes' |
                                  barcodes$V2 == 'Single Neuron near Nageotte' |
                                  barcodes$V2 == 'Single  neuron near nageotte' )]
    
    #Make ngn obj
    ngn_counts <- filtMat[ , which( colnames(filtMat) %in% ngnBC)]
    srat <- CreateSeuratObject(counts = ngn_counts, project = sampName)
    srat$sample <- sampName
    srat$origBarcode <- colnames(srat)
    srat[["percent.mt"]] <- PercentageFeatureSet(srat, assay="RNA", pattern = "^MT-")
    qc_violinsRNA(srat, "ngn")
    ngn <- append( ngn, srat)
    names(ngn)[sampNum] <- sampName
    
    #Make neu obj
    neu_counts <- filtMat[ , which( colnames(filtMat) %in% neuBC)]
    srat <- CreateSeuratObject(counts = neu_counts, project = sampName)
    srat$sample <- sampName
    srat$origBarcode <- colnames(srat)
    srat[["percent.mt"]] <- PercentageFeatureSet(srat, assay="RNA", pattern = "^MT-")
    qc_violinsRNA(srat, "ngn")
    neu <- append( neu, srat)
    names(neu)[sampNum] <- sampName

  }
  
  #Prepare output
  ngnObj <- first(ngn)
  if( length(samples) > 1 ) {
    ngn[1] <- NULL
    ngnObj <- merge( x = ngnObj, y = ngn )
  }
  
  neuObj <- first(neu)
  if( length(samples) > 1 ) {
    neu[1] <- NULL
    neuObj <- merge( x = neuObj, y = neu )
  }

  return( list(ngnObj, neuObj) )
  
}

