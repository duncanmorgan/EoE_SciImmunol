# these are general purpose functions. source this script at the beginning of every session
setwd("L:/Duncan/eoepaper_final")

reticulate::use_python('C:/Users/dmorgan/AppData/Local/Continuum/anaconda3/python.exe', required = TRUE)

library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(feather)
library(dplyr)
library(reshape2)
library(viridis)
library(tidyr)
library(pheatmap)

# import a dataset exported from python with the seuratExport function
pyImport = function(name, rawfile = 'cellsAll.feather') {
  
  # load normalized data
  normdata = read_feather(paste0(name, '.feather')) %>% as.data.frame()
  
  # load raw data 
  print('reading in raw data')
  rawdata = read_feather(rawfile) %>% as.data.frame()
    
  # transfer genes to rownames and drop columns  
  rownames(rawdata) = rawdata[,1]
  rawdata = rawdata[,-1]
    
  # subset raw data
  rawdata = rawdata[,colnames(rawdata) %in% colnames(normdata)]
  
  # read metadata
  metadata = read.csv(paste0(name, '_meta.txt'), row.names = 1, stringsAsFactors = FALSE)
  
  # create and return seurat object
  seurat = CreateSeuratObject(rawdata, min.cells = 5)
  seurat@meta.data = metadata
  seurat
}

# standard Seurat processing
seuratProcess = function(seurat) {
  seurat = NormalizeData(seurat) 
  seurat = FindVariableGenes(seurat, do.plot= FALSE)
  seurat = ScaleData(seurat, genes.use =seurat@var.genes, vars.to.regress = c('n_genes'), model.use = 'poisson')
  seurat = RunPCA(seurat, dims.use = seurat@var.genes, do.print = FALSE)
  seurat = RunUMAP(seurat, dims.use = 1:20)
  seurat@meta.data$UMAP1 = seurat@dr$umap@cell.embeddings[,1]
  seurat@meta.data$UMAP2 = seurat@dr$umap@cell.embeddings[,2]
  seurat
}

# take the files from the exportSeurat function, assemble them into a seurat object, and process
pyToSeurat = function(name, rawfile = 'cellsAll.feather') {
  print('reading in data')
  seurat = pyImport(name, rawfile)
  print('processing data')
  seurat = seuratProcess(seurat)
  print('converting to sparse format')
  seurat = MakeSparse(seurat)
  seurat
}

# add UMAP coordinates to seurat@meta.data
addUMAP = function(seurat) {
    seurat@meta.data$UMAP1 = seurat@dr$umap@cell.embeddings[,1]
    seurat@meta.data$UMAP2 = seurat@dr$umap@cell.embeddings[,2]
    seurat
} 

# randomly shuffle the rows in a dataframe (used primarily for plotting)
shuffle = function(data) {
    set.seed(1)
    data[sample(rownames(data), length(rownames(data))),]
}

# create FeaturePlots using the non-default Seurat color scheme
geneplot= function(seurat, genes) {
    plots = c()
    for (curr in genes){
        seurat@meta.data$gene = seurat@data[curr, rownames(seurat@meta.data)]
        plots[[curr]] =  ggplot(shuffle(seurat@meta.data), aes(x = UMAP1, y = UMAP2, color = gene)) + geom_point(size = .8) +
          scale_color_viridis_c() + labs(title = curr) + guides(color = FALSE) + theme(axis.title = element_blank(), axis.text = element_blank()) + remove_grid
    }
    gg = plot_grid(plotlist = plots)
    gg
}
set.seed(1)

pct = function(x) {
  sum(x >0)/length(x)
}
meanexp = function(x) {
    mean(x)
}


# import plotting elements
source('figure_parameters.R')