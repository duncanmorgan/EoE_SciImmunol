setwd("L:/Duncan/eoepaper")


library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(feather)
library(dplyr)

# update: 
# seurat is renormalizing the data (since python data is already normalized)
# i'm going to implement a brute-forcey work around, for now
# analysis-wise, this requires getting rid of poisson scaling

# read in files transferred from python
pyImport = function(name) {
  normdata = read_feather(paste0(name, '.feather')) %>% as.data.frame()
  #rownames(data) = data[,1]
  #data = data[,-1]
  
  print('reading in raw data')
  rawdata = read_feather('cellsAll.feather') %>% as.data.frame()
  rownames(rawdata) = rawdata[,1]
  rawdata = rawdata[,-1]
  rawdata = rawdata[,colnames(rawdata) %in% colnames(normdata)]
  
  
  metadata = read.csv(paste0(name, '_meta.txt'), row.names = 1, stringsAsFactors = FALSE)
  seurat = CreateSeuratObject(rawdata, min.cells = 5)
  seurat@meta.data = metadata
  seurat
}

# standard Seurat processing
pyProcess = function(seurat) {
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
pyToSeurat = function(name) {
  print('reading in data')
  seurat = pyImport(name)
  print('processing data')
  seurat = pyProcess(seurat)
  #print('converting to sparse format')
  #seurat = MakeSparse(seurat)
  seurat
}

# for plotting (in future should just source Andy's default plotting script)
umap_theme = theme_bw() + theme(panel.background = element_blank(),
                                panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(color = 'black'), 
                                panel.border=element_rect(color = 'black', size=1, fill = NA), 
                                text = element_text(family = "sans", size = 16))

# assemble the 500g object
seurat = pyToSeurat('all_raw')
seurat
DimPlot(seurat, 'umap')
seurat = MakeSparse(seurat)
seurat = FindClusters(seurat, dims.use = 1:20, resolution = .4)
#DotPlot(seurat, genes = c('SDC1', 'CPA3', 'KRT5', 'KRT13', 'IGKC', 'IGHM', 'IGHA1', 'PIGR'))

#FeaturePlot(seurat, c('IGHM', 'IGHA1' ,'SDC1', 'IGKC'), reduction.use = 'umap')
#FeaturePlot(seurat, c("S100A8", 'KRT13', 'CPA3', 'TRBC2'), reduction.use = 'umap')
seurat@meta.data$adjustment_pheno = 'Cell'
seurat@meta.data$adjustment_pheno[seurat@meta.data$res.0.4 == 4] = 'Plasma Cell'
seurat@meta.data$adjustment_pheno[seurat@meta.data$res.0.4 %in% c(13, 2, 11)] = "Esophagus Epithelium"


DimPlot(seurat, 'umap', group.by = 'adjustment_pheno')

# rename tissue in metadata
seurat@meta.data$tissue[seurat@meta.data$tissue == 'E'] = 'Esophagus'
seurat@meta.data$tissue[seurat@meta.data$tissue == 'D'] = 'Duodenum'

saveRDS(seurat, 'initial_assembly.RDS')
