# a quick script to read in data from Python and produce Seurat objects
source('../../functions.R')

seurat = pyToSeurat('Preprocessing/Blood/sort_raw', 'Preprocessing/Blood/sort_cellsAll.feather')

seurat = RunUMAP(seurat, dims.use = 1:10)
DimPlot(seurat, 'umap')

# remove one patient without EoE (processed as control) and cells with under 900 genes
seurat2 = SubsetData(seurat, seurat@cell.names[seurat@meta.data$Patient != "PNOIT" & seurat@meta.data$n_genes > 900])
seurat2 = seuratProcess(seurat2)

seurat2 = RunUMAP(seurat2, dims.use = 1:10)
DimPlot(seurat2, 'umap', group.by = 'Fraction')
seurat2 = FindClusters(seurat2, dims.use = 1:10, resolution = .15, print.output = FALSE, force.recalc = TRUE)
saveRDS(seurat2, 'Data/pb_no_tcr.RDS')