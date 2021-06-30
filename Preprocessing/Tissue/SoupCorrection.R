source("L:/Duncan/eoepaper/functions.R")
# Read in data 
# seurat is the seurat containing seurat raw data
seurat = readRDS('Data/initial_assembly.RDS')

library(SoupX)

seurat = SetAllIdent(seurat, 'orig.ident')
input.file = read.table("Preprocessing/input_file.txt", header = 1, sep = '\t', 
                        stringsAsFactors = FALSE, row.names = 1)

# for each individual seqwell chip
for(id in levels(seurat@ident)) {
  message(id)
  
  #subset data
  curr = SubsetData(seurat, WhichCells(seurat, id), subset.raw = TRUE)
  
  # extract counts, read in soup profile  
  mat = as.data.frame(as.matrix(curr@raw.data))
  profile = read.table(input.file[as.character(id), 'Soup'], stringsAsFactors = FALSE, skip = 1, row.names = 1)
  profile$gene = rownames(profile)
  
  # some unusual manipulation to get the "soup" profile incorporated into the soupX object 
  # could be avoided if you have a dge containing all cell barcodes
  profile = profile[profile$gene %in% rownames(mat),]
  mat = cbind(mat, rep(0, dim(mat)[1]))
  mat[rownames(profile), dim(mat)[2]] = profile[,1]
  profile$rat = profile$V2/sum(profile$V2)
  profile = profile[order(profile$V2, decreasing = TRUE),]
  print(head(profile, 20))
  
  # create soup, estimate nonexpressed genes for reference
  currsoup = SoupChannel(tod = mat, toc  = curr@raw.data, metadata = NULL,soupRange = c(max(curr@meta.data$n_counts) + 1, 1e99), keepDroplets =             TRUE)
  currsoup = estimateSoup(currsoup, soupRange = c(max(curr@meta.data$n_counts) + 1, 1e99))  
  
  # here is where to exclude cells from soup estimation
  # in the seuratphagus, i excluded epithelial cells (since these express keratin)
  # in the duodenum, I excluded plasma cells
  tissue = curr@meta.data$tissue[1]
  if (tissue == 'Esophagus'){
    boolmat = !curr@meta.data$adjustment_pheno %in% c("Esophagus Epithelium")
    boolmat = as.matrix(boolmat)
    
    # name of list has to match name of genes vector
    rownames(boolmat) =  colnames(currsoup$toc)
    genes.use = c("KRT4", "KRT5", "KRT13", 'KRT15', 'KRT6A')
    
  }
    
  if (tissue == 'Duodenum'){
    boolmat = !curr@meta.data$adjustment_pheno %in% c("Plasma Cell")
    boolmat = as.matrix(boolmat)
    
    # name of list has to match name of genes vector
    rownames(boolmat) =  colnames(currsoup$toc)
    genes.use = c("IGHA1", 'IGJ', "IGKC", "IGHA2", 'IGHM', 'IGKC')
    
  }
  
  # calculate contamination and create plot
  geneList = list()
  geneList[['genes.use']] = genes.use
  print(genes.use)
  currsoup = calculateContaminationFraction(currsoup,
                                            nonExpressedGeneList = geneList,
                                            useToEst = boolmat, 
                                            verbose = TRUE)
  
  print(head(currsoup$metaData, 10))
  
  # adjust data
  currsoup_mat= SoupX::adjustCounts(currsoup, method = 'subtraction')
  
  # change seurat@raw.data to new raw count data
  seurat@raw.data[rownames(curr@data), curr@cell.names] = as.matrix(ceiling(currsoup_mat))
  
   gOrig = apply(curr@raw.data, 1, sum)
   gNew = apply(ceiling(currsoup_mat), 1, sum) 
   rats = (gOrig - gNew)/gOrig
   rats = rats[!is.na(rats)]
   rats = rats[order(rats, decreasing = TRUE)]
   print(rats[1:50])
   message(sum(ceiling(currsoup_mat)) / sum(curr@raw.data))
  
}

seurat = MakeSparse(seurat)
saveRDS(seurat, 'Data/soup_output.RDS')
