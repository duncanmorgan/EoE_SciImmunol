import pandas as pd
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import matplotlib.backends
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import linear_model
import anndata as an
import scipy
import feather

# Input file - table containing sample names, locations of DGE files, locations of ReadCounts files, and locations of soup estimation files
inputFile = "input_file.txt"

# Input parameters
nGene = 300
nCell = 0

try:
    os.stat("Plots")
except:
    os.mkdir("Plots")
    
# Read in DGE Files    
fileNames = pd.read_csv(inputFile, sep = "\t")

# this code is essentially just reading in each file on the list, performining an initial filtering of low quality cells, and concatenating 
# these files to a matrix. the files are then saved as .feather files (compressed format that can be opened in R or Python)

tic = time.time()
first = True
for i in range(0, len(fileNames)):
        
    # read and print sample name
    sampName = fileNames.Sample.iloc[i]
    print(sampName)
        
    # read in DGE and readCounts files, calculate nUMI
    cells = pd.read_csv(fileNames.values[i,1], index_col = 0, header = 0, delim_whitespace = True)
    nUMIs = cells.sum(axis = 0)
    readCounts = pd.read_csv(fileNames.values[i,2], index_col = 0, header = 0, delim_whitespace = True)
    readCounts = readCounts.reindex(index = list(cells))
    cells.index = cells.index.str.upper()
    cells = cells.loc[~cells.index.duplicated(),:]    
    cells
        
    # plotting 
    spam = plt.figure()
    plt.rcParams['figure.figsize'] = [7,8]
    plt.plot(readCounts['ExonReads'], nUMIs, 'ko')
    plt.xlabel('Total Reads')
    plt.ylabel('Number of UMIs')
    plt.title(sampName + " Complexity")  
    regr = linear_model.LinearRegression()
    regr.fit(X = nUMIs.values.reshape(-1,1), y = readCounts['ExonReads'].values)
    plt.plot(X = nUMIs.values.reshape(-1,1), Y = regr.predict(nUMIs.values.reshape(-1,1)))
    c = np.array2string(regr.coef_[0])
    plt.annotate("Complexity = " + c, (0,plt.ylim()[1]*.95))
    plt.savefig(fname = "Plots/" + sampName + '.png')
    plt.close()
        
    # filter by nGene
    keepcells = (cells.values>0).sum(axis = 0) > nGene
    keepgenes = (cells.values > 0).sum(axis = 1) > nCell
    cells = cells.loc[keepgenes, keepcells]
        
    readCounts = readCounts[:][keepcells]
    n = cells.shape[1]

    # rename columns of data frame to prevent barcode collisions
    identsadd = [sampName]*n
    index = [identsadd[i] + "_" + str(i) for i in range(0, len(identsadd))]
    cells.columns = index
    bcs = readCounts.index.values
    print(len(bcs))
    
    if first:
        cellsAll = cells.copy()
        readCountsAll = readCounts.copy()
        bcsAll= bcs.copy()
        first = False
        idents_all = identsadd

    else:
        cellsAll = pd.DataFrame.join(cellsAll, cells, how = 'outer')
        readCountsAll = readCountsAll.append(readCounts)
        bcsAll = np.append(bcsAll, bcs)
        idents_all = np.append(idents_all, identsadd)
        
toc = time.time()
print(toc - tic)
cellsAll = cellsAll.fillna(0)
cellsAll = cellsAll.reset_index()
feather.write_dataframe(cellsAll, 'cellsAll.feather')  
feather.write_dataframe(readCountsAll, 'readCountsAll.feather')
np.savetxt('bcs.txt', bcsAll, fmt = "%s")

from functions import *
# read in and fix data
cells = pd.read_feather('cellsAll.feather')
cells.index = cells.loc[:, 'Gene']
cells = cells.drop('Gene', axis = 1)

# SCANPY
# this is an initial processing of the data
sc.settings.verbosity = 4
adata = an.AnnData(cells.values.transpose())
adata.var_names = cells.index
adata.obs_names = cells.columns
bcs = np.loadtxt('bcs.txt', dtype = str)
adata.obs['orig'] =[x.split('_')[0] for x in adata.obs_names]
adata.obs['bcs']= bcs
adata.obs['orig.ident'] = idents_all
adata.obs_names_make_unique()
adata.var_names_make_unique()
nGene =500
sc.pp.filter_genes(adata, min_cells = 3)
sc.pp.filter_cells(adata, min_genes = nGene)
adata

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata = process(adata)

# fill out metadata
samples = set(fileNames.Sample)

tissue = [x.split('-')[1][0] for x in samples]
patient = [x.split('-')[0] for x in samples]
metadata = pd.DataFrame({'Tissue' : tissue, 'Patient' : patient}, dtype = str)
metadata.index = samples
metadata['Diagnosis'] = 'Diseased'
metadata.Diagnosis.loc[metadata.Patient.isin(['392', '355', '468', '249'])] = 'Remission'
metadata

# add metadata to adata object
adata.obs['tissue'] = metadata.Tissue.loc[adata.obs['orig.ident']].values
adata.obs['patient']= metadata.Patient.loc[adata.obs['orig.ident']].values
adata.obs['diagnosis']= metadata.Diagnosis.loc[adata.obs['orig.ident']].values
adata.obs['ci'] = [adata.obs.diagnosis.loc[x] + '_' + adata.obs.tissue.loc[x] for x in adata.obs_names]
adata.obs.head()

sc.pl.umap(adata, color = ['leiden'])
umap(adata, 'leiden')

pickle.dump(adata, file = open('all_raw.pickle', 'wb'), protocol = 4)
seuratExport(adata, 'all_raw')





