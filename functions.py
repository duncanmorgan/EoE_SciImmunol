# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:41:49 2019

@author: dmorgan
"""
import os
import numpy as np
import pandas as pd
import feather
import anndata as an
import scanpy as sc
import igraph
import statsmodels.api as sm
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import sklearn.preprocessing
import pickle
from statsmodels.formula.api import ols
import statsmodels

def subset(origdata, cells):
    df = origdata.raw.X[origdata.obs_names.isin(cells),:]
    adata = an.AnnData(df)
    adata.obs_names = origdata.obs_names[origdata.obs_names.isin(cells)]
    adata.var_names = origdata.raw.var_names
    adata.obs = origdata.obs.loc[adata.obs_names,:]
    return adata

def process(adata, model = 'linear'):
    adata.raw = adata
    sc.pp.filter_genes_dispersion(adata, flavor = 'seurat')
    if (model == 'linear'):
        sc.pp.regress_out(adata, keys = 'n_counts')
        sc.pp.scale(adata)
        print('linear scaling')
    else:
        adata = poissonscale(adata, 'n_counts')
        print('poisson scaling')
    sc.tl.pca(adata, svd_solver = 'arpack')
    sc.pp.neighbors(adata, n_pcs = 20, n_neighbors = 30)
    sc.tl.umap(adata, min_dist = .3)
    sc.tl.leiden(adata, resolution = None)
    sc.pl.umap(adata,color ='leiden')
    return adata



def getgenes(adata, cluster, field = 'leiden', thresh = .25):
    data1 = adata.raw.X[adata.obs[field].isin(cluster),:]
    data2 = adata.raw.X[~adata.obs[field].isin(cluster),:]
    genes = adata.raw.var_names
    g1sums = (data1>0).sum(axis = 0)
    g2sums = (data2>0).sum(axis = 0)
    genes_use = genes[g1sums + g2sums > (data2.shape[1] + data1.shape[1])*.03]
    mean1 = np.log(np.mean(np.expm1(data1), axis = 0) + 1)
    mean2 = np.log(np.mean(np.expm1(data2), axis = 0) + 1)
    logfc = mean1 - mean2
    pval = scipy.stats.ttest_ind(data1, data2, axis = 0, equal_var = False)
    pct_1 = ((data1 > 0).sum(axis = 0)) / data1.shape[0]
    pct_2 = ((data2 > 0).sum(axis = 0)) / data2.shape[0]
    df = pd.DataFrame({'LogFC' : logfc,'pct1' : pct_1, 'pct2' : pct_2,  'pval' : pval.pvalue})
    df.index = genes
    df = df.sort_values(by = 'pval', ascending = True)
    df = df[(df.pct1 > .1) | (df.pct2 > .1)]
    df = df[abs(df.LogFC) > thresh]
    df['pval_adj'] = statsmodels.stats.multitest.multipletests(df.pval.values, method = 'bonferroni')[1]
    return df.sort_values(by = 'LogFC', ascending = False)


def diffgenes(adata, c1, c2, field = 'leiden', threshold = .25):
    data1 = adata.raw.X[adata.obs[field].isin(c1),:]
    data2 = adata.raw.X[adata.obs[field].isin(c2),:]
    genes = adata.raw.var_names
    g1sums = (data1>0).sum(axis = 0)
    g2sums = (data2>0).sum(axis = 0)
    genes_use = genes[g1sums + g2sums > (data2.shape[1] + data1.shape[1])*.03]
    mean1 = np.log(np.mean(np.expm1(data1), axis = 0) + 1)
    mean2 = np.log(np.mean(np.expm1(data2), axis = 0) + 1)
    logfc = mean1 - mean2
    pval = scipy.stats.ttest_ind(data1, data2, equal_var = False)
    pct_1 = ((data1 > 0).sum(axis = 0)) / data1.shape[0]
    pct_2 = ((data2 > 0).sum(axis = 0)) / data2.shape[0]
    df = pd.DataFrame({'LogFC' : logfc,'pct1' : pct_1, 'pct2' : pct_2,  'pval' : pval.pvalue})
    df.index = genes
    df = df.sort_values(by = 'pval', ascending = True)
    df = df[(df.pct1 > .1) | (df.pct2 > .1)]
    df = df[abs(df.LogFC) > threshold]
    return df.sort_values(by = 'LogFC', ascending = False)

def subprocess(adata, cells, model = 'linear'):
    adata = subset(adata, cells)
    adata = process(adata, model)
    return adata

def removecluster(adata, cluster, field = 'leiden', model = 'linear'):
    obs = adata.obs_names[~adata.obs[field].isin(cluster)]
    return subprocess(adata, obs,  model)

def whichcells(adata, cluster, field = 'leiden'):
     return adata.obs_names[adata.obs[field].isin(cluster)]   
    
  

def barplot(tab, index1 = 'size', index2 = 'leiden'):
    tab = tab.groupby([index1, index2]).size().unstack()
    tab = tab.fillna(0)
    tab = tab.div(tab.sum(axis = 1), axis = 0)
    if (index1 == 'size'):
        tab = tab.reindex(['S', 'M', "L", 'X'])
    elif (index1 == 'tumor'):
        inds = ['S1' ,"S3", 'S4', 'S6', 
               'M1', 'M2', 'M3', 'M4', 
               'L1', 'L2', 'L3', 'L4',
               'X1', 'X2', 'X3', 'X4']
        inds = [x for x in inds if (x in set(tab.index))]
        tab = tab.reindex(inds)
    tab.plot.bar(stacked = True).legend(bbox_to_anchor =(1,1))
    return tab

def scaleRow(row, regressOut):
    import statsmodels.api as sm
    poisson_model = sm.GLM(row, regressOut, family = sm.families.Poisson())
    results = poisson_model.fit()
    adjusted = results.resid_pearson - min(results.resid_pearson)
    scaled = sklearn.preprocessing.scale(adjusted, axis = 0)
    return scaled

def poissonscale(adata, regressOut):
    mat = adata.X
    regressOut = adata.obs[regressOut]
    for i in range(0, adata.X.shape[1]):
        row = mat[:,i]
        regressOut = sm.add_constant(regressOut)
        mat[:,i] = scaleRow(row, regressOut)
    adata.X = mat
    return adata
            
        
def seuratExport(adata, fname):
    df = pd.DataFrame(adata.raw.X.transpose())
    df.columns = adata.obs_names
    df.index = adata.raw.var_names
    df = df.reset_index()
    feather.write_dataframe(df, fname + '.feather') 
    adata.obs.to_csv(fname + '_meta.txt')
