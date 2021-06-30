# Preprocessing/Tissue

This folder contains the code that was used to process single-cell data collected from esophageal biopsies and duodenal biopsies of ten patients with EoE. In sum, four datasets are generated:
  * All cells collected from esophageal and duodeanl biopsies (used for Figure 1)
  * Cells from the granulocyte cluster (i.e. eosinophils; used for Figure 2)
  * Esophageal T cells (used for Figures 3, 4, and 5)
  * Duodenal T cells (used for Figure 3)
 
 The order in which these scripts were run is as follows:
 1. initial_python_assembly.py: This script takes as inputs the count matricies generated from Seq-Well (specifically, an in-house version of the Dropseq pipeline utilized by the Love lab), performs an initial filtering of cells and genes, and assembles all of the samples into a single matrix for export into R. 
 2. initial_R_assembly: Initial assembly of the data in R/Seurat.
 3. SoupCorrection.R: This script performs ambient RNA correction on the single-cell data using SoupX v0.3.0.
 4. TissueProcessing.ipynb: Reprocessing of the SoupX-corrected dataset and generation of the Seurat object used for Figure 1.
 5. EosinophilProcessing.ipynb: Sub-clustering of the granulocyte cluster (mostly eosinophils).
 6. TCellProcessing.ipynb: Sub-clustering of esophageal and duodenal T cells.
 7. tissue_tcr_correction.ipynb: Correction of the TCR sequencing data from esophageal and duodenal biopsies for mispairings due to barcode collision errors.
 8. TCR_integration: Integration of TCR data with whole-transcriptome data for T cell objects.
