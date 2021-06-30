# Preprocessing/Blood

This folder contains the code that was used to preprocess the data obtained from sorted subsets of memory CD4+ T cells obtained from the peripheral blood of patients with EoE. This data was used for Figures 6 and 7 of the paper.

The order in which these scripts were run was:
1. sort_assembly.ipynb: Initial processing of the Seq-Well output files in Python and export to R/Seurat.
2. Seurat_assembly.R: Initial processing of the data in Seurat and generation of the core Seurat object used in Figures 6 and 7.
3. umi_correction.ipynb: Repair of UMI errors in TCR data.
4. TCR_integration.ipynb: Integration of TCR data into Seurat object.
