# S-LDSC Pipelines for Nextflow

This repository contains three different Nextflow (NF) pipelines for performing S-LDSC (LD Score Regression). The main distinction among these pipelines lies in how S-LDSC annotations are generated.

All three pipelines utilize S-LDSC through a Singularity image. Instructions on how to operate Singularity with Nextflow can be found [here](https://github.com/roxyisat-rex/Nextflow_with_Singularity_on_PBS/tree/master).

## Pipelines Overview
- For all three pipelines, a table summarizing the results of the S-LDSC run is created alongside a plot in `.pdf` format for users to quickly review the results.

### 1. 0.2% Full Pipeline (`0.2%_Full_pipeline.nf`)

This pipeline creates S-LDSC annotations containing >= 0.2% of SNPs from the baseline model. It selects SNPs overlapping each peak by iterating through the input `.bed` files with the `.bim` files from the Plink files. The selection is based on peaks with the highest levels of cell-type specificity.

### 2. Deciles Full Pipeline (`Deciles_Full_pipeline.nf`)

This pipeline generates 10 S-LDSC annotations for each cell type. It selects all peaks within a cell type that are equal to or above a specified z-score threshold (or other customizable metrics) and splits them into deciles according to z-score levels. This is useful for examining relationships such as cis-regulatory element conservation and cell type specificity. Additionally, this pipeline includes an extra step that maps cCREs to cell-type-specific peaks. Various cCRE files are looped across different `.bed` files to identify cCREs residing in each peak, and the results are outputted in `.csv` format.

### 3. Non-Deciles Full Pipeline (`Non_Deciles_Full_pipeline.nf`)

This pipeline creates S-LDSC annotations without deciles and only selects for cell-type-specific peaks, determined by a z-score threshold.

## Configurations (`sample_NF.config`)
This file contains a sample config file for the pipelines. The HPC system here is PBS pro. This should be customized when using the pipeline as the HPC system is likely different, but hopefully this gives you a bit of a start. It is to be used with the (`.pbs`) file when submitting the NF pipelines to run. 

## Necessities

- It is advised to keep all `.py` and `.R` scripts used in the pipelines in the same directory once downloaded.

- UCSC's liftOver is required for this pipeline, and users should have a Conda R environment and a Python 3.7 environment.

- It is also advised that users go over the singularity with NF page at the start of the readme. 
