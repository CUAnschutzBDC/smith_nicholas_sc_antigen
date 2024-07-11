# Exploration of islet reactive B cells
Using single cell RNA-seq, single cell CITE-seq, and single cell LIBRA-seq

<!--ts-->
   * [Setting up the pipeline](#Setting-up-the-pipeline)
     * [Download required packages](#download-required-packages)
     * [`Snakemake` installation](#snakemake-installation)
       * [Setting up `docker` or `singularity`](#setting-up-docker-or-singularity)
     * [Updating the configfile](#updating-the-config-file)
     * [Running the pipeline](#running-the-pipeline)
     * [T1K](#t1k)
  * [R analysis](#r-analysis)
    * [Overall analysis steps](#overall-analysis-steps)
    * [Detailed description of each script](#detailed-description-of-each-script)
  * [Seurat object](#seurat-object)
    * [Assays](#assays)
    * [Most useful meta data columns](#most-useful-meta-data-columns)
    * [All meta data columns](#all-meta-data-columns)

<!--te-->


This repository contains a pipeline that can be used to fully replicate all analysis and figures associated with the manuscript [TODO add manuscript name].

In order to fully replicate analysis download the following docker images from dockerhub:
* [smith_2024_r_docker](https://hub.docker.com/r/kwellswrasman/smith_2024_r_docker)
* [smith_2024_dropkick](https://hub.docker.com/r/kwellswrasman/smith_2024_dropkick)
* [smith_2024_scar](https://hub.docker.com/r/kwellswrasman/smith_2024_scar)
* [smith_2024_t1k](https://hub.docker.com/r/kwellswrasman/smith_2024_t1k)

All aspects of this pipeline are contained within a `Snakemake` pipeline. This will be able to start with fastq files and generate a final directory containing all final figures numbered based on the figure in the manuscript.

This pipeline was written to interact with an lsf scheduler. Feel free to reach out if you need help modifying to run on your own server.

## Setting up the pipeline

### Download required packages

The majority of the packages required to run this are in docker images that are publicly available:
* [smith_2024_r_docker](https://hub.docker.com/r/kwellswrasman/smith_2024_r_docker)
* [smith_2024_dropkick](https://hub.docker.com/r/kwellswrasman/smith_2024_dropkick)
* [smith_2024_scar](https://hub.docker.com/r/kwellswrasman/smith_2024_scar)
* [smith_2024_t1k](https://hub.docker.com/r/kwellswrasman/smith_2024_t1k)

To load these packages as docker images:


Other required packages include
* `cellranger V 7.1.0` from [10x genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions)
* `singularity`, `apptainer`, or `docker`
  * [`singularity`/`apptainer` install instructions](https://apptainer.org/admin-docs/master/installation.html)
  * [`docker` install instructions](https://docs.docker.com/engine/install/)
* `snakemake` V 6.0.3
  * [`snakemake` install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

#### `Snakemake` installation
1. Download and install miniconda3: For Linux
```{bash}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install `Snakemake`:
```{bash}
conda install snakemake=6.3.0 -c bioconda -c conda-forge
```

#### Setting up `docker` or `singularity`
If you are using `docker`, you can download the images directly from dockerhub
[TODO] - Test how downloading works.
[TODO] - Test making sif from docker hub directly

### Updating the config file

Once everything is installed, you will need to update the config file to make sure all paths are correctly configured for your environment.

>* RAW_DATA: Path to the directory containing fastq files
>* SAMPLES: the list of samples you want to test. This is the name that will be in the output files. The order must be identical to the order of RNA_SAMPLES, ADT_SAMPLES, and VDJ_SAMPLES
>* AGGR_SAMPLES: Which samples should be aggregated together at the end of the pipeline. These should be samples that were split between multiple 10x runs. This was left blank for our experiment.
>* RNA_SAMPLES: The name of the RNA samples. Not the full fastq name, but the name that is the same for all samples. Likely ends in "GEX". This will likely need to be updated depending on how you downloaded files from GEO.
>* ADT_SAMPLES: *optional* The name of the samples from CITE-seq or hashtagging. If CITE-seq and hashtagging files are separate, include them both separated by a comma. Put samples in the same order as their RNA counterparts. If CITE-seq or hashtagging were not performed, leave this blank. 
>* VDJ_SAMPLES: *optional* The name of the samples from VDJ-seq. Put samples in the same order as their RNA counterparts. If VDJ-sequencing was not peformed, leave this blank.
>* RESULTS: Path to the output directory. This will be created if it doesn't already exist
>* GENOME: Path to the cellranger genome reference. We used refdata-gex-GRCh38-2020-A [downloaded from 10x genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads)
>* ADT_REF: *optional* Path to the ADT-reference. Should be a comma separated file with the columns described in the 10x tutorial: id, name, read, pattern, sequence, and feature_type. The feature_type will be Antibody Capture. The name will be the name in the final output matrix. Leave this blank if CITE-seq or hashtagging were not performed. This can be found [here](https://github.com/CUAnschutzBDC/smith_nicholas_sc_antigen/blob/main/files/antibodies.csv)
>* VDJ_REF: Path to the cellranger VDJ reference. If VDJ sequencing were not performed, leave this blank.
>* MAX_JOBS: The maximum number of jobs that can be submitted by cell ranger at a time
>* LSF_TEMPLATE: Path to an LSF template. One is included in this git repo.
>* CHEMISTRY: *optional* Arguments to the `--chemstiry` flag in cellranger count. If left blank, chemistry will be `auto`. Only use this if the pipeline failed because of a failure to detect the chemistry. Can be filled in for only some samples.
>* SCRIPT_PATH: The path to the scripts. Does not need to be updated if you downloaded this directory and haven't modified the paths.
>* SCRIPTS_RUN: Name of one-off scripts that need to be run
>* SAMPLE_SCRIPTS: A dictionary containing 1. What samples to run the scripts on. Either should match the sample names in SAMPLES or "all" for all samples. 2. What scripts should be run. These scripts will be found in `src/scripts/individual_analysis` and will be run on all samples. Any parameters that are sample-specific will be found in [the sample info file](https://github.com/CUAnschutzBDC/smith_nicholas_sc_antigen/blob/main/files/sample_info.tsv)
>* MERGE_SCRIPTS: A dictionary containing 1. What samples to merge, if set to "all", all samples will be merged. 2. What scripts should be run. These will be found in `src/scripts/integrated_analysis`. Any parameters such as the batch correction will be found in [the sample info file](https://github.com/CUAnschutzBDC/smith_nicholas_sc_antigen/blob/main/files/sample_info.tsv)
>* SAMPLE_METADATA: Any metadata to add to the sample objects. This file is found [here]()
>* SAMPLE_INFO: Path to the file containing sample info. If you are running analysis from scratch, you can fill this out as you get outputs to help make decisions. This file is [here](https://github.com/CUAnschutzBDC/smith_nicholas_sc_antigen/blob/main/files/sample_info.tsv)
>* RSRIPT_CONTAINER: Path to docker or singularity rscript file. See above for access.
>* DROPKICK_CONTAINER: Path to docker or singularity dropkick file. See above for access.
>* SCAR_CONTAINER: Path to docker or singularity scar file. See above for access.
>* T1K_CONTAINER: Path to docker or singularity T1K file. See above for access.

### Running the pipeline
The script to submit the job is [here](https://github.com/CUAnschutzBDC/smith_nicholas_sc_antigen/blob/main/snakecharmer.sh). You will need to update this command to work with your cluster or system.

This pipeline will run
* `Cellranger` data preprocessing
* `Dropkick` to determine cutoffs for cells [Dropkick](https://github.com/KenLauLab/dropkick)
* `Scar` to remove ambient protein reads [scar](https://github.com/Novartis/scar)
* `Immcantation` VDJ clone calling [Immcantation](https://immcantation.readthedocs.io/en/stable/index.html)
* `T1K` for HLA identification [T1K](https://github.com/mourisl/T1K)
* R analysis scripts to recreate analysis and figures (described below)

This runs through snakemake and will process the whole pipeline for you.

I highly recommend looking at the csv files that are generated and passed to cell ranger to ensure that the correct fastq files have been detected for each sample.

### T1K
Once you've run T1K, map to HLA types [here](https://www.ebi.ac.uk/ipd/imgt/hla/alleles/)

## R analysis

R scripts for the analysis are found in `src/scripts`. There are two directories here. One called `individual_analysis` includes scripts to process an individual sample. `intigrated_analysis` contains scripts to run analysis on multiple samples.

All packages needed to run this scripts are in the smith_r_docker described above.

### Overall analysis steps

Scripts within each directory are numbered based on the order in which they should be run.

They include:

1. Steps for intial processing and filtering of the data
2. Steps for doublet removal
3. Steps for normalizing ADTs (if they are present in your data)
4. Steps for dimensional reduction with PCA and UMAP either on RNA alone or with ADTs included
5. Steps to name clusters based on references with `clustifyr`
6. Steps for marker detection
7. Steps for gene ontology and pathway analysis
8. Steps for integrating multiple samples
9. Steps for differential expression when multiple samples are present

For many of these steps, I've included different methods of performing them (for example batch correction can be done with `harmony` or `fastMNN`) and ways of comparing the different approaches so you can make an informed decision about what method is best for your samples.

For all scripts in this analysis, please make sure you understand the processing steps so you know what are the appropriate parameters for each step. These scripts are not intended to be run blindly with no understanding of the underlying algorithims.

### Detailed description of each script

#### Individual analysis

* `01_dropkick.py` - Runs [`dropkick`](https://github.com/KenLauLab/dropkick) on each sample individually to identify cells to be filtered.
* `01b_run_scar.py` - Runs [`scar`](https://github.com/Novartis/scar) on each sample individually to estimate and remove background contamination based on ambient droplets
* `02_Initial_processing.R` - Sets up seurat objects for each sample individually. This runs the following steps.
  * Creates a seurat object that contains ADT, HTO, and RNA data
  * Determines the percentage of reads that map to the mitochondria.
  * Runs [`scuttle`](https://www.bioconductor.org/packages/release/bioc/html/scuttle.html) to estimate cutoffs based on quality metics and compares this to the output of `dropkick`. This comparison is output as a heatmap and barplots saved in `results/R_analysis/{sample_name}/images/dropkick_vs_cellqc.pdf`
  * Subsets the seurat object based on the `scuttle` cutoffs 
  * Reads in VDJ data using [`djvdj`](https://github.com/rnabioco/djvdj)
  * Processes tetramers
    * Pulls the tetramers out of the ADT slot
    * Normalizes the tetramers and ADT using the `CLR` method
    * Identifies the libra score from the raw tetramer data and adds the score as an assay
    ```R
    clr_matrix <- log2((all_tetramers + 1) / rowMeans(all_tetramers))

	  libra_score <- scale(clr_matrix, center = TRUE, scale = TRUE)
    ```
  * Reads in the output from `scar`
    * Pulls out the tetramers
    * Normalizes the tetramers and ADTs with the `CLR` method
    * Identifies the libra score form the scar corrected tetramer data and add the score as an assay
    ```R
    clr_matrix <- log2((all_tetramers + 1) / rowMeans(all_tetramers))

	  libra_score <- scale(clr_matrix, center = TRUE, scale = TRUE)
    ```
  * Runs an updated form of `HTOdemux` on the tetramers (both raw and scar corrected)
    * `HTOdemux` was updated because the previous form would take the highest counts for any feature even if it was not the feature with the highest value above the determined cutoff. Because of this, sometimes the actual feature that was higher than the cutoff was not returned while features that were below the  cutoff were returned. The updated function also returns a list of all features that were above the cutoff in the `full_hash_id` column.
  * Identifies positive values for each tetramer based on a libra cutoff of 1 for both the raw and `scar` corrected libra scores.
  * Runs SCT normalization on the RNA assay
* `03_remove_doublets.R` - Identifies and tags doublets using [`DoubletFinder`](https://github.com/chris-mcginnis-ucsf/DoubletFinder). Images related to this can be found in `results/R_analysis/{sample_name}/images/dropkick_vs_cellqc.pdf`.
* `04_adt_dsb_normalization.R` - Uses [`dsb`](https://github.com/niaid/dsb) to identify and remove background. This is similar to `scar` and was only used to compare to `scar` output.
* `05_PCA.R` - Runs PCA using methods from `Seurat`. Makes the output `results/R_analysis/{sample_name}/images/RNA_pca.pdf` to help determine the number of pcs to use for downstream processing.
* `06a_UMAP_find_resolution.R` - Uses the number of pcs found above (specified in the file `files/sample_info.tsv`) and uses `clustree` to visualize many clustering resolutions. The `clustree` and umaps of each resolution will be in a file `results/R_analysis/{sample_name}/images/clustering_resolution.pdf` to guide resolution selection for downstream analysis.
* `06b_UMAP.R` - Taking the number of pcs and resolution determined in the previous steps (specified in `files/sample_info.tsv`), generates the final umap and clustering for the object. The final images can be found in `results/R_analysis/{sample_name}/images/RNA_pca.pdf`
* `07_name_clusters.R` - Uses two references the [seurat reference](https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub) which can be downloaded [here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat) and a [BND reference previously published by Mia Smith](https://rupress.org/jem/article/220/8/e20221604/214110/Identification-of-an-anergic-BND-cell-derived) to name clusters. I name each cluster on individual samples to provide more confidence in the integration of samples.
* `08_find_markers.R` - Finds markers for each of the clusters and cell types using a one-vs-all approach with `FindAllMarkers` in Seurat. This outputs a file of markers that should be used to double check cell type determination.
* `09_demultiplex_tets.R` - More attempts to identify positive tetramers. Here, UMAPs are generated on the `ADT_CLR`, `TET_CLR`, `ADT_DSB` and the `TET_DSB` assays and clusters are found in the same. These clusters are compared to the tetramer deremination at other steps. Nothing from this script was used in the final manuscript.
* `10_compare_dsb_clr.R` - A script comparing the output of the dsb and clr normalization approaches. Nothing from this script was used in the final manuscript
* `11_remove_ambience.R` - A script that runs `removeAmbience` from `DropletUtils` to attempt to remove ambient contamination. This assay was used for all steps and compared to the non-abient removal downstream, but was not used in the final manuscript because ambience removal from a popultion of only B cells did not make much difference.
* `12_improve_cutoff.R` - Steps to try a tetramer identification that uses the non-b cells in our population to determine a cutoff. A cutoff of the 95th quartile of the non-b cells was used to identify a cutoff. A proportion above the cutoff was than calculated for each cell (and made into an assay) and tetramer calls were based on any tetramers with scores above 1. This scoring system ended up being used in the final manuscript.


*NOTE on tetramer labeling. While libra, HTODemux, scar, and raw data were all used, only the libra score on scar corrected values and the t-cell/myeloid cell cutoff on the scar corrected values were used to identify tetramer names the others were just used to test methods and compare results*

#### Integrated analysis
[TODO]

## Seurat object

Below is a detailed explanation of all parts of the seurat object

### Assays
* `RNA` - Raw and log normalized RNA-seq reads
* `ADT` - Raw and `CLR` normalized values for all ADTs
* `TET` - Raw and `CLR` normalized values for all tetramers
* `TET_LIBRA` - libra score on all raw tetramer values
* `SCAR_ADT` - `scar` corrected and `CLR` normalized `scar` corrected values for all ADTs
* `SCAR_ADT_LOG` - `scar` corrected and `log` normalized (with `log1p`) `scar` corrected values for all ADTs
* `SCAR_TET` - `scar` corrected and `CLR` normalized `scar` corrected values for all tetramers
* `SCAR_TET_LOG` - `scar` corrected and `log` normalized (with `log1p`) `scar` corrected values for all tetramers
* `SCAR_TET_LIBRA` - libra score on all `scar` corrected tetramer values
* `DSB_ADT` - ADTs with background corrected with `dsb`
* `SCT` - SCT normalized data
* `TET_PROPORTIONS` - The proportions above the cutoff based on the modified `HTOdemux` run on the raw tetramer data.
* `NEW_TET_PROPORTIONS` - The proportion above the cutoff based on the 95th quartile of the non-b cells included in the assay. The scores were determined based on the scar corrected tetramer data. This was added in `12_improve_cutoffs.R`

### Most useful meta data columns
These columns are all in the metadata available on GEO

* `Seurat` default columns
  * `orig.ident` - Sample name
  * `nCount_RNA` - Number of RNA molecules in the cell
  * `nFeature_RNA` - Number of genes in the cell
  * `nCount_ADT` - Number of ADT molecules in the cell
  * `nFeature_ADT` - Number of different ADTs in the cell
  * `percent.mt` - Percent of reads mapping to the mitochondria
* Meta data columns
  * `Sex` - Sex
  * `Collection.Date` - Date the blood was collected
  * `Age.at.Collection..years.` - Age of patient at date of collection in years
  * `Age.at.Collection..Months.` - Age of patient at date of collection in months
  * `Status` - Disease status (ND = non-diabetic, T1D = type 1 diabetes, AAB = autoantibody positive)
  * `Date.of.Diagnosis` - Date of T1D diagnosis
  * `Days.post.onset` - Number of days between diagnosis and collection
  * `millions.of.cells.frozen` - Number of cells collected
  * `HLA` - HLA type of the individual
  * `Autoantibodies` - What autoantibodies were detected
  * `date.processed.for.scSeq` - The date cells were thawed and prepped for scRNA-seq
  * `Notes..FDR.relationship.` - For first degree relative, how they are related
  * `sample` - Sample ID
  * `Ethnicity` - The ethnicity of the individual
* Cell cycle columns
  * `Phase` - Final cell cycle phase determination as determined by Seurat's `CellCycleScoring` function
* VDJ columns (Added by [`djvdj`](https://github.com/rnabioco/djvdj))
  * `chains` - What chains are present in the cell, each chain is separated by a semi colon
  * `n_chains` - How many chains are in the cell
  * `cdr3` - The amino acid sequence of the CDR3 region(s). Each chain is separated by a semi colon
  * `cdr3_nt`  - The nucleotide sequence of the CDR3 region(s). Each chain is separated by a semi colon
  * `cdr3_length` - The CDR3 amino acid length(s). Each chain is separated by a semi colon
  * `cdr3_nt_length` - The CDR3 nucleotide length(s). Each chain is separated by a semi colon
  * `v_gene` - The V gene(s) identified by cellranger. Each chain is separated by a semi colon
  * `d_gene` - The D gene(s) identified by cellranger. Each chain is separated by a semi colon
  * `j_gene` - The J gene(s) identified by cellranger. Each chain is separated by a semi colon
  * `c_gene` - The C gene(s) identified by cellranger. Each chain is separated by a semi colon
  * `isotype` - The isotype(s) identified by cellranger. Each chain is separated by a semi colon
  * `reads` - How many reads support each chain. Each chain is separated by a semi colon
  * `umis` - How many umis support each chain. Each chain is separated by a semi colon
  * `productive` - If each chain (separated by a simi colon) is productive
  * `full_length` - If each chain (separated by a simi colon) is full length
  * `paired` - If there were paired chains present
  * `all_mis_freq` - The number of mimatches in all region(s). Each chain is separated by a semi colon
* Columns relating to tetramer calling
  * `libra_tet_hash_id` - The classification of tetramers based on the libra score - the name of the tetramer if only one is above the cutoff or multi-reactive (Islet or Other) depending on what tetramers were above the cutoff. Defined in `02_Initial_processing.R`. This uses the libra score based on the raw tetramer data.
  * `libra_full_hash_id` - Full id based on libra scores. This gives all tetramers above the cutoff. Defined in `02_Initial_processing.R`. This uses the libra scores computed based on the raw tetramer data.
  * `scar_libra_tet_hash_id` - The classification of tetramers based on the libra score - the name of the tetramer if only one is above the cutoff or multi-reactive (Islet or Other) depending on what tetramers were above the cutoff. Defined in `02_Initial_processing.R`. This uses the libra score based on the scar corrected tetramer data.
  * `scar_libra_full_hash_id` - Full id based on libra scores. This gives all tetramers above the cutoff. Defined in `02_Initial_processing.R`. This uses the libra scores computed based on the scar corrected tetramer data.
  * `tet_name_cutoff` - Cutoff determined based on the non-b cells present for each sample. Here a cutoff was drawn ath the 95th quartile of the non-b cells. The value is the name of the tetramer if only one is above the cutoff or multi-reactive (Islet or Other) depending on what tetramers were above the cutoff. Defined in `12_improve_cutoff.R`. This uses the scar corrected tetramer data.
  * `full_tet_name_cutoff` - Cutoff determined based on the non-b cells present for each sample. Here a cutoff was drawn ath the 95th quartile of the non-b cells. This gives all tetramers above the cutoff. Defined in `12_improve_cutoff.R`. This uses the scar corrected tetramer data.
* Clustering and cell type columns
  * `RNA_cluster` - The final clusters used for downstream analysis
  * `cluster_celltype` - A combination of the final cluster and final cell type
  * `final_celltype` - Final cell types that were used for making the figures 
* immcantation columns
  * `final_clone` - Clone call by immcantation
  * `imcantation_isotype` - Isotype determined by immcantation

### All meta data columns
These columns will be generated when the seurat object is made following the pipeline in this repository. These are not all in the published meta data.

* `Seurat` default columns
  * `orig.ident` - Sample name
  * `nCount_RNA` - Number of RNA molecules in the cell
  * `nFeature_RNA` - Number of genes in the cell
  * `nCount_ADT` - Number of ADT molecules in the cell
  * `nFeature_ADT` - Number of different ADTs in the cell
  * `nCount_SCT` - Number of reads based on SCT normalization
  * `nFeature_SCT` - Number of genes based on SCT normalization 
  * `nCount_TET` - Number of reads from the tetramers
  * `nFeature_TET` - Number of features from the tetramers
  * `nCount_TET_LIBRA` - Number of reads associated with the libra score (ignore)
  * `nFeature_TET_LIBRA` - Number of genes associated with the libra score (ignore)
  * `nCount_SCAR_ADT_LOG` - Number of reads from the log normalized scar corrected ADTs (ignore)
  * `nFeature_SCAR_ADT_LOG` - Number of genes from the log normalized scar corrected ADTs (ignore)
  * `nCount_SCAR_ADT` - Number of reads from the scar corrected ADTs
  * `nFeature_SCAR_ADT` - Number of genes from the scar corrected ADTs
  * `nCount_SCAR_TET` - Number of reads from the scar corrected tetramers
  * `nFeature_SCAR_TET` - Number of genes from the scar corrected tetramers
  * `nCount_SCAR_TET_LOG` - Number of reads from the log normalized scar corrected tetramers (ignore)
  * `nFeature_SCAR_TET_LOG` - Number of genes from the log normalized scar corrected tetramers (ignore)
  * `nCount_SCAR_TET_LIBRA` - Number of reads associated with the tetramer libra score (ignore)
  * `nFeature_SCAR_TET_LIBRA` - Number of genes associated with the tetramer libra score (ignore)
  * `nCount_CLR_ADT` - Number of reads associated with the clr normalized adts (ignore)
  * `nFeature_CLR_ADT` - Number of genes associated with the clr normalized adts (ignore)
  * `nCount_CLR_TET` - Number of reads associated with the clr normalized tetramers (ignore)
  * `nFeature_CLR_TET` - Number of genes associated with the clr normalized tetramers (ignore)
  * `nCount_AMBRNA` - Number of reads from the ambient corrected rna
  * `nFeature_AMBRNA` - Number of genes from the ambient corrected rna
  * `nCount_TET_PROPORTIONS` - Number of reads from the hto determined proportions on the raw data (ignore)
  * `nFeature_TET_PROPORTIONS` - Number of features from the hto determined proportions on the raw data (ignore)
  * `nCount_SCAR_TET_PROPORTIONS` - Number of reads from the hto determined proportions on the scar corrected data (ignore)
  * `nFeature_SCAR_TET_PROPORTIONS` - Number of features from the quantile determined proportions on the scar corrected data (ignore)
  * `nCount_DSB_ADT` - Number of reads associated with the dsb normalized adts (ignore)
  * `nFeature_DSB_ADT` - Number of genes associated with the dsb normalized adts (ignore)
  * `nCount_DSB_TET` - Number of reads associated with the dsb normalized tetramers (ignore)
  * `nFeature_DSB_TET` - Number of genes associated with the dsb normalized tetramers (ignore)
  * `nCount_NEW_TET_PROPORTIONS` - Number of reads from the qantile determined proportions based on the scar corrected tetramers (ignore)
  * `nFeature_NEW_TET_PROPORTIONS` - Number of features from the quantile determined proportions based on the scar corrected tetramers (ignore)
  * `percent.mt` - Percent of reads mapping to the mitochondria
* Meta data columns
  * `ID` - Sample id
  * `Initials` - Individual initials
  * `Sample.Name` - Full sample name
  * `Sex` - Sex
  * `Collection.Date` - Date the blood was collected
  * `Age.at.Collection..years.` - Age of patient at date of collection in years
  * `Age.at.Collection..Months.` - Age of patient at date of collection in months
  * `Status` - Disease status (ND = non-diabetic, T1D = type 1 diabetes, AAB = autoantibody positive)
  * `Date.of.Diagnosis` - Date of T1D diagnosis
  * `Days.post.onset` - Number of days between diagnosis and collection
  * `millions.of.cells.frozen` - Number of cells collected
  * `HLA` - HLA type of the individual
  * `Autoantibodies` - What autoantibodies were detected
  * `date.processed.for.scSeq` - The date cells were thawed and prepped for scRNA-seq
  * `Notes..FDR.relationship.` - For first degree relative, how they are related
  * `sample` - Sample ID
  * `Ethnicity` - The ethnicity of the individual
* `scuttle` qc columns (Added by [`scuttle`](https://www.bioconductor.org/packages/release/bioc/html/scuttle.html))
  * `cell_qc_sum` - Same as `nCount_RNA`
  * `cell_qc_detected` - Same as `nFeature_RNA`
  * `cell_qc_subsets_Mito_sum` - Sum of mito counts per cell
  * `cell_qc_subsets_Mito_detected` - Number of mit genes per cell
  * `cell_qc_subsets_Mito_percent` - Percentmit per cell
  * `cell_qc_altexps_ADT_sum` - Same as `nCount_ADT`
  * `cell_qc_altexps_ADT_detected` - Same as `nFeature_ADT`
  * `cell_qc_altexps_ADT_percent` - Percent of ADT counts per cell
  * `cell_qc_total` - sum of counts for each cell across the main and alternative experiment
  * `cell_qc_low_lib_size` - If the cell passed QC based on library size
  * `cell_qc_low_n_features` - If the cell passed QC based on number of features
  * `cell_qc_high_subsets_Mito_percent` - If the cell passed QC based on mito percent
  * `cell_qc_discard` - If the cell is flagged to be removed by `perCellQCFilters`
* Cell cycle columns
  * `S.Score` - Cell cycle score for S phase as determined by Seurat's `CellCycleScoring` function
  * `G2M.Score` - Cell cycle score for G2M phase as determined by Seurat's `CellCycleScoring` function
  * `Phase` - Final cell cycle phase determination as determined by Seurat's `CellCycleScoring` function
* VDJ columns (Added by [`djvdj`](https://github.com/rnabioco/djvdj))
  * `clonotype_id` - Clonotype determined by `djvdj`
  * `exact_subclonotype_id` - Clonotype sub id determined by `djvdj`
  * `chains` - What chains are present in the cell, each chain is separated by a semi colon
  * `n_chains` - How many chains are in the cell
  * `cdr3` - The amino acid sequence of the CDR3 region(s). Each chain is separated by a semi colon
  * `cdr3_nt`  - The nucleotide sequence of the CDR3 region(s). Each chain is separated by a semi colon
  * `cdr3_length` - The CDR3 amino acid length(s). Each chain is separated by a semi colon
  * `cdr3_nt_length` - The CDR3 nucleotide length(s). Each chain is separated by a semi colon
  * `v_gene` - The V gene(s) identified by cellranger. Each chain is separated by a semi colon
  * `d_gene` - The D gene(s) identified by cellranger. Each chain is separated by a semi colon
  * `j_gene` - The J gene(s) identified by cellranger. Each chain is separated by a semi colon
  * `c_gene` - The C gene(s) identified by cellranger. Each chain is separated by a semi colon
  * `isotype` - The isotype(s) identified by cellranger. Each chain is separated by a semi colon
  * `reads` - How many reads support each chain. Each chain is separated by a semi colon
  * `umis` - How many umis support each chain. Each chain is separated by a semi colon
  * `productive` - If each chain (separated by a simi colon) is productive
  * `full_length` - If each chain (separated by a simi colon) is full length
  * `paired` - If there were paired chains present
  * `v_ins` - The number of insertions in the v region(s). Each chain is separated by a semi colon
  * `v_del` - The number of deletions in the v region(s). Each chain is separated by a semi colon
  * `v_mis` - The number of mismatches in the v region(s). Each chain is separated by a semi colon
  * `d_ins` - The number of insertions in the d region(s). Each chain is separated by a semi colon
  * `d_del` - The number of deletions in the d region(s). Each chain is separated by a semi colon
  * `d_mis` - The number of mismatches in the d region(s). Each chain is separated by a semi colon
  * `j_ins` - The number of insertions in the j region(s). Each chain is separated by a semi colon
  * `j_del` - The number of deletions in the j region(s). Each chain is separated by a semi colon
  * `j_mis` - The number of mismatches in the j region(s). Each chain is separated by a semi colon
  * `c_ins` - The number of insertions in the c region(s). Each chain is separated by a semi colon
  * `c_del` - The number of deletions in the c region(s). Each chain is separated by a semi colon
  * `c_mis` - The number of mismatches in the c region(s). Each chain is separated by a semi colon
  * `all_ins` - The number of insertions in all region(s). Each chain is separated by a semi colon
  * `all_del` - The number of deletions in all region(s). Each chain is separated by a semi colon
  * `all_mis` - The number of mismatches in all region(s). Each chain is separated by a semi colon
  * `vd_ins` - The number of insertions in the v and d region(s). Each chain is separated by a semi colon
  * `vd_del` - The number of deletions in the c and d region(s). Each chain is separated by a semi colon
  * `dj_ins` - The number of insertions in the d and j region(s). Each chain is separated by a semi colon
  * `dj_del` - The number of deletions in the d and j region(s). Each chain is separated by a semi colon
  * `v_mis_freq` - The number of mimatches in the v region(s). Each chain is separated by a semi colon
  * `d_mis_freq` - The number of mimatches in the d region(s). Each chain is separated by a semi colon
  * `j_mis_freq` - The number of mimatches in the j region(s). Each chain is separated by a semi colon
  * `c_mis_freq` - The number of mimatches in the c region(s). Each chain is separated by a semi colon
  * `all_mis_freq` - The number of mimatches in all region(s). Each chain is separated by a semi colon
* Columns relating to tetramer calling
  * `TET_maxID` - Output of running the new version of `HTOdemux` based on the `TET` assay (the `CLR` normalized raw tetramer counts). This gives the highest tet id. (Ignore)
  * `TET_secondID` - Output of running the new version of `HTOdemux` based on the `TET` assay (the `CLR` normalized raw tetramer counts). This gives the second highest tet id. (Ignore)
  * `TET_margin` - Output of running the new version of `HTOdemux` based on the `TET` assay (the `CLR` normalized raw tetramer counts). This gives the difference between the highest and second highest cutoffs (ignore)
  * `TET_classification` - Output of running the new version of `HTOdemux` based on the `TET` assay (the `CLR` normalized raw tetramer counts). This gives the two highest tet ids. (ignore)
  * `TET_classification.global` - Output of running the new version of `HTOdemux` based on the `TET` assay (the `CLR` normalized raw tetramer counts). This gives the values `Negative`, `Singlet` and `Doublet` depending on how many were positive (Ignore)
  * `hash.ID` - Ignore
  * `tet_hash_id` - Output of running the new version of `HTOdemux` based on the `TET` assay (the `CLR` normalized raw tetramer counts). This gives the final label - the name of the tetramer if only one is above the cutoff or multi-reactive (Islet or Other) depending on what tetramers were above the cutoff. Defined in `02_Initial_processing.R`
  * `full_hash_id` - Output of running the new version of `HTOdemux` based on the `TET` assay (the `CLR` normalized raw tetramer counts). This gives all tetramers above the cutoff. Defined in `02_Initial_processing.R`
  * `libra_tet_hash_id` - The classification of tetramers based on the libra score - the name of the tetramer if only one is above the cutoff or multi-reactive (Islet or Other) depending on what tetramers were above the cutoff. Defined in `02_Initial_processing.R`. This uses the libra score based on the raw tetramer data.
  * `libra_full_hash_id` - Full id based on libra scores. This gives all tetramers above the cutoff. Defined in `02_Initial_processing.R`. This uses the libra scores computed based on the raw tetramer data.
  * `old_hash_id` - ignore
  * `SCAR_TET_maxID` - Output of running the new version of `HTOdemux` based on the `SCAR_TET` assay (the `CLR` normalized scar corrected tetramer counts). This gives the highest tet id. (Ignore)
  * `SCAR_TET_secondID` - Output of running the new version of `HTOdemux` based on the `SCAR_TET` assay (the `CLR` normalized scar corrected tetramer counts). This gives the second highest tet id. (Ignore)
  * `SCAR_TET_margin` - Output of running the new version of `HTOdemux` based on the `SCAR_TET` assay (the `CLR` normalized scar corrected tetramer counts). This gives the difference between the highest and second highest cutoffs (ignore)
  * `SCAR_TET_classification` - Output of running the new version of `HTOdemux` based on the `SCAR_TET` assay (the `CLR` normalized scar corrected tetramer counts). This gives the two highest tet ids. (ignore)
  * `SCAR_TET_classification.global`
  * `scar_hash_id` - Output of running the new version of `HTOdemux` based on the `SCAR_TET` assay (the `CLR` normalized scar corrected tetramer counts). This gives the values `Negative`, `Singlet` and `Doublet` depending on how many were positive (Ignore)
  * `full_scar_hash_id` - Output of running the new version of `HTOdemux` based on the `SCAR_TET` assay (the `CLR` normalized scar corrected tetramer counts). This gives all tetramers above the cutoff. Defined in `02_Initial_processing.R`
  * `scar_libra_tet_hash_id` - The classification of tetramers based on the libra score - the name of the tetramer if only one is above the cutoff or multi-reactive (Islet or Other) depending on what tetramers were above the cutoff. Defined in `02_Initial_processing.R`. This uses the libra score based on the scar corrected tetramer data.
  * `scar_libra_full_hash_id` - Full id based on libra scores. This gives all tetramers above the cutoff. Defined in `02_Initial_processing.R`. This uses the libra scores computed based on the scar corrected tetramer data.
  * `old_scar_hash_id` - ignore
  * `tet_name_cutoff` - Cutoff determined based on the non-b cells present for each sample. Here a cutoff was drawn ath the 95th quartile of the non-b cells. The value is the name of the tetramer if only one is above the cutoff or multi-reactive (Islet or Other) depending on what tetramers were above the cutoff. Defined in `12_improve_cutoff.R`. This uses the scar corrected tetramer data.
  * `full_tet_name_cutoff` - Cutoff determined based on the non-b cells present for each sample. Here a cutoff was drawn ath the 95th quartile of the non-b cells. This gives all tetramers above the cutoff. Defined in `12_improve_cutoff.R`. This uses the scar corrected tetramer data.
* Clustering and cell type columns
  * `seurat_clusters` - Final clusters called by the clustering algorithm (ignore)
  * `RNA_cluster` - The final clusters used for downstream analysis
  * `RNA_celltype_seurat` - Cell types determined by mapping the RNA clusters from individual samples to the [seurat reference](https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub). Download of the reference is [here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat)
  * `RNA_celltype_bnd` - Cell types determined by mapping the RNA clusters from individual samples to a [BND reference previously published by Mia Smith](https://rupress.org/jem/article/220/8/e20221604/214110/Identification-of-an-anergic-BND-cell-derived)
  * `RNA_celltype` - Final cell type determined by a combined best score from the BND and Seurat references based on the RNA clusters from individual samples.
  * `cluster_celltype` - A combination of the final cluster and final cell type
  * `adtdsb_clusters` - Clusters determined by using the adt dsb normalized data (ignore)
  * `adtclr_clusters` - Clusters determined by using the adt clr normalized data (ignore)
  * `adtscar_clusters` - Clusters determined by using the adt scar normalized data (ignore)
  * `tetdsb_clusters` - Clusters determined by using the tet dsb normalized data (ignore)
  * `tetclr_clusters` - Clusters determined by using the tet clr normalized data (ignore)
  * `tetscar_clusters` - Clusters determined by using the tet scar normalized data (ignore)
  * `tetscarlog_clusters` - Clusters determined by using the tet scar log normalized data (ignore)
  * `tetdsb_nd_clusters`
  * `grouped_celltype`
  * `celltype_cluster`
  * `rna_uncorrected_cluster` - Clusters determined before RNA correction
  * `rna_harmony_clust` - Clusters determined using the `harmony` dimensionality reduction
  * `rna_mnn_clust` - Clusters determined using the `mnn` dimensionality reduction
  * `AMBRNA_cluster` - Clusters determined using the ambient corrected RNA data (ignore)
  * `ambience_uncorrected_cluster` - Clusters determined before RNA correction using the ambient RNA corrected data
  * `ambience_harmony_clust` - Clusters determined using the `harmony` dimensionality reduction using the ambient RNA corrected data
  * `ambience_mnn_clust` - Clusters determined using the `mnn` dimensionality reduction using the ambient RNA corrected data
  * `rna_corrected_cluster` - Final clustering based on the `mnn` reduction
  * `ambience_corrected_cluster` - Final clustering using the abmient RNA based on the `mnn` reduction
  * `RNA_comb_celltype_seurat`- Cell types determined by mapping the RNA clusters from the combined samples to the [seurat reference](https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub). Download of the reference is [here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat)
  * `RNA_comb_celltype_bnd`- Cell types determined by mapping the RNA clusters from the combined samples to a [BND reference previously published by Mia Smith](https://rupress.org/jem/article/220/8/e20221604/214110/Identification-of-an-anergic-BND-cell-derived)
  * `RNA_combined_celltype` - Final cell type determined by a combined best score from the BND and Seurat references based on the RNA clusters from individual samples.
  * `AMBRNA_comb_celltype_seurat` - Combined cell type based on the ambience removed assay using the [seurat reference](https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub). Download of the reference is [here](https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat)
  * `AMBRNA_comb_celltype_bnd` - Combined cell type based on the ambience removed assaying using a [BND reference previously published by Mia Smith](https://rupress.org/jem/article/220/8/e20221604/214110/Identification-of-an-anergic-BND-cell-derived)
  * `AMBRNA_combined_celltype` - Combined cell type based on the ambience removed assay using both references.
  * `final_celltype` - Final cell types that were used for making the figures 
* Doublet finder
  * `Doublet_finder` - Doublet finder results
* immcantation columns
  * `final_clone` - Clone call by immcantation
  * `imcantation_isotype` - Isotype determined by immcantation