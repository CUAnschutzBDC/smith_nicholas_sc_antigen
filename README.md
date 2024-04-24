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
    * [Meta data columns](#meta-data-columns)

<!--te-->

## TODO
* Clean up metadata/seurat object
  * Talk to Mia and Catherine to determine what needs to be removed from the metadata before we can publish
  * Remove all doublet finder columns (`PANN...`)
* Upload docker containers


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
  * Runs `HTOdemux` on the tetramers
  * Identifies positive values for each tetramer based on a libra cutoff of 1 for both the raw and `scar` corrected libra scores.
  * Runs SCT normalization on the RNA assay
* `03_remove_doublets.R` - Identifies and tags doublets using [`DoubletFinder`](https://github.com/chris-mcginnis-ucsf/DoubletFinder). Images related to this can be found in `results/R_analysis/{sample_name}/images/dropkick_vs_cellqc.pdf`.
* `04_adt_dsb_normalization.R` - Uses [`dsb`](https://github.com/niaid/dsb) to identify and remove background. This is similar to `scar` and was only used to compare to `scar` output.
* `05_PCA.R` - Runs PCA using methods from `Seurat`. Makes the output `results/R_analysis/{sample_name}/images/RNA_pca.pdf` to help determine the number of pcs to use for downstream processing.
* `06a_UMAP_find_resolution.R` - Uses the number of pcs found above (specified in the file `files/sample_info.tsv`) and uses `clustree` to visualize many clustering resolutions. The `clustree` and umaps of each resolution will be in a file `results/R_analysis/{sample_name}/images/clustering_resolution.pdf` to guide resolution selection for downstream analysis.
* `06b_UMAP.R` - Taking the number of pcs and resolution determined in the previous steps (specified in `files/sample_info.tsv`), generates the final umap and clustering for the object. The final images can be found in `results/R_analysis/{sample_name}/images/RNA_pca.pdf`
* 


*NOTE on tetramer labeling. While libra, HTODemux, scar, and raw data were all used, only the libra score on ______ and the t-cell/myeloid cell cutoff on ______ were used to identify tetramer names the others were just used to test methods and compare results*

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


### Meta data columns

* `Seurat` default columns
  * `orig.ident` - Sample name
  * `nCount_RNA` - Number of RNA molecules in the cell
  * `nFeature_RNA` - Number of genes in the cell
  * `nCount_ADT` - Number of ADT molecules in the cell
  * `nFeature_ADT` - Number of different ADTs in the cell
  * `nCount_SCT` - Number of reads based on SCT normalization
  * `nFeature_SCT` - Number of genes based on SCT normalization 
  * `nCount_TET`
  * `nFeature_TET`
  * `nCount_TET_LIBRA`
  * `nFeature_TET_LIBRA`
  * `nCount_SCAR_ADT_LOG`
  * `nFeature_SCAR_ADT_LOG`
  * `nCount_SCAR_ADT`
  * `nFeature_SCAR_ADT`
  * `nCount_SCAR_TET`
  * `nFeature_SCAR_TET`
  * `nCount_SCAR_TET_LOG`
  * `nFeature_SCAR_TET_LOG`
  * `nCount_SCAR_TET_LIBRA`
  * `nFeature_SCAR_TET_LIBRA`
  * `nCount_CLR_ADT`
  * `nFeature_CLR_ADT`
  * `nCount_CLR_TET`
  * `nFeature_CLR_TET`
  * `nCount_AMBRNA`
  * `nFeature_AMBRNA`
  * `nCount_TET_PROPORTIONS`
  * `nFeature_TET_PROPORTIONS`
  * `nCount_SCAR_TET_PROPORTIONS`
  * `nFeature_SCAR_TET_PROPORTIONS`
  * `nCount_DSB_ADT`
  * `nFeature_DSB_ADT`
  * `nCount_DSB_TET`
  * `nFeature_DSB_TET`
  * `nCount_NEW_TET_PROPORTIONS`
  * `nFeature_NEW_TET_PROPORTIONS`
  * `percent.mt` - Percent of reads mapping to the mitochondria
* Meta data columns
  * `ID`
  * `Initials`
  * `Sample.Name`
  * `Sex`
  * `Collection.Date`
  * `Age.at.Collection..years.`
  * `Age.at.Collection..Months.`
  * `Status`
  * `Date.of.Diagnosis`
  * `Days.post.onset`
  * `millions.of.cells.frozen`
  * `HLA.type`
  * `Autoantibodies`
  * `date.processed.for.scSeq`
  * `Notes..FDR.relationship.`
  * `old_status`
  * `sample`
* `scuttle` qc columns
  * `cell_qc_sum`
  * `cell_qc_detected`
  * `cell_qc_subsets_Mito_sum`
  * `cell_qc_subsets_Mito_detected`
  * `cell_qc_subsets_Mito_percent`
  * `cell_qc_altexps_ADT_sum`
  * `cell_qc_altexps_ADT_detected`
  * `cell_qc_altexps_ADT_percent`
  * `cell_qc_total`
  * `cell_qc_low_lib_size`
  * `cell_qc_low_n_features`
  * `cell_qc_high_subsets_Mito_percent`
  * `cell_qc_discard`
* Cell cycle columns
 * `S.Score`
 * `G2M.Score`
 * `Phase`
* VDJ columns
  * `clonotype_id`
  * `exact_subclonotype_id`
  * `chains`
  * `n_chains`
  * `cdr3`
  * `cdr3_nt`
  * `cdr3_length`
  * `cdr3_nt_length`
  * `v_gene`
  * `d_gene`
  * `j_gene`
  * `c_gene`
  * `isotype`
  * `reads`
  * `umis`
  * `productive`
  * `full_length`
  * `paired`
  * `v_ins`
  * `v_del`
  * `v_mis`
  * `d_ins`
  * `d_del`
  * `d_mis`
  * `j_ins`
  * `j_del`
  * `j_mis`
  * `c_ins`
  * `c_del`
  * `c_mis`
  * `all_ins`
  * `all_del`
  * `all_mis`
  * `vd_ins`
  * `vd_del`
  * `dj_ins`
  * `dj_del`
  * `v_mis_freq`
  * `d_mis_freq`
  * `j_mis_freq`
  * `c_mis_freq`
  * `all_mis_freq`
* Columns relating to tetramer calling
  * `tet_hash_id` - The tetramer id defined by `HTODemux` (highest score). Defined in `02_Initial_processing.R`
  * `full_hash_id` - All tetramers defined as positive by `HTODemux`. Defined in `02_Initial_processing.R`
  * `libra_tet_hash_id` - The tetramer id defined by the libra score. If more than 1 tetramer had a value higher than 1. Defined in `02_Initial_processing.R`
 * `TET_maxID`
 * `TET_secondID`
 * `TET_margin`
 * `TET_classification`
 * `TET_classification.global`
 * `hash.ID`
 * `tet_hash_id`
 * `full_hash_id`
 * `libra_tet_hash_id`
 * `libra_full_hash_id`
 * `old_hash_id`
 * `SCAR_TET_maxID`
 * `SCAR_TET_secondID`
 * `SCAR_TET_margin`
 * `SCAR_TET_classification`
 * `SCAR_TET_classification.global`
 * `scar_hash_id`
 * `full_scar_hash_id`
 * `scar_libra_tet_hash_id`
 * `scar_libra_full_hash_id`
 * `old_scar_hash_id`
* Clustering and cell type columns
  * `seurat_clusters`
  * `RNA_cluster`
  * `RNA_celltype_seurat`
  * `RNA_celltype_bnd`
  * `RNA_celltype`
  * `cluster_celltype`
  * `adtdsb_clusters`
  * `adtclr_clusters`
  * `adtscar_clusters`
  * `tetdsb_clusters`
  * `tetclr_clusters`
  * `tetscar_clusters`
  * `tetscarlog_clusters`
  * `tetdsb_nd_clusters`
  * `tet_name_cutoff`
  * `full_tet_name_cutoff`
  * `grouped_celltype`
  * `celltype_cluster`
  * `rna_uncorrected_cluster`
  * `rna_harmony_clust`
  * `rna_mnn_clust`
  * `AMBRNA_cluster`
  * `ambience_uncorrected_cluster`
  * `ambience_harmony_clust`
  * `ambience_mnn_clust`
  * `rna_corrected_cluster`
  * `ambience_corrected_cluster`
  * `RNA_comb_celltype_seurat`
  * `RNA_comb_celltype_bnd`
  * `RNA_combined_celltype`
  * `AMBRNA_comb_celltype_seurat`
  * `AMBRNA_comb_celltype_bnd`
  * `AMBRNA_combined_celltype`
  * `final_celltype`
  * `AMBRNA_snn_res.0.6`
* Doublet finder
  * `Doublet_finder`   
* immcantation columns
  * `final_clone`
  * `imcantation_isotype`            