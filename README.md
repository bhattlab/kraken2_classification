# Kraken2 classification
Short read classification with the kraken2 program

## Introduction
[Kraken2](http://ccb.jhu.edu/software/kraken/) is a short read classification system that is fast and memory efficient. It allows you to assign a taxonomic identification to each read from a sequencing run. Kraken assigns each read to the lowest commen ancestor (LCA) of all sequences it alignes to. Through the use of the Bracken package, you can also get accurate estimates of proportions of different species. This guide will cover some of the basics, but the full [manual](http://ccb.jhu.edu/software/kraken/MANUAL.html) is quite good and has more detail.

## NEW AS OF 2019-09-01!
The outputs of this pipeline have been vastly improtved! Both internally and saved data now use the GCTx data format, from the [CMapR](https://github.com/cmap/cmapR) package. Basically, a GCT object is a data matrix that has associated row and column metadata. This allows for consistent metadata to live with the classification data, for both the rows (taxonomy information) and columns (sample metadata). See section [8. GCTx data processing](manual/gctx.md) for more information and tools for working with the new implementation. 

Also as of this update, the taxonomy information used by Kraken is fitered and improved some before saving any data or figures. For example, there were previously many taxonomy levels simply labeled "environmental samples" that are now named with their pared taxa name to remove ambiguity. Also, levels without a proper rank designation (listed with an abreviation and a number in the kraken report) have been forced into a specific rank when nothing was below them. This makes the taxonomy technically incorrect, but much more practically useful in these cases. Contact me with any questions. 

## Table of contents
1. [Installation](manual/installation.md)
2. [Usage](manual/usage.md)
3. [Available databases](manual/databases.md)
4. [Downstream processing and plotting](manual/downstream_plotting.md)
5. [Additional considerations](manual/extra.md)
6. [Expanded database construction](manual/db_construction.md)
7. [Using metagenome assembled genomes as a database](manual/mag_db.md)
8. [GCTx data parsing](manual/gctx.md)

## Quickstart
#### Install
If you're in the Bhatt lab, most of this work will take place on the SCG cluster. Otherwise, set this up on your own cluster or local machine (with a small database only). You will have to [build a database](manual/db_construction.md), create the `taxonomy_array.tsv` file (instructions at previous link) and set the database options if you're not in the Bhatt lab. Thanks to Singularity, all you need to have installed is snakemake. See the instructions [here](https://github.com/bhattlab/bhattlab_workflows/) to set up snakemake and set up a profile to submit jobs to the cluster. 

Then clone the repo wherever is convenient for you. I use a directory in `~/projects`
```
cd ~/projects
git clone https://github.com/bhattlab/kraken2_classification.git
```
#### Run
Copy the `config.yaml` file into the working directory for your samples. Change the options to suit your projects and make sure you specify the right `samples.tsv` file. See [Usage](manual/usage.md) for more detail. You can then lauch the workflow with a snakemake command like so:
```
# Snakemake workflow - change options in config.yaml first
snakemake -s path/to/Snakefile --configfile config.yaml --use-singularity --singularity-args '--bind /oak/,/labs/,/home' --profile scg --jobs 99
```

## Parsing output reports
The Kraken reports `classification/sample.krak.report`, bracken reports `sample.krak_bracken.report`, and data matrices in the `processed_results` folder are the best for downstream analysis. See [Downstream processing and plotting](manual/downstream_plotting.md) for details on using the data in R. 

## Example output
The pipeline outputs data, results and figures in the structure below. 
```
 | classification 
   - sample.krak                    Kraken results - classification of each read. These files
                                    can get very large and are unnecessary if you only want the reports. 
   - sample.krak.report             Kraken report - reads and percentages at each taxonomic level.
   - sample.krak.report.bracken     Standard bracken report at species level. Not useful, use the one below.
   - sample.krak_bracken.report     Most useful format of the the Bracken results.

 | processed_results
   - diversity.txt                  Diversity of each sample at each taxonomic level

   | ALDEX2_differential_abuncance  Compositional data analysis done with the ALDEx2 package.
                                    Only done if you have 2 groups in the sample_groups file.
     - aldex_result_[].tsv          Differential abundance at the given taxonomic level. 
     - aldex_scatter_[].pdf         Scatterplot of effect vs dispersion with significant hits highlighted
     - aldex_significant_boxplots_[].pdf  Boxplot of any significant hits

   | braycurtis_matrices
     - bravcurtis_distance_[].tsv   Matrix of braycurtis distance between samples at each taxonomic level
   
   | plots                          Lots of plots!
      - classified_taxonomy_barplot_[].pdf   Barplot at each taxonomic level.
      - compositional_PCA_plot.pdf  PCA done on clr values
      - diversity_allsamples.pdf    Diversity barplot
      - diversity_by_group.pdf      Diversity barplot, stratified by sample group
      - PCoA_2D_plot_labels.pdf     Principal coordinates analysis, calculated on braycurtis distances
      - PCoA_2D_plot_nolabels.pdf   Same as above without the point labels
      - rarefaction_curve.pdf       Rarefaction "collectors" curve plot

   | taxonomy_gctx_classified_only  
      - bracken_[]_reads.gctx       GCTx file with matrix containing reads classified at each level 
      - bracken_[]_percentage.gctx  Same, but percentage of classified reads

   | taxonomy_matrices_classified_only
      - bracken_[]_reads.txt        Matrix of taxa by sample classified reads
      - bracken_[]_percentage.txt   Same, but percentage of classified reads
      - clr_values_[].txt           From compositional data analysis, centered log-ratio values at each level

 | processed_results_krakenonly
   | Same as above, but using the results without Bracken. Also includes taxonomy matrices
   | that have unclassified reads in them (as bracken no longer reports unclassified reads)
 
 | unmapped reads
   - sample_unmapped_1.fq           Only present if selected in the config file; reads that are not 
   - sample_unmapped_2.fq           classified, as paired end fastq.
```

_Taxonomic barplot_
![default_barplot](images/taxonomic_barplot.png "Default taxonomic barplot output")

_PCoA plot example_
![pcoa_plot](images/pcoa_plot.png "PCoA plot")
