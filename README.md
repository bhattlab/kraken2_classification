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
*Install*
If you're in the Bhatt lab, most of this work will take place on the SCG cluster. Otherwise, set this up on your own cluster or local machine (with a small database only). You will have to [build a database](manual/db_construction.md) and set the database options if you're not in the Bhatt lab. Thanks to Singularity, all you need to have installed is snakemake. See the instructions [here](https://github.com/bhattlab/bhattlab_workflows/) to set up snakemake and set up a profile to submit jobs to the cluster. 

Then clone the repo wherever is convenient for you. I use a directory in `~/projects`
```
cd ~/projects
git clone https://github.com/bhattlab/kraken2_classification.git
```
*Run*
Copy the `config.yaml` file into the working directory for your samples. Change the options to suit your projects and make sure you specify the right `samples.tsv` file. See [Usage](manual/usage.md) for more detail. You can then lauch the workflow with a snakemake command like so:
```
# Snakemake workflow - change options in config.yaml first
snakemake -s path/to/Snakefile --configfile config.yaml --use-singularity --singularity-args '--bind /labs/ --bind /scratch --bind /home/' --profile scg --jobs 99
```

## Parsing output reports
The Kraken reports `classification/sample.krak.report`, bracken reports `sample.krak_bracken.report`, and data matrices in the `processed_results` folder are the best for downstream analysis. See [Downstream processing and plotting](manual/downstream_plotting.md) for details on using the data in R. 
