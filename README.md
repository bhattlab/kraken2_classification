# Kraken2 classification
Short read classification with the kraken2 program

## Introduction
[Kraken2](http://ccb.jhu.edu/software/kraken/) is a short read classification system that is fast and memory efficient. It allows you to assign a taxonomic identification to each read from a sequencing run. Kraken assigns each read to the lowest commen ancestor (LCA) of all sequences it alignes to. Through the use of the Bracken package, you can also get accurate estimates of proportions of different species. This guide will cover some of the basics, but the full [manual](http://ccb.jhu.edu/software/kraken/MANUAL.html) is quite good and has more detail.

## Table of contents
1. [Installation](manual/installation.md)
2. [Usage](manual/usage.md)
3. [Available databases](manual/databases.md)
4. [Downstream processing and plotting](manual/downstream_plotting.md)
5. [Additional considerations](manual/extra.md)
6. [Custom database construction](manual/db_construction.md)

## Quickstart
*Install*
Most of this work will take place on SCG. Thanks to Singularity, all you need to have installed is snakemake. See the instructions [here](https://github.com/bhattlab/bhattlab_workflows/) to set up snakemake and set up a profile to submit jobs to the cluster. 

Then clone the repo wherever is convenient for you. I use a directory in `~/projects`
```
cd ~/projects
git clone git@github.com:bhattlab/kraken2_classification.git
```
*Run*
Copy the `config.yaml` file into the working directory for your samples. Change the options to suit your projects and make sure you specify the right `samples.tsv` file. See [Usage](manual/usage.md) for more detail. You can then lauch the workflow with a snakemake command like so:
```
# Snakemake workflow - change options in config.yaml first
snakemake -s path/to/Snakefile --configfile config.yaml --use-singularity --singularity-args '--bind /labs/ --bind /scratch --bind /home/' --profile scg --jobs 99
```

## Parsing output reports
The Kraken reports `classification/sample.krak.report`, bracken reports `sample.krak_bracken.report`, and data matrices in the `processed_results` folder are the best for downstream analysis. See [Downstream processing and plotting](manual/downstream_plotting.md) for details on using the data in R. 
