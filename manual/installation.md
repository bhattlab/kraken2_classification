## Installation 
To install kraken2 on SCG (or any other machine), you should create a conda environment with all the requirements. See the [scg_tools](https://github.com/bhattlab/scg_tools#setting-up-your-environment-with-conda) repo for more information. 

Clone this repository, then use the commands below to install. You must first `cd` to the location where you cloned the repository.

It is also tasteful to copy datasets.tsv and config.yaml to a working directory for modification and workflow execution.
```
# crete an environment classification2 with all the requirements
conda env create -f classification2.yaml
# update taxonomy used with Krona plotting software
# not necessary for just classification
mkdir $(which kraken2 | sed 's/envs\/classification2.*$//g')/envs/classification2/bin/taxonomy
ktUpdateTaxonomy.sh 
```
