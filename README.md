# kraken2_classification
Short read classification with the kraken2 program

## Introduction
[Kraken2](http://ccb.jhu.edu/software/kraken/) is a short read classification system that is fast and memory efficient. It allows you to assign a taxonomic identification to each read from a sequencing run. Kraken assigns each read to the lowest commen ancestor (LCA) of all sequences it alignes to. Through the use of the Bracken package, you can also get accurate estimates of proportions of different species. This guide will cover some of the basics, but the full [manual](http://ccb.jhu.edu/software/kraken/MANUAL.html) is quite good and has more detail.

## Installation
```
conda env create -f classification2.yaml
```

## Usage
Run the following in an interactive session, or submit it as a job.  Do not add `--profile scg` as is typically done with other workflows.
```
source activate classification2
snakemake -s path/to/Snakefile --configfile config.yaml
```
The database stays in memory after loading, so it's very quick to run many samples consecutively. This is why it is faster to run the workflow on a single node than parallelizing.

## Datasets table
The input data are specified in the form of a tab-delimited table taking the below form:
```
Sample	Timepoint	Condition	Reads1.fq[.gz][,Reads2.fq[.gz]]
a	1	control	a_1.fq,a_2.fq
b	2	case	b_1.fq,b_2.fq
c	1	case	c_1.fq.gz,c_2.fq.gz
d	2	case	d_1.fq.gz,d_2.fq.gz
```

This yields a barplot that looks like this:

![example barplot](images/taxonomic_composition.png "barplot!")

### Snakemake workflow
A workflow and configuration file are provided in this repo. By default, samples are classified using the high quality genbank database. Edit the config.yaml file to change this behavior. Currently this depends on a list of sample names in the configuration file, and read files to be named following the sample name convention. 


### Memory usage
Memory usage depends on the database. For a rough estimate of memory requirements for a given db, check the size of the `hash.k2d` file, which has to fit completely in memory. The manual recommends at least 175GB for the standard database. In my experiene, usage maxes out at the following values:

- Genbank high quality: 64g
- Genbank all: 256g

You should request at least this much memory on a node you plan to do classification on. 

### Time usage
Kraken2 is very fast compared to kraken1 or other classification tools. After loading the database in memory the first time (about 5 minutes for genbank hq), subsequet runs proceed very quickly. Using 8 cores on SCG I had the following results with the genbank hq database: 
`10399748 sequences (2356.32 Mbp) processed in 85.562s (7292.8 Kseq/m, 1652.36 Mbp/m).`
That's 85 seconds for a dataset with about 10 million read pairs. Obviously more cores can make things faster for larger datasets.

## Classification percentages
I've evauluated some of the databases below on sets different metagenomic datasets, and compared the unclassified percentages to our old Kraken1 custom database. See [this spreadsheet](https://docs.google.com/spreadsheets/d/15nVMno4w4Q-DVO9tdp1DBpwKslOkLb2TVqnSnuEofTY/edit?usp=sharing) for a summary of results. Overall, Kraken2 with the high quality database classifies less reads than Kraken1 but much quicker. The full genbank database has increased classification percentages and speed.

## Available databases
Specify the right database for your classification needs. To see if your organism of interest is present in a database (and therefore is able to be classified in your reads), search the `inspect.out` file in the database folder. If this has the name or taxonomic identifier you're interested in, you're good to go! In this example we see that crassphage is indeed present in the DB:
```
$ grep -i crassphage /labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria/inspect.out 
>  0.00  13159   13159   S       1211417         uncultured crAssphage

```
#### Genbank bacteria (high quality assemblies)
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria/`

Includes all bacterial sequences from genbank that were assembled to "Complete Genome" or "Chromosome" status as of October 2018. Will miss organisms that don't have high-quality assemblies in genbank. Very fast and low memory for simple bacterial classification. 64g memory is sufficient for this database.

Has had Prevotella copri and crAssphage sequences manually added. 

**statistics**
- 13653 sequences went into the construction of this database
- 4053 species are present
- 8383 species, subspecies or strains are present

Check the `inspect.out` file for a closer look. 

#### Genbank bacteria and archaea (ALL assemblies)
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria_complete/`

Includes all bacterial sequences from genbank assembled to any quality (includes "Scaffold" and "Contig" level assemblies) as of October 2018. 256g memory is necessary for this database, and runtimes will be longer. I've found results from this database to be more noisy, but with higher classification percentages. There is a long tail of very low abundance species present. Running Bracken is essential when using this database.

Has had Prevotella copri and crAssphage sequences manually added. 

**statistics**
- 167229 sequences went into the construction of this database
- 38789 species are present
- 63871 species, subspecies or strains are present

Check the `inspect.out` file for a closer look. 

#### Standard (high quality refseq assemblies)
`/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/standard/`

Contains high-quality refseq assemblies of archaea, bacteria, human, UniVec_Core and viral sequences.

**statistics**
- 20934 sequences went into the construction of this database
- 11578 species are present
- 17506 species, subspecies or strains are present

Check the `inspect.out` file for a closer look. 

#### Standard protein (high quality refseq assemblies)
`/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/standard_protein/`

Protein database of high-quality refseq assemblies of archaea, bacteria, human, UniVec_Core and viral sequences.

#### nt, env_nt
The equivalent of the blast nt and env_nt databases can be found at 
`/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/nt`

`/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/env_nt`

I'm not sure how useful these are and haven't tried them in any real classification tasks.

#### Custom sequences
Don't see your favorite bug in the `inspect.out` file? Have a newly assembled organism you want added to the database? It's an easy process to put custom sequences into the database and can be done in under a day of processing time. Contact Ben for requests (I'll batch these to once a month or something if there are lots of requests.)


## Downstream processing
### Abundance estimation with Bracken
I highly recommend estimating abundances with Bracken! It makes the reports much more understandable. 

[Bracken publication](https://peerj.com/articles/cs-104/)

Due to Kraken's LCA reporting, clades with many similar species will only have species-level assignments for unique regions, leaving most reads "stranded" above the species level. The number of reads classified directly to a species may be far lower than the actual number present. Therefore, any assumption that Kraken’s raw read assignments can be directly translated into species- or strain-level abundance estimates is flawed, as ignoring reads at higher levels of the taxonomy will grossly underestimate some species, and creates the erroneous impression that Kraken’s assignments themselves were incorrect.Bracken (Bayesian Reestimation of Abundance after Classification with KrakEN) estimates species abundances in metagenomics samples by probabilistically re-distributing reads in the taxonomic tree. 
The [manual](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual) has an explanation of the command line options.

#### Databases
I've built a database for 100 & 150bp reads in each of the db folders mentioned above. Just specify the same database you used for classification in the bracken command. If you have different read lengths, contact Ben to build a DB or use the `bracken-build` command (use 32 cores for speed).

### Visualization with Krona
I like Krona for visualization of individual samples. It's a great way to explore the taxonomic identifications and proportions. Install the KronaTools package on your local machine from the [github repo](https://github.com/marbl/Krona/wiki/KronaTools). I then call the _ImportTaxonomy_ script on the Kraken report. Multiple reports can be specified on the command line to generate an html file with a tab for each input. Bracken reports can be used as an input here as well. 
```
~/software/KronaTools-2.7/scripts/ImportTaxonomy.pl -m 3 -s 0 -q 0 -t 5 -i krak.report -o out.html
```
![Krona_classification](images/krona_classification.png "Using Krona is a great way to explore classification at many levels.")

### Sample comparison with Pavian
[Pavian](http://ccb.jhu.edu/software/pavian/) is a web application fo comparing kraken classification results. I find it useful to visualize results across many samples. This application should be installed on your local machine in it's own conda environment. The [github repo](https://github.com/fbreitwieser/pavian) has installation instructions. I first have to run this option in R to work, then start the app: 
```
options(browser="google-chrome")
pavian::runApp()
```
Bracken reports can be used here as an input as well. Pavian also makes nice classification flowcharts, see below for an example:

![Pavian_classification](images/pavian_classification_flow.png "Using pavian can make nice classification flowcharts")

### Metaphlan2 style heatmaps
If you want to make a heatmap like in those made in metaphlan2, you can use the `--use-mpa-style` option to get a compatable report. Normalize to the total number of classified reads at the domain level:
```
sum=$(grep -vP "\|" out.mpareport | cut -f 2 | awk '{sum += $1} END {printf ("%.2f\n", sum/100)}')
awk -v sum="$sum" 'BEGIN {FS="\t"} {OFS="\t"} {print $1,$2/sum}' out.mpareport > out.mpareport.norm

```
Then reports can be merged and plotted just like in metaphlan2 (you might need to remove the -c viridis option here):
```
merge_metaphlan_tables.py *.norm >  merge_metaphlan.txt
grep -E "(s__)|(^ID)"  merge_metaphlan.txt | grep -v "t__" | sed 's/^.*s__//g' >  merge_metaphlan_species.txt
metaphlan_hclust_heatmap.py --in merge_metaphlan.txt  --top 25 --minv 0.1 -s log --out merge_metaphlan_heatmap.png -f braycurtis -d braycurtis -c viridis 
metaphlan_hclust_heatmap.py --in merge_metaphlan.txt  --top 150 --minv 0.1 -s log --out merge_metaphlan_heatmap_big -f braycurtis -d braycurtis -c viridis
```
![metaphlan2_heatmap](images/metaphlan2_heatmap.png "Using metaphlan2-style heatmaps from Kraken2 data.")
