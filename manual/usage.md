## Usage
Modify the config.yaml file to suit your needs - you can change the database, output file and samples.tsv file. `samples.tsv` must be a tab-delimited file containing the sample name and paths to forward and reverse reads. If you're using the preprocessing workflow, this file is outputted at `preprocessing/01_processing/classification_input.txt ` Example:

```
sample_1_name	/path/to/s1_read_1.fq	/path/to/s1_read_2.fq
sample_2_name	/path/to/s2_read_1.fq	/path/to/s2_read_2.fq
sample_3_name   /path/to/s3_read_1.fq   /path/to/s3_read_2.fq
sample_4_name	/path/to/s4_read_1.fq	/path/to/s4_read_2.fq
```

Run the following in an interactive session, or submit it as a job.  Do not add `--profile scg` as is typically done with other workflows. By default, samples are classified using the high quality genbank database (see Databases section below).

```
source activate classification2
snakemake -s path/to/Snakefile --configfile config.yaml
```
The database stays in memory after loading, so it's very quick to run many samples consecutively. This is why it is faster to run the workflow on a single node than parallelizing.

_Refresher_ To get an interactive session, run something like this (change cores and time to suit you needs):
```
srun -t 1:00:00 -p interactive  -n 1 -c 8 --mem=64000 --pty bash
```

### Downstream processing
Kraken/Bracken reports are processed into matrices containing read counts classified at the genus and species level, as well as relative proportions. Several plots are produced with taxonomic barplots and diversity metrics. You must provide the location of the `scripts` directory from this repository, or a symlinked folder for downstream processing to work. Samples can be split into separate groups for a patient, treatment group, etc. To do so, specify a tab delimited file with sample names and groups. Sample names must match the `samples.tsv` above. By default this is called `sample_groups.tsv` Example:
```
sample_1_name   group_1
sample_2_name   group_1
sample_3_name   group_2
sample_4_name   group_2
```

If no groups file is specified, all samples will be treated as a single group.

For more advanced options, see [Downstream processing and plotting](manual/downstream_plotting.md).