## Usage
First, modify the config.yaml file to suit your needs - you can change the database, output file and `samples.tsv` file. `samples.tsv` must be a tab-delimited file containing the sample name and paths to forward and reverse reads. If you're using the Bhatt lab preprocessing workflow, this file is outputted at `preprocessing/01_processing/classification_input.txt ` Example:

```
sample_1_name	/path/to/s1_read_1.fq	/path/to/s1_read_2.fq
sample_2_name	/path/to/s2_read_1.fq	/path/to/s2_read_2.fq
sample_3_name   /path/to/s3_read_1.fq   /path/to/s3_read_2.fq
sample_4_name	/path/to/s4_read_1.fq	/path/to/s4_read_2.fq
```


### In the Bhatt lab
Usage will be slightly different if you're in the Bhatt lab or not. In the Bhatt lab, this pipeline should be run with jobs submitted to the SCG SLURM cluster. Instructions for configuring tools for this can be found at our [bhattlab_workflows repository](https://github.com/bhattlab/bhattlab_workflows/blob/master/manual/setup.md). Then, a command like this can be used to run the workflow with up to 99 concurrent jobs submitted to the SLURM scheduler. Note the use of bind arguments to ensure all the different filesystems play together nicely. 
```
snakemake -s path/to/Snakefile --configfile config.yaml --use-singularity --singularity-args '--bind /oak/,/labs/,/home' --profile scg --jobs 99
```

After running the workflow and you're satisfied the results, run the cleanup command to remove temporary files that are not needed anymore. 
```
snakemake cleanup -s path/to/Snakefile --configfile config.yaml
```
### In other settings
You should run this pipeline in a setting with enough RAM for your database of choice. More CPU cores will also speedup processing large sequencing datasets. You could setup a [snakemake profile](https://github.com/Snakemake-Profiles/slurm) for submission to a SLURM cluster if desired. Running the pipeline is similar to the above, but you might need to add singularity bind arguments or a profile for SLURM job submission depending on your configuration. This example uses 8 cores, but that can be changed to reflect available resources.
```
snakemake -s path/to/Snakefile --configfile config.yaml --use-singularity --jobs 8 --cores 8
```

The cleanup rule removes extra `.krak` files which take up lots of space (one line per sequencing read). 
```
snakemake cleanup -s path/to/Snakefile --configfile config.yaml
```

### Downstream processing
Kraken/Bracken reports are processed into matrices containing read counts classified at the genus and species level, as well as relative proportions. Several plots are produced with taxonomic barplots and diversity metrics. Samples can be split into separate groups for a patient, treatment group, etc. To do so, specify a tab delimited file with sample names and groups. Sample names must match the `samples.tsv` above. By default this is called `sample_groups.tsv` Example:
```
sample_1_name   group_1
sample_2_name   group_1
sample_3_name   group_2
sample_4_name   group_2
```

If no groups file is specified, all samples will be treated as a single group.

For more advanced options, see [Downstream processing and plotting](manual/downstream_plotting.md).