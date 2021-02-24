## Usage
Modify the config.yaml file to suit your needs - you can change the database, output file and samples.tsv file. `samples.tsv` must be a tab-delimited file containing the sample name and paths to forward and reverse reads. If you're using the preprocessing workflow, this file is outputted at `preprocessing/01_processing/classification_input.txt ` Example:

```
sample_1_name	/path/to/s1_read_1.fq	/path/to/s1_read_2.fq
sample_2_name	/path/to/s2_read_1.fq	/path/to/s2_read_2.fq
sample_3_name   /path/to/s3_read_1.fq   /path/to/s3_read_2.fq
sample_4_name	/path/to/s4_read_1.fq	/path/to/s4_read_2.fq
```

To run Kraken2, either get an interactive session on the cluster, or submit batch jobs by adding the `--profile scg  --jobs 999` flags to the commands below. 

```
# Snakemake workflow - change options in config.yaml first
source activate classification2
snakemake -s path/to/Snakefile --configfile config.yaml --use-singularity --singularity-args '--bind /labs/ --bind /scratch --bind /home/' --profile scg --jobs 99

# Or if you want to run on the command line and have kraken/bracken installed
# get an interactive session first. 256Gb mem necessary for this database
DB=/labs/asbhatt/data/program_indices/kraken2/kraken_custom_feb2019/genbank_genome_chromosome_scaffold
threads=8
read_length=150   # must be 150 or 100
kraken2 --db "$DB" --threads "$threads" --output out.krak --report out.krak.report --paired r1.fq r2.fq
bracken -d "$DB" -t "$threads" -i out.krak.report -o out out.krak.report.bracken -r "$read_length"
```

_Refresher_ To get an interactive session, run something like this (change cores and time to suit you needs):
```
srun -t 1:00:00 -p interactive  -n 1 -c 8 --mem=256000 --pty bash
```

After running the workflow and you're satisfied the results, run the cleanup command to remove temporary files that are not needed anymore. 
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