# Pipeline configuration
## If config file not specified, use the default in the current working directory.
if "outdir" not in config.keys():
    if not exists("config.yaml"):
        sys.exit("Specify configfile on the command line: --configfile config.yaml")
    else:
        configfile: "config.yaml"

# If running this whole pipeline, sample_reports_file cannot be specified.
if config['sample_reports_file'] != '' and config['sample_reports_file'] is not None:
    sys.exit("sample_reports_file cannot be specified in the config if running the whole pipeline.")

# output base directory
outdir = config['outdir']

# Set backwards compatability and default options
## Set remove_chordata option
if not "remove_chordata" in config:
    config['remove_chordata'] = 'FALSE'
## Set desired confidence threshold for Kraken
if not 'confidence_threshold' in config:
    config['confidence_threshold'] = 0.0
confidence_threshold = config['confidence_threshold']

# Read in sample names and sequencing files from sample_reads_file
# Set options depending on if the input is paired-end or not
sample_reads, paired_end = get_sample_reads(config['sample_reads_file'])
if paired_end:
    paired_string = '--paired'
else:
    paired_string = ''
sample_names = sample_reads.keys()
downstream_processing_input_kraken = expand(join(outdir, "classification/{samp}.krak.report"), samp=sample_names)
downstream_processing_input_bracken = []

# Determine extra output files if certain steps are defined in the config
extra_run_list =[]
## Add bracken outputs
if config['run_bracken']:
    extra_run_list.append('bracken')
    extra_run_list.append('bracken_processed')
    downstream_processing_input_bracken = expand(join(outdir, "classification/{samp}.krak_bracken_species.report"), samp=sample_names)

## Add unmapped read extraction
if config['extract_unmapped']:
    if paired_end:
        extra_run_list.append('unmapped_paired')
    else:
        extra_run_list.append('unmapped_single')
## Define the actual outputs. Some options unused currently.
extra_files = {
    "bracken": expand(join(outdir, "classification/{samp}.krak_bracken_species.report"), samp=sample_names),
    "bracken_processed": join(outdir, 'processed_results_bracken/plots/classified_taxonomy_barplot_species.pdf'),
    "unmapped_paired": expand(join(outdir, "unmapped_reads/{samp}_unmapped_1.fq"), samp=sample_names),
    "unmapped_single": expand(join(outdir, "unmapped_reads/{samp}_unmapped.fq"), samp=sample_names),
    "barplot": join(outdir, "plots/taxonomic_composition.pdf"),
    "krona": expand(join(outdir, "krona/{samp}.html"), samp = sample_names),
    "mpa_heatmap": join(outdir, "mpa_reports/merge_metaphlan_heatmap.png"),
    "biom_file": join(outdir, "table.biom"),
}
run_extra_all_outputs = [extra_files[f] for f in extra_run_list]
# print("run Extra files: " + str(run_extra_all_outputs))

## set some resource requirements
if config['database'] in ['/labs/asbhatt/data/program_indices/kraken2/kraken_custom_feb2019/genbank_genome_chromosome_scaffold',
                          '/labs/asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/genbank_genome_chromosome_scaffold',
                          '/labs/asbhatt/data/program_indices/kraken2/kraken_custom_dec2021/genbank_genome_chromosome_scaffold',
                          '/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_feb2019/genbank_genome_chromosome_scaffold',
                          '/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/genbank_genome_chromosome_scaffold',
                          '/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_dec2021/genbank_genome_chromosome_scaffold']:
    kraken_memory = 256
    kraken_threads = 8
    bracken_memory = 64
    bracken_threads = 1
elif config['database'] == 'test_data/db/':
    kraken_memory = 8
    kraken_threads = 4
    bracken_memory = 8
    bracken_threads = 1
else:
    kraken_memory = 64
    kraken_threads = 4
    bracken_memory = 64
    bracken_threads = 1

# Taxonomic level can only be species right now. A future fix could look at the 
# output file name of Bracken and adjust based on taxonomic level.
if config['taxonomic_level'] != 'S':
    sys.exit('taxonomic_level setting can only be S')

# Set config['codata_min_reads'] to lower than the default if working with test data
if config['codata_min_reads'] > 1000 \
    and "M0018C_2_month" in sample_reads \
    and sample_reads["M0018C_2_month"][0] == "test_data/SRR7281046_sub_1.fastq.gz":
    print("Changing config['codata_min_reads'] to 1000 based on test data ")
    config['codata_min_reads'] = 1000