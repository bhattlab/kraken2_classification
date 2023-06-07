# Pipeline configuration
## If config file not specified, use the default in the current working directory.
if "outdir" not in config.keys():
    if not exists("config.yaml"):
        sys.exit("Specify configfile on the command line: --configfile config.yaml")
    else:
        configfile: "config.yaml"

# If running this downstream analysis only pipeline, sample_reads_file cannot be specified.
if config['sample_reads_file'] != '' and config['sample_reads_file'] is not None:
    sys.exit("sample_reads_file cannot be specified in the config if doing downstream processing only.")

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
# Set options depending on if the input has bracken reports or not
sample_reports, bracken_reports = get_sample_reports(config['sample_reports_file'])
sample_names = sample_reports.keys()
downstream_processing_input_kraken = [a[0] for a in sample_reports.values()]
downstream_processing_input_bracken = []
if bracken_reports != config['run_bracken']:
    sys.exit("Conflicting inputs: run_bracken in config file and bracken report listed in sample_reports")

# Determine extra output files if certain steps are defined in the config
extra_run_list =[]
## Add bracken outputs
if config['run_bracken']:
    downstream_processing_input_bracken = [a[1] for a in sample_reports.values()]
    extra_run_list.append('bracken_processed')

## Define the actual outputs. Some options unused currently.
extra_files = {
    "bracken_processed": join(outdir, 'processed_results_bracken/plots/classified_taxonomy_barplot_species.pdf'),
    "barplot": join(outdir, "plots/taxonomic_composition.pdf"),
}
run_extra_all_outputs = [extra_files[f] for f in extra_run_list]

# Set some resource requirements
bracken_memory = 8
bracken_threads = 1

# Taxonomic level can only be species right now. A future fix could look at the 
# output file name of Bracken and adjust based on taxonomic level.
if config['taxonomic_level'] != 'S':
    sys.exit('taxonomic_level setting can only be S')
