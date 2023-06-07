# Run this script from the base of the repository, ie
#  bash tests/run_tests.sh

# This runs four tests of basic functionality
## 1. Classification and downstream processing, paired end data
## 2. Classification and downstream processing, single end data
## 3. Downstream processing only, Kraken and Bracken reports
## 4. Downstream processing only, Kraken reports only


snakemake -s Snakefile --configfile tests/test_config/config_pe.yaml -j1  --use-singularity  --rerun-incomplete
snakemake -s Snakefile --configfile tests/test_config/config_se.yaml -j1  --use-singularity  --rerun-incomplete
snakemake -s Snakefile_downstream_only --configfile tests/test_config/config_downstream_only_bracken.yaml -j1  --use-singularity  --rerun-incomplete
snakemake -s Snakefile_downstream_only --configfile tests/test_config/config_downstream_only_kraken.yaml -j1  --use-singularity  --rerun-incomplete