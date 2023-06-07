# Snakemake pipeline for kraken2 classification of metagenomic samples
# Developed by Ben Siranosian 2018-2023
# Originally developed in the Bhatt Lab, Stanford Genetics
# MIT Licensed. https://github.com/bhattlab/kraken2_classification/

from os.path import join, exists
import sys
import snakemake
import time
localrules: downstream_processing_kraken, downstream_processing_bracken, bracken, copy_files_processing, create_taxonomy_array

# Include code from other files
# scripts/setup.smk interprets config file and sets pipeline options
include: "scripts/functions.smk"
include: "scripts/setup.smk"

rule all:
    input:
        expand(join(outdir, "classification/{samp}.krak.report"), samp=sample_names),
        join(outdir, 'processed_results_kraken/plots/classified_taxonomy_barplot_species.pdf'),
        run_extra_all_outputs,
        join(outdir, "kraken2_processing_completed.txt")

rule create_taxonomy_array:
    input:
        join(config['database'], 'taxo.k2d')
    output:
        join(config['database'], 'taxonomy_array.tsv')
    params: 
        db = config['database'],
        improve_taxonomy_script = join(workflow.basedir, 'scripts', 'improve_taxonomy.py')
    conda: "envs/anytree/anytree.yaml"
    container: "docker://bsiranosian/anytree:2.8"
    shell: """
        python {params.improve_taxonomy_script} {params.db}
    """

# Running the pipeline with singularity had a strange bug, where these files had 
# to be present in the output directory. This rule accomplishes that.
rule copy_files_processing:
    input: 
        tax_array = join(config['database'], 'taxonomy_array.tsv')
    output:
        join(outdir,  'taxonomy_array.tsv'),
        join(outdir,  'scripts/post_classification_workflow.R')
    params:
        scriptdir = join(workflow.basedir, 'scripts')
    shell: """
    cp -r {params.scriptdir} {outdir}
    cp {input.tax_array} {outdir}
    """

# Run classification with Kraken2
rule kraken:
    input:
        reads = lambda wildcards: sample_reads[wildcards.samp],
    output:
        krak = join(outdir, "classification/{samp}.krak"),
        krak_report = join(outdir, "classification/{samp}.krak.report")
    params:
        db = config['database'],
        paired_string = paired_string,
        confidence_threshold = confidence_threshold
    threads: kraken_threads
    resources:
        mem=kraken_memory,
        time=6
    singularity: "docker://quay.io/biocontainers/kraken2:2.1.2--pl5262h7d875b9_0"
    shell: """
        time kraken2 --db {params.db} --threads {threads} --output {output.krak} \
        --report {output.krak_report} {params.paired_string} {input.reads} \
        --confidence {params.confidence_threshold} --use-names
        """

# Run Bracken, if specified in the config file
rule bracken:
    input:
        krak_report = join(outdir, "classification/{samp}.krak.report"),
        krak = join(outdir, "classification/{samp}.krak")
    output:
        join(outdir, "classification/{samp}.krak_bracken_species.report"),
    params:
        db = config['database'],
        readlen = config['read_length'],
        level = config['taxonomic_level'],
        threshold = 10,
        outspec = join(outdir, "classification/{samp}.krak.report.bracken"),
    threads: bracken_threads
    resources:
        mem = bracken_memory,
        time = 1
    singularity: "docker://quay.io/biocontainers/bracken:2.8--py310h0dbaff4_1"
    shell: """
        bracken -d {params.db} -i {input.krak_report} -o {params.outspec} -r {params.readlen} \
        -l {params.level} -t {params.threshold}
        """

# Downstream processing with R
## Run for Kraken, and also Bracken if the tool was run
rule downstream_processing_kraken:
    input:
        downstream_processing_input_kraken,
        tax_array = join(outdir, 'taxonomy_array.tsv'),
        script_test = join(outdir, 'scripts/post_classification_workflow.R')
    params:
        sample_reads_file = config["sample_reads_file"],
        sample_reports_file = config["sample_reports_file"],
        sample_groups_file = config["sample_groups_file"],
        workflow_outdir = outdir,
        result_dir = join(outdir, 'processed_results_kraken'),
        use_bracken_report = False,
        remove_chordata = config['remove_chordata']
    singularity: "shub://bhattlab/kraken2_classification:kraken2_processing"
    output:
        join(outdir, 'processed_results_kraken/plots/classified_taxonomy_barplot_species.pdf')
    script:
        'scripts/post_classification_workflow.R'

rule downstream_processing_bracken:
    input:
        downstream_processing_input_bracken,
        tax_array = join(outdir, 'taxonomy_array.tsv'),
        script_test = join(outdir, 'scripts/post_classification_workflow.R')
    params:
        sample_reads_file = config["sample_reads_file"],
        sample_reports_file = config["sample_reports_file"],
        sample_groups_file = config["sample_groups_file"],
        workflow_outdir = outdir,
        result_dir = join(outdir, 'processed_results_bracken'),
        use_bracken_report = config['run_bracken'],
        remove_chordata = config['remove_chordata']
    singularity: "shub://bhattlab/kraken2_classification:kraken2_processing"
    output:
        join(outdir, 'processed_results_bracken/plots/classified_taxonomy_barplot_species.pdf')
    script:
        'scripts/post_classification_workflow.R'

# Remove file copied files during setup
rule remove_files_processing:
    input: 
        rules.downstream_processing_kraken.output,
        run_extra_all_outputs
    output:
        join(outdir, "kraken2_processing_completed.txt")
    params:
        workflow_outdir = outdir
    shell: """
    rm -rf {params.workflow_outdir}/scripts
    rm -f {params.workflow_outdir}/taxonomy_array.tsv
    touch {output}
    """

rule krona:
    input: rules.kraken.output.krak_report
    output: join(outdir, "krona/{samp}.html")
    shell: """
        ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i {input} -o {output} \
        -tax $(which kraken2 | sed 's/envs\/classification2.*$//g')/envs/classification2/bin/taxonomy
        """

# Optional rule to extract unmapped reads from the Kraken2 output
## Two versions: paired and single-end
rule extract_unmapped_paired:
    input:
        krak = join(outdir, "classification/{samp}.krak"),
        r1 = lambda wildcards: sample_reads[wildcards.samp][0],
        r2 = lambda wildcards: sample_reads[wildcards.samp][1],
    output:
        r1 = join(outdir, "unmapped_reads/{samp}_unmapped_1.fq"),
        r2 = join(outdir, "unmapped_reads/{samp}_unmapped_2.fq")
    params:
        taxid = str(0),
        tempfile = "{samp}_" + str(0) + "_reads.txt"
    resources:
        mem = 64
    singularity: "quay.io/biocontainers/bbmap:38.86--h1296035_0"
    shell: """
        awk '$1=="U" {{ print }}' {input.krak} | cut -f 2 > {params.tempfile}
        filterbyname.sh in={input.r1} in2={input.r2} names={params.tempfile} include=true out={output.r1} out2={output.r2}
        rm {params.tempfile}
    """

rule extract_unmapped_single:
    input:
        krak = join(outdir, "classification/{samp}.krak"),
        r1 = lambda wildcards: sample_reads[wildcards.samp],
    output:
        r1 = join(outdir, "unmapped_reads/{samp}_unmapped.fq"),
    params:
        taxid = str(0),
        tempfile = "{samp}_" + str(0) + "_reads.txt"
    singularity: "quay.io/biocontainers/bbmap:38.86--h1296035_0"
    shell: """
        awk '$1=="U" {{ print }}' {input.krak} | cut -f 2 > {params.tempfile}
        filterbyname.sh in={input.r1} names={params.tempfile} include=true out={output.r1}
        # rm {params.tempfile}
    """


################################################################################
# Cleanup rule - can be run after everything is done. Removes *.krak files,
# which contain information on every single read and can therefore be quite large
rule cleanup:
    input: join(outdir, 'processed_results/plots/classified_taxonomy_barplot_species.pdf')
    output: join(outdir, "cleaned")
    params:
        rmdir_1 = join(outdir, 'classification'),
    shell: """
        rm -f {params.rmdir_1}/*.krak
        touch {output}
    """

# Older pipeline rules... not used currently but could be re-enabled if needed
'''
# convert bracken to mpa syle report if desired
rule convert_bracken_mpa:
    input:
        rules.bracken.output
    output:
        "outputs/mpa_reports/{samp}.krak.report.bracken.mpa"
    script:
        "scripts/convert_report_mpa_style.py"


rule norm_mpa:
    input:
        rules.convert_bracken_mpa.output
    output:
        "outputs/mpa_reports/{samp}.krak.report.bracken.mpa.norm"
    shell:
        """
        sum=$(grep -vP "\\|" {input} | cut -f 2 | awk '{{sum += $1}} END {{printf ("%.2f\\n", sum/100)}}')
        awk -v sum="$sum" 'BEGIN {{FS="\\t"}} {{OFS="\\t"}} {{print $1,$2/sum}}' {input} > {output}
        """

rule merge_mpa:
    input:
        expand("outputs/mpa_reports/{samp}.krak.report.bracken.mpa.norm", samp=sample_names)
    output:
        merge = "outputs/mpa_reports/merge_metaphlan.txt",
        merge_species = "outputs/mpa_reports/merge_metaphlan_species.txt"
    shell:
        """
        source activate biobakery2
        merge_metaphlan_tables.py {input} >  {output.merge}
        grep -E "(s__)|(^ID)"  {output.merge} | grep -v "t__" | sed 's/^.*s__//g' >  {output.merge_species}
        """

rule hclust_mpa:
    input:
        merge = "outputs/mpa_reports/merge_metaphlan.txt"
    output:
        heamap1 = "outputs/mpa_reports/merge_metaphlan_heatmap.png",
        heamap2 = "outputs/mpa_reports/merge_metaphlan_heatmap_big.png"
    shell:
        """
        source activate biobakery2
        metaphlan_hclust_heatmap.py --in {input} --top 25 --minv 0.1 -s log --out {output.heatmap1} -f braycurtis -d braycurtis -c viridis
        metaphlan_hclust_heatmap.py --in {input} --top 150 --minv 0.1 -s log --out {output.heatmap2} -f braycurtis -d braycurtis -c viridis
        """

# make biom formatted tables for use with Qiime2
rule make_biom:
    input:
        expand("outputs/{samp}.krak.report.bracken", samp=sample_names)
    output:
        "outputs/table.biom"
    shell:
        """
        kraken-biom outputs/*_bracken.report -o {output}
        """
'''
