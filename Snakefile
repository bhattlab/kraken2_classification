# Snakefile for kraken2 classification
# Does classification, plotting, etc
from os.path import join
import sys
import snakemake
import time
# output base directory
outdir = config['outdir']
localrules: downstream_processing, downstream_processing_krakenonly, bracken

#perform a check on the Lathe git repo and warn if not up to date
onstart:
    print("Checking for updates or modifications to workflow")
    import git
    repo_dir = os.path.dirname(workflow.snakefile)
    repo = git.Repo(repo_dir)
    assert not repo.bare
    repo_git = repo.git
    stat = repo_git.diff('origin/master')
    if stat != "":
        print()
        print("#")
        print("##")
        print("###")
        print("####")
        print("#####")
        print("######")
        print()
        print('WARNING: Differences to latest version detected. Please reset changes and/or pull repo.')
        print()
        print("######")
        print("#####")
        print("####")
        print("###")
        print("##")
        print("#")
        # time.sleep(5)
    else:
        print("No updates or modifications found")

print('WORKFLOW DIR: ')
print(workflow.basedir)

def get_sample_reads(sample_file):
    sample_reads = {}
    paired_end = ''
    with open(sample_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            # paired end specified
            if (len(s)==3):
                reads = [s[1],s[2]]
                if paired_end != '' and not paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = True
            # single end specified
            elif len(s)==2:
                reads=s[1]
                if paired_end != '' and paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = False
            if sample in sample_reads:
                raise ValueError("Non-unique sample encountered!")
            sample_reads[sample] = reads
    return (sample_reads, paired_end)


# read in sample info and reads from the sample_file
sample_reads, paired_end = get_sample_reads(config['sample_file'])
if paired_end:
    paired_string = '--paired'
else:
    paired_string = ''
sample_names = sample_reads.keys()

# extra specified files to generate from the config file
extra_run_list =[]
# add bracken to extra files if running it
if config['run_bracken']:
    extra_run_list.append('bracken')
    extra_run_list.append('krakenonly_processed')
    downstream_processing_input = expand(join(outdir, "classification/{samp}.krak_bracken.report"), samp=sample_names)
else:
    downstream_processing_input = expand(join(outdir, "classification/{samp}.krak.report"), samp=sample_names)

# do we want to extract unmapped reads?
if config['extract_unmapped']:
    if paired_end:
        extra_run_list.append('unmapped_paired')
    else:
        extra_run_list.append('unmapped_single')

# additional outputs determined by whats specified in the readme
extra_files = {
    "bracken": expand(join(outdir, "classification/{samp}.krak_bracken.report"), samp=sample_names),
    "krakenonly_processed": join(outdir, 'processed_results_krakenonly/plots/classified_taxonomy_barplot_species.pdf'),
    "unmapped_paired": expand(join(outdir, "unmapped_reads/{samp}_unmapped_1.fq"), samp=sample_names),
    "unmapped_single": expand(join(outdir, "unmapped_reads/{samp}_unmapped.fq"), samp=sample_names),
    "barplot": join(outdir, "plots/taxonomic_composition.pdf"),
    "krona": expand(join(outdir, "krona/{samp}.html"), samp = sample_names),
    "mpa_heatmap": join(outdir, "mpa_reports/merge_metaphlan_heatmap.png"),
    "biom_file": join(outdir, "table.biom"),
}
run_extra_all_outputs = [extra_files[f] for f in extra_run_list]
# print("run Extra files: " + str(run_extra_all_outputs))

# set some resource requirements
if config['database'] in ['/labs/asbhatt/data/program_indices/kraken2/kraken_custom_feb2019/genbank_genome_chromosome_scaffold',
                          '/labs/asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/genbank_genome_chromosome_scaffold',
                          '/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_feb2019/genbank_genome_chromosome_scaffold',
                          '/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/genbank_genome_chromosome_scaffold']:
    kraken_memory = 256
    kraken_threads = 8
else:
    kraken_memory = 64
    kraken_threads = 4


rule all:
    input:
        expand(join(outdir, "classification/{samp}.krak.report"), samp=sample_names),
        expand(join(outdir, "classification/{samp}.krak"), samp=sample_names),
        join(outdir, 'processed_results/plots/classified_taxonomy_barplot_species.pdf'),
        run_extra_all_outputs,
        join(outdir, "kraken2_processing_completed.txt")
        # expand(join(outdir, "krona/{samp}.html"), samp = sample_names)

rule kraken:
    input:
        reads = lambda wildcards: sample_reads[wildcards.samp],
    output:
        krak = join(outdir, "classification/{samp}.krak"),
        krak_report = join(outdir, "classification/{samp}.krak.report")
    params:
        db = config['database'],
        paired_string = paired_string
    threads: kraken_threads
    resources:
        mem=kraken_memory,
        time=6
    singularity: "docker://quay.io/biocontainers/kraken2:2.0.9beta--pl526hc9558a2_0"
    shell: """
        time kraken2 --db {params.db} --threads {threads} --output {output.krak} \
        --report {output.krak_report} {params.paired_string} {input.reads} --use-names
        """

rule bracken:
    input:
        krak_report = join(outdir, "classification/{samp}.krak.report"),
        krak = join(outdir, "classification/{samp}.krak")
    output:
        join(outdir, "classification/{samp}.krak_bracken.report"),
    params:
        db = config['database'],
        readlen = config['read_length'],
        level = config['taxonomic_level'],
        threshold = 10,
        outspec = join(outdir, "classification/{samp}.krak.report.bracken"),
    threads: 1
    resources:
        mem = 64,
        time = 1
    singularity: "docker://quay.io/biocontainers/bracken:2.2--py27h2d50403_1"
    shell: """
        bracken -d {params.db} -i {input.krak_report} -o {params.outspec} -r {params.readlen} \
        -l {params.level} -t {params.threshold}
        """

# copy over scripts and taxonomy_array
# so that the R script works properly
# darn you, singularity for making me do this dumb thing
rule copy_files_processing:
    input: 
        tax_array = join(config['database'], 'taxonomy_array.tsv')
    output:
        'taxonomy_array.tsv',
        script_test = 'scripts/post_classification_workflow.R'
    params:
        scriptdir = join(workflow.basedir, 'scripts')
    shell: """
    cp -r {params.scriptdir} .
    cp {input.tax_array} .
    """

# must also have the processed taxonomy file generated manually 
rule downstream_processing:
    input:
        downstream_processing_input,
        tax_array = 'taxonomy_array.tsv',
        script_test = 'scripts/post_classification_workflow.R'
    params:
        sample_reads = config["sample_file"],
        sample_groups = config["sample_groups_file"],
        workflow_outdir = outdir,
        result_dir = join(outdir, 'processed_results'),
        use_bracken_report = config['run_bracken']
    singularity: "shub://bhattlab/kraken2_classification:kraken2_processing"
    output:
        join(outdir, 'processed_results/plots/classified_taxonomy_barplot_species.pdf')
    script:
        'scripts/post_classification_workflow.R'

rule downstream_processing_krakenonly:
    input:
        downstream_processing_input,
        tax_array = 'taxonomy_array.tsv',
        script_test = 'scripts/post_classification_workflow.R'
    params:
        sample_reads = config["sample_file"],
        sample_groups = config["sample_groups_file"],
        workflow_outdir = outdir,
        result_dir = join(outdir, 'processed_results_krakenonly'),
        use_bracken_report = False
    singularity: "shub://bhattlab/kraken2_classification:kraken2_processing"
    output:
        join(outdir, 'processed_results_krakenonly/plots/classified_taxonomy_barplot_species.pdf')
    script:
        'scripts/post_classification_workflow.R'

# remove these copied files now
rule remove_files_processing:
    input: 
        rules.downstream_processing.output
    output:
        join(outdir, "kraken2_processing_completed.txt")
    params:
        scriptdir = join(workflow.basedir, 'scripts')
    shell: """
    rm -rf scripts
    rm -f taxonomy_array.tsv
    touch {output}
    """


rule krona:
    input: rules.kraken.output.krak_report
    output: join(outdir, "krona/{samp}.html")
    shell: """
        ktImportTaxonomy -m 3 -s 0 -q 0 -t 5 -i {input} -o {output} \
        -tax $(which kraken2 | sed 's/envs\/classification2.*$//g')/envs/classification2/bin/taxonomy
        """

# optional rule to extract unmapped reads
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
