import glob
import os 
from os.path import join, basename, splitext
#import glob 

read_basedir = config['read_basedir']
sample_names = config['sample_names']
read_suffix = config["read_specification"]
extension = config["extension"]

# remove markfiles prior to running workflow
#[os.remove(a) for a in glob.glob('outputs/*.mark')]

rule all:
    input:
        #"outputs/mpa_reports/merge_metaphlan.txt" 
        #expand("outputs/{samp}.krak.mark", samp=sample_names),
        #expand("outputs/{samp}.krak.report.mark", samp=sample_names),
        expand("outputs/{samp}.krak.report.bracken", samp=sample_names),
        # "outputs/mpa_reports/merge_metaphlan_heatmap.png"
        # expand("kraken2_genbank_all/{samp}.krak.mark", samp=sample_names),
        # expand("kraken2_genbank_all/{samp}.krak.report.mark", samp=sample_names),
        # expand("kraken2_genbank_all/{samp}.krak.report.bracken", samp=sample_names),

rule kraken:
    input: 
        r1 = join(read_basedir, "{samp}_" + read_suffix[0] + extension),
        r2 = join(read_basedir, "{samp}_" + read_suffix[1] + extension),
    output:
        krak = "outputs/{samp}.krak.mark",
        krak_report = "outputs/{samp}.krak.report.mark"
    params: 
        db = config['database']
    resources:
        mem=48,
        time=6
    shell:
        "kraken2 --db {params.db} --threads {threads} --output {output.krak} --report {output.krak_report} " + \
        "--paired {input.r1} {input.r2}"

rule bracken: 
    input:
        rules.kraken.output
    output:
        "outputs/{samp}.krak.report.bracken"
    params: 
        db = config['database'],
        readlen = config['read_length']
    threads: 1
    resources:
        mem = 64,
        time = 1
    shell:
        "bracken -d {params.db} -i {input[1]} -o {output} -r {params.readlen}"

# convert bracken to mpa syle report if desired 
rule convert_bracken_mpa:
    input:
        rules.bracken.output
    output:
        "outputs/mpa_reports/{samp}.krak.report.bracken.mpa"
    script:
        "convert_report_mpa_style.py"


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
