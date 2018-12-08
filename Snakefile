sample_reads = {}
with open(config['datasets']) as asmf:
    basedir = os.path.dirname(config['datasets'])
    for l in asmf.readlines():
        s = l.strip().split("\t")
        if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample':
            continue

        sample = s[0]
        reads = s[3].split(',')
        if sample in sample_reads:
            raise ValueError("Non-unique sample encountered!")
        sample_reads[sample] = reads

sample_names = sample_reads.keys()

rule all:
    input:
        "outputs/plot/taxonomic_composition.pdf"

rule kraken:
    input: 
        r1 = lambda wildcards: sample_reads[wildcards.samp][0],
        r2 = lambda wildcards: sample_reads[wildcards.samp][1]
    output:
        krak = "outputs/{samp}.krak",
        krak_report = "outputs/{samp}.krak.report"
    params: 
        db = config['database']
    resources:
        mem=48,
        time=6
    shell:
        "time kraken2 --db {params.db} --threads {threads} --output {output.krak} --report {output.krak_report} " + \
        "--paired {input.r1} {input.r2}" 

rule bracken: 
    input:
        rules.kraken.output
    output:
        "outputs/{samp}.krak.report.bracken"
    params: 
        db = config['database'],
        readlen = config['read_length'],
        level = config['taxonomic_level']
    threads: 1
    resources:
        mem = 64,
        time = 1
    shell:
        "bracken -d {params.db} -i {input[1]} -o {output} -r {params.readlen} -l {params.level}"

rule collect_results:
    input: expand("outputs/{samp}.krak.report.bracken", samp = sample_names)
    output: "outputs/class_long.tsv"
    script: "scripts/collate_results.py"

rule barplot:
    input: rules.collect_results.output
    output: "outputs/plot/taxonomic_composition.pdf"
    params: taxlevel='G'
    script: "scripts/composition_barplot.R"



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