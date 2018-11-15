import glob
import os 
from os.path import join, basename, splitext

configfile: "config.yaml"
read_basedir = config['read_basedir']
sample_names = config['sample_names']
read_suffix = config["read_specification"]
extension = config["extension"]

rule all:
    input: 
        expand("kraken2_genbank_hq/{samp}.krak.mark", samp=sample_names),
        expand("kraken2_genbank_hq/{samp}.krak.report.mark", samp=sample_names),
        expand("kraken2_genbank_hq/{samp}.krak.report.bracken", samp=sample_names),
        # expand("kraken2_genbank_all/{samp}.krak.mark", samp=sample_names),
        # expand("kraken2_genbank_all/{samp}.krak.report.mark", samp=sample_names),
        # expand("kraken2_genbank_all/{samp}.krak.report.bracken", samp=sample_names),

rule kraken_classify_hq:
    input:
        r1 = expand(join(read_basedir, "{samp}_" + read_suffix[0] + extension), samp=sample_names),
        r2 = expand(join(read_basedir, "{samp}_" + read_suffix[1] + extension), samp=sample_names)
    output:
        krak = expand("kraken2_genbank_hq/{samp}.krak.mark", samp=sample_names),
        krak_report = expand("kraken2_genbank_hq/{samp}.krak.report.mark", samp=sample_names)
    params: 
        db = "/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria"
    threads: 16
    resources:
        mem = 64, 
        time = 24
    run:
        # print(params)
        # have any files already been run? If so dont run them again
        files_krak = [splitext(basename(x))[0] for x in glob.glob('kraken2_genbank_hq/*.krak')]
        files_report = [splitext(splitext(basename(x))[0])[0] for x in glob.glob('kraken2_genbank_hq/*.krak.report')]
        files_both = list(set(files_krak) & set(files_report))
        # generate markfiles for those that have already been run
        for f in files_both:
            t1 = join('kraken2_genbank_hq', f + '.krak.mark')
            t2 = join('kraken2_genbank_hq', f + '.krak.report.mark')
            open(t1, 'a').close()
            open(t2, 'a').close()

        do_samples = list(set(sample_names) - set(files_both))
        r1_new = expand(join(read_basedir, "{samp}_" + read_suffix[0] + extension), samp=do_samples)
        r2_new = expand(join(read_basedir, "{samp}_" + read_suffix[1] + extension), samp=do_samples)
        krak_new = expand("kraken2_genbank_hq/{samp}.krak", samp=do_samples)
        krak_report_new = expand("kraken2_genbank_hq/{samp}.krak.report", samp=do_samples)
        print('running on ' + str(len(do_samples)) + ' out of ' + str(len(sample_names)))
        # print(len(r1_new))
        # print(len(r2_new))
        # print(len(krak_new))
        # print(len(krak_report_new))
        for r1, r2, outf, outfr in zip(r1_new, r2_new, krak_new, krak_report_new):
            sub = "kraken2 --db {db} --threads {threads} --output {outf} \
                   --report {outfr} --paired {r1} {r2}".format(\
                    db=params, threads=threads, outf=outf, outfr=outfr,
                    r1=r1, r2=r2)
            # print(sub)
            os.system(sub)
            open(outf + '.mark','a').close()
            open(outfr + '.mark','a').close()

rule kraken_classify_all:
    input:
        r1 = expand(join(read_basedir, "{samp}_" + read_suffix[0] + extension), samp=sample_names),
        r2 = expand(join(read_basedir, "{samp}_" + read_suffix[1] + extension), samp=sample_names)
    output:
        krak = expand("kraken2_genbank_all/{samp}.krak.mark", samp=sample_names),
        krak_report = expand("kraken2_genbank_all/{samp}.krak.report.mark", samp=sample_names)
    params: 
        db = "/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria"
    threads: 16
    resources:
        mem = 256, 
        time = 48
    run:
        # print(params)
        # have any files already been run? If so dont run them again
        files_krak = [splitext(basename(x))[0] for x in glob.glob('kraken2_genbank_all/*.krak')]
        files_report = [splitext(splitext(basename(x))[0])[0] for x in glob.glob('kraken2_genbank_all/*.krak.report')]
        files_both = list(set(files_krak) & set(files_report))
        # generate markfiles for those that have already been run
        for f in files_both:
            t1 = join('kraken2_genbank_all', f + '.krak.mark')
            t2 = join('kraken2_genbank_all', f + '.krak.report.mark')
            open(t1, 'a').close()
            open(t2, 'a').close()

        do_samples = list(set(sample_names) - set(files_both))
        r1_new = expand(join(read_basedir, "{samp}_" + read_suffix[0] + extension), samp=do_samples)
        r2_new = expand(join(read_basedir, "{samp}_" + read_suffix[1] + extension), samp=do_samples)
        krak_new = expand("kraken2_genbank_all/{samp}.krak", samp=do_samples)
        krak_report_new = expand("kraken2_genbank_all/{samp}.krak.report", samp=do_samples)
        print('running on ' + str(len(do_samples)) + ' out of ' + str(len(sample_names)))
        # print(len(r1_new))
        # print(len(r2_new))
        # print(len(krak_new))
        # print(len(krak_report_new))
        for r1, r2, outf, outfr in zip(r1_new, r2_new, krak_new, krak_report_new):
            sub = "kraken2 --db {db} --threads {threads} --output {outf} \
                   --report {outfr} --paired {r1} {r2}".format(\
                    db=params, threads=threads, outf=outf, outfr=outfr, r1=r1, r2=r2)
            # print(sub)
            os.system(sub)
            open(outf + '.mark','a').close()
            open(outfr + '.mark','a').close()

rule bracken_hq: 
    input:
        krak_report = "kraken2_genbank_hq/{samp}.krak.report.mark"
    output:
        bracken_report = "kraken2_genbank_hq/{samp}.krak.report.bracken"
    params: 
        db = "/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria",
        readlen = 150,
        actual_input = "kraken2_genbank_hq/{samp}.krak.report"
    threads: 1
    resources:
        mem = 64,
        time = 1
    shell:
        """
        bracken -d {params.db} -i {params.actual_input} -o {output.bracken_report} -r {params.readlen}
        """

rule bracken_all: 
    input:
        krak_report = "kraken2_genbank_all/{samp}.krak.report.mark"
    output:
        bracken_report = "kraken2_genbank_all/{samp}.krak.report.bracken"
    params: 
        db = "/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria_all",
        readlen = 150,
        actual_input = "kraken2_genbank_all/{samp}.krak.report"
    threads: 1
    resources:
        mem = 64,
        time = 1
    shell:
        """
        bracken -d {params.db} -i {params.actual_input} -o {output.bracken_report} -r {params.readlen}
        """
