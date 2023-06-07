# Extra functions for the kraken2_classification pipeline

def get_sample_reads(sample_reads_file):
    sample_reads = {}
    paired_end = ''
    with open(sample_reads_file) as sf:
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
                reads=[s[1]]
                if paired_end != '' and paired_end:
                    sys.exit('All samples must be paired or single ended.')
                paired_end = False
            if sample in sample_reads:
                raise ValueError("Non-unique sample encountered!")
            sample_reads[sample] = reads
    return (sample_reads, paired_end)

def get_sample_reports(sample_reads_file):
    sample_reports = {}
    bracken_reports = ''
    with open(sample_reads_file) as sf:
        for l in sf.readlines():
            s = l.strip().split("\t")
            if len(s) == 1 or s[0] == 'Sample' or s[0] == '#Sample' or s[0].startswith('#'):
                continue
            sample = s[0]
            # bracken report specified
            if (len(s)==3):
                reports = [s[1],s[2]]
                if bracken_reports != '' and not bracken_reports:
                    sys.exit('All samples must be have the same set of reports.')
                bracken_reports = True
            # bracken report not specified
            elif len(s)==2:
                reports=[s[1]]
                if bracken_reports != '' and bracken_reports:
                    sys.exit('All samples must have the same set of reports.')
                bracken_reports = False
            if sample in sample_reports:
                raise ValueError("Non-unique sample encountered!")
            sample_reports[sample] = reports
    return (sample_reports, bracken_reports)
