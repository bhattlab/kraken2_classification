# Extra functions for the kraken2_classification pipeline

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
