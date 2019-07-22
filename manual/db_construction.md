#Custom database construction
To get the most use out of Kraken, I recommend creating a custom database for classification. The defautlt bacteria database is incomplete, and excludes genomes that aren't assembled to "complete genome" or "chromosome" quality in NCBI. This misses key commensals like _Bacteroides intestinalis_. Here, we'll switch the database to use Genbank instead of Refseq and include genomes of "scaffold" quality as well. This requires modifying some of the Kraken2 code. 

###Set ftp server to genbank
First, find where the supplementary scripts to kraken2 are installed
```
$ grep "my \$KRAKEN2_DIR =" $(which kraken2)
my $KRAKEN2_DIR = "/data/bsiranos/miniconda3/envs/kraken2/libexec";
```
This libexec folder is what you want. Backup and modify the script, changing the libexec variable here to the output of the above command.
```
libexec=/data/bsiranos/miniconda3/envs/kraken2/libexec
cp "$libexec/download_genomic_library.sh" "$libexec/download_genomic_library.sh.bak" 
sed -i "s/refseq/genbank/g" "$libexec/download_genomic_library.sh"
```

###Less stringent genome quality filters
```
libexec=/data/bsiranos/miniconda3/envs/kraken2/libexec
cp "$libexec/rsync_from_ncbi.pl " "$libexec/rsync_from_ncbi.pl .bak" 
sed -i "s/\"Chromosome\"/\"Chromosome\", "Scaffold"/g" "$libexec/rsync_from_ncbi.pl"
```

###Build the database
You need to download the taxonomy and respective database sets (bacteria, fungi, etc) in the same way you usually would. Then build the kraken and brakcen databses, follwing the [standard instructions](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases). This db build will take a long time - it took ~24h with 32 cores. 