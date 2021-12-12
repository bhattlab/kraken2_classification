# Setting up your own database 
If you're not in the Bhatt lab or are using a new database, build it with instructions from the [Kraken2 github](https://github.com/DerrickWood/kraken2/) or the custom instructions below. The only other thing you need to do is create the `taxonomy_array.tsv` file, which contains an improved version of the ncbi taxonomy for downstream analysis. To do this, run `scripts/improve_taxonomy.py` with the database directory as the only argument. You need to have the python packages _numpy_ and _anytree_ installed for this script. 

# Custom database construction
To get the most use out of Kraken, I recommend creating a custom database for classification. The defautlt bacteria database is incomplete, and excludes genomes that aren't assembled to "complete genome" or "chromosome" quality in NCBI. This misses key commensals like _Bacteroides intestinalis_. Here, we'll switch the database to use Genbank instead of Refseq and include genomes of "scaffold" quality as well. This requires modifying some of the Kraken2 code. 

### Set ftp server to genbank
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

### Less stringent genome quality filters
```
libexec=/data/bsiranos/miniconda3/envs/kraken2/libexec
cp "$libexec/rsync_from_ncbi.pl " "$libexec/rsync_from_ncbi.pl .bak" 
sed -i "s/\"Chromosome\"/\"Chromosome\", "Scaffold"/g" "$libexec/rsync_from_ncbi.pl"
```

### Build the database
You need to download the taxonomy and respective database sets (bacteria, fungi, etc) in the same way you usually would. Then build the kraken and brakcen databses, follwing the [standard instructions](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases). This db build will take a long time - it took ~24h with 32 cores.


### Instructions for building Dec2021 databases
Here are the commands and steps I used to build the december 2021 databases. This was done on a 96-core GCP machine. The genbank kraken2 and bracken database construction took 48 hours to finish roughly. 
```
# need to modify some code in kraken because https://github.com/DerrickWood/kraken2/issues/518

# download taxonomy
kraken2-build --download-taxonomy --db refseq
kraken2-build --download-library archaea --db refseq --threads 8
kraken2-build --download-library bacteria --db refseq --threads 8 &
kraken2-build --download-library plasmid --db refseq --threads 8 &
kraken2-build --download-library viral --db refseq --threads 8 &
kraken2-build --download-library fungi --db refseq --threads 8 &
kraken2-build --download-library human --db refseq --threads 8 &

# additional genomes to add to the database are here: gs://gbsc-gcp-lab-bhatt_user-bsiranos/kraken_db_construction_dec2021/add_genomes
# add the custom genomes that go into refseq
# refseq has mouse and crassphages
# mouse is easy
kraken2-build --add-to-library add_genomes/mm10.fa --db refseq --threads 8

# crassphage much more complicated
# taxid of crAss-like phages clade: 1978007
# we can remove all of the existing entries that are below this node because we replace with the guerin taxonomy
# if we just filter linking up to 1978007 with this fancy grep expression it should stop those from appearing in the database
mv refseq/taxonomy/nodes.dmp refseq/taxonomy/nodes_bak.dmp 
grep -v "[0-9][[:blank:]]|[[:blank:]]1978007" refseq/taxonomy/nodes_bak.dmp > refseq/taxonomy/nodes.dmp
# then add the crasslike nodes and names to the bottom of these files
cat add_genomes/add_guerin_names.txt >> refseq/taxonomy/names.dmp
cat add_genomes/add_guerin_nodes.txt >> refseq/taxonomy/nodes.dmp

# do the actual genome adding
ls add_genomes/Guerin_Data_S1/split_genomes/modified_headers/* | xargs -I {} -P 64 sh -c "echo{}; kraken2-build --add-to-library {}  --db refseq"

# do viral construction after that download finishes
mkdir -p refseq_viral/library
ln -s $PWD/refseq/taxonomy refseq_viral/
# doesnt work as a symlink for some reason...
cp -r $PWD/refseq/library/viral refseq_viral/library
# want all the crassphage genome, but not mouse
# this gives files to exclude
grep 10090 refseq/library/added/*
cp -r refseq/library/added refseq_viral/library
rm refseq_viral /library/added/JG8Uo9uGN6*
rm refseq_viral /library/added/prelim_map_J_5GO1nzcM.txt

# BUILD THE REFSEQ VIRAL DB
kraken2-build --build --db refseq_viral --threads 48
# For bracken, need to pull code from github and use fix described here: https://github.com/jenniferlu717/Bracken/issues/54#issuecomment-968787387
~/Bracken/bracken-build -d refseq_viral -t 48 -k 35 -l 150
~/Bracken/bracken-build -d refseq_viral -t 48 -k 35 -l 100
kraken2-inspect --db refseq_viral > refseq_viral/inspect.out


# BUILD THE REFSEQ COMPLETE DB
kraken2-build --build --db refseq --threads 48
# For bracken, need to pull code from github and use fix described here: https://github.com/jenniferlu717/Bracken/issues/54#issuecomment-968787387
~/Bracken/bracken-build -d refseq -t 48 -k 35 -l 150
~/Bracken/bracken-build -d refseq -t 48 -k 35 -l 100
kraken2-inspect --db refseq > refseq/inspect.out

# prepare for genbank construction
mkdir -p genbank/library
mkdir -p genbank_viral/library
ln -s $PWD/refseq/taxonomy genbank
ln -s $PWD/refseq/taxonomy genbank_viral
# copy added genomed from refseq
cp -r refseq/library/added genbank/library/
cp -r refseq_viral/library/added genbank_viral/library/

# modify code to download from genbank
grep "my \$KRAKEN2_DIR =" $(which kraken2)
libexec=/home/bsiranos/miniconda3/libexec
cp "$libexec/download_genomic_library.sh" "$libexec/download_genomic_library.sh.bak" 
sed -i "s/refseq/genbank/g" "$libexec/download_genomic_library.sh"
sed -i "s/genbank_old/genbank/g" "$libexec/download_genomic_library.sh"
#Less stringent genome quality filters
cp "$libexec/rsync_from_ncbi.pl" "$libexec/rsync_from_ncbi.pl.bak" 
sed -i "s/\"Chromosome\"/\"Chromosome\", \"Scaffold\"/g" "$libexec/rsync_from_ncbi.pl"

# download refseq
kraken2-build --download-library archaea --db genbank --threads 8 & 
kraken2-build --download-library bacteria --db genbank --threads 8 &
kraken2-build --download-library viral --db genbank --threads 8 &
kraken2-build --download-library fungi --db genbank --threads 8 &
# human is the same as before
cp -r refseq/library/human genbank/library
# fungi and plasmid had issues downloading
cp -r refseq/library/plasmid genbank/library
cp -r refseq/library/fungi genbank/library

# BUILD THE GENBANK COMPLETE DB
kraken2-build --build --db genbank --threads 48
# For bracken, need to pull code from github and use fix described here: https://github.com/jenniferlu717/Bracken/issues/54#issuecomment-968787387
~/Bracken/bracken-build -d genbank -t 48 -k 35 -l 150
~/Bracken/bracken-build -d genbank -t 48 -k 35 -l 100
kraken2-inspect --db genbank > genbank/inspect.out

# genbank viral
cp -r genbank/library/viral genbank_viral/library
# BUILD THE GENBANK VIRAL DB
kraken2-build --build --db genbank_viral --threads 48
# For bracken, need to pull code from github and use fix described here: https://github.com/jenniferlu717/Bracken/issues/54#issuecomment-968787387
~/Bracken/bracken-build -d genbank_viral -t 48 -k 35 -l 150
~/Bracken/bracken-build -d genbank_viral -t 48 -k 35 -l 100
kraken2-inspect --db genbank_viral > genbank_viral/inspect.out
```