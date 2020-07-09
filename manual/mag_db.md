# Creating Kraken2 databses from Metagenome Assembled Genomes
Metagenome Assembled Genomes (MAGs) help with the fact that reference databases for microbiome research are incomplete. By assembling and binning genomes from many studies, MAG databases have been used to classify new organisms and find new microbial associations with disease. For a more detailed summary, check out my [blog post on the topic](https://www.bsiranosian.com/bioinformatics/metagenome-assembled-genomes-enhance-short-read-classification/).

In the Bhatt lab, we were interested in using MAG databases for short read metagenomic classification with Kraken2. Here's how you set up a Kraken2 database from a collection of MAGs, using Almeida et al. (2019) as an example. This uses the Unclassified MAGs (UMGS) that represent new candidate species. The species are set up in a flat taxonomy - every bin is a species that links directly to root. In the future, you could use some of the taxonomy information provided in the paper to create a better taxonomy for the databse.

```
mkdir ~/almeida_genomes/
cd ~/almeida_genomes/
# download umgs from ftp
wget -r ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/umgs_analyses/umgs/genomes/
mv ftp.ebi.ac.uk/pub/databases/metagenomics/umgs_analyses/umgs/genomes/ .
rm -r ftp.ebi.ac.uk
cd genomes

# Unzip all UMGS genomes
gunzip *.gz

# add kraken_taxid to each
for i in $(seq 1 2079); do
    if [[ -f UMGS$i.fa ]]; then
        echo $i
        sed -i "/>/ c\>Almeida_UMGS_$i|kraken:taxid|1$i" UMGS$i.fa
    fi
done

# add to library - set this directory
dbdir=~/almeida_genomes/kraken2_db
find -name 'UMGS*.fa' -print0 | xargs -0 -I{} -P 8 -n1 kraken2-build --add-to-library {} --db $dbdir

mkdir ~/almeida_genomes/kraken2_db/taxonomy/
# construct nodes.dmp and names.dmp taxonomy files
echo -e "1\t|\t1\t|\tno rank\t|\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|" > ~/almeida_genomes/kraken2_db/taxonomy/nodes.dmp
while read line; do
    UMGSID=1"$line"
    echo -e "$UMGSID\t|\t1\t|\tspecies\t|\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|" >> ~/almeida_genomes/kraken2_db/taxonomy/nodes.dmp
done < ~/almeida_genomes/umgsids.txt

echo -e "1\t|\tall\t|\t\t|\tsynonym\t|" > ~/almeida_genomes/kraken2_db/taxonomy/names.dmp
echo -e "1\t|\troot\t|\t\t|\tscientific name\t|" >> ~/almeida_genomes/kraken2_db/taxonomy/names.dmp
while read line; do
    taxid=1"$line"
    UMGSID="$line"
    echo -e "$taxid\t|\tAlmeida_UMGSID_$UMGSID\t|\t\t|\tAlmeida name\t|" >> ~/almeida_genomes/kraken2_db/taxonomy/names.dmp
done < ~/almeida_genomes/umgsids.txt

# then build in parallel
kraken2-build --build --threads 16 --db ~/almeida_genomes/kraken2_db/
                               
```

## Using a dereplicated genome set
Another way to create a custom MAG database is to assemble a set of genomes from your own samples, run them through [dRep](https://drep.readthedocs.io/en/latest/) to get a non-redundant genome set, and then create a database out of those genomes. Here, I had metagenomic samples from mouse stool in different experimental conditions and sampled at different days. Because there were 5 replicates at each condition, the samples were co-assembled and co-binned to get the highest quality genomes possible (acknowledging that they may be "franken-genomes" and not truly representative of the species present). dRep was used to create the non-redundant genome set. I then made a kraken database as follows (this is kind of a hack that could be made much better by using a mapping of genome>taxid):

```
# form the root of the dRep output
mkdir -p kraken_db/genomes

# add kraken_taxid to each dereplicated genome
# start with taxid 11 and increment by 1
cd dereplicated_genomes
taxid=11
outdir=../kraken_db/genomes

for fa in *.fa; do
    echo $taxid
    fa_name=$(echo $fa | sed 's/.fa//g')
    sed "/>/ c\>$fa_name|kraken:taxid|$taxid" $fa > "$outdir"/"$fa"
    taxid=$(expr $taxid + 1)
done
cd ../kraken_db

# add to library - set this directory
dbdir=.
find -name '*.fa' -print0 | xargs -0 -I{} -P 8 -n1 kraken2-build --add-to-library {} --db $dbdir


mkdir taxonomy
# construct nodes.dmp and names.dmp taxonomy files
echo -e "1\t|\t1\t|\tno rank\t|\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|" > taxonomy/nodes.dmp
echo -e "1\t|\tall\t|\t\t|\tsynonym\t|" > taxonomy/names.dmp
echo -e "1\t|\troot\t|\t\t|\tscientific name\t|" >> taxonomy/names.dmp

taxid=11
for fa in genomes/*.fa; do
    echo $taxid
    fa_name=$(echo $fa | sed 's/.fa//g')
    fa_name=$(basename $fa_name)
    echo -e "$taxid\t|\t1\t|\tspecies\t|\t|\t0\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|" >> taxonomy/nodes.dmp
    echo -e "$taxid\t|\t$fa_name\t|\t\t|\tdrep bin name\t|" >> taxonomy/names.dmp
    taxid=$(expr $taxid + 1)
done


# then build in parallel
kraken2-build --build --threads 16 --db $dbdir

# you can't build a bracken database because there's no real taxonomy information here. Make sure to set run_braken to false in the configfile when using this database

# to get the classification pipeline to work correctly, need to create a taxonomy_array.tsv file 
# this is a tab delimited file that provides the taxonomic information of each taxid 
# something like this, with a line for each mag, should do the trick
unclassified\t0\t\t\t\t\t\t\t\t\t
root\t1\troot\t\t\t\t\t\t\t\t
MAG_NAME\t11\troot\t\t\t\t\t\t\t\t
.
.
... fill in the rest


```


1. Almeida, A. et al. A new genomic blueprint of the human gut microbiota. Nature 1 (2019). doi:10.1038/s41586-019-0965-1
