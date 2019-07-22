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

1.Almeida, A. et al. A new genomic blueprint of the human gut microbiota. Nature 1 (2019). doi:10.1038/s41586-019-0965-1
