## Available databases
Specify the right database for your classification needs. To see if your organism of interest is present in a database (and therefore is able to be classified in your reads), search the `inspect.out` file in the database folder. If this has the name or taxonomic identifier you're interested in, you're good to go! In this example we see that crassphage is indeed present in the DB:
```
$ grep -i crassphage /labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria/inspect.out 
>  0.00  13159   13159   S       1211417         uncultured crAssphage

```

### Recommended database
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_jan2019/genbank_custom/`

**NOTE** the bracken database for this has not been built yet. You will not be able to use bracken. 

This should suit most classification needs. Includes all sequences from genbank that were assembled to "Complete Genome" or "Chromosome" status (basically hiqh quality assemblies) as of January 2019 from the following classes 
- Bacteria 
- Viruses
- Archea
- Fungi
- Human
- Mouse
- Plasmid

#### Statistics

There are the following number of genomes from each domain of life in the database

Domain | Genome count
-------|-------------
d__Archaea | 264
d__Bacteria | 4282
d__Eukaryota | 106
d__Viroids | 45
d__Viruses | 9230

#### Genbank bacteria (high quality assemblies)
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria/`

Includes all bacterial sequences from genbank that were assembled to "Complete Genome" or "Chromosome" status as of October 2018. Will miss organisms that don't have high-quality assemblies in genbank. Very fast and low memory for simple bacterial classification. 64g memory is sufficient for this database.

Has had Prevotella copri and crAssphage sequences manually added. 

**statistics**
- 13653 sequences went into the construction of this database
- 4053 species are present
- 8383 species, subspecies or strains are present

Check the `inspect.out` file for a closer look. 

#### Genbank bacteria and archaea (ALL assemblies)
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria_complete/`

Includes all bacterial sequences from genbank assembled to any quality (includes "Scaffold" and "Contig" level assemblies) as of October 2018. 256g memory is necessary for this database, and runtimes will be longer. I've found results from this database to be more noisy, but with higher classification percentages. There is a long tail of very low abundance species present. Running Bracken is essential when using this database.

Has had Prevotella copri and crAssphage sequences manually added. 

**statistics**
- 167229 sequences went into the construction of this database
- 38789 species are present
- 63871 species, subspecies or strains are present

Check the `inspect.out` file for a closer look. 

#### Standard (high quality refseq assemblies)
`/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/standard/`

Contains high-quality refseq assemblies of archaea, bacteria, human, UniVec_Core and viral sequences.

**statistics**
- 20934 sequences went into the construction of this database
- 11578 species are present
- 17506 species, subspecies or strains are present

Check the `inspect.out` file for a closer look. 

#### Standard protein (high quality refseq assemblies)
`/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/standard_protein/`

Protein database of high-quality refseq assemblies of archaea, bacteria, human, UniVec_Core and viral sequences.

#### nt, env_nt
The equivalent of the blast nt and env_nt databases can be found at 
`/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/nt`

`/labs/asbhatt/data/program_indices/kraken2/kraken_unmod/env_nt`

I'm not sure how useful these are and haven't tried them in any real classification tasks.

#### Custom sequences
Don't see your favorite bug in the `inspect.out` file? Have a newly assembled organism you want added to the database? It's an easy process to put custom sequences into the database and can be done in under a day of processing time. Contact Ben for requests (I'll batch these to once a month or something if there are lots of requests.)

