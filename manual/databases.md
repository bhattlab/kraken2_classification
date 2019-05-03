## Available databases
Specify the right database for your classification needs. To see if your organism of interest is present in a database (and therefore is able to be classified in your reads), search the `inspect.out` file in the database folder. If this has the name or taxonomic identifier you're interested in, you're good to go! In this example we see that crassphage is indeed present in the DB:
```
$ grep -i crassphage /labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria/inspect.out 
>  0.00  13159   13159   S       1211417         uncultured crAssphage

```

### Recommended database
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_feb2019/genbank_genome_chromosome_scaffold`

This should suit most classification needs. Includes all sequences from genbank that were assembled to Complete Genome, Chromosome or Scaffold status as of February 2019 from the following classes:
- Bacteria 
- Archea
- Fungi
- Human
- Mouse
- Plasmid

High memory useage and higher classification time, but gives the highest classification percentages. 256Gb memory needed. This is accounted for when running via batch jobs on SCG. 

#### Viruses
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_feb2019/viral`
Only interested in viral sequences? This is the quick and easy classification database for you! Includes all sequences from genbank that were assembled to Complete Genome, Chromosome or Scaffold status as of February 2019 from the following classes:
- Viruses

#### Viruses with Guerin crAss-like phages 
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_feb2019/viral`
Same as the above database, but has had the "crAss-like phages" clade replaced with the genomes from Guerin et al. 2018, where they assembled a bunch of crAss-like phage genomes. Classification is done at the family and genus level reported in that paper. 

_Guerin E, Shkoporov A, Stockdale SR, Clooney AG, Ryan FJ, Sutton TDS, et al. Biology and Taxonomy of crAss-like Bacteriophages, the Most Abundant Virus in the Human Gut. Cell Host & Microbe [Internet]. 2018 Oct 25 [cited 2018 Oct 25];0(0). Available from: http://www.cell.com/cell-host-microbe/abstract/S1931-3128(18)30524-9_


#### Genbank bacteria (high quality assemblies only)
`/labs/asbhatt/data/program_indices/kraken2/kraken_custom_oct2018/genbank_bacteria/`

Includes all bacterial sequences from genbank that were assembled to "Complete Genome" or "Chromosome" status as of October 2018. Will miss organisms that don't have high-quality assemblies in genbank. Very fast and low memory for simple bacterial classification. 64g memory is sufficient for this database.

Has had Prevotella copri and crAssphage sequences manually added. 

**statistics**
- 13653 sequences went into the construction of this database
- 4053 species are present
- 8383 species, subspecies or strains are present

Check the `inspect.out` file for a closer look. 

#### Custom sequences
Don't see your favorite bug in the `inspect.out` file? Have a newly assembled organism you want added to the database? It's an easy process to put custom sequences into the database and can be done in under a day of processing time. Contact Ben for requests (I'll batch these to once a month or something if there are lots of requests.)

