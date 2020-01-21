## New databases as of 2020-01-12! 
There are new databases available for 2020. Chose the right one for your needs:
- Refseq standard: Refseq genomes assembled to "Complete Genome" or "Scaffold" quality from the following categories.
    - `/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/refseq`
    - Bacteria 
    - Archaea
    - Fungi
    - Human
    - Mouse
    - Plasmid
- Refseq viral: same as above, but only viral sequences
    - `/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/refseq_viral`
- Genbank expanded: Genbank genomes assembled to "Complete Genome" "Scaffold" or "Chromosome" quality from the following categories. The largest and most permissive database. Note that somewhere around 256G memory needed for this database.
    - `/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/genbank_genome_chromosome_scaffold`
    - Bacteria 
    - Archaea
    - Fungi
    - Human
    - Mouse
    - Plasmid
- Genbank viral: same as above, but only viral sequences
    - `/oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/genbank_viral`

Be aware that Genbank contains sequences that were not high quality enough to be added to Refseq. This includes things like MAGs, environmental samples, etc. Some of those genomes are of interest to us, like surveillance projects. See [this page](https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/) for information on Genbank vs Refseq.

Additionally, each database has had the "crAss-like phages" clade replaced with the genomes from Guerin et al. 2018, following the clustering hierarchy from that paper. The rank and taxid of each cluster:
 
 - Guerin_crAss_Alpha       Genus      4000100
 - Guerin_crAss_Alpha_01    Species    4000110
 - Guerin_crAss_Alpha_03    Species    4000120
 - Guerin_crAss_Alpha_04    Species    4000130
 - Guerin_crAss_Alpha_09    Species    4000140
 - Guerin_crAss_Beta        Genus      4000200
 - Guerin_crAss_Beta_06     Species    4000210
 - Guerin_crAss_Gamma       Genus      4000300
 - Guerin_crAss_Gamma_02    Species    4000310
 - Guerin_crAss_Gamma_05    Species    4000320
 - Guerin_crAss_Delta       Genus      4000400
 - Guerin_crAss_Delta_07    Species    4000410
 - Guerin_crAss_Delta_08    Species    4000420
 - Guerin_crAss_Delta_10    Species    4000430

The following genomes have also been added to the Genbank expanded database.
- _Prevotella copri_ complete genome from Eli's Nanopore assembly
- Mouse microbiome Nanopore MAGs from Jessy's project
- Anything else

To see if your organism of interest is present in a database (and therefore is able to be classified in your reads), search the `inspect.out` file in the database folder. If this has the name or taxonomic identifier you're interested in, you're good to go! In this example we see that a crAss-like phage cluster is indeed present in the Refseq viral database:
```
$ grep -i Guerin_crAss_Alpha_01 /oak/stanford/scg/lab_asbhatt/data/program_indices/kraken2/kraken_custom_jan2020/refseq_viral/inspect.out 
>    0.31  227983  227983  S       4000110               Guerin_crAss_Alpha_01
```

*If something breaks, the databases are actually still building at the moment*

#### Custom sequences
Don't see your favorite bug in the `inspect.out` file? Have a newly assembled organism you want added to the database? It's an easy process to put custom sequences into the database and can be done in under a day of processing time. Contact Ben for requests (I'll batch these to once a month or something if there are lots of requests.)

_Guerin E, Shkoporov A, Stockdale SR, Clooney AG, Ryan FJ, Sutton TDS, et al. Biology and Taxonomy of crAss-like Bacteriophages, the Most Abundant Virus in the Human Gut. Cell Host & Microbe [Internet]. 2018 Oct 25 [cited 2018 Oct 25];0(0). Available from: http://www.cell.com/cell-host-microbe/abstract/S1931-3128(18)30524-9_
