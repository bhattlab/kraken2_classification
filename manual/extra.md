## Additional considerations
### Memory usage
Memory usage depends on the database. For a rough estimate of memory requirements for a given db, check the size of the `hash.k2d` file, which has to fit completely in memory. The manual recommends at least 175GB for the standard database. In my experiene, usage maxes out at the following values:

- Genbank high quality: 64g
- Genbank all: 256g

You should request at least this much memory on a node you plan to do classification on. 

### Time usage
Kraken2 is very fast compared to kraken1 or other classification tools. After loading the database in memory the first time (about 5 minutes for genbank hq), subsequet runs proceed very quickly. Using 8 cores on SCG I had the following results with the genbank hq database: 
`10399748 sequences (2356.32 Mbp) processed in 85.562s (7292.8 Kseq/m, 1652.36 Mbp/m).`
That's 85 seconds for a dataset with about 10 million read pairs. Obviously more cores can make things faster for larger datasets.

### Classification percentages
I've evaluated some of the databases below on sets different metagenomic datasets, and compared the unclassified percentages to our old Kraken1 custom database. See [this spreadsheet](https://docs.google.com/spreadsheets/d/15nVMno4w4Q-DVO9tdp1DBpwKslOkLb2TVqnSnuEofTY/edit?usp=sharing) for a summary of results. Overall, Kraken2 with the high quality database classifies less reads than Kraken1 but much quicker. The full genbank database has increased classification percentages and speed.


### Taxonomy modifications
When modifying the NCBI taxonomy with the [improve_taxonomy.py](../scripts/improve_taxonomy.py) script (which creates the `taxonomy_array.tsv` file necessary for downstream processing steps), the following changes are made to make the taxonomy more human interpretable: 

 - Anything below a species gets set as rank subspecies
 - Prune the taxonomy to a specified set of levels.
 - To fix cases of "unclassified XXX" or "environmental samples" being at weird places in the tree and not having a useful name, use the heuristic:
   - if the name is "unclassified X" or "environmental samples" and everything below it is a species, and parent is not a genus or species than mark that as a genus
 - A few custom changes
   - Set crass-like viruses to genus
   - Set unclassified bacterial viruses to family
   - Set all superkingdoms to kingdom
   - Set Fungi and metazoa to be direct descendants of root