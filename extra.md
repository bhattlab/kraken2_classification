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
I've evauluated some of the databases below on sets different metagenomic datasets, and compared the unclassified percentages to our old Kraken1 custom database. See [this spreadsheet](https://docs.google.com/spreadsheets/d/15nVMno4w4Q-DVO9tdp1DBpwKslOkLb2TVqnSnuEofTY/edit?usp=sharing) for a summary of results. Overall, Kraken2 with the high quality database classifies less reads than Kraken1 but much quicker. The full genbank database has increased classification percentages and speed.
