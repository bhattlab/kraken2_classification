## How to work with GCT implementation and the CMapR package.

A GCT object is simply a data matrix with embedded row and column metadata. Using this package has many benefits, namely that the metadata for the matrix lives WITH the data at all times. Any operation that subsets the data also subsets the metadata, and so on. More documentaion is available at the CMapR [github](https://github.com/cmap/cmapR) and readthedocs.

In R, a GCT object has 5 slots. These are accesed with the "@" notation: 
 - gct@mat: matrix of data values
 - gct@rdesc: dataframe of of row metadata. In this implementation - taxonomic information
 - gct@rid: vector of row identifiers. Should be the same as rownames(gct@mat) and rownames(gct@rdesc)
 - gct@cdesc: dataframe of of column metadata. In this implementation - sample information
 - gct@cid: vector of column identifiers. Should be the same as colnames(gct@mat) and colnames(gct@cdesc)

There are also a few helpful functions to easily parse and subset data:
 - subset.gct(gct, RID or CID): from the CMapR package, use this to subset a gct object from either row ids, column ids, or vectors of indices.
 - subset_kgct_to_level(kgct, level): from scripts/process_classification_gctx.R, subset a taxonomy GCT to a specific taxonomic level
 - normalize_kgct(kgct): from scripts/process_classification_gctx.R, normalize a GCT to percentages
 
To load and save GCTx objects (which are compressed binary objects, you can also use the GCT format for a text based file) use the `parse.gctx()` and `save.gctx()` methods.