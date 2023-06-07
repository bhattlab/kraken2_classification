# processing classification results - new version that uses cmapR
suppressMessages(library(cmapR, quietly = TRUE, warn.conflicts = FALSE))
# need this taxonomy file defined from processing the database
# should be loaded in workflow from snakemake arg
# tax.array.file <- '~/bhatt_local/kraken2_testing/taxonomy_parsing/tax_array.tsv'
# tax.array <- read.table(tax.array.file, sep='\t', quote='', header=F, comment.char = '', colClasses = 'character')
# colnames(tax.array) <- c('id', 'taxid', 'root', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies')
# rownames(tax.array) <- tax.array$taxid

# simple reader for a kraken report file
kraken_file_to_df <- function(fname){
    if (!file.exists(fname)){
        stop(paste('File does not exist:', fname))
    }
    df <- read.table(fname, sep='\t', quote='', header = F,
                     col.names=c('pct', 'reads.below','reads.direct','tax.level','taxid','name'))
    # trim whitespace, ensure names.
    # if no name, set it to the taxid
    for (i in 1:nrow(df)){
        if (trimws(df[i, 'name']) == ''){
            df[i, 'name'] <- df[i, 'taxid']
        }
    }
    return(df)
}

# workflow for this
# read files into a list of dataframes
# simplify df to just taxid, reads.below
# merge into a single matrix by taxid
# set to gtcx
# subset to each tax level

# combine a list of dataframes from reading Kraken reports and
# merge into a large matrix of read counts
# names will be taken from the names of the list
# output will have taxonmic name in rows, sample in columns
merge_kraken_df_list <- function(df.list){
    # if(length(df.list) <2){
    #     stop('Must be a list longer than 1 - TODO implement this with one element')}
    if(is.null(names(df.list))){
        warning('Names of list is null, colnames of result will be empty')}

    keep.cols <- c('taxid', 'reads.below')
    df.list.simple <- lapply(df.list, function(x) x[, keep.cols])

    # merge on the taxid column
    merge.colname <- 'taxid'
    # Reduce to merge list of data frames into one
    merge.temp <- suppressWarnings(Reduce(function(x,y) merge(x, y, all=TRUE, by=merge.colname, sort=F), df.list.simple))
    rownames(merge.temp) <- merge.temp[,1]
    merge.mat <- as.matrix(merge.temp[,c(-1),drop=FALSE])
    colnames(merge.mat) <- names(df.list)
    merge.mat[is.na(merge.mat)] <- 0
    return(merge.mat)
}

# construct new GCT from three pieces of data
# matrix of classified reads, sample metadata and taxonomy metadata
make_gct_from_kraken <- function(merge.mat, sample.metadata, tax.array){
    # subset tax array for making the gctx
    # TODO: option to keep all taxids?
    common.taxids <- intersect(tax.array$taxid, rownames(merge.mat))
    tax.array.subset <- tax.array[common.taxids,]
    kgct <- new('GCT')
    kgct@rdesc <- tax.array.subset
    rownames(kgct@rdesc) <- tax.array.subset$id
    kgct@rid <- tax.array.subset$id
    kgct@mat <- merge.mat[common.taxids,,drop=FALSE]
    rownames(kgct@mat) <- tax.array.subset$id
    kgct@cid <- colnames(merge.mat)
    kgct@cdesc <- sample.metadata
    return(kgct)
}

# and now a new suite of functions to filter this to a specific tax level
# TODO: include classified above this level (as an option?)
subset_kgct_to_level <- function(kgct, level){
    # must do this on taxonomy with root
    if (!('root' %in% kgct@rid)){
        stop('subset must be done on gct that contains root ')}
    # print(level)
    valid.levels <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies')
    if (!(level %in% valid.levels)){
        stop(paste(level, 'not in valid tax levels:', paste(valid.levels, collapse= ', ')))
    }
    if (level == 'subspecies'){
        below.levels <- ''
    } else {
        below.levels <- valid.levels[(which(valid.levels == level) +1):length(valid.levels)]
        if ( ("subspecies" %in% below.levels ) & (!'subspecies' %in% colnames(kgct@rdesc))){
            below.levels <- below.levels[below.levels != 'subspecies']
        }
    }
    # to filter to a level, tax entry at that level must be nonempty
    # and all levels BELOW musst be empty
    # print(level)
    # print(below.levels)
    keep.inds <- which((kgct@rdesc[, level] != '') &
                       (apply(kgct@rdesc[, below.levels, drop=FALSE], 1, function(x) all(x == ''))))
    # always keep unclassified
    if ('unclassified' %in% kgct@rid){
        keep.inds <- c(which(kgct@rid == 'unclassified'), keep.inds)
    }
    if(length(keep.inds)==0){
        stop('No valid ids at this level in the taxonomy')
    }
    # names(keep.inds) <- NULL
    # print(keep.inds)
    # subset to thse indices
    kgct.subset <- subset_gct(kgct, rid=keep.inds)

    # something to deal with reads unclassified, and classified above this level
    # convention: reads classified above this level have taxid -1
    if ('unclassified' %in% kgct.subset@rid){
        total.reads <- apply(kgct@mat[c('unclassified', 'root'),,drop=FALSE],2,sum)
        total.reads.subset <- apply(kgct.subset@mat,2,sum)
        classified.higher <- total.reads - total.reads.subset
        # add this to kgct.subset
        kgct.subset@mat <- rbind(classified.higher, kgct.subset@mat)
        rownames(kgct.subset@mat)[1] <- 'classified at a higher level'
        kgct.subset@rid <- c('classified at a higher level',  kgct.subset@rid)
        kgct.subset@rdesc <- rbind(c('classified at a higher level', '-1', '', '','','','','','','',''),  kgct.subset@rdesc)
    }
    return(kgct.subset)
}

# check if the gct has been filtered to a particular taxonomic level
# return T/F
check_gct_filtered <- function(kgct){
    if (class(kgct) != 'GCT'){
        stop('Call on GCT object')
    }
    # work on the rdesc, look at each column from the back
    # and see if all are filled.
    # can't do this on unclassified rows!
    valid.levels <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies')
    test.colnames <- rev(colnames(kgct@rdesc))
    test.colnames <- test.colnames[test.colnames %in% valid.levels]
    for (n in test.colnames){
        remove.ids <- c('unclassified', 'classified at a higher level')
        this.vec <- kgct@rdesc[!(kgct@rdesc$id %in% remove.ids), n]
        if (all(this.vec=='')){
            next
        } else if (all(this.vec != '')){
            return(T)
        } else {
            return(F)
        }
    }
    # if get to this point without returning, probably a mag database
    return(T)
}

# normalize a gct matrix from reads to percentages
# and filter rows with minimal abundance
normalize_kgct <- function(kgct, min_otu_percentage = 0.001, check_filtered=TRUE){
    # this can only be done at a specific taxonomic level
    if (!check_gct_filtered(kgct) & check_filtered) {
        stop('GCT must be filtered to a single tax level before normalizing')
    }

    # catch edge case with one row
    if (nrow(kgct@mat) >1 ) {
        kgct@mat <- apply(kgct@mat, 2, function(x) x/sum(x) * 100)
    } else {
        kgct@mat[kgct@mat>0] <- 100
    }

    # round to x decimal places
    round.places <- ceiling(-log10(min_otu_percentage))
    kgct@mat <- round(kgct@mat, round.places)
    # filter to only rows with > filter_thresh abundance in one sample
    # catch edge case with one column
    if (ncol(kgct@mat)==1) {
        keep.rows <- kgct@rid[which(kgct@mat > min_otu_percentage)]
    } else{
        keep.rows <- kgct@rid[apply(kgct@mat, 1, function(x) sum(x>min_otu_percentage)>=1)]
    }
    # must always have unclassified rows if starting with them
    unclassified.rownames <- c('classified at a higher level', 'unclassified')
    keep.rows <- unique(c(unclassified.rownames[unclassified.rownames %in% kgct@rid], keep.rows))
    kgct <- subset_gct(kgct, keep.rows)

    # Normalize again so that columns sum to 100 after we removed some rows
    # catch edge case with one row
    if (nrow(kgct@mat) >1 ) {
        kgct@mat <- apply(kgct@mat, 2, function(x) x/sum(x) * 100)
    } else {
        kgct@mat[kgct@mat>0] <- 100
    }
    return(kgct)
}
