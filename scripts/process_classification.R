# toolkit for processing classification data


# read in a bracken/kraken table and return a dataframe
# at the specified taxonomic classification level
# due to some cases of Species not having a reported genus
# (or other higher level taxonomic classification), when reporting
# not at the species level, reads will be merged to the next highest 
# level available. 
# include.unmapped option will include unmapped read counts
# report.taxid option allows df to report taxid instead of name
read_kraken_report <- function(fname, filter.tax.level="S", include.unmapped=F, report.taxid=F){
    if(report.taxid){
        report.column <- 'taxid'
    } else {
        report.column <- 'name'
    }
    
    valid.tax.levels <- c('D','P','C','O','F','G','S')
    if (!file.exists(fname)){
        stop('File does not exist!')
    } else if(!(filter.tax.level %in% valid.tax.levels)) {
         stop(paste('filter tax.level must be in', paste(valid.tax.levels, collapse = ', ')))
    }
    df <- read.table(fname, sep='\t', quote='', header = F,
                     col.names=c('pct', 'reads.below','reads.direct','tax.level','taxid','name'))
    
    # take out everything that maps to this level
    df.filter <- df[df$tax.level==filter.tax.level,]
    
    # check for mismapped
    found.level <- F
    add.df <- data.frame()
    below.levels <- valid.tax.levels[(which(valid.tax.levels==filter.tax.level)+1):length(valid.tax.levels)]
    below.sublevels <- paste(below.levels, rep(1:9, each=length(below.levels)), sep='')
    ilter.tax.sublevels <- paste(filter.tax.level, 1:9, sep='')
    for (i in 1:nrow(df)){
        this.level <- df[i, 'tax.level']
        if (this.level==filter.tax.level){
            found.level <- T
        } else if (!(this.level %in% c(ilter.tax.sublevels, below.levels, below.sublevels))){
            found.level <- F
            last.line <- df[i,]
        }
        
        if (!(found.level) & (this.level  %in% c(below.levels, below.sublevels))){
            # print(df[i,'name'])
            # print(last.line)
            # when this happens, add counts to last line
            # and include it in our exported genus list
            if (!(last.line$taxid %in% add.df$taxid)){
                add.df <- rbind(add.df, last.line)
            }
        }
    }
    
    # append to df.filter
    df.filter <- rbind(df.filter, add.df)
    # add unmpped if desired
    if (include.unmapped){ 
        df.filter <- rbind(df.filter, df[df$taxid==0,])}
    df.filter <- df.filter[,c(report.column, 'reads.below')]
    # sort decreasing
    df.filter <- df.filter[order(df.filter$reads.below, decreasing = T),]
    # strip whitespace from name
    df.filter[,1] <- trimws(df.filter[,1])
    return(df.filter)
}

# combine a list of dataframes from reading Kraken reports and 
# merge into a large matrix of read counts
# names will be taken from the names of the list
# output will have taxonmic name in rows, sample in columns
merge_kraken_df_list <- function(df.list){
    if((class(df.list) != "list") | length(df.list) <2){
        stop('Must be a list longer than 1')
    }
    if(is.null(names(df.list))){
        warning('Names of list is null, colnames of result will be empty')
    }
    
    # merge on the first column name
    merge.colname <- colnames(df.list[[1]])[1]
    # ensure they're all the same
    if (!all(sapply(df.list, function(x) colnames(x)[1]==merge.colname))){
        stop('First column name must be all identical to merge')
    }
    
    # first merge
    merge.temp <- merge(df.list[[1]], df.list[[2]], by = merge.colname, ALL=TRUE, sort=F)
    # subsequent merges
    if(length(df.list) > 2){
        for (i in 3:length(df.list)){
            merge.temp <- suppressWarnings(merge(merge.temp, df.list[[i]], by = merge.colname, all.x=TRUE, all.y=TRUE, sort=F))
        }
    }

    rownames(merge.temp) <- merge.temp[,1]
    merge.mat <- as.matrix(merge.temp[,c(-1)])
    colnames(merge.mat) <- names(df.list)
    merge.mat[is.na(merge.mat)] <- 0
    return(merge.mat)
}

# convert a merged table of reads to a merged table of percentages
reads_matrix_to_percentages <- function(reads.matrix){
    return(round(apply(reads.matrix, 2, function(x) (x / sum(x)*100)), 3))
}

# given many file names, process them all at the specified level
# and combine into a large matrix
# takes in a named vector of files
many_files_to_matrix <- function(files, filter.tax.level="S", include.unmapped=F, report.taxid=F, percentages=F){
    report.list <- lapply(files, function(x) read_kraken_report(x, filter.tax.level, include.unmapped, report.taxid))
    names(report.list) <- names(files)
    merge.mat <- merge_kraken_df_list(report.list)
    if (percentages){
        merge.mat <- reads_matrix_to_percentages(merge.mat)
    }
    return(merge.mat)
}