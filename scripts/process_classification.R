# toolkit for processing classification data
options(stringsAsFactors = F)

# simple reader for a kraken file
kraken_file_to_df <- function(fname){
  if (!file.exists(fname)){
    stop(paste('File does not exist:', fname))
  }
  df <- read.table(fname, sep='\t', quote='', header = F,
                   col.names=c('pct', 'reads.below','reads.direct','tax.level','taxid','name'))
  # ensure names, if not add taxid
  for (i in 1:nrow(df)){
    if (trimws(df[i, 'name']) == ''){
      df[i, 'name'] <- df[i, 'taxid']
    }
  }
  return(df)
}

# take a dataframe from a bracken/kraken file and return
# a dataframe filtered at the specified taxonomic classification level
# due to some cases of Species not having a reported genus
# (or other higher level taxonomic classification), when reporting
# not at the species level, reads will be merged to the next highest
# level available.
# include.unclassified option will include unmapped read counts
# report.taxid option allows df to report taxid instead of name
parse_kraken_report <- function(df, filter.tax.level="S", include.unclassified=F, report.taxid=F){
  valid.tax.levels <- c('D','P','C','O','F','G','S')
  if(!(filter.tax.level %in% valid.tax.levels)) {
    stop(paste('filter tax.level must be in', paste(valid.tax.levels, collapse = ', ')))
  }
  if(report.taxid){
    report.column <- 'taxid'
  } else {
    report.column <- 'name'
  }

  # if domain, add in Kingdom too becuase theyre all classified that way...
  if (filter.tax.level == 'D'){
    filter.tax.levels <- c('D','K')
  } else {
    filter.tax.levels <- c(filter.tax.level)
  }

  # take out everything that maps to this level
  df.filter <- df[df$tax.level %in% filter.tax.levels,]

  # check for mismapped
  found.level <- F
  add.df <- data.frame()
  below.levels <- valid.tax.levels[(which(valid.tax.levels %in% filter.tax.levels)+1):length(valid.tax.levels)]
  below.sublevels <- paste(below.levels, rep(1:9, each=length(below.levels)), sep='')
  filter.tax.sublevels <- paste(filter.tax.level, 1:9, sep='')
  for (i in 1:nrow(df)){
    this.level <- df[i, 'tax.level']
    if (this.level %in% filter.tax.levels){
      found.level <- T
    } else if (!(this.level %in% c(filter.tax.sublevels, below.levels, below.sublevels))){
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
  # dont do this if looking at domain level, get weird results
  if (!('D' %in% filter.tax.levels)){
    df.filter <- rbind(df.filter, add.df)
  }

  # add unmpped if desired
  # first is unmapped totally - unclassified reads
  # second is reads Classified at a higher level
  # that's (all reads - (unclassified + classified at given level))
  if (include.unclassified){
    # print(df[df$taxid==0,])
    unclassified.completely <- df[df$taxid==0, 'reads.direct']
    unclassified.thislevel <- max(0, sum(df$reads.direct) - (sum(df.filter$reads.below) + unclassified.completely))
    unclassified.df <- data.frame(pct=NA, reads.below=c(unclassified.completely, unclassified.thislevel),
                                  reads.direct=c(unclassified.completely, unclassified.thislevel),
                                  taxid=0, tax.level=NA, name=c('Unclassified', 'Classified at a higher level'))
    df.filter <- rbind(df.filter, unclassified.df)}
  df.filter <- df.filter[,c(report.column, 'reads.below')]
  # strip whitespace from name
  df.filter[,1] <- trimws(df.filter[,1])
  # aggregate any duplicate rows - sometimes happens with 'environmental samples'
  df.filter <- aggregate(reads.below ~ name, data=df.filter, FUN=sum)
  # sort decreasing
  df.filter <- df.filter[order(df.filter$reads.below, decreasing = T),]
  return(df.filter)
}

# convenience funnction to read and filter all in one
read_kraken_report <- function(fname, ...){
  df <- kraken_file_to_df(fname)
  df.filter <- parse_kraken_report(df, ...)
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

  # Reduce to merge list of data frames into one
  merge.temp <- suppressWarnings(Reduce(function(x,y) merge(x, y, all=TRUE, by=merge.colname, sort=F), df.list))

  rownames(merge.temp) <- merge.temp[,1]
  merge.mat <- as.matrix(merge.temp[,c(-1)])
  colnames(merge.mat) <- names(df.list)
  merge.mat[is.na(merge.mat)] <- 0
  return(merge.mat)
}

# convert a merged table of reads to a merged table of percentages
reads_matrix_to_percentages <- function(reads.matrix){
  # catch edge case with one row
  if (nrow(reads.matrix) >1 ) {
    reads.matrix.pct <- apply(reads.matrix, 2, function(x) x/sum(x) * 100)
  } else {
    reads.matrix.pct <- reads.matrix
    reads.matrix.pct[reads.matrix.pct>0] <- 100
  }
  return(reads.matrix.pct)
}

# takes in a named vecor of files
# given many file names, process them all at the specified levels
# combines into a list of matrices, one for each taxonomic level
many_files_to_matrix_list <- function(files, filter.tax.levels=c("S"), include.unclassified=F, report.taxid=F, percentages=F){
  df.list <- lapply(files, function(x) kraken_file_to_df(x))
  merge.mat.list <- list()
  for (filter.tax.level in filter.tax.levels){
    report.list <- lapply(df.list, function(x) parse_kraken_report(x, filter.tax.level, include.unclassified, report.taxid))
    names(report.list) <- names(files)
    # only merge if >1 sample
    if (length(report.list) >1){
      merge.mat <- merge_kraken_df_list(report.list)
    } else {
      merge.mat <- matrix(report.list[[1]][,2], ncol=1)
      rownames(merge.mat) <- report.list[[1]][,1]
      colnames(merge.mat) <- names(report.list)
      }
      if (percentages){
          merge.mat <- reads_matrix_to_percentages(merge.mat)
      }
      merge.mat.list[[filter.tax.level]] <- merge.mat
    }
    return(merge.mat.list)
}

# convenience function for a single tax level
many_files_to_matrix <- function(files, filter.tax.level="S"){
 return(many_files_to_matrix_list(files, filter.tax.levels = c(filter.tax.level),
        include.unclassified=F, report.taxid=F, percentages=F)[[1]])
}

# load up all matching bracken/kraken files in a folder and return a list
# if reads=F, return percentage matrices
# search.string can be specified for kraken/bracken matrices
read_processed_kraken_matrices <- function(dir, reads=TRUE, search.string='bracken'){
    tax.level.names <- tolower(c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
    matrix.list <- list()
    found.files <- list.files(dir, pattern=search.string)
    if(reads){
        reads.str <- 'reads'
    } else {
        reads.str <- 'percentage'
    }
    for (tax.level in tax.level.names){
        file.str <- paste(tax.level, reads.str, sep='_')
        read.file <- file.path(dir, grep(file.str, found.files, value=T))
        mat <- as.matrix(read.table(read.file, sep='\t', quote='', header = T, row.names = 1))
        matrix.list[[tax.level]] <- mat
    }
    return(matrix.list)
}
