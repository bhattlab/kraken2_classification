# toolkit for processing classification data


# read in a bracken/kraken table and return a dataframe
# at the specified taxonomic classification level
# due to some cases of Species not having a reported genus
# (or other higher level taxonomic classification), when reporting
# not at the species level, reads will be merged to the next highest 
# level available. 
read.kraken.report <- function(fname, filter.tax.level="S"){
    valid.tax.levels <- c('D','P','C','O','F','G','S')
    if (!file.exists(fname)){
        stop('File does not exist!')
    } else if(!(filter.tax.level %in% valid.tax.levels)) {
         stop(paste('filter tax.level must be in', paste(valid.tax.levels, collapse = ', ')))
    }
    df <- read.table(fname, sep='\t', quote='', header = F,
                     col.names=c('pct', 'reads.below','reads.direct','tax.level','taxid','name'))
    
    # take out everything that maps to this level
    df.filter <- df[df$tax.level==filter.tax.level, c('name', 'reads.below')]
    
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
    df.filter <- rbind(df.filter, add.df[,c('name', 'reads.below')])
    return(df.filter)
}
