# plotting tools for kraken classification results
# uses input matrices from the process_classification.R script
library(RColorBrewer)
library(ggpubr)

# testing color pal stuff
if(F){
    library(rafalib)
    mypar(3,2)
    pie(rep(1,8), col=brewer.pal(8,'Accent'), main='Accent')
    pie(rep(1,8), col=brewer.pal(8,'Dark2'), main='Dark2')
    pie(rep(1,12), col=brewer.pal(12,'Paired'), main='Paired')
    pie(rep(1,9), col=brewer.pal(9,'Set1'), main='Set1')
    pie(rep(1,8), col=brewer.pal(8,'Set2'), main='Set2')
    pie(rep(1,12), col=brewer.pal(12,'Set3'), main='Set3')
    
    mypar(1,1)
    n <- 12
    pie(rep(1,n), col=colorRampPalette(brewer.pal(12,'Paired'))(n))
    pie(rep(1,4), col=colorRampPalette(brewer.pal(12,'Paired')[11:12])(4))
    
    sapply(1:12, function(i) colorRampPalette(brewer.pal(12)[i,i+1])(3))
    
    # Getting all the R color brew pal options
    n <- 20
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    pie(rep(1,n), col=col_vector[1:n])
    
    # old cols
    cols <- c(brewer.pal(12, 'Set3'), 'grey80')
    cols[9] = '#CD6155'

}

# plot many samples in a barplot horizontally
# samples must be in normalized fraction first
# include.unclassified: add unmpped reads into the mix barplot
# expects two rows in the matrix named 'unclassified completely' and 
#  'unclassified at this level'
# plots them in shades of grey and changes y axis name
plot_many_samples <- function(kraken.mat, n.colors=12, scale.name = 'Species', include.unclassified =T){
    if (any(kraken.mat > 100)){
        stop('Must give normalized percent matrix')
    }

    # ncolors must be equal or less to nrow kraken.mat
    n.colors <- min(n.colors, nrow(kraken.mat))
    
    nmax <- 12
    cols <- c(colorRampPalette(brewer.pal(12,'Paired'))(nmax)[1:n.colors], 'grey80')
    # replace 11th color
    cols[11] <-  colorRampPalette(brewer.pal(12,'Paired')[11:12])(4)[2]
    # limit to ncolors 
    cols <- cols[1:(n.colors+1)]
    
    # proceed with classified reads here    
    unclassified.rownames <- c('unclassified completely', 'unclassified at this level')
    kraken.mat.classified <- kraken.mat[!(rownames(kraken.mat) %in% unclassified.rownames),]
    keep.rows <- names(sort(rowSums(kraken.mat.classified), decreasing = T))[1:n.colors]
    kraken.mat.plot <- as.matrix(kraken.mat.classified[keep.rows,])

    # add unclassified taxa if desired
    if (include.unclassified){
        if (!(all(unclassified.rownames %in% rownames(kraken.mat)))){
            stop('Must have unclassified rows in matrix')
        }
        kraken.mat.plot <- rbind(kraken.mat.plot, kraken.mat[unclassified.rownames,])
        cols <- c(cols, 'grey60', 'grey40')
        # kraken.mat.plot['Other',] <- 100 - colSums(kraken.mat.plot)
        
    } else{ 
        # renorm percentages
        kraken.mat.plot <- rbind(kraken.mat.plot, 100 - colSums(kraken.mat.plot))
        kraken.mat.plot <- (kraken.mat.plot / (colSums(kraken.mat.plot)) * 100
    }
    rownames(kraken.mat.plot)[nrow(kraken.mat.plot)] <- 'Other'
    
    
    toplot.df <- melt(kraken.mat.plot, varnames = c('taxa', 'sample'))

    taxplot <- ggplot(data=toplot.df, aes(x = sample, y = value, fill = taxa)) + 
        geom_bar(stat = 'identity') +
        scale_fill_manual(values = cols, name = scale.name) +
        # theme_bw() + 
        # theme(panel.border = element_rect(), axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
              # plot.margin = unit(c(0,0,0,0),'lines')) +
        labs(x='Sample', y='Percent of classified reads') 
        # guides(fill=FALSE)
    # taxplot
    return(taxplot)
}


plot_many_samples_with_diversity <- function(kraken.mat, diversity.df, y.title='Shannon Diversity', ...){
    if (ncol(diversity.df) != 2){
        stop('Must be two column dataframe: sample names in 1, diversity in 2')
    }
    taxplot <- plot_many_samples(kraken.mat, ...)
    diversity.df[,1] <- factor(diversity.df[,1], levels=diversity.df[,1])
    divplot <- ggplot(diversity.df, aes(x=sample, group=1, y=shannon)) +
        geom_point(stat='summary', fun.y=sum, size=5) +
        stat_summary(fun.y=sum, geom="line") +
        # theme_bw() + 
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
              # plot.margin = unit(c(0,0,0,0),'lines')) +
        ylim(0, max(diversity.df$shannon)*1.20) +
        labs(y=y.title)
    # divplot
    
    figure <- ggarrange(divplot, taxplot,ncol=1, align = 'v', nrow=2, heights = c(1,4),legend = 'right', common.legend = T)
    return(figure)
    
}

# figure <- plot_many_samples_with_diversity(bracken.species.fraction, diversity.df)
# annotate_figure(figure, top = text_grob('TEST', size=14))
# figure
# 
# taxplot <- plot_many_samples(bracken.species.fraction, scale.name = 'Species') #+ labs(title = 'TEST')
# diversity.df <- data.frame(sample=names(div.list.species$shannon), shannon=div.list.species$shannon)


plot_many_samples_with_diversity_barplot <- function(kraken.mat, diversity.df, y.title='Shannon Diversity', plot.title='', ...){
    if (ncol(diversity.df) != 2){
        stop('Must be two column dataframe: sample names in 1, diversity in 2')
    }
    colnames(diversity.df) <- c('sample', 'shannon')
    taxplot <- plot_many_samples(kraken.mat, ...)
    diversity.df[,1] <- factor(diversity.df[,1], levels=diversity.df[,1])
    divplot <- ggplot(diversity.df, aes(x=sample, group=1, y=shannon)) +
        # as a bar plot instead?
        geom_bar(stat='identity', width = 0.5) + 
        # theme_bw() + 
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
              # plot.margin = unit(c(0,0,0,0),'lines')) +
        ylim(0, max(diversity.df$shannon)*1.20) +
        labs(y=y.title)
    
    print(taxplot)
    print(divplot)
    
    figure <- ggarrange(divplot, taxplot,ncol=1, align = 'v', nrow=2, heights = c(1,4),legend = 'right', common.legend = T)
    if (plot.title!=''){
        figure <- annotate_figure(figure, top = text_grob(plot.title, size=14), fig.lab.pos='top.left')
        # figure <- annotate_figure(figure, fig.lab=plot.title, fig.lab.pos='top.left')
    }
    return(figure)
    
}