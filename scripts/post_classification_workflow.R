# Processing kraken results into matrices and plots
# Ben Siranosian - Bhatt lab - Stanford Genetics 
# bsiranosian@gmail.com
# January 2019 - July 2019

suppressMessages(library(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(vegan, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(reshape2, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(RColorBrewer, quietly = TRUE, warn.conflicts = FALSE))

# options we need from snakemake
scripts.folder <- snakemake@params[['scripts_folder']]
sample.reads.f <- snakemake@params[['sample_reads']]
sample.groups.f <- snakemake@params[['sample_groups']]
workflow.outdir <- snakemake@params[['workflow_outdir']]
result.dir <- snakemake@params[['result_dir']]
use.bracken.report <- snakemake@params[['use_bracken_report']]
scripts.folder <- snakemake@scriptdir

# # Testing Args
# scripts.folder <- '~/scg/projects/kraken2_classification/scripts/'
# sample.reads.f <- '~/scg_scratch/ssrc_conventional/kraken2_classify/samples.tsv'
# sample.groups.f <- '~/scg_scratch/ssrc_conventional/kraken2_classify/sample_groups.tsv'
# classification.folder <- '~/scg_scratch/ssrc_conventional/kraken2_classify/kraken2_classification/'
# use.bracken.report <- T

# # for segata debug
# scripts.folder <- '~/scg/projects/kraken2_classification/scripts/'
# sample.reads.f <- '~/scg_lab/transmit_crass/kraken2_classification/samples_2018_unmapped.tsv'
# sample.groups.f <- ''
# classification.folder <- '~/scg_lab/transmit_crass/kraken2_classification/kraken2_classification_2018_unmapped_segata/classification/'
# use.bracken.report <- F
# result.dir <- '~/scg_lab/transmit_crass/kraken2_classification/kraken2_classification_2018_unmapped_segata/processed_results/'
# outfolder.matrices.taxonomy <- file.path(result.dir, 'taxonomy_matrices')
# outfolder.matrices.bray <- file.path(result.dir, 'braycurtis_matrices')
# outfolder.plots <- file.path(result.dir, 'plots')
# # # for  debug
# scripts.folder <- '~/scg/projects/kraken2_classification/scripts/'
# scripts.folder <- '~/projects/kraken2_classification/scripts/'
# sample.reads.f <- '~/scg_lab/transmit_crass/kraken2_classification/samples.tsv'
# sample.groups.f <- ''
# classification.folder <- '~/scg_lab/transmit_crass/kraken2_classification/kraken2_classification_genbank2019/classification/'
# use.bracken.report <- F
# workflow.outdir <- '~/scg_lab/transmit_crass/kraken2_classification/kraken2_classification_genbank2019/'
# result.dir <- '~/scg_lab/transmit_crass/kraken2_classification/kraken2_classification_genbank2019//processed_results_TESTING/'

# # FOR FIONA TESTING
# sample.reads.f <- '~/scg/fiona_debug/3.kraken2/samples.tsv'
# sample.groups.f <- ''
# classification.folder <- '~/scg/fiona_debug/3.kraken2/kraken2_classification/classification/'
# use.bracken.report <- T
# result.dir <- '~/scg/fiona_debug/3.kraken2/kraken2_classification/processed_results_BAS'

# set up directories and make those that don't exist
classification.folder <- file.path(workflow.outdir, 'classification')
# result.dir <- file.path(workflow.outdir, 'processed_results')
outfolder.matrices.taxonomy <- file.path(result.dir, 'taxonomy_matrices')
outfolder.matrices.taxonomy.classified <- file.path(result.dir, 'taxonomy_matrices_classified_only')
outfolder.matrices.bray <- file.path(result.dir, 'braycurtis_matrices')
outfolder.plots <- file.path(result.dir, 'plots')
for (f in c(result.dir, outfolder.matrices.bray, outfolder.matrices.taxonomy,
            outfolder.matrices.taxonomy.classified, outfolder.plots)){
    if (!dir.exists(f)){ dir.create(f, recursive = T)}
}

# load other data processing and plotting scripts
source.script.1 <- file.path(scripts.folder, 'process_classification.R')
source.script.2 <- file.path(scripts.folder, 'plotting_classification.R')
if (!(file.exists(source.script.1) & file.exists(source.script.2))){
    stop('Specify right source script dir')
}
suppressMessages(source(source.script.1))
suppressMessages(source(source.script.2))

# read sample groups file
sample.reads <- read.table(sample.reads.f, sep='\t', quote='', header=F, comment.char = "#", colClasses = 'character')
colnames(sample.reads) <- c('sample', 'r1', 'r2')[1:ncol(sample.reads)]
# ensure we skip first row as 'sample'
if(tolower(sample.reads[1,'sample'] == 'sample')){
    sample.reads <- sample.reads[2:nrow(sample.reads), ]
}

# if a groups file is specified, read it. Otherwise assign everything to one group
if (sample.groups.f != '') {
    sample.groups <- read.table(sample.groups.f, sep='\t', quote='', header=F, comment.char = "#", colClasses = 'character')
    colnames(sample.groups) <- c('sample', 'group')
    # if first sample is sample, discard the line
    if(tolower(sample.groups[1,'sample']) == 'sample'){
        sample.groups <- sample.groups[2:nrow(sample.groups),]
    }
} else {
    sample.groups <- data.frame(sample=sample.reads$sample, group='All')
}


# ensure reads and groups have the same data
if (!(all(sample.groups$sample %in% sample.reads$sample) & all(sample.reads$sample %in% sample.groups$sample))){
    message('sample.reads:')
    print(sample.reads)
    message('sample.groups:')
    print(sample.groups)
    message('sample.groups$sample[!(sample.groups$sample %in% sample.reads$sample)]')
    print(sample.groups$sample[!(sample.groups$sample %in% sample.reads$sample)])
    message('sample.reads$sample[!(sample.reads$sample %in% sample.groups$sample)]')
    print(sample.reads$sample[!(sample.reads$sample %in% sample.groups$sample)])
    stop('Sample reads and sample groups dont contain the same samples... check inputs')
}

# get sample names
if (use.bracken.report){f.ext <- '.krak_bracken.report'} else {f.ext <- '.krak.report'}
flist <- sapply(sample.groups$sample, function(x) file.path(classification.folder, paste(x, f.ext, sep='')))
names(flist) <- sample.groups$sample
if (!(all(file.exists(flist)))){
    # print which files don't exist 
    message('The follwing files do not exist:')
    print(flist[!sapply(flist, file.exists)])
    stop("Some classification files do not exist!")
}

# read in all as a list of matrices instead - much quicker and can do all tax levels
tax.level.names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
tax.level.abbrev <- c('D','P','C','O','F','G','S')
# need something to limit to Segata database
test.df <- kraken_file_to_df(flist[1])
# like number of genus classifications is zero or something
segata <- F
if (sum(sum(test.df$tax.level=='G')) ==0){
    segata <- T
    tax.level.names <- c('Species')
    tax.level.abbrev <- c('S')
}

# load into list of matrices, one for each tax level
message(paste('Loading data from ', length(flist), ' kraken/bracken results.', sep=''))
message(paste('Processing at taxonomic levels: ', paste(tax.level.names, collapse=', '), sep=''))
bracken.reads.matrix.list <- many_files_to_matrix_list(flist, filter.tax.levels = tax.level.abbrev, include.unclassified = T)
# separate matrices including classified and unclassified
# god this naming scheme is getting ugly
unclassified.rownames <- c('Unclassified', 'Classified at a higher level')
# special case for one sample
if (ncol(bracken.reads.matrix.list[[1]]) ==1){
    bracken.reads.matrix.list.classified <- lapply(bracken.reads.matrix.list, function(x) {
        y <- as.matrix(x[!(rownames(x) %in% unclassified.rownames), ])
        colnames(y) <- colnames(x)
        y})
    } else {
    bracken.reads.matrix.list.classified <- lapply(bracken.reads.matrix.list, function(x) x[!(rownames(x) %in% unclassified.rownames),,drop=FALSE])
}

# normalize to percentages
bracken.fraction.matrix.list <- lapply(bracken.reads.matrix.list, reads_matrix_to_percentages)
names(bracken.reads.matrix.list) <- tax.level.names
names(bracken.fraction.matrix.list) <- tax.level.names
# for classified only
bracken.fraction.matrix.list.classified <- lapply(bracken.reads.matrix.list.classified, reads_matrix_to_percentages)
names(bracken.reads.matrix.list.classified) <- tax.level.names
names(bracken.fraction.matrix.list.classified) <- tax.level.names

message('Saving classification matrices...')
# save matrices
if (use.bracken.report){mat.name <- 'bracken'} else {mat.name <- 'kraken'}
for (tn in tax.level.names){
    out.mat.reads <- file.path(outfolder.matrices.taxonomy, paste(mat.name, tolower(tn), 'reads.txt', sep='_'))
    out.mat.fraction <- file.path(outfolder.matrices.taxonomy, paste(mat.name, tolower(tn), 'percentage.txt', sep='_'))
    write.table(bracken.reads.matrix.list[[tn]], out.mat.reads, sep='\t', quote=F, row.names = T, col.names = T)
    # now only write rows of pct matrix that are nonzero
    # that is, at least one sample has 0.001 fraction (due to rounding)
    write.table(bracken.fraction.matrix.list[[tn]][rowSums(bracken.fraction.matrix.list[[tn]])>0,,drop=FALSE], out.mat.fraction, sep='\t', quote=F, row.names = T, col.names = T)
}
# save classified only also
for (tn in tax.level.names){
    out.mat.reads <- file.path(outfolder.matrices.taxonomy.classified, paste(mat.name, tolower(tn), 'reads.txt', sep='_'))
    out.mat.fraction <- file.path(outfolder.matrices.taxonomy.classified, paste(mat.name, tolower(tn), 'percentage.txt', sep='_'))
    write.table(bracken.reads.matrix.list.classified[[tn]], out.mat.reads, sep='\t', quote=F, row.names = T, col.names = T)
    # now only write rows of pct matrix that are nonzero
    # that is, at least one sample has 0.001 fraction (due to rounding)
    keep.idx <- rowSums(bracken.fraction.matrix.list.classified[[tn]])>0
    write.table(bracken.fraction.matrix.list.classified[[tn]][keep.idx,,drop=FALSE], out.mat.fraction, sep='\t', quote=F, row.names = T, col.names = T)
}


#################################################################################
## Diversity calculation and plots ##############################################
#################################################################################
message('Doing diversity calculations and saving figures...')

#plot a rarefaction curve.
plot_rarefaction_curve(bracken.reads.matrix.list, file.path(outfolder.plots, 'rarefaction_curve.pdf'))

# diversity calculations
div.methods <- c('shannon', 'simpson')
div.tax.levels <- tax.level.names
# list of lists, oh my!
div.level.method <- lapply(div.tax.levels, function(x) {
    # if not enough valid rows
    use.matrix <- bracken.reads.matrix.list.classified[[x]]
    if (nrow(use.matrix) <3){
        dl <- lapply(div.methods, function(x) rep(0, times=ncol(use.matrix)))
    } else {
        dl <- lapply(div.methods, function(y) diversity(use.matrix, index=y, MARGIN = 2))
    }
    names(dl) <- div.methods
    dl
})
names(div.level.method) <- div.tax.levels

# there's definitely a better way to do this...
div.df <- as.data.frame(div.level.method)
rownames(div.df) <- sample.reads$sample
div.df <- cbind(data.frame(sample=rownames(div.df)), div.df)
div.df <- melt(div.df, id.vars = 'sample')
div.df$tax.level <- sapply(div.df$variable, function(x) strsplit(as.character(x), split="\\.")[[1]][1])
div.df$method <- sapply(div.df$variable, function(x) strsplit(as.character(x), split="\\.")[[1]][2])
div.df <- div.df[,c('sample', 'tax.level','method','value')]
# round to a sensible number
div.df$value <- round(div.df$value, 3)
# write out dataframe
out.div <- file.path(result.dir, 'diversity.txt')
write.table(div.df, out.div, sep='\t', quote=F, row.names=F, col.names=T)


# barplot of diversity under different methods
# one with everything, one separated by sample group
div.df$sample <- factor(div.df$sample, levels = sample.groups$sample)
div.plot.list <- list()
for (g in unique(sample.groups$group)){
    plot.samples <- sample.groups[sample.groups$group==g, "sample"]
    plot.df <- div.df[div.df$sample %in% plot.samples, ]
    p <- ggplot(plot.df, aes(x=sample, y=value)) +
        geom_bar(stat='identity') +
        facet_grid(rows = vars(tax.level), cols= vars(method), scales = 'fixed') +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(title=paste('Diversity for group:', g),
             y='Diversity', x='Sample')
    div.plot.list[[g]] <- p
}


# save one page for each group
div.group.pdf <- file.path(outfolder.plots, 'diversity_by_group.pdf')
pdf(div.group.pdf, height=8, width = 10)
for(p in div.plot.list){print(p)}
trash <- dev.off()

# also a big figure with one page for each method, all samples
div.plot.all.list <- list()
for (m in unique(div.df$method)){
    plot.df <- div.df[div.df$method == m, ]
    p <- ggplot(plot.df, aes(x=sample, y=value)) +
            geom_bar(stat='identity') +
            facet_grid(rows = vars(tax.level), cols= vars(method), scales = 'fixed') +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            labs(title=paste('Diversity for all samples:', m),
                 y='Diversity', x='Sample')
    div.plot.all.list[[m]] <- p
}

# width of plot = 9 + quarter inch for each sample over 10
max.samps <- max(table(sample.groups$group))
pdf.width <- max(9, (9 + (0.25 * (max.samps-10))))
div.all.pdf <- file.path(outfolder.plots, 'diversity_allsamples.pdf')
pdf(div.all.pdf, height=5, width = pdf.width)
for(p in div.plot.all.list){print(p)}
trash <- dev.off()


#################################################################################
## Taxonomic Barplots ###########################################################
#################################################################################
message('Generating taxonomic barplots...')
# page in the pdf for each sample group
taxlevel.plots <- list()
for (tn in tax.level.names){
    group.plots <- list()
    for (g in unique(sample.groups$group)){
        plot.samples <- sample.groups[sample.groups$group==g, "sample"]
        div.df.sub <- div.df[div.df$tax.level==tn & div.df$method=='shannon' & div.df$sample %in% plot.samples,]
        # print(div.df.sub)
        div.df.plot <- div.df.sub[, c('sample', 'value')]
        rownames(div.df.plot) <- div.df.plot$sample
        # div.df.sub <- div.df.sub[plot.samples, ]
        plot.title <- paste('Taxonomy and diversity: ', g, ', ', tn, sep='')
        plot.mat <- bracken.fraction.matrix.list[[tn]][,plot.samples, drop=FALSE]
        # colnames(plot.mat) <- plot.samples
        # print(head(plot.mat))
        group.plots[[g]] <- plot_many_samples_with_diversity_barplot(plot.mat,
                                                                  div.df.plot, plot.title = plot.title, include.unclassified=T, tax.level.name = tn)
        # plot_many_samples_with_diversity_barplot(bracken.fraction.matrix.list[[tn]][,plot.samples],
        #                                          div.df.plot, plot.title = plot.title, include.unclassified=T)
    }
    taxlevel.plots[[tn]] <- group.plots
}

# save pdfs
# width of plot = 9 + quarter inch for each sample over 10
max.samps <- max(table(sample.groups$group))
pdf.width <- max(9, (9 + (0.25 * (max.samps-10))))
for (tn in tax.level.names){
    tax.pdf <- file.path(outfolder.plots, paste('taxonomy_barplot_', tolower(tn), '.pdf', sep=''))
    pdf(tax.pdf, height=6, width=pdf.width)
    for (p in taxlevel.plots[[tn]]){print(p)}
    trash <- dev.off()
}

# REDO for classified only
taxlevel.plots.classified <- list()
for (tn in tax.level.names){
    group.plots <- list()
    # print(tn)
    for (g in unique(sample.groups$group)){
        plot.samples <- sample.groups[sample.groups$group==g, "sample"]
        div.df.sub <- div.df[div.df$tax.level==tn & div.df$method=='shannon' & div.df$sample %in% plot.samples,]
        # print(div.df.sub)
        div.df.plot <- div.df.sub[, c('sample', 'value')]
        rownames(div.df.plot) <- div.df.plot$sample
        # div.df.sub <- div.df.sub[plot.samples, ]
        plot.title <- paste('Taxonomy and diversity: ', g, ', ', tn, sep='')
        plot.mat <- bracken.fraction.matrix.list.classified[[tn]][,plot.samples, drop=FALSE]

        group.plots[[g]] <- plot_many_samples_with_diversity_barplot(plot.mat,
                                                                  div.df.plot, plot.title = plot.title, include.unclassified=F, tax.level.name = tn)
    }
    taxlevel.plots.classified[[tn]] <- group.plots
}

# save pdfs
# width of plot = 9 + quarter inch for each sample over 10
max.samps <- max(table(sample.groups$group))
pdf.width <- max(9, (9 + (0.25 * (max.samps-10))))
for (tn in tax.level.names){
    tax.pdf <- file.path(outfolder.plots, paste('classified_taxonomy_barplot_', tolower(tn), '.pdf', sep=''))
    pdf(tax.pdf, height=6, width=pdf.width)
    for (p in taxlevel.plots.classified[[tn]]){print(p)}
    trash <- dev.off()
}


#################################################################################
## PCoA plots ###################################################################
#################################################################################
message('Doing PCoA calculations...')
# only do this if we have >=3 samples
if (nrow(sample.reads) >=3){
    # do for each tax level
    plotlist.nolabels <- list()
    plotlist.labels <- list()
    # don't do for domain
    do.tn <- tax.level.names[tax.level.names!='Domain']
    for (tn in do.tn){
        # print(paste('    for:', tn))
        fraction.mat <- bracken.fraction.matrix.list.classified[[tn]]
        pcoa.res <- capscale(t(fraction.mat)~1, distance='bray')
        pcoa.df <- data.frame(sample.groups, scores(pcoa.res)$sites)
        pcoa.df.melt <- melt(pcoa.df, id.vars = c('sample','group'))
        pcoa.variance <- summary(pcoa.res)$cont$importance[2,1:2]
        # set colors based on number of groups
        ng <- length(unique(sample.groups$group))
        if (ng <10){cols <- brewer.pal(max(ng,3), 'Set1')} else {cols <- colorRampPalette(brewer.pal(9,'Set1'))(ng)}

        plotlist.nolabels[[tn]] <- ggplot(pcoa.df, aes(x=MDS1, y=MDS2, color=group, label=sample)) +
            geom_point(size=3) +
            # scale_color_brewer(palette='Set1') +
            scale_color_manual(values = cols) +
            theme_bw() +
            labs(title=paste('Microbiome beta diversity, Bray-Curtis, ', tn, sep=''),
                 subtitle='Principal Coordinates Analysis plot',
                 x = paste('PC1 (', round(pcoa.variance[1], 3) * 100, '% of variance)', sep=''),
                 y = paste('PC2 (', round(pcoa.variance[2], 3) * 100, '% of variance)', sep='')) +
            xlim(c(min(pcoa.df$MDS1), max(pcoa.df$MDS1) * 1.25)) +
            ylim(c(min(pcoa.df$MDS2), max(pcoa.df$MDS2) * 1.15)) +
            guides(color=guide_legend(title='Group'))

        # add a version with labels above the points
        nudge.x <- sum(abs(range(pcoa.df$MDS1))) / 12
        nudge.y <- sum(abs(range(pcoa.df$MDS2))) / 30
        plotlist.labels[[tn]] <- plotlist.nolabels[[tn]] + geom_text(nudge_x = nudge.x, nudge_y = nudge.y)
    }

    # save plot
    pcoa.pdf.nolabels <- file.path(outfolder.plots, 'PCoA_2D_plot_nolabels.pdf')
    pdf(pcoa.pdf.nolabels, height=6.5, width=9)
    for (tn in rev(do.tn)){
        print(plotlist.nolabels[[tn]])
    }
    trash <- dev.off()
    pcoa.pdf.labels <- file.path(outfolder.plots, 'PCoA_2D_plot_labels.pdf')
    pdf(pcoa.pdf.labels, height=6.5, width=9)
    for (tn in rev(do.tn)){
        print(plotlist.labels[[tn]])
    }
    trash <- dev.off()
} else {
    # warn the user and make fake plots
    warning('Less than 3 samples, not doing PCoA plots')
    pcoa.pdf.nolabels <- file.path(outfolder.plots, 'PCoA_2D_plot_nolabels.pdf')
    pdf(pcoa.pdf.nolabels, height=6.5, width=9)
    plot(0,0, col='white')
    text(0,0, 'Not enough samples for PCoA plot', col='firebrick')
    trash <- dev.off()
    pcoa.pdf.labels <- file.path(outfolder.plots, 'PCoA_2D_plot_labels.pdf')
    pdf(pcoa.pdf.labels, height=6.5, width=9)
    plot(0,0, col='white')
    text(0,0, 'Not enough samples for PCoA plot', col='firebrick')
    trash <- dev.off()
}

# simple bray curtis distance metrics
for (tn in tax.level.names){
    bray.dist <- as.matrix(vegdist(t(bracken.fraction.matrix.list.classified[[tn]])))
    out.bray <- file.path(outfolder.matrices.bray, paste('braycurtis_distance_', tolower(tn), '.txt', sep=''))
    write.table(bray.dist, out.bray, sep='\t', quote=F, row.names = T, col.names = T)
}

message('Done! :D')
