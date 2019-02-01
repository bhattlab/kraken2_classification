# what to do with kraken results once they're done any set of samples
library(ggplot2)
library(vegan)
library(reshape2)

# 1) parse into matrix at species and genus level
# 2) diversity calculations at species and genus level
# 3) taxonomic barplot with diversity at species and genus level
# 4) pcoa plot?

# options we need from snakemake
scripts.folder <- snakemake@params[['scripts_folder']]
sample.reads.f <- snakemake@params[['sample_reads']]
sample.groups.f <- snakemake@params[['sample_groups']]
workflow.outdir <- snakemake@params[['outdir']]
bracken.report <- snakemake@params[['bracken_report']]
classification.folder <- file.path(workflow.outdir, 'classification')
outfolder.matrices <- file.path(workflow.outdir, 'processed_results')
outfolder.plots <- file.path(outfolder.matrices, 'plots')

# # Testing Args
# scripts.folder <- '~/projects/kraken2_classification/scripts/'
# sample.reads.f <- '~/scg_scratch/ssrc_conventional/kraken2_classify/samples.tsv' 
# sample.groups.f <- '~/scg_scratch/ssrc_conventional/kraken2_classify/sample_groups.tsv' 
# classification.folder <- '~/scg_scratch/ssrc_conventional/kraken2_classify/kraken2_classification/'
# bracken.report <- T
# outfolder.matrices <- '~/scg_scratch/ssrc_conventional/kraken2_classify/processed_results'
# outfolder.plots <- file.path(outfolder.matrices, 'plots')

# source other scripts for processing and plotting
source.script.1 <- file.path(scripts.folder, 'process_classification.R')
source.script.2 <- file.path(scripts.folder, 'plotting_classification.R')
if (!(file.exists(source.script.1) & file.exists(source.script.2))){
    stop('Specify right source script dir')
}
source(source.script.1)
source(source.script.2)
if (!dir.exists(outfolder.plots)){ dir.create(outfolder.plots, recursive = T)}

# read sample groups file
sample.reads <- read.table(sample.reads.f, sep='\t', quote='', header=F, comment.char = "#", colClasses = 'character')
colnames(sample.reads) <- c('sample', 'r1', 'r2')[1:ncol(sample.reads)]
# if a groups file is specified, read it. Otherwise assign everything to one group
if (sample.groups.f != '') {
    sample.groups <- read.table(sample.groups.f, sep='\t', quote='', header=F, comment.char = "#", colClasses = 'character')
    colnames(sample.groups) <- c('sample', 'group')
} else {
    sample.groups <- data.frame(sample=sample.reads$sample, group='All')
}
# ensure reads and groups have the same data
if (!(all(sample.groups$sample %in% sample.reads$sample) & all(sample.reads$sample %in% sample.groups$sample))){
    stop('Sample reads and sample groups dont contain the same samples... check inputs')
}

# get sample names
if (bracken.report){f.ext <- '.krak_bracken.report'} else {f.ext <- '.krak.report'}
flist <- sapply(sample.groups$sample, function(x) file.path(classification.folder, paste(x, f.ext, sep='')))
names(flist) <- sample.groups$sample
if (!(all(file.exists(flist)))){
    stop("Some classification files do not exist!")
}

# read in to matrices
# print(flist)
bracken.species.reads <- many_files_to_matrix(flist)
bracken.genus.reads <- many_files_to_matrix(flist,filter.tax.level = 'G')
bracken.species.fraction <- reads_matrix_to_percentages(bracken.species.reads)
bracken.genus.fraction <- reads_matrix_to_percentages(bracken.genus.reads)
# save matrices
if (bracken.report){mat.name <- 'bracken'} else {mat.name <- 'kraken'}
out.mat.1 <- file.path(outfolder.matrices, paste(mat.name, 'species_reads.txt', sep='_'))
out.mat.2 <- file.path(outfolder.matrices, paste(mat.name, 'genus_reads.txt', sep='_'))
out.mat.3 <- file.path(outfolder.matrices, paste(mat.name, 'species_percentage.txt', sep='_'))
out.mat.4 <- file.path(outfolder.matrices, paste(mat.name, 'genus_percentage.txt', sep='_'))
write.table(bracken.species.reads, out.mat.1, sep='\t', quote=F, row.names = T, col.names = T)
write.table(bracken.genus.reads, out.mat.2, sep='\t', quote=F, row.names = T, col.names = T)
write.table(bracken.species.fraction, out.mat.3, sep='\t', quote=F, row.names = T, col.names = T)
write.table(bracken.genus.fraction, out.mat.4, sep='\t', quote=F, row.names = T, col.names = T)


#################################################################################
## Diversity calculation and plots ##############################################
#################################################################################
# diversity calculations
div.methods <- c('shannon', 'simpson')
div.list.species <- lapply(div.methods, function(x) diversity(bracken.species.reads, index=x, MARGIN = 2))
div.list.genus <- lapply(div.methods, function(x) diversity(bracken.genus.reads, index=x, MARGIN = 2))
div.df <- cbind(do.call(cbind, div.list.species), do.call(cbind, div.list.genus))
colnames(div.df) <- paste(div.methods, rep(c('species', 'genus'), each=length(div.methods)), sep='_')
# write out dataframe
out.div <- file.path(outfolder.matrices, 'diversity.txt')
write.table(round(div.df,3), out.div, sep='\t', quote=F, row.names=T, col.names=T)

# barplot of diversity under different methods
# one with everything, one separated by sample group
div.df.melt <- melt(div.df)
colnames(div.df.melt) <- c('sample','method','value')
div.df.melt$sample <- factor(div.df.melt$sample, levels = sample.groups$sample)
div.plot.list <- list()
# print(sample.groups)
# print(div.df.melt)
for (g in unique(sample.groups$group)){
    plot.samples <- sample.groups[sample.groups$group==g, "sample"]
    plot.df <- div.df.melt[div.df.melt$sample %in% plot.samples, ]
    p <- ggplot(plot.df, aes(x=sample, y=value)) + 
            geom_bar(stat='identity') + 
            facet_wrap(plot.df$method,dir = 'v', ncol=2, scales = 'fixed') +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            labs(title=paste('Diversity for group:', g),
                 y='Diversity', x='Sample')
    div.plot.list[[g]] <- p
}
# save one page for each group
div.group.pdf <- file.path(outfolder.plots, 'diversity_by_group.pdf')
pdf(div.group.pdf, height=8, width = 10)
for(p in div.plot.list){print(p)}
dev.off()

# also a big figure with one page for each method, all samples
div.plot.all.list <- list()
for (m in unique(div.df.melt$method)){
    plot.df <- div.df.melt[div.df.melt$method == m, ]
    p <- ggplot(plot.df, aes(x=sample, y=value)) + 
            geom_bar(stat='identity') + 
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
dev.off()

#################################################################################
## Taxonomic Barplots ###########################################################
#################################################################################
# page in the pdf for each sample group
tax.pdf.species <- file.path(outfolder.plots, 'taxonomy_barplot_species.pdf')
tax.pdf.genus <- file.path(outfolder.plots, 'taxonomy_barplot_genus.pdf')
plotlist.species <- list()
plotlist.genus <- list()
for (g in unique(sample.groups$group)){
    plot.samples <- sample.groups[sample.groups$group==g, "sample"]
    # species level
    div.df.g.species <- data.frame(plot.samples, div.df[plot.samples, 'shannon_species'])
    plot.title.species <- paste('Taxonomy and diversity:', g, 'species')
    plotlist.species[[g]] <- plot_many_samples_with_diversity_barplot(bracken.species.fraction[,plot.samples], div.df.g.species, plot.title = plot.title.species)
    # genus level
    div.df.g.genus <- data.frame(plot.samples, div.df[plot.samples, 'shannon_genus'])
    plot.title.genus <- paste('Taxonomy and diversity:', g, 'genus')
    plotlist.genus[[g]] <- plot_many_samples_with_diversity_barplot(bracken.genus.fraction[,plot.samples], div.df.g.genus, plot.title = plot.title.genus, scale.name='Genus')
}

# width of plot = 9 + quarter inch for each sample over 10
max.samps <- max(table(sample.groups$group))
pdf.width <- max(9, (9 + (0.25 * (max.samps-10))))
pdf(tax.pdf.species, height=6, width=pdf.width)
for (p in plotlist.species){print(p)}
dev.off()
pdf(tax.pdf.genus, height=6, width=pdf.width)
for (p in plotlist.genus){print(p)}
dev.off()

#################################################################################
## PCoA plots ###################################################################
#################################################################################
pcoa.res <- capscale(t(bracken.species.fraction)~1, distance='bray')
pcoa.df <- data.frame(sample.groups, scores(pcoa.res)$sites)
pcoa.df.melt <- melt(pcoa.df)
pcoa.variance <- summary(pcoa.res)$cont$importance[2,1:2]
# set colors based on number of groups
ng <- length(unique(sample.groups$group))
if (ng <10){cols <- brewer.pal(max(ng,3), 'Set1')} else {cols <- colorRampPalette(brewer.pal(9,'Set1'))(ng)}

g1 <- ggplot(pcoa.df, aes(x=MDS1, y=MDS2, color=group, label=sample)) +
    geom_point(size=3) +
    # scale_color_brewer(palette='Set1') +
    scale_color_manual(values = cols) +
    theme_bw() +
    labs(title='Microbiome beta diversity (Bray-Curtis metric)',
         subtitle='Principal Coordinates Analysis plot',
         x = paste('PC1 (', round(pcoa.variance[1], 3) * 100, '% of variance)', sep=''),
         y = paste('PC2 (', round(pcoa.variance[2], 3) * 100, '% of variance)', sep='')) + 
    xlim(c(min(pcoa.df$MDS1), max(pcoa.df$MDS1) * 1.25)) +
    ylim(c(min(pcoa.df$MDS2), max(pcoa.df$MDS2) * 1.1)) +
    guides(color=guide_legend(title='Group'))

# add a version with labels above the points
nudge.x <- sum(abs(range(pcoa.df$MDS1))) / 12
nudge.y <- sum(abs(range(pcoa.df$MDS2))) / 30
g2 <- g1 + geom_text(nudge_x = nudge.x, nudge_y = nudge.y)

# save plot
pcoa.pdf <- file.path(outfolder.plots, 'PCoA_2D_plot.pdf')
pdf(pcoa.pdf, height=6.5, width=9)
print(g1)
print(g2)
dev.off()

