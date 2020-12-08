# Processing kraken results into matrices and plots
# Ben Siranosian - Bhatt lab - Stanford Genetics
# bsiranosian@gmail.com
# January 2019 - November 2019

suppressMessages(library(ggplot2, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(rafalib, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(vegan, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(reshape2, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(RColorBrewer, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(cmapR, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(compositions, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(zCompositions, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(ALDEx2, quietly = TRUE, warn.conflicts = FALSE))
suppressMessages(library(ggpubr, quietly = TRUE, warn.conflicts = FALSE))
options(stringsAsFactors = F)
# suppressMessages(library(CoDaSeq, quietly = TRUE, warn.conflicts = FALSE))

# options we need from snakemake
# scripts.folder <- snakemake@params[['scripts_folder']]
sample.reads.f <- snakemake@params[['sample_reads']]
sample.groups.f <- snakemake@params[['sample_groups']]
workflow.outdir <- snakemake@params[['workflow_outdir']]
result.dir <- snakemake@params[['result_dir']]
use.bracken.report <- snakemake@params[['use_bracken_report']]
# now just use locations in the working directory
scripts.folder <- "scripts"
# scripts.folder <- snakemake@scriptdir
tax.array.file <- 'taxonomy_array.tsv'
# tax.array.file <- snakemake@input[['tax_array']]
print(paste('scriptdir:', scripts.folder))
############ Testing Args ############################################################################
# sample.reads.f <- '~/bhatt_local/kraken2_testing/small_hct_dataset/samples.tsv'
# sample.groups.f <- '~/bhatt_local/kraken2_testing/small_hct_dataset/sample_groups.tsv'
# # sample.reads.f <- '~/bhatt_local/kraken2_testing/small_hct_dataset/samples_1.tsv'
# # sample.groups.f <- ''
# workflow.outdir <- '~/bhatt_local/kraken2_testing/small_hct_dataset/kraken2_classification_feb2019/'
# result.dir <- '~/bhatt_local/kraken2_testing/small_hct_dataset/test_processing_out/'
# # result.dir <- '~/bhatt_local/kraken2_testing/small_hct_dataset/test_processing_out_1/'
# use.bracken.report <- T
# scripts.folder <- '~/projects/kraken2_classification/scripts/'
# tax.array.file <- '~/bhatt_local/kraken2_testing/taxonomy_parsing/tax_array.tsv'

# # chris testing
# sample.reads.f <- '~/scg/Ctrl_preNov2019/controls_new.txt'
# sample.groups.f <- ''
# workflow.outdir <- '~/scg/Ctrl_preNov2019/'
# result.dir <- '~/bhatt_local/kraken2_testing/small_hct_dataset/test_processing_out2/'
# use.bracken.report <- T
# scripts.folder <- '~/projects/kraken2_classification/scripts/'
# tax.array.file <- '~/bhatt_local/kraken2_testing/taxonomy_parsing/tax_array.tsv'

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
######################################################################################################

# set up directories and make those that don't exist
classification.folder <- file.path(workflow.outdir, 'classification')
# result.dir <- file.path(workflow.outdir, 'processed_results')
outfolder.matrices.taxonomy <- file.path(result.dir, 'taxonomy_matrices')
outfolder.matrices.taxonomy.classified <- file.path(result.dir, 'taxonomy_matrices_classified_only')
outfolder.gctx.taxonomy <- file.path(result.dir, 'taxonomy_gctx')
outfolder.gctx.taxonomy.classified <- file.path(result.dir, 'taxonomy_gctx_classified_only')
outfolder.matrices.bray <- file.path(result.dir, 'braycurtis_matrices')
outfolder.plots <- file.path(result.dir, 'plots')
outfolder.aldex <- file.path(result.dir, 'ALDEx2_differential_abundance')

for (f in c(result.dir, outfolder.matrices.bray, outfolder.matrices.taxonomy.classified,
            outfolder.gctx.taxonomy.classified, outfolder.plots, outfolder.aldex)){
    if (!dir.exists(f)){ dir.create(f, recursive = T)}
}

# dont create dirs we dont need from bracken
if (!(use.bracken.report)){
    for (f in c(outfolder.matrices.taxonomy, outfolder.gctx.taxonomy)){
        if (!dir.exists(f)){ dir.create(f, recursive = T)}
    }
}

# load other data processing and plotting scripts
source.script.process <- file.path(scripts.folder, 'process_classification_gctx.R')
source.script.plot <- file.path(scripts.folder, 'plotting_classification.R')
source.script.codaseq <- file.path(scripts.folder, 'CoDaSeq_functions.R')
if (!(file.exists(source.script.process) & file.exists(source.script.plot) & file.exists(source.script.codaseq))) {
    # if these don't exist it could be due to a singularity error.
    # Try loading from a backup location on scg
    warning('Cannot find processing scripts in the scripts dir relative to this snakefile. Attempting to load processing scripts from backup directory.... (/oak/stanford/scg/lab_asbhatt/tools/kraken2_classification/scripts)')
    scripts.folder <- '/oak/stanford/scg/lab_asbhatt/tools/kraken2_classification/scripts'
    source.script.process <- file.path(scripts.folder, 'process_classification_gctx.R')
    source.script.plot <- file.path(scripts.folder, 'plotting_classification.R')
    source.script.codaseq <- file.path(scripts.folder, 'CoDaSeq_functions.R')
    if (!(file.exists(source.script.process) & file.exists(source.script.plot) & file.exists(source.script.codaseq))){
            stop('processing scripts do not exist at any of the attempted directories. Exiting. ')
    }
}
suppressMessages(source(source.script.process))
suppressMessages(source(source.script.plot))
suppressMessages(source(source.script.codaseq))

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

# get taxonomy array
# classification at each level for each entry in the database.
# simplified after processing krakens default taxonomy
tax.array <- read.table(tax.array.file, sep='\t', quote='', header=F, comment.char = '', colClasses = 'character')
colnames(tax.array) <- c('id', 'taxid', 'root', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'subspecies')
# Bug in generation code gave muliple zero taxids. Eliminate that here now
tax.array <- tax.array[!duplicated(tax.array$taxid),]
rownames(tax.array) <- tax.array$taxid
# some duplicate names in this. if so, change them to include the tax level
dup.ids <- tax.array$id[duplicated(tax.array$id)]
dup.inds <- which(tax.array$id %in% dup.ids)
# fix these by adding taxid to the name
tax.array[tax.array$id %in% dup.ids, "id"] <-
    paste(tax.array[tax.array$id %in% dup.ids, "id"], ' (', tax.array[tax.array$id %in% dup.ids, "taxid"], ')', sep='')

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

# special case for MAG databases.
# need something to limit to Segata database
test.df <- kraken_file_to_df(flist[1])
# like number of genus classifications is zero or something
segata <- sum(sum(test.df$tax.level=='G')) == 0

# load classification results from each sample and process into gct format
message(paste('Loading data from ', length(flist), ' kraken/bracken results.', sep=''))
df.list <- lapply(flist, function(x) kraken_file_to_df(x))
# merge classification matrix
merge.mat <- merge_kraken_df_list(df.list)
# sample metadata just has groups for now
sample.metadata <- data.frame(id=sample.groups$sample, group=sample.groups$group)
# construct gct object from reads matrix, sample metadata, row metadata
kgct <- make_gct_from_kraken(merge.mat, sample.metadata, tax.array)

# filter to each of the taxonomy levels
if (segata){
    filter.levels <- c('species')
    kgct.filtered.list <- list(species=subset.gct(kgct, rid=kgct@rid[kgct@rid != 'root']))
} else {
    filter.levels <-  c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    kgct.filtered.list <- lapply(filter.levels, function(level) subset_kgct_to_level(kgct, level))
    names(kgct.filtered.list) <- filter.levels
}
# print("kgct.filtered.list[['species']]@mat")
# print(kgct.filtered.list[['species']]@mat)

# subset to classified taxa only
unclassified.rownames <- c('unclassified', 'classified at a higher level')
kgct.filtered.classified.list <- lapply(kgct.filtered.list, function(x) subset.gct(x, rid=x@rid[!(x@rid %in% unclassified.rownames)]))

# normalize to percentages
min.frac <- 0.001
kgct.filtered.percentage.list <- lapply(kgct.filtered.list, function(x) normalize_kgct(x, min.frac=min.frac))
kgct.filtered.classified.percentage.list <- lapply(kgct.filtered.classified.list, function(x) normalize_kgct(x, min.frac=min.frac))

message('Saving classification matrices...')
# save matrices and gctx objects
if (use.bracken.report){mat.name <- 'bracken'} else {mat.name <- 'kraken'}
for (tn in filter.levels){
    outf.mat.reads <- file.path(outfolder.matrices.taxonomy, paste(mat.name, tolower(tn), 'reads.txt', sep='_'))
    outf.mat.percentage <- file.path(outfolder.matrices.taxonomy, paste(mat.name, tolower(tn), 'percentage.txt', sep='_'))
    outf.mat.reads.classified <- file.path(outfolder.matrices.taxonomy.classified, paste(mat.name, tolower(tn), 'reads.txt', sep='_'))
    outf.mat.percentage.classified <- file.path(outfolder.matrices.taxonomy.classified, paste(mat.name, tolower(tn), 'percentage.txt', sep='_'))
    # gctx names
    outf.gctx.reads <- file.path(outfolder.gctx.taxonomy, paste(mat.name, tolower(tn), 'reads.gctx', sep='_'))
    outf.gctx.percentage <- file.path(outfolder.gctx.taxonomy, paste(mat.name, tolower(tn), 'percentage.gctx', sep='_'))
    outf.gctx.reads.classified <- file.path(outfolder.gctx.taxonomy.classified, paste(mat.name, tolower(tn), 'reads.gctx', sep='_'))
    outf.gctx.percentage.classified <- file.path(outfolder.gctx.taxonomy.classified, paste(mat.name, tolower(tn), 'percentage.gctx', sep='_'))

    # save matrices
    if (!use.bracken.report){
        write.table(kgct.filtered.list[[tn]]@mat, outf.mat.reads, sep='\t', quote=F, row.names = T, col.names = T)
    }
    write.table(kgct.filtered.classified.list[[tn]]@mat, outf.mat.reads.classified, sep='\t', quote=F, row.names = T, col.names = T)

    mat.percentage <- kgct.filtered.percentage.list[[tn]]@mat
    # order rows by mean abundance across samples
    row.order <- rownames(mat.percentage)[order(rowMeans(mat.percentage), decreasing = T)]
    row.order <- row.order[!(row.order %in% unclassified.rownames)]
    row.order <- c(unclassified.rownames[unclassified.rownames %in% rownames(mat.percentage)], row.order)
    mat.percentage <- mat.percentage[row.order,]
    # classified only matrix
    mat.percentage.classified <- kgct.filtered.classified.percentage.list[[tn]]@mat
    row.order.classified <- order(rowMeans(mat.percentage.classified), decreasing = T)
    mat.percentage.classified <- mat.percentage.classified[row.order.classified,]

    if (!use.bracken.report){
        write.table(mat.percentage, outf.mat.percentage, sep='\t', quote=F, row.names = T, col.names = T)
    }
    write.table(mat.percentage.classified, outf.mat.percentage.classified, sep='\t', quote=F, row.names = T, col.names = T)

    # save gctx
    # only save with unclassified if not using Bracken
    if (!use.bracken.report){
        suppressMessages(write.gctx(kgct.filtered.list[[tn]], outf.gctx.reads, appenddim = F))
        suppressMessages(write.gctx(kgct.filtered.percentage.list[[tn]], outf.gctx.percentage, appenddim = F))
    }
    suppressMessages(write.gctx(kgct.filtered.classified.list[[tn]], outf.gctx.reads.classified, appenddim = F))
    suppressMessages(write.gctx(kgct.filtered.classified.percentage.list[[tn]], outf.gctx.percentage.classified, appenddim = F))
}


#################################################################################
## Diversity calculation and plots ##############################################
#################################################################################
message('Doing diversity calculations and saving figures...')

# DISABLE rarefaction curve as it's useless with so many species
# Wait actually people want it so Im going to re-enable it
plot_rarefaction_curve(kgct.filtered.classified.list$species@mat,
                       file.path(outfolder.plots, 'rarefaction_curve.pdf'))

# diversity calculations
div.methods <- c('shannon', 'simpson')
# calculate for each tax level
div.level.method <- lapply(filter.levels, function(x) {
    use.matrix <- kgct.filtered.classified.list[[x]]@mat
    # if not enough valid rows, skip it
    if (nrow(use.matrix) <3){
        dl <- lapply(div.methods, function(x) rep(0, times=ncol(use.matrix)))
    } else {
        dl <- lapply(div.methods, function(y) diversity(use.matrix, index=y, MARGIN = 2))
    }
    names(dl) <- div.methods
    dl
})
names(div.level.method) <- filter.levels

# coerce forcefully into a dataframe
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
    plot.df$sample <- factor(plot.df$sample, levels = sample.groups$sample)
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
taxlevel.plots.classified <- list()
for (tn in filter.levels){
    group.plots <- list()
    group.plots.classified <- list()
    for (g in unique(sample.groups$group)){
        plot.samples <- sample.groups[sample.groups$group==g, "sample"]
        div.df.sub <- div.df[div.df$tax.level==tn & div.df$method=='shannon' & div.df$sample %in% plot.samples,]
        div.df.plot <- div.df.sub[, c('sample', 'value')]
        rownames(div.df.plot) <- div.df.plot$sample
        plot.title <- paste('Taxonomy and diversity: ', g, ', ', tn, sep='')
        plot.mat <- kgct.filtered.percentage.list[[tn]]@mat[,plot.samples, drop=FALSE]

        # no longer doing include unclassifed with Bracken report
        # because new version of report doesn't include this data
        if (!use.bracken.report){
        group.plots[[g]] <- plot_many_samples_with_diversity_barplot(plot.mat,
                            div.df.plot, plot.title = plot.title, include.unclassified=T, tax.level.name = tn)
        }
        group.plots.classified[[g]] <- plot_many_samples_with_diversity_barplot(plot.mat,
                            div.df.plot, plot.title = plot.title, include.unclassified=F, tax.level.name = tn)
    }
    taxlevel.plots[[tn]] <- group.plots
    taxlevel.plots.classified[[tn]] <- group.plots.classified
}

# save pdfs
# width of plot = 9 + quarter inch for each sample over 10
max.samps <- max(table(sample.groups$group))
pdf.width <- max(9, (9 + (0.25 * (max.samps-10))))
for (tn in filter.levels){
    if (!use.bracken.report){
        tax.pdf <- file.path(outfolder.plots, paste('taxonomy_barplot_', tolower(tn), '.pdf', sep=''))
        pdf(tax.pdf, height=6, width=pdf.width)
        for (p in taxlevel.plots[[tn]]){print(p)}
        trash <- dev.off()
    }
    tax.pdf.classified <- file.path(outfolder.plots, paste('classified_taxonomy_barplot_', tolower(tn), '.pdf', sep=''))
    pdf(tax.pdf.classified, height=6, width=pdf.width)
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
    # don't do for kingdom
    do.tn <- filter.levels[filter.levels!='kingdom']
    for (tn in do.tn){
        print(paste('    for:', tn))
        fraction.mat <- kgct.filtered.classified.percentage.list[[tn]]@mat
        if(nrow(fraction.mat) > 2){
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
        } else{
            plotlist.nolabels[[tn]] <- plot(1,1)
            plotlist.labels[[tn]] <- plot(1,1)
        }
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
for (tn in filter.levels){
    bray.dist <- as.matrix(vegdist(t(kgct.filtered.classified.percentage.list[[tn]]@mat)))
    out.bray <- file.path(outfolder.matrices.bray, paste('braycurtis_distance_', tolower(tn), '.txt', sep=''))
    write.table(bray.dist, out.bray, sep='\t', quote=F, row.names = T, col.names = T)
}

#################################################################################
## Compositional data analysis and plots ########################################
#################################################################################
if (length(unique(sample.groups$group)) != 2){
    warning('Will not do ALDEx2 differential abundance with !=2 groups')}

if (nrow(sample.groups) < 3){
    warning('Will not do compositional data analysis with < 3 samples')
} else {
    pca.plot.list <- list()
    if(segata){
        do.tax.levels <- c('species')
    } else {
        do.tax.levels <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    }

    print('Compositional data analysis....')
    for (tax.level in do.tax.levels[do.tax.levels != 'kingdom']){
        print(paste('.....', tax.level))
        # AT A SPECIFIC TAX LEVEL: filtering
        # keep only those samples with > min.reads
        min.reads <- 5000
        # keep only OTUs with an abundance of at least 0.01
        min.prop = 0.01
        # keep OTUs that are found in at least 30% of samples
        cutoff = .3
        use.mat <- kgct.filtered.classified.list[[tax.level]]@mat
        if(nrow(use.mat) > 2){
            reads.filtered <- tryCatch(codaSeq.filter(use.mat,
                min.reads=min.reads, min.occurrence=cutoff, min.prop=min.prop, samples.by.row=FALSE),
            error=function(e) matrix(0))
            } else {
                reads.filtered <- matrix(0)
            }
        if(nrow(reads.filtered) > 2){
            # replace 0 values with an estimate of the probability that the zero is not 0
            # only if necessary
            # at this point samples are by row, and variables are by column
            # we can use the GBM or CZM methods. Function from zCompositions.
            if (sum(reads.filtered==0)>0){
                rfz <- t(cmultRepl(t(reads.filtered),  label=0, method="CZM", output="p-counts"))
                # round to integers and round anything low to 1
                rfz <- round(rfz)
                rfz[rfz==0] <- 1
            } else {
                rfz <- reads.filtered
            }

            # convert to CLR
            clr.aldex <- aldex.clr(rfz, conds=NA, mc.samples=128, denom = 'all', verbose = F)
            # get clr from mean of the MC instances
            rfz.clr <- t(sapply(clr.aldex@analysisData , function(x) rowMeans(x)))
            # save CLR matrix
            outf <- file.path(outfolder.matrices.taxonomy.classified, paste('clr_values_', tax.level, '.tsv', sep=''))
            write.table(round(rfz.clr, 4), outf, sep='\t', quote=F, row.names = T, col.names = T)

            rfz.clr.pca <- prcomp(rfz.clr)
            rfz.clr.mvar <- mvar(rfz.clr)

            # PCA biplot
            plot.df <- data.frame(rfz.clr.pca$x[,1:2], group=sample.groups[sample.groups$sample %in% rownames(rfz.clr), 'group'])
            pc1.var <- round(sum(rfz.clr.pca$sdev[1]^2)/rfz.clr.mvar * 100, 1)
            pc2.var <- round(sum(rfz.clr.pca$sdev[2]^2)/rfz.clr.mvar * 100, 1)
            pca.plot <- ggplot(plot.df, aes(x=PC1, y=PC2, col=group)) +
                geom_point(size=3) +
                geom_text(label=rownames(plot.df),
                          nudge_x = sum(abs(range(plot.df$PC1)))/ 50,
                          nudge_y = sum(abs(range(plot.df$PC2)))/ 50) +
                theme_bw() +
                labs(x=paste("PC1: ", pc1.var, "% of variace", sep=""),
                     y=paste("PC2: ", pc2.var, "% of variace", sep=""),
                     title=paste('Compositional PCA,', tax.level)) +
                scale_color_brewer(palette='Set1')
            pca.plot.list[[tax.level]] <- pca.plot


            # need to limit to two groups for aldex
            if (length(unique(sample.groups$group)) == 2){
                # get two groups
                group.pos <- unique(sample.groups$group)[1]
                group.neg <- unique(sample.groups$group)[2]
                s.pos <- sample.groups$sample[sample.groups$group==group.pos]
                s.neg <- sample.groups$sample[sample.groups$group==group.neg]

                # aldex differential expression between groups
                aldex.res <- aldex(rfz[,c(s.pos, s.neg)],
                    conditions = c(rep(group.pos, length(s.pos)), rep(group.neg, length(s.neg))),
                    mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE,
                    denom="all", verbose=FALSE)

                # scatterplot
                outf.scatter <- file.path(outfolder.aldex, paste('aldex_scatter_', tax.level, '.pdf', sep=''))
                pdf(outf.scatter, width = 7, height = 4)
                mypar(1,2, mar = c(3.5, 3.5, 2.5, 1.1))
                aldex.plot(aldex.res, type="MW", test="wilcox", xlab="Dispersion",ylab="Difference")
                mtext(paste(tax.level, ', wilcox', sep=''), line=0.25, cex=1.25)
                aldex.plot(aldex.res, type="MW", test="welch", xlab="Dispersion",ylab="Difference")
                mtext(paste(tax.level, ', welch', sep=''), line=0.25, cex=1.25)
                dev.off()

                # improve dataframe
                aldex.res <- signif(aldex.res, 4)
                # add sample numbers
                aldex.res$n.pos <- length(s.pos)
                aldex.res$n.neg <- length(s.neg)
                # add in name
                aldex.res$taxa <- rownames(aldex.res)
                # reorder
                aldex.res <- aldex.res[, c('taxa', 'n.pos', 'n.neg', colnames(aldex.res)[1:11])]
                aldex.res <- aldex.res[order(aldex.res$we.eBH),]
                # add abs effect
                aldex.res$abs.effect <- abs(aldex.res$effect)

                # save dataframe of differential results
                outf.res <- file.path(outfolder.aldex, paste('aldex_result_', tax.level, '.tsv', sep=''))
                write.table(aldex.res, outf.res, sep='\t', quote=F, row.names = F, col.names = T)

                # make boxplots
                # make plots from everything with this value or less
                plot.thresh <- 0.1
                plot.thresh.col <- 'we.eBH'
                aldex.res.plot <- aldex.res[aldex.res[, plot.thresh.col] <= plot.thresh, ]
                outf.boxplot <- file.path(outfolder.aldex, paste('aldex_significant_boxplots_', tax.level, '.pdf', sep=''))
                if (nrow(aldex.res.plot) >0 ){
                    # make clr values
                    m.prop <-  apply(rfz, 2, function(x) x/sum(x))
                    # rfz.clr <- t(codaSeq.clr(rfz, samples.by.row=FALSE))
                    clr.aldex <- aldex.clr(rfz[,c(s.pos, s.neg)],
                        conds = c(rep(group.pos, length(s.pos)), rep(group.neg, length(s.neg))),
                        mc.samples=128, denom = 'all', verbose = F)
                    # get clr from mean of the MC instances
                    rfz.clr <- sapply(clr.aldex@analysisData , function(x) rowMeans(x))

                    # save all significant in one boxplot
                    pdf(outf.boxplot, height=5, width = 6, onefile=TRUE)
                    toplot <- min(nrow(aldex.res.plot), 50)
                    for ( i in 1:toplot){
                        g <- aldex.res.plot[i, 'taxa']
                        test.df <- data.frame(group = c(rep(group.pos, length(s.pos)), rep(group.neg, length(s.neg))),
                          sample = c(s.pos, s.neg),
                          proportion = c(m.prop[g, s.pos], m.prop[g, s.neg]),
                          clr = c(rfz.clr[g, s.pos], rfz.clr[g, s.neg]))


                        g1 <- ggplot(test.df, aes(x=group, y=proportion, fill=group)) +
                        geom_boxplot() +
                        theme_bw() +
                        labs(subtitle=g) +
                        ylim(c(0,max(test.df$proportion)*1.2)) +
                        stat_compare_means(method = 'wilcox') +
                        stat_compare_means(method = 't.test', label.y.npc = 0.95) + guides(fill=FALSE)

                        g2 <- ggplot(test.df, aes(x=group, y=clr, fill=group)) +
                        geom_boxplot() +
                        theme_bw() +
                        labs(subtitle=g) +
                        ylim(c(min(test.df$clr), max(test.df$clr)*1.2)) +
                        stat_compare_means(method = 'wilcox') +
                        stat_compare_means(method = 't.test', label.y.npc = 0.95)

                        print(ggarrange(g1,g2,widths = c(1,1), common.legend = T))
                    }
                    dev.off()
                }
            }
        }
    }

    # save pca plots from each level
    pdf(file.path(outfolder.plots, 'compositional_PCA_plot.pdf'), height=6.5, width=9)
    for (tn in rev(do.tax.levels)){
        print(pca.plot.list[[tn]])
    }
    trash <- dev.off()
}



message('Done! :D')
