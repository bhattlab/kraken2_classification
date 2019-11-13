indir <- '~/scg/crassphage_public_data/shao_2019/kraken2_classification_viral_guerin/classification/'
f.names <- sapply(list.files(indir, pattern='.krak.report$', full.names =F), function(x) strsplit(x, '\\.')[[1]][1])
fs.full <- list.files(indir, pattern='.krak.report$', full.names =T)

keep.taxids <- c(19780071, 19780072, 19780073, 19780074, 197800711, 197800712, 197800713, 197800714, 197800721, 197800731, 197800732, 197800741, 197800742, 197800743)
keep.names <- c("Guerin_Alphacrassvirinae","Guerin_Betacrassvirinae","Guerin_Gammacrassvirinae","Guerin_Deltacrassvirinae","Guerin_crAss_alpha_cluster_01","Guerin_crAss_alpha_cluster_03","Guerin_crAss_alpha_cluster_04","Guerin_crAss_alpha_cluster_09","Guerin_crAss_beta_cluster_06","Guerin_crAss_gamma_cluster_02","Guerin_crAss_gamma_cluster_05","Guerin_crAss_delta_cluster_07","Guerin_crAss_delta_cluster_08","Guerin_crAss_delta_cluster_10")

class.mat <- matrix(0, nrow=length(keep.taxids), ncol=length(f.names), dimnames=list(keep.names, f.names))

for (i in 1:length(f.names)){
    print(i)
    this.name <- f.names[i]
    this.file <- fs.full[i]
    df <- read.table(this.file, sep='\t', header=F, quote='')
    df$V6 <- trimws(df$V6)
    df <- df[df$V5 %in% keep.taxids,]
    for (j in 1:nrow(df)){
        class.mat[df[j, 6], this.name] <- df[j,2]
    }
}

write.table(class.mat, '~/scg/crassphage_public_data/shao_2019/crasslike_classified_reads.tsv', sep='\t', quote=F, row.names=T, col.names=T)
write.table(t(class.mat), '~/scg/crassphage_public_data/shao_2019/crasslike_classified_reads_transpose.tsv', sep='\t', quote=F, row.names=T, col.names=T)
