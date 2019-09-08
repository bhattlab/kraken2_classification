# functions from the CoDaSeq package (https://github.com/ggloor/CoDaSeq)
# by Greg Gloor
# included this way because the package wouldn't install in my singularity container!

codaSeq.clr <- function (x, IQLR = FALSE, aitch = FALSE, samples.by.row = TRUE)
{
    if (min(x) <= 0)
        stop("only positive real values permitted")
    if (samples.by.row == TRUE)
        margin = 1
    if (samples.by.row == FALSE)
        margin = 2
    if (aitch == FALSE) {
        if (IQLR == FALSE) {
            return(t(apply(x, margin, function(x) {
                log(x) - mean(log(x))
            })))
        }
        else if (IQLR == TRUE) {
            reads.clr <- t(apply(x, margin, function(x) {
                log(x) - mean(log(x))
            }))
            reads.var <- apply(reads.clr, 2, var)
            reads.qtl <- quantile(unlist(reads.var))
            mid.set <- which(reads.var < (reads.qtl[4]) & reads.var >
                                 (reads.qtl[2]))
            return(t(apply(x, margin, function(x) log(x) - mean(log(x[mid.set])))))
        }
    }
    if (aitch == TRUE) {
        aitchison.mean <- function(n, log = TRUE) {
            n <- round(as.vector(n, mode = "numeric"))
            if (any(n < 0))
                stop("counts cannot be negative")
            a <- n + 0.5
            sa <- sum(a)
            log.p <- digamma(a) - digamma(sa)
            log.p <- log.p - mean(log.p)
            if (log)
                return(log.p)
            p <- exp(log.p - max(log.p))
            p <- p/sum(p)
            return(p)
        }
        if (samples.by.row == TRUE)
            x <- t(x)
        return(aitchison.mean(x, log = TRUE))
    }
}


codaSeq.filter <- function (x, min.reads = 5000, min.prop = 0.001, max.prop = 1,
          min.occurrence = 0, var.filt = FALSE, min.count = 0, samples.by.row = TRUE)
{
    if (samples.by.row == FALSE)
        data <- x
    if (samples.by.row == TRUE)
        data <- t(x)
    if (length(rownames(data)) == 0)
        stop("rownames cannot be empty")
    if (length(colnames(data)) == 0)
        stop("colnames cannot be empty")
    if (any(round(data) != data))
        stop("not all values are integers")
    if (any(data < 0))
        stop("one or more values are negative")
    if (var.filt == FALSE & min.count == 0) {
        data.0 <- data[, which(apply(data, 2, sum) > min.reads)]
        d.frac <- apply(data.0, 2, function(x) {
            x/sum(x)
        })
        data.1 <- data.0[(which((apply(d.frac, 1, max) > min.prop) &
                                    (apply(d.frac, 1, max) < max.prop))), ]
        rm(d.frac)
        data.2 <- data.frame(data.1[which(apply(data.1, 1, function(x) {
            length(which(x != 0))/length(x)
        }) > min.occurrence), ], stringsAsFactors = FALSE, check.names = FALSE)
    }
    else if (var.filt == TRUE) {
        if (min.count > 0)
            warning("filtering on variance will not filter on read count")
        warning("filtering only on sample read count and feature variance")
        data.0 <- data[, which(apply(data, 2, sum) > min.reads)]
        data.1 <- data.0[which(apply(data.0, 1, max) >= min.count),
                         ]
        d.n0 <- cmultRepl(t(data.1), method = "CZM", label = 0)
        d.clr <- codaSeq.clr(d.n0, samples.by.row = TRUE)
        var.clr <- apply(d.clr, 2, var)
        names.hvar <- names(var.clr)[which(var.clr > median(var.clr))]
        data.2 <- data.frame(data.1[which(apply(data.1[names.hvar,
                                                       ], 1, max) > 0), ], stringsAsFactors = FALSE, check.names = FALSE)
    }
    else if (min.count > 0) {
        warning("filtering on sample read count and minimum feature read count only")
        data.0 <- data[, which(apply(data, 2, sum) > min.reads)]
        data.2 <- data.frame(data.0[which(apply(data.0, 1, max) >=
                                              min.count), ], stringsAsFactors = FALSE, check.names = FALSE)
    }
    return(data.2)
}
