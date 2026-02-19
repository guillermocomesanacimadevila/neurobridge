#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(susieR) })

args <- commandArgs(trailingOnly=TRUE)
sumstats_path <- args[1]
ld_path <- args[2]
out_dir <- args[3]

res_dir <- file.path(out_dir, "res")
dir.create(res_dir, recursive=TRUE, showWarnings=FALSE)
out_prefix <- file.path(res_dir, basename(out_dir))

d <- read.table(sumstats_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)

if (!("Z" %in% names(d))) d$Z <- d$BETA / d$SE
if (!("N" %in% names(d))) stop("missing N")
if (!all(c("SNP","A1","A2","BETA","SE","Z") %in% names(d))) stop("missing cols")

R <- as.matrix(read.table(gzfile(ld_path), header=FALSE))
snps <- readLines(sub("ld\\.ld\\.gz$", "ld.snplist", ld_path))

d <- d[d$SNP %in% snps, ]
d <- d[match(snps, d$SNP), ]
if (nrow(d) != nrow(R)) stop(paste("LD dim mismatch:", nrow(R), "vs", nrow(d)))

ref <- read.table(file.path(out_dir, "ld", "ref_alleles.tsv"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(ref) <- c("SNP","REF_A1","REF_A2")
ref <- ref[match(snps, ref$SNP), ]

okref <- !is.na(ref$SNP)
d <- d[okref, ]
ref <- ref[okref, ]
R <- R[okref, okref, drop=FALSE]

flip <- (d$A1 == ref$REF_A2 & d$A2 == ref$REF_A1)
keep <- (d$A1 == ref$REF_A1 & d$A2 == ref$REF_A2) | flip

d <- d[keep, ]
ref <- ref[keep, ]
R <- R[keep, keep, drop=FALSE]

d$Z[flip[keep]] <- -d$Z[flip[keep]]
d$BETA[flip[keep]] <- -d$BETA[flip[keep]]

n <- as.integer(median(d$N, na.rm=TRUE))
fit <- susie_rss(z=d$Z, R=R, n=n, L=10, estimate_residual_variance=FALSE, max_iter=500)

write.table(data.frame(SNP=d$SNP, PIP=fit$pip),
            paste0(out_prefix, ".pip.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

write.table(summary(fit)$cs,
            paste0(out_prefix, ".credible_sets.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

saveRDS(fit, paste0(out_prefix, ".susie.rds"))
