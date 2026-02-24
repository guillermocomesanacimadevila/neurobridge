#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript best_L_pair.R <trait1> <trait2> <locus>")
}

trait1 <- args[1]
trait2 <- args[2]
locus <- args[3]

base_dir <- "../outputs/susie"
pick_best_L <- function(trait) {
  trait_dir <- file.path(base_dir, trait, paste0("locus_", locus))
  Ldirs <- list.dirs(trait_dir, full.names=TRUE, recursive=FALSE)
  Ldirs <- Ldirs[grepl("/L[0-9]+$", Ldirs)]
  if (length(Ldirs) == 0) return(NA)
  getL <- function(d) as.numeric(sub(".*?/L([0-9]+)$", "\\1", d))
  Lvals <- sort(sapply(Ldirs, getL))
  best_L <- NA
  for (L in Lvals) {
    rds_path <- file.path(trait_dir, paste0("L", L), paste0("res_", trait, "_L", L, ".rds"))
    if (!file.exists(rds_path)) next
    fit <- readRDS(rds_path)
    if (!is.null(fit$converged) && fit$converged) {
      best_L <- L
    }
  }
  return(best_L)
}

best1 <- pick_best_L(trait1)
best2 <- pick_best_L(trait2)

cat(trait1, "best L =", best1, "\n")
cat(trait2, "best L =", best2, "\n")