library("RColorBrewer")
library("ggplot2")
library("plyr")

setwd("/PATH/05_aDNA/00_eager_TE_mappedtoprobefasta/results/mapping/bwa/00_coverage")

cov <- read.delim("all.coverage_across_reference", header=TRUE, row.names=1, sep="\t")
cols <- c("DSZ007_TE_PE", "DVT014_TE_PE", "DVT017_TE_PE", "DVT022_TE_PE", "KCZ003_TE_PE", "NCP001_TE_PE", "TDN002_TE_PE", "TPC001_TE_PE", "TPC019_TE_PE", "WMP006_TE_PE")
cov[cols] <- log(cov[cols])
tcov <- t(cov)

tiff("cov.tiff")
heatmap(tcov, Colv = NA, Rowv = NA, scale = "none")
dev.off()

cov <- read.delim("outlier_genic.coverage_across_reference", header=TRUE, row.names=1, sep="\t")
cols <- c("DSZ007_TE_PE", "DVT014_TE_PE", "DVT017_TE_PE", "DVT022_TE_PE", "KCZ003_TE_PE", "NCP001_TE_PE", "TDN002_TE_PE", "TPC001_TE_PE", "TPC019_TE_PE", "WMP006_TE_PE")
cov[cols] <- log(cov[cols])
tcov <- t(cov)

tiff("outlier.cov.tiff")
heatmap(tcov, Colv = NA, Rowv = NA, scale = "none")
dev.off()

