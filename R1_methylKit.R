rm(list = ls())
library(methylKit)


output_dir <- "D:/test2/test/plots2"
dir.create(output_dir, showWarnings = FALSE)


setwd("D:/test2/test/data")

## 需要解压
file.list <- list("D:/test2/test/data/data/XM-ZX-HiME-4_L3_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov",
                  "D:/test2/test/data/data/XM-ZX-HiME-5_L2_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov",
                  "D:/test2/test/data/data/XM-ZX-HiME-6_L2_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov")

sample.ids <- list("HiME-4", "HiME-5", "HiME-6")


myobj <- methRead(file.list,
                  sample.id = sample.ids,
                  assembly = "hg38",
                  treatment = c(0, 0, 0),
                  context = "CpG",
                  pipeline = "bismarkCoverage",
                  mincov = 10)


pdf(file.path(output_dir, "methylation_stats_all_samples.pdf"), width = 14, height = 5) # 调整宽度
par(mfrow = c(1, 3))  
for (i in seq_along(myobj)) {
  getMethylationStats(myobj[[i]], plot = TRUE, both.strands = FALSE)

}
dev.off()


pdf(file.path(output_dir, "coverage_stats_all_samples.pdf"), width = 14, height = 5) 
par(mfrow = c(1, 3))  
for (i in seq_along(myobj)) {
  getCoverageStats(myobj[[i]], plot = TRUE, both.strands = FALSE)

}
dev.off()


meth = unite(myobj, destrand = FALSE)

pdf(file.path(output_dir, "correlation_plot.pdf"))
getCorrelation(meth, plot = TRUE)
dev.off()

pdf(file.path(output_dir, "cluster_plot.pdf"))
clusterSamples(meth, dist = "correlation", method = "ward", plot = TRUE)
dev.off()

pdf(file.path(output_dir, "pca_plot.pdf"))
PCASamples(meth)
dev.off()
