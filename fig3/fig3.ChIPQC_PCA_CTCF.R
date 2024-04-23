library(ChIPQC)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BiocParallel)
register(SerialParam())

samples <- read.csv("fig3.ChIPQC_PCA_CTCF.csv")

exampleExp <- ChIPQC(samples, annotation = "hg38", consensus = TRUE,
                     chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                     "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
                                     "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))

result = plotPrincomp(exampleExp, attributes = "Condition", label = "Replicate", dotSize = 1, labelSize = 0, components = 1:3,
                      ylim = c(-10, 10), xlim = c(-10, 10), alpha = 0.9, vColors = c("#342D80", "#C4671B", "#769A52", "#D45575"),
                      labelCols = c("#342D80", "#C4671B", "#769A52", "#D45575"), box.ratio = 1)
result$main <- NULL
result$aspect.ratio <- 1
result$ylab <- sub("Principal Component #", "PC", result$ylab)
result$xlab <- sub("Principal Component #", "PC", result$xlab)
plotData <- environment(result[["panel"]])[["plotData"]]

plotData$sample <- rownames(plotData)
plotData$condition <- environment(result[["panel"]])[["classvec"]]
plotData$condition <- factor(plotData$condition, levels = c("75mM", "150mM", "225mM", "X-ChIP"))

ggplot(plotData) +
  geom_point(aes(x = PC1, y = PC2, color = condition), alpha = 0.8, size = 4, shape = 16) +
  labs(title = NULL, x = result$xlab, y = result$ylab) +
  theme_bw() +
  theme(legend.box = "horizontal",
        legend.key.size = unit(18, "pt"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"),
        legend.text = element_text(colour = "black", size = 15),
        panel.border = element_rect(fill = NA, size = 1, color = 'black'),
        axis.title.x = element_text(colour = "black", size = 15, vjust = 0),
        axis.title.y = element_text(colour = "black", size = 15, vjust = 2),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 13, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 13, colour = "black"),
        plot.title = element_text(colour = "black", size = 15, vjust = 0.5, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        text = element_text(size = 15, family = "sans", colour = "black")) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_color_manual(NULL,
                     labels = c("N-ChIP 75 mM", "N-ChIP 150 mM", "N-ChIP 225 mM", "X-ChIP"),
                     values = c("#386CAF", "#9f79d3", "#7FC87F", "#fc9533")) +
  scale_x_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 5)) +
  scale_y_continuous(limits = c(-11, 11), breaks = seq(-10, 10, 5)) +
  coord_fixed(ratio = 1)

ggsave(filename = "fig3.ChIPQC_PCA_CTCF.pdf", width = 4.7, height = 3)