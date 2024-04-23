library(DiffBind)
library(ggplot2)

samples <- read.csv("fig6.Volcano_plot_of_CTCF_N-ChIP_input.csv")

dbObj <- dba(sampleSheet = samples)

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)

dbObj <- dba.contrast(dbObj, reorderMeta = list(Condition = "DMSO"), categories = DBA_CONDITION, minMembers = 2)

dbObj <- dba.analyze(dbObj, method = DBA_ALL_METHODS)

res_deseq <- dba.report(dbObj, method = DBA_DESEQ2, contrast = 1, th = 1)

result <- as.data.frame(res_deseq)

result$threshold <- factor(ifelse(result$FDR < 0.05 & abs(result$Fold) >= 1, ifelse(result$Fold <= -1, 'TPA_loss', 'TPA_gain'), 'stable'), levels = c('TPA_loss', 'stable', 'TPA_gain'))

ggplot(result, aes(x = Fold, y = -log(FDR, 10), color = threshold)) +
  geom_point(size = 0.1, alpha = 1, shape = 16) +
  scale_color_manual(values = c("#000000", "#AEACAB", "#903529")) +
  theme_bw() +
  theme_classic() +
  annotate("text", x = 5.5, y = 13.5, hjust = 1, label = paste0(format(nrow(result[result$threshold == "TPA_loss",]), big.mark = ","), " sites"), size = 6.5, color = "black") +
  annotate("text", x = 5.5, y = 12.3, hjust = 1, label = paste0(format(nrow(result[result$threshold == "stable",]), big.mark = ","), " sites"), size = 6.5, color = "#AEACAB") +
  annotate("text", x = 5.5, y = 11.1, hjust = 1, label = paste0(format(nrow(result[result$threshold == "TPA_gain",]), big.mark = ","), " sites"), size = 6.5, color = "#903529") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.line.x = element_line(linetype = 1, color = "black", size = 0),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0),
        axis.ticks.x = element_line(color = "black", size = 0.5),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.title.x = element_text(size = 20, vjust = 0),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, colour = 'black'),
        axis.text.y = element_text(size = 20, colour = 'black'),
        plot.title = element_text(colour = "black", size = 20, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.3, 0.3, 0.3), "lines"),
        legend.position = "none") +
  labs(title = "CTCF N-ChIP", x = expression("log"[2] * " (TPA / DMSO)"), y = expression("-log"[10] * " (FDR)"), size = 20) +
  scale_x_continuous(breaks = seq(-6, 6, 3), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 20, 5), expand = c(0, 0)) +
  coord_fixed(ratio = 14 / 15 * 1.11, xlim = c(-7, 7), ylim = c(0, 15))

ggsave("fig6.Volcano_plot_of_CTCF_N-ChIP.pdf", width = 5, height = 5)