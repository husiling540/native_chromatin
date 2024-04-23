library(ggplot2)
library(ggpubr)
require(grid)

input_name <- read.table(file = paste0("bw_K562_public_data_to_plot_profile.sim.txt"))
samples <- c(input_name$V1)

for (sample in samples) {
  cat(sample, "============", "\n")
  enrichment <- read.table(file = paste0("computeMatrix_outFileMatrix/K562_Xu_TF_K562_ChIP-seq_", sample, ".hg38.matrix"), sep = "\t", header = F, skip = 3)
  sample_n <- nrow(enrichment) / 4
  N1 <- data.frame(experiment = "N1", pos = c(1:200), sig = colMeans(enrichment[c(1:sample_n),]))
  N2 <- data.frame(experiment = "N2", pos = c(1:200), sig = colMeans(enrichment[c((sample_n + 1):(sample_n * 2)),]))
  N3 <- data.frame(experiment = "N3", pos = c(1:200), sig = colMeans(enrichment[c((sample_n * 2 + 1):(sample_n * 3)),]))
  XChIP <- data.frame(experiment = "XChIP", pos = c(1:200), sig = colMeans(enrichment[c((sample_n * 3 + 1):(sample_n * 4)),]))

  enrichment_merge <- rbind(N1, N2, N3, XChIP)
  enrichment_merge$experiment <- factor(enrichment_merge$experiment, levels = c("N1", "N2", "N3", "XChIP"))
  max_signal <- max(enrichment_merge$sig)
  min_signal <- min(enrichment_merge$sig)
  max_s <- max_signal + (max_signal - min_signal) / 7
  min_s <- min_signal - (max_signal - min_signal) / 7
  break_y <- seq(round(min_signal, 2), (round(min_signal, 2) + round((max_signal - min_signal) / 2, 2) * 2), round((max_signal - min_signal) / 2, 2))
  break_n <- length(break_y)
  p = ggplot(enrichment_merge, aes(x = pos, y = sig)) +
    ggtitle(NULL) +
    annotate("text", x = 100, y = max_signal, label = sample, vjust = -0.5, size = 4.5) +
    geom_line(lwd = 1, aes(color = experiment)) +
    labs(x = NULL, y = NULL) +
    guides(color = guide_legend(ncol = 1)) +
    theme_bw() +
    theme(plot.background = element_rect(fill = NA, color = NA),
          plot.margin = unit(c(.0, .0, 0.1, .1), 'inches'),
          plot.title = element_text(size = 14, hjust = .5),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black"),
          axis.ticks.y = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          axis.title.x = element_text(size = 14, colour = 'black'),
          axis.title.y = element_text(size = 14, colour = 'black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, colour = 'black'),
          legend.text = element_text(size = 14),
          legend.position = c(0.85, 0.9), legend.background = element_rect(fill = "transparent")) +
    scale_colour_manual(name = NULL, values = c("#7FC87F", "#9f79d3", "#386CAF", "#fc9533"), labels = c("N1", "N2", "N3", "X-ChIP")) +
    scale_x_continuous(limits = c(1, 200), breaks = c(1, 100, 200), labels = c(" -1kb", "center", "1kb ")) +
    scale_y_continuous(limits = c(round(min_s, 2), round(max_s, 2)),
                       breaks = break_y,
                       labels = c(round(min_signal, 2), rep("", (break_n - 2)), (round(min_signal, 2) + round((max_signal - min_signal) / 2, 2) * 2)))
  assign(paste0("p_", sample), p)

}

for (sample in c("GRO-seq")) {
  cat(sample, "============", "\n")
  enrichment <- read.table(file = paste0("computeMatrix_outFileMatrix/K562_transcription_", sample, ".matrix"), sep = "\t", header = F, skip = 3)
  sample_n <- nrow(enrichment) / 4
  N1 <- data.frame(experiment = "N1", pos = c(1:200), sig = colMeans(enrichment[c(1:sample_n),]))
  N2 <- data.frame(experiment = "N2", pos = c(1:200), sig = colMeans(enrichment[c((sample_n + 1):(sample_n * 2)),]))
  N3 <- data.frame(experiment = "N3", pos = c(1:200), sig = colMeans(enrichment[c((sample_n * 2 + 1):(sample_n * 3)),]))
  XChIP <- data.frame(experiment = "XChIP", pos = c(1:200), sig = colMeans(enrichment[c((sample_n * 3 + 1):(sample_n * 4)),]))

  enrichment_merge <- rbind(N1, N2, N3, XChIP)
  enrichment_merge$experiment <- factor(enrichment_merge$experiment, levels = c("N1", "N2", "N3", "XChIP"))
  max_signal <- max(enrichment_merge$sig)
  min_signal <- min(enrichment_merge$sig)
  max_s <- max_signal + (max_signal - min_signal) / 7
  min_s <- min_signal - (max_signal - min_signal) / 7
  break_y <- seq(round(min_signal, 2), (round(min_signal, 2) + round((max_signal - min_signal) / 2, 2) * 2), round((max_signal - min_signal) / 2, 2))
  break_n <- length(break_y)
  p = ggplot(enrichment_merge, aes(x = pos, y = sig)) +
    ggtitle(NULL) +
    annotate("text", x = 100, y = max_signal, label = sample, vjust = -0.5, size = 4.5) +
    geom_line(lwd = 1, aes(color = experiment)) +
    labs(x = NULL, y = NULL) +
    guides(color = guide_legend(ncol = 1)) +
    theme_bw() +
    theme(plot.background = element_rect(fill = NA, color = NA),
          plot.margin = unit(c(.0, .0, 0.1, .1), 'inches'),
          plot.title = element_text(size = 14, hjust = .5),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black"),
          axis.ticks.y = element_line(color = "black"),
          axis.ticks.length.x = unit(-0.12, "cm"),
          axis.line = element_blank(),
          axis.title.x = element_text(size = 14, colour = 'black'),
          axis.title.y = element_text(size = 14, colour = 'black'),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, colour = 'black'),
          legend.text = element_text(size = 14),
          legend.position = c(0.85, 0.9), legend.background = element_rect(fill = "transparent")) +
    scale_colour_manual(name = NULL, values = c("#7FC87F", "#9f79d3", "#386CAF", "#fc9533"), labels = c("N1", "N2", "N3", "X-ChIP")) +
    scale_x_continuous(limits = c(1, 200), breaks = c(1, 100, 200), labels = c(" -1kb", "center", "1kb ")) +
    scale_y_continuous(limits = c(round(min_s, 2), 0.052),
                       breaks = break_y,
                       labels = c(round(min_signal, 2), rep("", (break_n - 2)), (round(min_signal, 2) + round((max_signal - min_signal) / 2, 2) * 2)))
  assign(paste0("p_", sample), p)
}

p_POLR2A = p_POLR2A + theme(axis.ticks.x = element_line(color = "black"), axis.ticks.length.x = unit(-0.12, "cm"))
p_H3K36me3 = p_H3K36me3 + theme(axis.ticks.x = element_line(color = "black"), axis.ticks.length.x = unit(-0.12, "cm"))
p_SUMO23 = p_SUMO23 + theme(axis.ticks.x = element_line(color = "black"), axis.ticks.length.x = unit(-0.12, "cm"))

figure <- ggarrange(p_POLR2A, p_POLR2AphosphoS5, p_POLR2AphosphoS2, ncol = 3, nrow = 1, widths = c(1:1:1), align = "hv", common.legend = TRUE, legend = "right") +
  theme(plot.background = element_rect(fill = 'white', color = 'white'))
annotate_figure(figure, left = textGrob("ChIP signal", rot = 90, vjust = 1, gp = gpar(cex = 1.2))) +
  theme(plot.margin = unit(c(.1, .1, .1, .1), 'inches'))
ggsave("fig4.Average_profile_Pol2_N1N2N3XChIP.pdf", width = 8.27, height = 2.2, dpi = 600)

figure <- ggarrange(p_H2AFZ, p_H3K4me1, p_H3K4me2, p_H3K4me3, p_H3K9ac, p_H3K27ac, p_H3K36me3, p_H3K79me2, p_H3K9me3_PE, p_H3K27me3_PE,
                    ncol = 4, nrow = 3, widths = c(1:1:1:1:1), align = "hv", common.legend = TRUE, legend = "right") +
  theme(plot.background = element_rect(fill = 'white', color = 'white'))
annotate_figure(figure, left = textGrob("ChIP signal", rot = 90, vjust = 1, gp = gpar(cex = 1.2))) +
  theme(plot.margin = unit(c(.1, .1, .1, .1), 'inches'))
ggsave(paste0("fig4.Average_profile_histone_N1N2N3XChIP.pdf"), width = 10.5, height = 6.2, dpi = 600)

figure <- ggarrange(`p_GRO-seq`, ncol = 1, nrow = 1, widths = c(1:1:1:1), align = "hv", common.legend = TRUE, legend = "right") +
  theme(plot.background = element_rect(fill = 'white', color = 'white'))
annotate_figure(figure, left = textGrob("GRO-seq signal", rot = 90, vjust = 1, gp = gpar(cex = 1.2))) +
  theme(plot.margin = unit(c(.1, .1, .1, .1), 'inches'))
ggsave("figEV4.Average_profile_GRO-seq_N1N2N3XChIP.pdf", width = 3.815, height = 2.2, dpi = 600)

figure <- ggarrange(p_SUMO23, ncol = 1, nrow = 1, widths = c(1:1:1:1), align = "hv", common.legend = TRUE, legend = "right") +
  theme(plot.background = element_rect(fill = 'white', color = 'white'))
annotate_figure(figure, left = textGrob("SUMO2/3 signal", rot = 90, vjust = 1, gp = gpar(cex = 1.2))) +
  theme(plot.margin = unit(c(.1, .1, .1, .1), 'inches'))
ggsave("figEV4.Average_profile_SUMO_N1N2N3XChIP.pdf", width = 3.815, height = 2.2, dpi = 600)