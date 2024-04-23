library(ggseqlogo)
library(ggplot2)
library(ggpubr)

N1_motif <- read.table(file = "N1_peak_motif_top100.flank30bp.sequence.txt")
N1_motif_seqs <- substr(N1_motif$V1, 1, 50)
p1 <- ggseqlogo(N1_motif_seqs) +
  ggtitle(NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, seq(5, 60, 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = c(0, 1, 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0.5),
        axis.text.x = element_text(size = 15, colour = "black", angle = 0),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 12, color = "black", hjust = 0.5),
        axis.title.y = element_text(size = 15, color = "black", hjust = 0.5),
        plot.title = element_text(colour = "black", face = "bold",
                                  size = 12, vjust = 0, hjust = 0.5)) +
  theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
  annotate('segment', x = 6, xend = 21, y = 1.7, yend = 1.7, size = 1, colour = "black") +
  annotate('text', x = 13.5, y = 1.9, label = '11-8', size = 6, colour = "black") +
  annotate('segment', x = 22, xend = 35, y = 1.7, yend = 1.7, size = 1, colour = "black") +
  annotate('text', x = 28.5, y = 1.9, label = '7-4', size = 6, colour = "black") +
  annotate('segment', x = 36, xend = 44, y = 1.7, yend = 1.7, size = 1, colour = "black") +
  annotate('text', x = 40, y = 1.9, label = '3-1', size = 6, colour = "black") +
  annotate('text', x = 3, y = 1.9, label = 'ZF', size = 6, colour = "black") +
  annotate('text', x = 3, y = 1.5, label = 'N1', size = 6, colour = "black") +
  annotate('rect', xmin = 5.5, xmax = 14.5, ymin = 0, ymax = 1.9, alpha = .1, fill = '#EBAC3C', size = 0) +
  annotate('rect', xmin = 38.5, xmax = 47.5, ymin = 0, ymax = 1.9, alpha = .1, fill = "#79BFC1", size = 0) +
  coord_fixed(ratio = 6)

N2_motif <- read.table(file = "N2_peak_motif_top100.flank30bp.sequence.txt")
N2_motif_seqs <- substr(N2_motif$V1, 1, 50)
p2 <- ggseqlogo(N2_motif_seqs) +
  ggtitle(NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, seq(5, 60, 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = c(0, 1, 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0.5),
        axis.text.x = element_text(size = 15, colour = "black", angle = 0),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 12, color = "black", hjust = 0.5),
        axis.title.y = element_text(size = 15, color = "black", hjust = 0.5),
        plot.title = element_text(colour = "black", face = "bold",
                                  size = 12, vjust = 0, hjust = 0.5)) +
  theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
  annotate('text', x = 3, y = 1.5, label = 'N2', size = 6, colour = "black") +
  annotate('rect', xmin = 5.5, xmax = 14.5, ymin = 0, ymax = 1.9, alpha = .1, fill = '#EBAC3C', size = 0) +
  annotate('rect', xmin = 38.5, xmax = 47.5, ymin = 0, ymax = 1.9, alpha = .1, fill = "#79BFC1", size = 0) +
  coord_fixed(ratio = 6)

N3_motif <- read.table(file = "N3_peak_motif_top100.flank30bp.sequence.txt")
N3_motif_seqs <- substr(N3_motif$V1, 1, 50)
p3 <- ggseqlogo(N3_motif_seqs) +
  ggtitle(NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, seq(5, 60, 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = c(0, 1, 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0.5),
        axis.text.x = element_text(size = 15, colour = "black", angle = 0),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 12, color = "black", hjust = 0.5),
        axis.title.y = element_text(size = 15, color = "black", hjust = 0.5),
        plot.title = element_text(colour = "black", face = "bold",
                                  size = 12, vjust = 0, hjust = 0.5)) +
  theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
  annotate('text', x = 3, y = 1.5, label = 'N3', size = 6, colour = "black") +
  annotate('rect', xmin = 5.5, xmax = 14.5, ymin = 0, ymax = 1.9, alpha = .1, fill = '#EBAC3C', size = 0) +
  annotate('rect', xmin = 38.5, xmax = 47.5, ymin = 0, ymax = 1.9, alpha = .1, fill = "#79BFC1", size = 0) +
  coord_fixed(ratio = 6)

XChIPspecific_motif <- read.table(file = "XChIPspecific_peak_motif_top100.flank30bp.sequence.txt")
XChIPspecific_motif_seqs <- substr(XChIPspecific_motif$V1, 1, 50)
p4 <- ggseqlogo(XChIPspecific_motif_seqs) +
  ggtitle(NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, seq(5, 60, 5))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = c(0, 1, 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0.5),
        axis.text.x = element_text(size = 15, colour = "black", angle = 0),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.x = element_text(size = 12, color = "black", hjust = 0.5),
        axis.title.y = element_text(size = 15, color = "black", hjust = 0.5),
        plot.title = element_text(colour = "black", face = "bold",
                                  size = 12, vjust = 0, hjust = 0.5)) +
  theme(plot.background = element_rect(fill = 'white', colour = 'white')) +
  annotate('text', x = 5, y = 1.5, label = 'XChIP', size = 6, colour = "black") +
  annotate('rect', xmin = 5.5, xmax = 14.5, ymin = 0, ymax = 1.9, alpha = .1, fill = '#EBAC3C', size = 0) +
  annotate('rect', xmin = 38.5, xmax = 47.5, ymin = 0, ymax = 1.9, alpha = .1, fill = "#79BFC1", size = 0) +
  coord_fixed(ratio = 6)

ggarrange(p1, p2, p3, p4, ncol = 1, nrow = 4)
ggsave("fig3.Motif_logo_representation_of_CTCF_motifs.pdf", width = 7, height = 7)