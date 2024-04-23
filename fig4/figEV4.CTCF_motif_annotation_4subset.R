library(ggplot2)
library(ggpubr)
library(dplyr)

base <- read.table(file = "figEV4.CTCF_motif_annotation_4subset.txt", header = TRUE)

by_XM <- base %>% group_by(xchip_signal)
base$group_index <- by_XM %>% group_indices()
base$dup <- duplicated(base$group_index)
base_group <- base[base$group_index %in% base$group_index[base$dup == TRUE],]
base_group_index <- unique(sort(base_group$group_index))
subset <- data.frame(matrix(ncol = ncol(base_group), nrow = 0))
colnames(subset) <- colnames(base_group)

set.seed(10086)
for (index in base_group_index) {
  m <- base_group[base_group$group_index == index,]
  if (length(unique(m$peak_group)) >= 4) {
    length_m <- as.data.frame(table(m$peak_group))
    min_l <- min(length_m$Freq)
    m_a1 <- m[m$peak_group == "N1",]
    m_a2 <- m[m$peak_group == "N2",]
    m_a3 <- m[m$peak_group == "N3",]
    m_a4 <- m[m$peak_group == "XChIP",]

    m_a1 <- m_a1[sample(nrow(m_a1), min_l),]
    m_a2 <- m_a2[sample(nrow(m_a2), min_l),]
    m_a3 <- m_a3[sample(nrow(m_a3), min_l),]
    m_a4 <- m_a4[sample(nrow(m_a4), min_l),]

    subset <- rbind(subset, m_a1, m_a2, m_a3, m_a4)
  }
}

subset <- subset[subset$xchip_signal >= 0.5,]
subset_n <- nrow(subset) / 4

N1 <- subset[subset$peak_group == "N1",]
N2 <- subset[subset$peak_group == "N2",]
N3 <- subset[subset$peak_group == "N3",]
XChIP <- subset[subset$peak_group == "XChIP",]

cat(mean(N1$xchip_signal), "\t", mean(N2$xchip_signal), "\t", mean(N3$xchip_signal), "\t", mean(XChIP$xchip_signal), "\n")

write.table(N1[, c(1:8)], file = "figEV4.sameXChIP_differNChIP_N1.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(N2[, c(1:8)], file = "figEV4.sameXChIP_differNChIP_N2.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(N3[, c(1:8)], file = "figEV4.sameXChIP_differNChIP_N3.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(XChIP[, c(1:8)], file = "figEV4.sameXChIP_differNChIP_XChIP.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

all_group <- rbind(N1, N2, N3, XChIP)
my_comparisons <- list(c("N1", "N2"), c("N2", "N3"), c("N3", "XChIP"))

ggplot(all_group, aes(x = peak_group, y = salt_tolerance, fill = peak_group, color = peak_group)) +
  geom_boxplot(width = 0.9, outlier.shape = NA, size = 0.5, alpha = 0.9) +
  scale_fill_manual(values = c("#7FC87F", "#9f79d3", "#386CAF", "#fc9533")) +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 3, colour = "black", size = 21, angle = 0),
        axis.text.y = element_text(size = 21, face = "plain", colour = "black", angle = 0),
        axis.title.x = element_text(size = 21, face = "plain", colour = "black", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 21, face = "plain", colour = "black", vjust = 1, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "lines"),
        plot.title = element_text(colour = "black", face = "plain", size = 18, vjust = 8, hjust = 1)) +
  coord_cartesian(ylim = c(-1, 8.5)) +
  scale_x_discrete(breaks = c("N1", "N2", "N3", "XChIP"), labels = c("N1", "N2", "N3", "X-ChIP")) +
  scale_y_continuous(breaks = seq(0, 8, 2), labels = seq(0, 8, 2)) +
  labs(title = NULL, x = "Subset of CTCF motif", y = "S-score", size = 21) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(7.5, 6, 4.5, 3, 1.5, 0.2), method = "wilcox.test", label = "p.signif")
ggsave("figEV4.Boxplot_S-score.pdf", width = 4, height = 4, dpi = 600)

ggplot(all_group, aes(x = peak_group, y = xchip_signal, fill = peak_group, color = peak_group)) +
  geom_boxplot(width = 0.9, outlier.shape = NA, size = 0.5, alpha = 0.9) +
  scale_fill_manual(values = c("#7FC87F", "#9f79d3", "#386CAF", "#fc9533")) +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 3, colour = "black", size = 21, angle = 0),
        axis.text.y = element_text(size = 21, face = "plain", colour = "black", angle = 0),
        axis.title.x = element_text(size = 21, face = "plain", colour = "black", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 21, face = "plain", colour = "black", vjust = 1, hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "lines"),
        plot.title = element_text(colour = "black", face = "plain", size = 18, vjust = 8, hjust = 1)) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_x_discrete(breaks = c("N1", "N2", "N3", "XChIP"), labels = c("N1", "N2", "N3", "X-ChIP")) +
  scale_y_continuous(breaks = seq(0, 30, 15), labels = seq(0, 30, 15)) +
  labs(title = NULL, x = "Subset of CTCF motif", y = "X-ChIP signal", size = 21)
ggsave("figEV4.Boxplot_XChIP_signal.pdf", width = 4, height = 4, dpi = 600)