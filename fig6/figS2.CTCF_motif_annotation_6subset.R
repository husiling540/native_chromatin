#Filter motifs with equal XChIP and loMNase signals but different S-score
library(ggplot2)
library(ggpubr)
library(dplyr)

all_motif <- read.table(file = "figS2.CTCF_motifs_within_CTCF_binding_sites.txt", header = T)

all_motif <- all_motif[nchar(all_motif$motif_chr) <= 5 &
                         all_motif$XChIP >= 0.5 &
                         all_motif$loMNase >= 0.5,]
all_motif$XChIP_ <- round(all_motif$XChIP, 1)
all_motif$loMNase_ <- round(all_motif$loMNase, 1)


by_XM <- all_motif %>% group_by(XChIP_, loMNase_)
all_motif$group_index <- by_XM %>% group_indices()
all_motif$dup <- duplicated(all_motif$group_index)
all_motif_group <- all_motif[all_motif$group_index %in% all_motif$group_index[all_motif$dup == "TRUE"],]
all_motif_group_index <- unique(sort(all_motif_group$group_index))

num_data_frames <- 6

for (i in 1:num_data_frames) {
  name <- paste0("s", i)
  code <- paste0(name, " <- data.frame(matrix(ncol = ncol(all_motif_group), nrow = 0)); colnames(", name, ") <- colnames(all_motif_group)")
  eval(parse(text = code))
}

all_motif_group$site_type <- factor(all_motif_group$site_type, levels = c("nonnative", "native"))

set.seed(10086)

for (index in all_motif_group_index) {

  m <- all_motif_group[all_motif_group$group_index == index,]
  m_n <- nrow(m)

  if (m_n >= 6) {
    sample_count <- m_n %/% 6
    m <- m[order(m$site_type, m$salt_tolerance, decreasing = T),]
    m_s1 <- m[c(1:sample_count),]
    m_s2 <- m[c((sample_count + 1):(sample_count * 2)),]
    m_s3 <- m[c((sample_count * 2 + 1):(sample_count * 3)),]
    m_s4 <- m[c((sample_count * 3 + 1):(sample_count * 4)),]
    m_s5 <- m[c((sample_count * 4 + 1):(sample_count * 5)),]
    m_s6 <- m[c((m_n - sample_count + 1):(m_n)),]

    if (all(unique(m_s1$site_type) == "native" &&
              unique(m_s6$site_type) == "nonnative" &&
              max(m_s4$salt_tolerance) <= 2)) {
      s1 <- rbind(s1, m_s1)
      s2 <- rbind(s2, m_s2)
      s3 <- rbind(s3, m_s3)
      s4 <- rbind(s4, m_s4)
      s5 <- rbind(s5, m_s5)
      s6 <- rbind(s6, m_s6)
    }
  }
}

write.table(s1[, c(1:9)], file = "figS2.sameXChIP_differNChIP_S1.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(s2[, c(1:9)], file = "figS2.sameXChIP_differNChIP_S2.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(s3[, c(1:9)], file = "figS2.sameXChIP_differNChIP_S3.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(s4[, c(1:9)], file = "figS2.sameXChIP_differNChIP_S4.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(s5[, c(1:9)], file = "figS2.sameXChIP_differNChIP_S5.txt", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(s6[, c(1:9)], file = "figS2.sameXChIP_differNChIP_S6.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#boxplot
s1$group_type <- "S1"
s2$group_type <- "S2"
s3$group_type <- "S3"
s4$group_type <- "S4"
s5$group_type <- "S5"
s6$group_type <- "S6"

all_group <- rbind(s1, s2, s3, s4, s5, s6)

create_boxplot <- function(data, y_var, y_lim, y_breaks, y_labels, output_pdf, y_title) {
  p <- ggplot(data, aes(x = group_type, y = { { y_var } }, fill = group_type, color = group_type)) +
    geom_boxplot(width = 0.9, outlier.shape = NA, size = 0.5, alpha = 0.9) +
    scale_fill_manual(values = c("#A67BB0", "#86A6CE", "#AAD5A8", "#D3D3D3", "#AFA476", "#D3AFAC")) +
    scale_color_manual(values = c("black", "black", "black", "black", "black", "black")) +
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
    coord_cartesian(ylim = y_lim) +
    scale_y_continuous(breaks = y_breaks, labels = y_labels) +
    labs(title = NULL, x = "Subset of CTCF motif", y = y_title, size = 21)
  ggsave(output_pdf, plot = p, width = 4.3, height = 4.5)
}

create_boxplot(all_group, XChIP_, c(0, 28), seq(0, 30, 14), seq(0, 30, 14), "figS2.Boxplot_S1-S6_XChIP_signal.pdf", "X-ChIP signal")

create_boxplot(all_group, loMNase_, c(0.5, 1.3), seq(0.5, 1.5, 0.4), seq(0.5, 1.5, 0.4), "figS2.Boxplot_S1-S6_loMNase_signal.pdf", "loMNase-seq signal")

my_comparisons <- list(c("S1", "S2"), c("S2", "S3"), c("S3", "S4"), c("S4", "S5"), c("S5", "S6"))

ggplot(all_group, aes(x = group_type, y = salt_tolerance, fill = group_type, color = group_type)) +
  geom_boxplot(width = 0.9, outlier.shape = NA, size = 0.5, alpha = 0.9) +
  scale_fill_manual(values = c("#A67BB0", "#86A6CE", "#AAD5A8", "#D3D3D3", "#AFA476", "#D3AFAC")) +
  scale_color_manual(values = c("black", "black", "black", "black", "black", "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 3, colour = "black", size = 21, angle = 0),
        axis.text.y = element_text(size = 21, face = "plain", colour = "black", angle = 0),
        axis.title.x = element_text(size = 21, face = "plain", colour = "black", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 21, face = "plain", colour = "black", vjust = 1, hjust = 0.5),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "lines"),
        plot.title = element_text(colour = "black", face = "plain", size = 18, vjust = 8, hjust = 1)) +
  coord_cartesian(ylim = c(-0.5, 8.5)) +
  scale_y_continuous(breaks = seq(0, 8, 2), labels = seq(0, 8, 2)) +
  labs(title = NULL, x = "Subset of CTCF motif", y = "S-score", size = 21) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(7.5, 6, 4.5, 3, 1.5, 0.2), method = "t.test", label = "p.signif")

ggsave("figS2.Boxplot_S1-S6_S-score.pdf", width = 4.3, height = 4.5, dpi = 600)