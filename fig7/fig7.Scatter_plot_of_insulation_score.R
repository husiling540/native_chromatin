library(ggplot2)
library(ggpubr)
library(reshape2)

#Violin plot----
IS_motif <- read.table("K562_MicroC_TAD_cooltools_IS_10kb.DMSO_TPA_boundaries.CTCF_motif.bed", header = FALSE, sep = "\t")
colnames(IS_motif) <- c("DMSO_chr", "DMSO_start", "DMSO_end", "DMSO_IS", "TPA_chr", "TPA_start", "TPA_end", "TPA_IS")
IS_motif$motif <- "CTCF_motif"
t.test(IS_motif$DMSO_IS, IS_motif$TPA_IS)

IS_nomotif <- read.table("K562_MicroC_TAD_cooltools_IS_10kb.DMSO_TPA_boundaries.nonCTCF_motif.bed", header = FALSE, sep = "\t")
colnames(IS_nomotif) <- c("DMSO_chr", "DMSO_start", "DMSO_end", "DMSO_IS", "TPA_chr", "TPA_start", "TPA_end", "TPA_IS")
IS_nomotif$motif <- "no_CTCF_motif"
t.test(IS_nomotif$DMSO_IS, IS_nomotif$TPA_IS)

IS <- rbind(IS_motif, IS_nomotif)

IS_m <- melt(IS[, c("DMSO_IS", "TPA_IS", "motif")], variable.name = "type", value.name = "score")
IS_m$border <- with(IS_m, paste0(motif, "_", type))
my_comparisons <- list(c("CTCF_motif_DMSO_IS", "CTCF_motif_TPA_IS"), c("no_CTCF_motif_DMSO_IS", "no_CTCF_motif_TPA_IS"))

ggplot(IS_m, aes(x = border, y = score)) +
  geom_violin(trim = FALSE, color = NA, aes(fill = type), size = 0.2, width = 0.9, scale = "area", alpha = 0.8) +
  geom_boxplot(aes(color = motif), width = 0.3, position = position_dodge(0.9), size = 0.5, alpha = 1, outlier.shape = NA, fill = NA) +
  scale_fill_manual(values = c("#386CAF", "#DB9052")) +
  scale_color_manual(values = c("black", "black", "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, colour = "black", size = 15, angle = 0),
        axis.text.y = element_text(size = 16, face = "plain", colour = "black"),
        axis.title.x = element_text(size = 16, face = "plain", colour = "black", vjust = -1.5, hjust = 0.5),
        axis.title.y = element_text(size = 16, face = "plain", colour = "black", vjust = 2, hjust = 0.5),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), "lines"),
        plot.title = element_text(colour = "black", face = "plain", size = 15, vjust = 8, hjust = 1)) +
  labs(title = "", x = "", y = "Insulation score") +
  scale_x_discrete(breaks = c("CTCF_motif_DMSO_IS", "CTCF_motif_TPA_IS", "no_CTCF_motif_DMSO_IS", "no_CTCF_motif_TPA_IS"), labels = c("DMSO\n+", "TPA\n+", "DMSO\n-", "TPA\n-")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(-2, 2, 1), labels = seq(-2, 2, 1)) +
  coord_cartesian(ylim = c(-1.7, 1.5)) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(1, 1), method = "t.test", label = "p.signif", paired = FALSE, alternative = "two.sided")

ggsave('figS3.Violin_plot_of_insulation_score.pdf', width = 4, height = 3)

#Scatter plot----
t.test(IS$DMSO_IS, IS$TPA_IS, alternative = c("less"), paired = TRUE)
n_Y <- nrow(IS[IS$motif == "CTCF_motif",])
n_N <- nrow(IS[IS$motif == "no_CTCF_motif",])
IS$motif <- factor(IS$motif, levels = c("CTCF_motif", "no_CTCF_motif"))

p = ggplot(IS, aes(x = DMSO_IS, y = TPA_IS)) +
  geom_point(aes(color = motif), size = 0.1, alpha = 1, shape = 16) +
  scale_color_manual("CTCF motif", values = c("#DB9052", "#386CAF"), labels = c(paste0("Y (n = ", format(n_Y, big.mark = ","), ")"), paste0("N (n = ", format(n_N, big.mark = ","), ")"))) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#C2212B") +
  stat_cor(data = IS, method = "pearson", color = "#7f4f24", size = 10) +
  theme_bw() +
  theme_classic() +
  theme(panel.border = element_blank(),
        legend.position = c(0.7, 0.15),
        legend.text = element_text(colour = "black", size = 15),
        legend.title = element_text(colour = "black", size = 15),
        legend.background = element_rect(fill = "transparent"),
        axis.line.x = element_line(linetype = 1, color = "black", size = 0.8),
        axis.line.y = element_line(linetype = 1, color = "black", size = 0.8),
        axis.ticks.x = element_line(color = "black", size = 0.8),
        axis.ticks.y = element_line(color = "black", size = 0.8),
        axis.title.x = element_text(size = 20, vjust = 0, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        plot.title = element_text(colour = "black", size = 20, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.3, 0.3, 0.3), "lines")) +
  labs(x = "Insulation score (DMSO)", y = "Insulation score (TPA)", size = 20) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(breaks = seq(-2, 2, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-2, 2, 1), expand = c(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(-2, 1.5), ylim = c(-2, 1.5))

IS1 <- IS[IS$motif == "CTCF_motif",]
p + geom_point(data = IS1, aes(x = DMSO_IS, y = TPA_IS, color = motif), size = 0.1, alpha = 1, shape = 16)

ggsave('fig7.Scatter_plot_of_insulation_score.pdf.pdf', width = 4.5, height = 4.5, dpi = 600)