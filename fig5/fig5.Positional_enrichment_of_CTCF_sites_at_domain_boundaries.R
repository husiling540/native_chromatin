library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

for (boundary_type in 1:10) {
  native_bins <- read.table(file = paste0("native_CTCF_sites_at_domain_boundaries_", boundary_type, "cell.bed"), sep = "\t", header = F)
  colnames(native_bins) <- c("chr", "start", "end", "boundary_bin", "overlap")
  native_bins_sep <- separate(native_bins, boundary_bin, c("boundary", "bin"), "_")
  native_bins_sep_group <- native_bins_sep %>%
    group_by(bin) %>%
    summarise(enrichment = sum(overlap) / (nrow(native_bins) / 100) * 100)
  native_bins_sep_group$peak_type <- "native"

  non_native_bins <- read.table(file = paste0("nonnative_CTCF_sites_at_domain_boundaries_", boundary_type, "cell.bed"), sep = "\t", header = F)
  colnames(non_native_bins) <- c("chr", "start", "end", "boundary_bin", "overlap")
  non_native_bins_sep <- separate(non_native_bins, boundary_bin, c("boundary", "bin_non"), "_")
  non_native_bins_sep_group <- non_native_bins_sep %>%
    group_by(bin_non) %>%
    summarise(enrichment_non = sum(overlap) / (nrow(non_native_bins) / 100) * 100)
  non_native_bins_sep_group$peak_type <- "non_native"

  profile <- cbind(native_bins_sep_group[, c("bin", "enrichment")], non_native_bins_sep_group[, c("bin_non", "enrichment_non")])
  profile$bin <- as.integer(profile$bin)
  profile$bin_non <- as.integer(profile$bin_non)
  max_signal <- max(c(profile$enrichment, profile$enrichment_non))
  min_signal <- min(c(profile$enrichment, profile$enrichment_non))

  p <- ggplot(profile) +
    ggtitle(paste0(boundary_type, " cell types")) +
    geom_line(aes(x = bin, y = enrichment, color = "native"), lwd = 0.8, alpha = 1) +
    geom_line(aes(x = bin_non, y = enrichment_non, color = "non-native"), lwd = 0.8, alpha = 1) +
    theme_bw() +
    labs(x = '', y = '') +
    scale_x_continuous(limits = c(1, 100), breaks = c(1, 50.5, 100), labels = c("-500", "0", "500")) +
    scale_y_continuous(limits = c(min_signal, max((round(min_signal, 1) + round((max_signal - min_signal) / 2, 1) * 2), max_signal)), breaks = seq(round(min_signal, 1), (round(min_signal, 1) + round((max_signal - min_signal) / 2, 1) * 2), round((max_signal - min_signal) / 2, 1))) +
    guides(color = guide_legend(ncol = 2)) +
    theme(plot.margin = unit(c(0, 0.05, 0.05, 0), 'inches'),
          plot.title = element_text(size = 15, hjust = 0.5),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, size = 1, color = 'black'),
          axis.title.x = element_text(size = 15, vjust = -1),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 14, colour = 'black'),
          axis.text.y = element_text(size = 14, colour = 'black'),
          legend.position = c(0.25, 0.8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.text = element_text(size = 15, colour = 'black')) +
    scale_color_manual("", values = c(native = "#941416", `non-native` = "#386CAF"))

  assign(paste0("p_", boundary_type), p)
}

figure <- ggarrange(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10,
                    ncol = 5, nrow = 2, widths = c(1, 1, 1, 1, 1), align = "hv", common.legend = TRUE) +
  theme(plot.background = element_rect(fill = 'white', color = 'white'))

annotate_figure(figure, left = textGrob("Normalized overlap (%)", rot = 90, vjust = 1, gp = gpar(cex = 1.5)),
                bottom = textGrob("Distance from arrowhead boundaries (kb)", rot = 0, vjust = -0.5, gp = gpar(cex = 1.5))) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inches'))

ggsave("fig5.Positional_enrichment_of_CTCF_sites_at_domain_boundaries.pdf", width = 16, height = 6)