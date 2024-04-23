library(ggplot2)
library(dplyr)
library(tidyr)

native_bins <- read.table(file = "K562_loops.120bins.native_sites.bed", sep = "\t", header = FALSE)
colnames(native_bins) <- c("chr", "start", "end", "loop_bin", "bool")
native_bins_sep <- separate(native_bins, loop_bin, into = c("loop", "bin"), sep = "_")
native_bins_sep_group <- native_bins_sep %>%
  group_by(bin) %>%
  summarise(enrichment = sum(bool) / (nrow(native_bins) / 120) * 100) %>%
  mutate(peak_type = "native")

non_native_bins <- read.table(file = "K562_loops.120bins.nonnative_sites.bed", sep = "\t", header = FALSE)
colnames(non_native_bins) <- c("chr", "start", "end", "loop_bin", "bool")
non_native_bins_sep <- separate(non_native_bins, loop_bin, into = c("loop", "bin_non"), sep = "_")
non_native_bins_sep_group <- non_native_bins_sep %>%
  group_by(bin_non) %>%
  summarise(enrichment_non = sum(bool) / (nrow(non_native_bins) / 120) * 100) %>%
  mutate(peak_type = "non_native")

profile <- cbind(native_bins_sep_group[, c("bin", "enrichment")], non_native_bins_sep_group[, c("bin_non", "enrichment_non")])
profile$bin <- as.integer(profile$bin)
profile$bin_non <- as.integer(profile$bin_non)

ggplot(profile) +
  geom_line(aes(x = bin, y = enrichment, color = "native"), lwd = 1, alpha = 1) +
  geom_line(aes(x = bin_non, y = enrichment_non, color = "non-native"), lwd = 1, alpha = 1) +
  geom_vline(xintercept = c(11, 110), color = '#708090', linetype = 'dashed', lwd = 0.5) +
  theme_bw() +
  labs(x = NULL, y = 'Normalized overlap (%)') +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 120), breaks = c(seq(1, 11), 60, seq(110, 120)), labels = c(rep("", 10), "left", "loop domain", "right", rep("", 10))) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.4)) +
  guides(color = guide_legend(ncol = 1)) +
  theme(plot.margin = unit(c(.1, .1, .1, .1), 'inches'),
        plot.title = element_text(size = 15, hjust = .5),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.2),
        axis.line = element_line(colour = "black", size = 0.2),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 15, vjust = -1),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14, colour = 'black'),
        axis.text.y = element_text(size = 14, colour = 'black'),
        legend.position = c(0.5, 0.6),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = 'black'),
        axis.line.y = element_blank()) +
  scale_color_manual("", values = c(native = "#C4671B", `non-native` = "#342D80"))

ggsave("fig5.Profile_of_CTCF_sites_at_loop_domains.pdf", width = 3.8, height = 2.5)