library(ggplot2)

base <- read.table(file = "fig3.Cumulative_percentage_plot.txt", header = T)

base_s <- base[order(base$salt_tolerance, decreasing = T),]
base_s$rank <- seq(1, nrow(base_s))
base_s <- base_s[base_s$U == 1,]
base_s$rank_type <- "S-score"

base_x <- base[order(base$xchip_signal, decreasing = T),]
base_x$rank <- seq(1, nrow(base_x))
base_x <- base_x[base_x$U == 1,]
base_x$rank_type <- "X-ChIP"

base_run <- base[order(base$CUTRUN_signal, decreasing = T),]
base_run$rank <- seq(1, nrow(base_run))
base_run <- base_run[base_run$U == 1,]
base_run$rank_type <- "CUTRUN"

base_tag <- base[order(base$CUTTag_signal, decreasing = T),]
base_tag$rank <- seq(1, nrow(base_tag))
base_tag <- base_tag[base_tag$U == 1,]
base_tag$rank_type <- "CUTTag"

all_rank <- rbind(base_s[, c("rank_type", "rank")], base_x[, c("rank_type", "rank")], base_run[, c("rank_type", "rank")], base_tag[, c("rank_type", "rank")])

all_rank$rank <- all_rank$rank / nrow(base)

my_comparisons <- list(c("X-ChIP", "S-score"))

all_rank$rank_type <- factor(all_rank$rank_type, levels = c("S-score", "CUTRUN", "CUTTag", "X-ChIP"))

ggplot(all_rank, aes(x = rank, color = rank_type)) +
  stat_ecdf(size = 1) +
  scale_color_manual(values = c("#BC1419", "grey", "#083B7B", "#9f79d3"),
                     labels = c("S-score", "CUT&RUN", "CUT&Tag", "X-ChIP")) +
  labs(title = NULL, x = "Rank of motifs (%)", y = "Cumulative % with U") +
  theme_bw() +
  theme(plot.margin = unit(c(.1, .3, .1, .1), 'inches'),
        plot.title = element_text(size = 15, hjust = .5),
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.line = element_line(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 15, vjust = -1),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, colour = 'black'),
        axis.text.y = element_text(size = 15, colour = 'black'),
        legend.position = c(0.7, 0.3),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = 'black')) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c(0, 50, 100), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = c(0, format(seq(20, 100, 20), big.mark = ",")), expand = c(0, 0))

ggsave("fig3.Cumulative_percentage_plot.pdf", width = 3.2, height = 3)