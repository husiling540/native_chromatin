library(ggplot2)
library(scales)

Matrix <- read.table(file = "motifcenter.K562_footprint.Matrix", sep = "\t", header = F, skip = 3)

XChIP_3border <- data.frame(experiment = "XChIP_border3", pos = c(1:200), sig = colMeans(Matrix[, c(1:200)]) * -1)
XChIP_5border <- data.frame(experiment = "XChIP_border5", pos = c(1:200), sig = colMeans(Matrix[, c(201:400)]))
NChIP_3border <- data.frame(experiment = "NChIP_border3", pos = c(1:200), sig = colMeans(Matrix[, c(401:600)]) * -1)
NChIP_5border <- data.frame(experiment = "NChIP_border5", pos = c(1:200), sig = colMeans(Matrix[, c(601:800)]))
CUTRUN_3border <- data.frame(experiment = "CUT-RUN_border3", pos = c(1:200), sig = colMeans(Matrix[, c(801:1000)]) * -1)
CUTRUN_5border <- data.frame(experiment = "CUT-RUN_border5", pos = c(1:200), sig = colMeans(Matrix[, c(1001:1200)]))
CUTTag_3border <- data.frame(experiment = "CUT-Tag_border3", pos = c(1:200), sig = colMeans(Matrix[, c(1201:1400)]) * -1)
CUTTag_5border <- data.frame(experiment = "CUT-Tag_border5", pos = c(1:200), sig = colMeans(Matrix[, c(1401:1600)]))
ChIPexo_3border <- data.frame(experiment = "ChIPexo_border3", pos = c(1:200), sig = colMeans(Matrix[, c(1601:1800)]) * -1)
ChIPexo_5border <- data.frame(experiment = "ChIPexo_border5", pos = c(1:200), sig = colMeans(Matrix[, c(1801:2000)]))

ChIPexo_max <- max(abs(ChIPexo_3border$sig), abs(ChIPexo_5border$sig))
NChIP_max <- max(abs(NChIP_3border$sig), abs(NChIP_5border$sig))
CUTRUN_max <- max(abs(CUTRUN_3border$sig), abs(CUTRUN_5border$sig))
CUTTag_max <- max(abs(CUTTag_3border$sig), abs(CUTTag_5border$sig))
XChIP_max <- max(abs(XChIP_3border$sig), abs(XChIP_5border$sig))

ChIPexo_3border_sum <- sum(abs(ChIPexo_3border$sig))
ChIPexo_5border_sum <- sum(abs(ChIPexo_5border$sig))
NChIP_3border_sum <- sum(abs(NChIP_3border$sig))
NChIP_5border_sum <- sum(abs(NChIP_5border$sig))
CUTRUN_3border_sum <- sum(abs(CUTRUN_3border$sig))
CUTRUN_5border_sum <- sum(abs(CUTRUN_5border$sig))
CUTTag_3border_sum <- sum(abs(CUTTag_3border$sig))
CUTTag_5border_sum <- sum(abs(CUTTag_5border$sig))
XChIP_3border_sum <- sum(abs(XChIP_3border$sig))
XChIP_5border_sum <- sum(abs(XChIP_5border$sig))

Matrix_all_profile <- data.frame(pos = c(1:200),
                                 NChIP_5border = NChIP_5border$sig, NChIP_3border = NChIP_3border$sig,
                                 CUTRUN_5border = CUTRUN_5border$sig, CUTRUN_3border = CUTRUN_3border$sig,
                                 CUTTag_5border = CUTTag_5border$sig, CUTTag_3border = CUTTag_3border$sig,
                                 ChIPexo_5border = ChIPexo_5border$sig, ChIPexo_3border = ChIPexo_3border$sig,
                                 XChIP_5border = XChIP_5border$sig, XChIP_3border = XChIP_3border$sig)

annotation <- data.frame(x = c(30, 30), y = c(0.027, -0.027), label = c("up", "down"))

ggplot(Matrix_all_profile, aes(x = pos)) +
  geom_area(aes(y = NChIP_5border / NChIP_5border_sum, fill = "NChIP"), linetype = 1, lwd = 0.8, alpha = 0.5) +
  geom_area(aes(y = NChIP_3border / NChIP_3border_sum, fill = "NChIP"), linetype = 1, lwd = 0.8, alpha = 0.5) +
  geom_line(aes(y = CUTRUN_5border / CUTRUN_5border_sum, colour = "CUTRUN"), lwd = 0.8) +
  geom_line(aes(y = CUTRUN_3border / CUTRUN_3border_sum, colour = "CUTRUN"), lwd = 0.8) +
  geom_line(aes(y = CUTTag_5border / CUTTag_5border_sum, colour = "CUTTag"), lwd = 0.8) +
  geom_line(aes(y = CUTTag_3border / CUTTag_3border_sum, colour = "CUTTag"), lwd = 0.8) +
  geom_line(aes(y = XChIP_5border / XChIP_5border_sum, colour = "XChIP"), lwd = 0.8) +
  geom_line(aes(y = XChIP_3border / XChIP_3border_sum, colour = "XChIP"), lwd = 0.8) +
  geom_line(aes(y = ChIPexo_5border / ChIPexo_5border_sum, color = "ChIPexo"), lwd = 0.8) +
  geom_line(aes(y = ChIPexo_3border / ChIPexo_3border_sum, color = "ChIPexo"), lwd = 0.8) +
  geom_text(data = annotation, aes(x = x, y = y, label = label), size = 7) +
  scale_fill_manual(values = c(NChIP = "#4A68A0")) +
  scale_colour_manual(values = c(ChIPexo = "#9f79d3", CUTRUN = "#8FBA48", CUTTag = "#DC9338", XChIP = "#2D68B1")) +
  scale_linetype_manual(values = c(ChIPexo = "solid", CUTRUN = "solid", CUTTag = "solid", XChIP = "solid")) +
  theme_bw() +
  labs(x = "Distance from CTCF motif (bp)", y = "Scaled signal") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 200), breaks = c(0, 50, 100, 150, 200), labels = c("-100", "-50", "0", "50", "100")) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.035, 0.035), breaks = seq(-0.04, 0.04, 0.01)) +
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) +
  theme(plot.margin = unit(c(.2, .2, 0.2, .2), 'inches'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1, colour = "black"),
        axis.line = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 20))

ggsave("fig2.Scaled_pileup_of_fragment_ends.pdf", width = 7.5, height = 5.1)