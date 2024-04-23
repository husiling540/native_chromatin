library(ChIPpeakAnno)

SUMO_R1 <- read.table("K562_NChIP_CTCF_SUMO_R1_macs2/K562_NChIP_CTCF_SUMO_R1_peaks.rmBlacklist.narrowPeak", header = F, sep = "\t")
names(SUMO_R1) <- c("chr", "start", "end", "name")
gr_SUMO_R1 <- GRanges(SUMO_R1, format = "BED", header = FALSE)

SUMO_R2 <- read.table("K562_NChIP_CTCF_SUMO_R2_macs2/K562_NChIP_CTCF_SUMO_R2_peaks.rmBlacklist.narrowPeak", header = F, sep = "\t")
names(SUMO_R2) <- c("chr", "start", "end", "name")
gr_SUMO_R2 <- GRanges(SUMO_R2, format = "BED", header = FALSE)

noSUMO_R1 <- read.table("K562_NChIP_CTCF_noSUMO_R1_macs2/K562_NChIP_CTCF_noSUMO_R1_peaks.rmBlacklist.narrowPeak", header = F, sep = "\t")
names(noSUMO_R1) <- c("chr", "start", "end", "name")
gr_noSUMO_R1 <- GRanges(noSUMO_R1, format = "BED", header = FALSE)

noSUMO_R2 <- read.table("K562_NChIP_CTCF_noSUMO_R2_macs2/K562_NChIP_CTCF_noSUMO_R2_peaks.rmBlacklist.narrowPeak", header = F, sep = "\t")
names(noSUMO_R2) <- c("chr", "start", "end", "name")
gr_noSUMO_R2 <- GRanges(noSUMO_R2, format = "BED", header = FALSE)

gr_SUMO_R1 <- unique(gr_SUMO_R1)
gr_SUMO_R2 <- unique(gr_SUMO_R2)
gr_noSUMO_R1 <- unique(gr_noSUMO_R1)
gr_noSUMO_R2 <- unique(gr_noSUMO_R2)

ol_SUMO <- findOverlapsOfPeaks(gr_SUMO_R1, gr_SUMO_R2)
ol_noSUMO <- findOverlapsOfPeaks(gr_noSUMO_R1, gr_noSUMO_R2)

SUMO_2rep <- ol_SUMO$peaklist$`gr_SUMO_R1///gr_SUMO_R2`
names(SUMO_2rep) <- "SUMO"
noSUMO_2rep <- ol_noSUMO$peaklist$`gr_noSUMO_R1///gr_noSUMO_R2`
names(noSUMO_2rep) <- "noSUMO"
ol <- findOverlapsOfPeaks(SUMO_2rep, noSUMO_2rep)

pdf("fig4.VennDiagram_of_CTCF-SUMO23_sequential_N-ChIP_datasets.pdf", width = 8, height = 8)
makeVennDiagram(
  ol,
  fill = c("#C36518", "#2A5CA3"),
  col = NA,
  cat.col = 'black',
  cat.cex = 4,
  cex = 4,
  lwd = 3,
  NameOfPeaks = c("CTCF(+)SUMO2/3(+)", "CTCF(+)SUMO2/3(-)"),
  margin = 0.1,
  fontfamily = "sans"
)
dev.off()

SUMO_noSUMO <- data.frame(c(ol$peaklist$`SUMO_2rep///noSUMO_2rep`))
noSUMO_unique <- data.frame(c(ol$peaklist$noSUMO_2rep))

write.table(SUMO_noSUMO, file = "K562_CTCF_SUMO_sites.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(noSUMO_unique, file = "K562_CTCF_non-SUMO_sites.bed", sep = "\t", quote = F, row.names = F, col.names = F)
