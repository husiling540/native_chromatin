library(ChIPpeakAnno)

peak_XChIP <- read.table("K562_xchip_sites.bed", header = F, sep = "\t")
names(peak_XChIP)[1:4] <- c("chr", "start", "end", "name")
gr_peak_XChIP <- GRanges(peak_XChIP, format = "BED", header = FALSE)

peak_NChIP <- read.table("K562_native_sites.bed", header = F, sep = "\t")
names(peak_NChIP)[1:4] <- c("chr", "start", "end", "name")
gr_peak_NChIP <- GRanges(peak_NChIP, format = "BED", header = FALSE)

peak_LNCaP_persistent <- read.table("LNCaP_ChIP-seq_CTCF_siRNA_persistent_peaks.rmBlacklist.narrowPeak", header = F, sep = "\t")
names(peak_LNCaP_persistent)[1:4] <- c("chr", "start", "end", "name")
gr_peak_LNCaP_persistent <- GRanges(peak_LNCaP_persistent, format = "BED", header = FALSE)

gr_peak_XChIP <- unique(gr_peak_XChIP)
gr_peak_NChIP <- unique(gr_peak_NChIP)
gr_peak_LNCaP_persistent <- unique(gr_peak_LNCaP_persistent)

ol_XChIP <- findOverlapsOfPeaks(gr_peak_XChIP, gr_peak_NChIP, gr_peak_LNCaP_persistent)

pdf("fig5.VennDiagram.pdf", width = 9, height = 9)
makeVennDiagram(
  ol_XChIP,
  fill = c("#A085BD", "#82B2E2", "#E0853F"),
  col = NA,
  cat.col = 'black',
  cat.cex = 4,
  cex = 4,
  lwd = 0,
  NameOfPeaks = c("XChIPsites", "NChIPsites", "LNCaP_persistent"),
  margin = 0.1,
  fontfamily = "sans"
)
dev.off()