## Placenta ASD DMRs Enrichment test

#################################################
######### Load the packages #####################
#################################################
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(reshape2)
library(LOLA)
library(simpleCache)
library(GenomicRanges)
library(qvalue)

#################################################
######### Use LOLA to check the enrichment ######
#################################################
# Load Files
cat("\nReading DMRs\n")
DMRs <- readBed(file = "DMR.bed")
Background <- readBed(file = "Background.bed")
DMRs_List <- GRangesList(DMRs)

# chromHMM
cat("\nLoading Regions\n")
regionDB <- loadRegionDB(dbLocation = "./", useCache = TRUE, limit = NULL, collections = "encode_placenta_chipseq")
cat("\nRegions finished loading\n")

cat("\nRunning LOLA\n")
Results <- runLOLA(userSets = DMRs_List, userUniverse = Background, regionDB = regionDB, minOverlap = 1, cores=2, redefineUserSets = FALSE)
cat("\nLOLA Finished Running\n")

cat("\nPrinting Results\n")
writeCombinedEnrichment(combinedResults = Results, outFolder = "Hyper_DMRs_roadmap_placenta_new", includeSplits=FALSE)

cat("\nDone!\n")

#################################################
######### Plot the enrichment results ###########
#################################################
######### CpG location #####################
CpG_location =read.csv("CpG_location.csv", header=T)
# Create column of significance labels
CpG_location$stars <- cut(CpG_location$qValue, 
                          breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
                          label=c("***", "**", "*", ""))  

head(CpG_location)

pdf(file = "CpG_location.PDF", width = 6, height = 7, family = "Helvetica")
gg <- ggplot(data = CpG_location)
gg +
  geom_tile(aes(y = CpG_location, x = Category, fill = oddsRatio)) +
  scale_fill_gradientn("logOR", colors = c("Black", "#FF0000")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1),
        axis.ticks = element_line(size = 1), legend.key = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.background = element_blank(),
        plot.margin = unit(c(1,7,1,1), "lines"), axis.text.y = element_text(size = 14, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
        axis.title = element_blank(), 
        axis.text.x = element_text(size = 14, color = "Black"),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 2))+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_discrete(expand=c(0,0))
dev.off()

