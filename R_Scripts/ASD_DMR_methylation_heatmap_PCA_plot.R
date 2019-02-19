## ASD DMR methylation heatmap and PCA plot
# Packages ####
library(ggdendro)
library(ggplot2)
library(reshape2)
library(grid)
library(scales)
library(plyr)
library(reshape)
library(ggbiplot)
library(ggplot2)
library(extrafont)

################################################
# Functions ####
################################################
combineFiles <- function(chroms, prefix, suffix){
# Combines chromosome-level DMR methylation files into one data.frame
        DMRs <- NULL
        for(i in 1:length(chroms)){
                temp <- NULL
                if(file.exists(paste(prefix,chroms[i],suffix, sep=""))){
                        temp <- read.delim(paste(prefix,chroms[i],suffix, sep=""), header = TRUE, sep = "\t", stringsAsFactors=FALSE)
                        DMRs <- rbind(DMRs, temp)
                }
        }
        DMRs
}

# Heatmap with Pheno Data Functions
mydplot_pheno <- function(ddata, row=!col, col=!row, labels=col) {
# plots a dendrogram
        yrange <- range(ddata$segments$y)
        yd <- yrange[2] - yrange[1]
        nc <- max(nchar(as.character(ddata$labels$label)))
        tangle <- if(row) { 0 } else { 90 }
        tshow <- col
        p <- ggplot() +
                geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend), lwd = 0.45) +
                labs(x = NULL, y = NULL) + theme_dendro()
        if(row) {
                p <- p +
                        scale_x_continuous(expand=c(0.5/length(ddata$labels$x),0)) +
                        coord_flip()
        } else {
                p <- p +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = "black"))
        }
        return(p)
}

g_legend_pheno <-function(a.gplot){
# plots a legend
# from http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

ggheatmap.show_pheno <- function(L, widths=c(0.02,0.79,0.17,0.02), heights=c(0.02,0.17,0.04,0.75,0.02)){
# plots the heatmap
        grid.newpage()
        top.layout <- grid.layout(nrow = 5, ncol = 4,
                                  widths = unit(widths, "null"),
                                  heights = unit(heights, "null"))
        pushViewport(viewport(layout=top.layout))
        print(L$col, vp=viewport(layout.pos.col=2, layout.pos.row=2))
        print(L$row, vp=viewport(layout.pos.col=3, layout.pos.row=4))
        ## print centre without legend
        print(L$centre +
                      theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            axis.title.x=element_blank(),axis.title.y=element_blank(),
                            legend.position="none",
                            panel.background=element_blank(),
                            panel.border=element_blank(),panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),plot.background=element_blank()),
              vp=viewport(layout.pos.col=2, layout.pos.row=4))
        
        print(L$phenoData +
                      theme_bw(base_size = 24) +
                      theme(panel.grid.major = element_blank(), panel.border = element_blank(),
                            legend.key = element_blank(), legend.key.size = unit(1, "lines"),
                            panel.grid.minor = element_blank(), legend.position = "none", 
                            legend.background = element_blank(), legend.text = element_text(size = 12, color = "Black"),
                            plot.margin = unit(c(0,0,0,-0.45), "lines"), 
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.title = element_blank(), 
                            legend.title = element_blank(),
                            plot.title = element_blank()), vp=viewport(layout.pos.col=2, layout.pos.row=3))
        
        ## add heatmap legend
        legend <- g_legend_pheno(L$centre +
                                         theme(legend.title = element_blank(), 
                                               legend.text = element_text(size = 15),
                                               legend.background = element_blank(),
                                               legend.position = c(0.54, -0.2)))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
        grid.draw(legend)
        upViewport(0)
        
        ## add pheno legend
        phenoLegend <- g_legend_pheno(L$phenoData +
                                              theme(legend.title = element_text(size=16), 
                                                    legend.text = element_text(size = 15),
                                                    legend.direction = "vertical",
                                                    legend.position = c(0.917, 0.92),
                                                    legend.background = element_blank()))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=3))
        grid.draw(phenoLegend)
        upViewport(0)
}

ggheatmap2_pheno <- function(x, phenoData, name="Diagnosis", breaks=c("TD", "ASD"), values=c("TD"="#3366CC", "ASD"="#FF3366"), hm.colours=c("#0072ff", "black", "#FF0000"), 
                             my.values=c(0,0.5,1), low=min(x), high=max(x)) {
# makes the heatmap
        if(is.null(colnames(x)))
                colnames(x) <- sprintf("col%s",1:ncol(x))
        if(is.null(rownames(x)))
                rownames(x) <- sprintf("row%s",1:nrow(x))
        ## plot a heatmap
        ## x is an expression matrix
        row.hc <- hclust(dist(x), "ward.D")
        col.hc <- hclust(dist(t(x)), "ward.D")
        row.dendro <- dendro_data(as.dendrogram(row.hc),type="rectangle")
        col.dendro <- dendro_data(as.dendrogram(col.hc),type="rectangle")
        
        ## dendro plots
        col.plot <- mydplot_pheno(col.dendro, col=TRUE, labels=FALSE) +
                theme(plot.margin = unit(c(0,-1.8,0,-1.9), "lines"),
                      axis.text.x = element_blank())
        row.plot <- mydplot_pheno(row.dendro, row=TRUE, labels=FALSE) +
                theme(plot.margin = unit(c(0,2,0,0), "lines"))
        
        ## order of the dendros
        col.ord <- match(col.dendro$labels$label, colnames(x))
        row.ord <- match(row.dendro$labels$label, rownames(x))
        xx <- x[row.ord,col.ord]
        dimnames(xx) <- NULL
        xx <- melt(xx)
        
        # Heatmap
        centre.plot <- ggplot(xx, aes(X2,X1)) + 
                geom_tile(aes(fill=value, colour=value)) +
                scale_fill_gradientn(colours = hm.colours, values = my.values, limits = c(low, high), na.value = "black") +
                scale_colour_gradientn(colours = hm.colours, values = my.values, limits = c(low, high), na.value = "black") +
                labs(x = NULL, y = NULL) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0),breaks = NULL) +
                theme(plot.margin = unit(rep(0, 4), "lines"))
        
        # phenoData
        sample.ord <- match(col.dendro$labels$label, as.character(phenoData$Sample))
        phenoData$Sample <- factor(as.character(phenoData$Sample), levels = as.character(phenoData$Sample)[sample.ord], ordered = TRUE)
        phenoData_m <- melt(phenoData, id.vars = "Sample")
        phenoData_m$variable <- factor(phenoData_m$variable, levels = rev(unique(phenoData_m$variable)), ordered = TRUE)
        phenoData.plot <- ggplot(phenoData_m, aes(Sample, variable)) +
                geom_tile(aes(fill=value, color=value)) +
                scale_x_discrete(expand=c(0,0)) +
                scale_y_discrete(expand=c(0,0)) +
                scale_color_manual(name=name, breaks = breaks, values = values) +
                scale_fill_manual(name=name, breaks = breaks, values = values)
        ret <- list(col=col.plot,row=row.plot,centre=centre.plot, phenoData=phenoData.plot)
        invisible(ret)
}

################################################
# Load Methylation and outcome data ####
################################################
# Methylation
meth <- combineFiles(chroms=paste("chr",1:22, sep=""), prefix="DMR",
                     suffix=".txt")
meth <- meth[order(meth$chr, meth$start),]
dim(meth)
write.csv(meth, "meth.csv")

# import DMRs information in 
DMRs <- read.csv("Diagnosis_DMR_infor.csv", header=T)
DMRs <- DMRs[order(DMRs$chr, DMRs$start),]
DMRs <- DMRs[!duplicated(DMRs[,1:3]),]

# check on whether all TRUE
table(DMRs$start == meth$start) # All TRUE, must be the same order

# Samples information with outcome
samples <- read.csv("Sample_Info.csv", header=TRUE, stringsAsFactors = FALSE)

################################################
# Plots on Heatmap and PCA ####
################################################
# Heatmap with Pheno Data ####
meth <- meth[,4:ncol(meth)]
phenoData <- samples[,c("Sequencing.ID", "Diagnosis")]
colnames(phenoData)[1] <- "Sample"
phenoData <- phenoData[match(colnames(meth), phenoData$Sample),]
table(colnames(meth) == phenoData$Sample) # All TRUE must be the same order
phenoData$Sample <- factor(phenoData$Sample, levels=unique(phenoData$Sample), ordered=TRUE)
phenoData$Diagnosis[phenoData$Diagnosis == "TD"] <- "TD"
phenoData$Diagnosis[phenoData$Diagnosis == "ASD"] <- "ASD"
meth <- meth[,phenoData$Diagnosis %in% c("TD", "ASD")]
phenoData <- subset(phenoData, Diagnosis %in% c("TD", "ASD"))
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels=c("TD", "ASD"), ordered=TRUE)

methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg # subtract mean methylation for each DMR across all samples
methdiff <- as.matrix(methdiff)*100 # Transform to 0-100% scale

methplot <- ggheatmap2_pheno(x=methdiff, phenoData=phenoData, name="Diagnosis", breaks=c("TD", "ASD"), values=c("TD"="#3366CC", "ASD"="#FF3366")) 
pdf(file="Heatmap.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show_pheno(methplot)
dev.off()



# PCA Plot ####
data <- t(as.matrix(meth))
diagnosis <- c(rep("TD", 21), rep("ASD", 20))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)
# PC1 
# PC2 

g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, 
              var.axes = FALSE, varname.abbrev = FALSE, choices = 1:2,ellipse.prob = 0.95) # Ellipses are 95% confidence intervals
g + 
  theme_bw(base_size = 25) +
  theme(legend.direction = 'vertical', legend.position = c(0.85, 0.9), panel.grid.major = element_blank(), 
        panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), axis.title=element_text(size=18),
        legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size=18), legend.key.size=unit(1.5, "line"),
        axis.text = element_text(color = "black", size=18), legend.background = element_blank(), legend.text=element_text(size=18)) +
  coord_cartesian(xlim = c(-20, 20), ylim = c(-20,20)) +
  xlab("PC1 (18.29% of Variance)") +
  ylab("PC2 (4.41% of Variance)") +
  scale_color_manual(name="Diagnosis", breaks = c("TD", "ASD"), 
                     values = c("TD"="#3366CC", "ASD"="#FF3366")) +
  scale_x_continuous(breaks=pretty_breaks(n=4)) +
  scale_y_continuous(breaks=pretty_breaks(n=4)) +
  geom_point(aes(color = diagnosis), size=3)+
  theme(text=element_text(family="Arial", size=14))

ggsave("PCA plot.png", dpi = 600, width = 8, height = 6, units = "in")
## Script helped by Charles Mordaunt
