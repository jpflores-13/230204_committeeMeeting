## Draft figure for Fig. 1B

# Load data ---------------------------------------------------------------
library(plotgardener)
library(InteractionSet)
library(glue)
library(dbscan)
library(data.table)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(mariner)
library(plyranges)
library(stringr)

diff_loopCounts <- readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions()
  
# Setting static and lost loops -------------------------------------------

## gained loops with a p-adj. value of < 0.1 and a (+) log2FC
gained_adj <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange > 0)
lost_loops <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange < 0)

## filter for the 100 best gained loops
bestGained <- head(gained_adj[order(gained_adj$log2FoldChange, decreasing = T)],100)
bestGained <- head(gained_adj[order(gained_adj$padj, decreasing = F)], 100)

# top 100 gained loops
## If the below function doesn't work, might need to use swapAchors()

loopRegions_gained <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = bestGained, "first"))),
          ranges = IRanges(start = start(anchors(bestGained, "first")),
                           end = end(anchors(bestGained, "second"))))

## Expand regions by buffer
# buffer <- 200e3
loopRegions_gained <- 
  GRanges(seqnames = "chr1",
          ranges = IRanges(start = 64860000,
                           end = 65860000))

seqlevelsStyle(loopRegions_gained) <- "ENSEMBL"

# loopRegions_gained <- GenomicRanges::shift(loopRegions_gained + buffer, 100e3)
loopRegions_gained <- as.data.frame(loopRegions_gained)

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/gainedLoop_region.pdf",
    width = 5.75,
    height = 5.5)

## Loop through each region

## Define parameters
p <- pgParams(assembly = "hg38",
              resolution = 10e3,
              chrom = paste0("chr",loopRegions_gained$seqnames),
              chromstart = loopRegions_gained$start,
              chromend = loopRegions_gained$end,
              zrange = c(0,100),
              norm = "SCALE",
              fontfamily = "Menlo",
              x = 0.25,
              width = 5,
              length = 5,
              height = 2,
              fill = "#37a7db",
              linecolor = "#37a7db")


# Begin Visualization -----------------------------------------------------
## Make page
pageCreate(width = 5.75, height = 5.5,
           xgrid = 0, ygrid = 0, showGuides = F)

## Plot middle Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls 

control <- plotHicRectangle(data = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                            params = p,
                            y = 0.5)

annoHeatmapLegend(control, orientation = "v",
                  fontsize = 8,
                  fontcolor = "black",
                  digits = 2,
                  x = 5.5,
                  y = 0.5,
                  width = 0.1,
                  height = 1.5,
                  just = c("left", "top"),
                  default.units = "inches")

annoPixels(control,
           data = gained_adj,
           shift = 0.5,
           type = "arrow",
           col = "#DC3220")

annoPixels(control,
           data = lost_loops,
           shift = 0.5,
           type = "arrow",
           col = "#005AB5")


## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls

sorb <- 
  plotHicRectangle(data = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                   params = p,
                   y = 2.6)

annoHeatmapLegend(sorb, orientation = "v",
                  fontsize = 8,
                  fontcolor = "black",
                  digits = 2,
                  x = 5.5,
                  y = 2.6,
                  width = 0.1,
                  height = 1.5,
                  just = c("left", "top"),
                  default.units = "inches")

annoPixels(sorb,
           data = gained_adj,
           shift = 0.5,
           type = "arrow",
           col = "#DC3220")

annoPixels(sorb,
           data = lost_loops,
           shift = 0.5,
           type = "arrow",
           col = "#005AB5")

## Plot genes
plotGenes(param = p,
          chrom = p$chrom,
          x = 0.25,
          y = 4.6,
          height = 0.5)

## Plot genome label
plotGenomeLabel(params = p,
                x = 0.25,
                y = 5.1)

# Annotate Hi-C rectangles by treatment ------------------------------------

plotText(label = "untreated",
         x = 0.25,
         y = 0.5,
         just = c("top", "left"),
         fontfamily = "Menlo")

plotText(label = "+ sorbitol",
         x = 0.25,
         y = 2.6,
         just = c("top", "left"),
         fontfamily = "Menlo")

plotText(label = "gained (n = 355)",
         x = 1.35,
         y = 3.15,
         just = c("top", "left"),
         col = "#DC3220",
         fontsize = 10,
         fontfamily = "Menlo")

plotText(label = "lost (n = 16,077)",
         x = 2.8,
         y = 1.05,
         just = c("top", "left"),
         col = "#005AB5",
         fontsize = 10,
         fontfamily = "Menlo")

dev.off()
