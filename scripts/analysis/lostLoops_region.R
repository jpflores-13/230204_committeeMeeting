## Tile chromosome 1 and observe differences between untreated & sorbitol-treated cells

## Load libraries
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Set up tiles ------------------------------------------------------------

tile <- 
  tileGenome(seqlengths = seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene),
             tilewidth = 25e5,
             cut.last.tile.in.chrom = T) |> 
  keepStandardChromosomes(species = "Homo_sapiens", pruning.mode = "coarse")

tile <- tile[seqnames(tile) == "chr1"]

# Set-up params -----------------------------------------------------------
pdf(file = "plots/lostLoops_region.pdf",
    width = 5.75,
    height = 5.5)


p <- pgParams(assembly = "hg38",
              resolution = 10e3,
              chrom = as.character(seqnames(tile)[31]),
              chromstart = start(tile)[31],
              chromend = end(tile)[31],
              zrange = c(0,50),
              norm = "SCALE",
              x = 0.25,
              width = 5,
              length = 5,
              height = 2,
              fill = "#37a7db",
              linecolor = "#37a7db")


# Begin visualization -----------------------------------------------------

## Make page
pageCreate(width = 5.75, height = 5.5,
           xgrid = 0, ygrid = 0, showGuides = F)


cont <- plotHicRectangle(data = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                         params = p,
                         y = 0.5)

annoHeatmapLegend(cont,
                  orientation = "v",
                  fontsize = 8,
                  fontcolor = "black",
                  digits = 2,
                  x = 5.5,
                  y = 0.5,
                  width = 0.1,
                  height = 1.5,
                  just = c("left", "top"),
                  default.units = "inches")
# 
# annoPixels(cont,
#            data = readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/processed/hic/noDroso_cont_loops_10kb.rds"),
#            shift = 0.5,
#            type = "arrow",
#            col = "#005AB5")

sorb <- plotHicRectangle(data = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic",
                         params = p,
                         y = 2.6)

annoHeatmapLegend(sorb,
                  orientation = "v",
                  fontsize = 8,
                  fontcolor = "black",
                  digits = 2,
                  x = 5.5,
                  y = 2.5,
                  width = 0.1,
                  height = 1.5,
                  just = c("left", "top"),
                  default.units = "inches")

# annoPixels(sorb,
#            data = readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/processed/hic/noDroso_cont_loops_10kb.rds"),
#            shift = 0.5,
#            type = "arrow",
#            col = "#005AB5")


## Plot Gene Track
plotGenes(params = p,
          chrom = p$chrom,
          height = 0.5,
          y = 4.6)

## Plot Gene Region
plotGenomeLabel(params = p,
                y = 5.1)

plotText(label = "control",
         x = 0.25,
         y = 0.5,
         just = c("top", "left"))

plotText(label = "+ sorbitol",
         x = 0.25,
         y = 2.6,
         just = c("top", "left"))



dev.off()
