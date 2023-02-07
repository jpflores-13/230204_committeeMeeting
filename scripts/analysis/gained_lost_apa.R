## create APA plots for Hi-C data

## load packages
library(InteractionSet)
library(mariner)
library(plotgardener)
library(glue)
library(hictoolsr)
library(RColorBrewer)

# load data ---------------------------------------------------------------
cond <- c("cont", "sorb")
hicFiles <- list.files(glue("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/{cond}"),
                          full.names = T)

diffLoops <- readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  mariner::as_ginteractions() |> 
  filterBedpe(res = 10e3, buffer = 10)

mcols(diffLoops)$loop_size <- pairdist(diffLoops)

# Subset diffLoops --------------------------------------------------------

gainedLoops <- diffLoops |> 
  subset(padj < 0.1 & log2FoldChange > 0)

lostLoops <- diffLoops |> 
  subset(padj < 0.1 & log2FoldChange < 0)

# calculate APA matrices --------------------------------------------------
gained_contAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = hicFiles[1], ##cont
                  half = "both",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

gained_sorbAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = hicFiles[2], ##sorb
                  half = "both",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

lost_contAPA <- lostLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = hicFiles[1], ##cont
                  half = "both",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

lost_sorbAPA <- lostLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = hicFiles[2], ##sorb
                  half = "both",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

## Divide each condition by nLoops
gainedN <- length(gainedLoops)
lostN <- length(lostLoops)

gained_contAPA_norm <- (gained_contAPA/gainedN)
gained_sorbAPA_norm <- (gained_sorbAPA/gainedN)

lost_contAPA_norm <- (lost_contAPA/lostN)
lost_sorbAPA_norm <- (lost_sorbAPA/lostN)

## combine APA matrices to pull out the max value for zrange max
gained_mats_combined <- c(gained_contAPA_norm, gained_sorbAPA_norm)
lost_mats_combined <- c(lost_contAPA_norm, lost_sorbAPA_norm)

# Visualize APAs ----------------------------------------------------------
pdf(file = "plots/gained_lost_apa.pdf",
    width = 5.5, height = 5.5)
pageCreate(width = 5.5, height = 5.5, showGuides = F)

p <- pgParams(assembly = "hg38",
              height = 2,
              width = 2,
              palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu")))


plotApa(gained_contAPA_norm,
        x = 0.5,
        y = 0.5,
        zrange = c(0, max(gained_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 2.6,
                    y = 1,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(gained_sorbAPA_norm,
        x = 0.5,
        y = 2.6,
        zrange = c(0, max(gained_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 2.6,
                    y = 3.5,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(lost_contAPA_norm,
        x = 3,
        y = 0.5,
        zrange = c(0, max(lost_mats_combined)),
        params = p)  |> 
  annoHeatmapLegend(x = 5.1,
                    y = 1,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotApa(lost_sorbAPA_norm,
        x = 3,
        y = 2.6,
        zrange = c(0, max(lost_mats_combined)),
        params = p) |> 
  annoHeatmapLegend(x = 5.1,
                    y = 3.5,
                    width = 0.1,
                    height = 1,
                    fontcolor = "black")

plotText(label = "Gained",
         fontcolor = "#0076BA",
         fontfamily = "Menlo",
         x = 1.5,
         y = 0.35,
         just = c("center", "top"))

plotText(label = "Lost/Existing",
         fontcolor = "#0076BA",
         fontfamily = "Menlo",
         x = 4,
         y = 0.35,
         just = c("center", "top"))

plotText(label = "untreated",
         fontcolor = "#0076BA",
         fontfamily = "Menlo",
         x = 0.35,
         y = 1.5,
         rot = 90)

plotText(label = "+ sorbitol",
         fontcolor = "#0076BA",
         fontfamily = "Menlo",
         x = 0.35,
         y = 3.5,
         rot = 90)
dev.off()
