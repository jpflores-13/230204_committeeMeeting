library(InteractionSet)
library(plotgardener)
library(mariner)
library(tidyverse)

loops <- readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions()

## create an mcol for loop size
mcols(loops)$loop_size <- pairdist(loops)

## create an mcol for loop type
mcols(loops)$loop_type <- case_when(
  mcols(loops)$padj < 0.1 & mcols(loops)$log2FoldChange > 1 & 
    mcols(loops)$loop_size >= 150000 ~ "truegained",
  mcols(loops)$padj < 0.1 & mcols(loops)$log2FoldChange > 0 ~ "gained",
  mcols(loops)$padj < 0.1 & mcols(loops)$log2FoldChange < 0 ~ "lost",
  mcols(loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

loops <- loops[mcols(loops)$loop_type == c("lost") &
                 mcols(loops)$loop_size >= 300e3]
loops <- loops[order(mcols(loops)$loop_size, 150, decreasing = T)]

loops_gr <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = loops, "first"))),
          ranges = IRanges(start = start(anchors(loops, "first")),
                           end = end(anchors(loops, "second"))))

# saveRDS(loops, "data/processed/hic/hg38/diffLoops/noDroso/fish_loopCandidates.rds")

# Create survey plots to visualize FISH candidates ------------------------

## Expand regions by buffer
buffer <- 250e3
loops_buffed <- loops_gr + buffer

##make pdf
pdf(file = "plots/findFISHcandidates_lostLoops.pdf",
    width = 5.75,
    height = 5.6)

## Loop through each region
for(i in seq_along(loops_buffed)){
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = as.character(seqnames(loops_buffed))[i],
                chromstart = start(loops_buffed)[i],
                chromend = end(loops_buffed)[i],
                zrange = c(0,100),
                norm = "SCALE",
                x = 0.25,
                width = 5,
                height = 2,
                length = 5,
                fill = "#37a7db",
                linecolor = "#37a7db")
  
  
  # Begin Visualization -----------------------------------------------------
  ## Make page
  pageCreate(width = 5.75, height = 5.6,
             xgrid = 0, ygrid = 0, showGuides = F)
  
  ## Plot middle Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls 
  
  control <- plotHicRectangle(data = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic",
                              params = p,
                              y = 0.5)
  
  annoPixels(plot = control,
             data = loops,
             shift = 0.5,
             type = "arrow",
             col = "#005AB5")
  
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
  
  ## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls
  
  sorb <- 
    plotHicRectangle(data = "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic", 
                     params = p,
                     y = 2.6) 
  
  annoPixels(plot = sorb,
             data = loops,
             shift = 0.5,
             type = "arrow",
             col = "#005AB5")
  
  annoHeatmapLegend(sorb, orientation = "v",
                    fontsize = 8,
                    fontcolor = "black",
                    digits = 2,
                    x = 5.5,
                    y = 2.6,
                    width = 0.1,
                    height = 1.5,
                    length = 5,
                    just = c("left", "top"),
                    default.units = "inches")
  
  ## Plot genes
  plotGenes(param = p,
            chrom = p$chrom,
            x = 0.25,
            y = 4.65,
            height = 0.5)
  
  ## Plot genome label
  plotGenomeLabel(params = p,
                  x = 0.25,
                  y = 5.25)
  
  ## Plot statistics
  plotText(label = paste0("Size:", mcols(loops)$loop_size[i]),
           params = p, 
           x = 4.25,
           y = 2.65,
           just = c("top", "left"),
           fontsize = 10, 
           fontcolor = "black")
  
  plotText(label = paste0("log2FC : ", signif(mcols(loops)$log2FoldChange[i], digits = 3)),
           params = p, 
           x = 4.25,
           y = 2.8,
           just = c("top", "left"),
           fontsize = 10,
           fontcolor = "black")
  
  plotText(label = paste0("padj < ", signif(mcols(loops)$pvalue[i], digits = 3)),
           params = p, 
           x = 4.25,
           y = 2.95,
           just = c("top", "left"),
           fontsize = 10,
           fontcolor = "black")
  
  plotText(label = "untreated",
           x = 0.25,
           y = 0.5,
           just = c("top", "left"),
           fontcolor = "#0076BA")
  
  plotText(label = "+ sorbitol",
           x = 0.25,
           y = 2.6,
           just = c("top", "left"),
           fontcolor = "#0076BA")
  
}
dev.off()