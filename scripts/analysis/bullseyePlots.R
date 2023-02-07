## create bullseye plots for Hi-C data

## load packages
library(InteractionSet)
library(mariner)
library(plotgardener)
library(glue)
library(hictoolsr)
library(RColorBrewer)

source("scripts/utils/plotBullseye.R")

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

## create an mcol for loop type
mcols(diffLoops)$loop_type <- case_when(
  mcols(diffLoops)$padj < 0.05 & mcols(diffLoops)$log2FoldChange > 1 ~ "truegained",
  mcols(diffLoops)$padj < 0.1 & mcols(diffLoops)$log2FoldChange > 0 ~ "gained",
  mcols(diffLoops)$padj < 0.1 & mcols(diffLoops)$log2FoldChange < 0 ~ "lost",
  mcols(diffLoops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

## create aggregated sorb counts mcol
mcols(diffLoops)$sorb_contacts <- mcols(diffLoops)$YAPP_HEK_sorbitol_4_2_inter_30.hic +
  mcols(diffLoops)$YAPP_HEK_sorbitol_5_2_inter_30.hic + mcols(diffLoops)$YAPP_HEK_sorbitol_6_2_inter_30.hic

## create aggregated cont counts mcol
mcols(diffLoops)$cont_contacts <- mcols(diffLoops)$YAPP_HEK_control_1_2_inter_30.hic +
  mcols(diffLoops)$YAPP_HEK_control_2_2_inter_30.hic + mcols(diffLoops)$YAPP_HEK_control_3_2_inter_30.hic

## create aggregated sorb+cont counts mcol
mcols(diffLoops)$agg_contacts <- mcols(diffLoops)$sorb_contacts + mcols(diffLoops)$cont_contacts

## log transform loop_size to get a normal distribution for matchRanges
mcols(diffLoops)$loop_size <- log(mcols(diffLoops)$loop_size)
hist(mcols(diffLoops)$loop_size) # check for normality

## convert all `0` agg_contact values to NAs and log transform for matchRanges
mcols(diffLoops)$agg_contacts[mcols(diffLoops)$agg_contacts == 0] <- NA
mcols(diffLoops)$agg_contacts <- log((mcols(diffLoops)$agg_contacts + 1))
hist(mcols(diffLoops)$agg_contacts) # check for normality

# Matched set for gained YAPP loops ---------------------------------------
## use matchRanges to select a null set of control sample loops that is matched for size & contact frequency
focal <- diffLoops[!diffLoops$loop_type %in% c("static", "lost","other")] 
pool <- diffLoops[diffLoops$loop_type %in% c("static", "lost","other")] 

nullSet <- matchRanges(focal = focal,
                       pool = pool,
                       covar = ~ agg_contacts + loop_size, 
                       method = 'stratified',
                       replace = F)

## QC
# plotCovariate(nullSet)
# plotCovariate(nullSet, covar = "loop_size")
# plotPropensity(nullSet, sets = c('f', 'p', 'm'), log = 'x')

# Calculate loop decay for all loop sets (S/O Manjari!)
gainedLoops <- focal

# calculate APA matrices --------------------------------------------------

gained_sorbAPA <- gainedLoops |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = hicFiles[2], ##sorb
                  half = "both",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

gained_sorbAPA_norm <- (gained_sorbAPA/length(gainedLoops))

lost_contAPA <- nullSet |> 
  pixelsToMatrices(buffer = 10) |> 
  pullHicMatrices(binSize = 10e3,
                  files = hicFiles[1], ##cont
                  half = "both",
                  norm = "NONE",
                  matrix = "observed") |>
  aggHicMatrices(FUN = sum) |> 
  as.matrix()

lost_contAPA_norm <- (lost_contAPA/length(nullSet))

# Visualize Bullseye ------------------------------------------------------

plotBullseye(x = gained_sorbAPA,
             title = "",
             nl = length(gainedLoops),
             zrange=c(0,max(gained_sorbAPA_norm)),
             cols = brewer.pal(n = 9, name = "YlGnBu"))

plotBullseye(x = lost_contAPA,
             title = "",
             nl = length(nullSet),
             zrange=c(0,max(lost_contAPA_norm)),
             cols = brewer.pal(n = 9, name = "YlGnBu"))

