library(InteractionSet)
library(strawr)
library(tidyverse)
library(glue)
library(nullranges)
library(data.table)
library(mariner)
library(raster)
library(reshape2)

source("scripts/utils/aggregateLoops.R")
source("scripts/utils/aggregateTAD.R")
source("scripts/utils/plotAggTAD.R")

#------------------------- Read in and filter data ---------------------#

## Load in each HiC replicate 
cond <- c("cont", "sorb")
hicFiles <- list.files(glue("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220722_dietJuicerCore/{cond}"),
                       full.names = T)

## Load in merged HiC files for human 
merged_hicFiles <- list.files(glue("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/{cond}"),
                              full.names = T)

## Load in loops
noDroso_loops <- readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |> 
  interactions() |> 
  BiocGenerics::as.data.frame() |> 
  mariner::as_ginteractions() |>
  filterBedpe(res = 10e3, buffer = 10) |> 
  pullHicPixels(binSize = 10e3,
                files = hicFiles,
                half = "both",
                norm = "NONE",
                matrix = "observed")

## replace 10kb extracted counts from original .rds with 1kb extracted counts
# mcols(noDroso_loops) <- cbind(as.matrix(assay(noDroso_loops)), mcols(noDroso_loops)[-c(1:6)])

## create an mcol for loop size
mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

## create an mcol for loop type
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 & 
    mcols(noDroso_loops)$loop_size >= 150000 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 &
    mcols(noDroso_loops)$loop_size >= 150000 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

## create aggregated sorb counts mcol
mcols(noDroso_loops)$sorb_contacts <- mcols(noDroso_loops)$YAPP_HEK_sorbitol_4_2_inter_30.hic +
  mcols(noDroso_loops)$YAPP_HEK_sorbitol_5_2_inter_30.hic + mcols(noDroso_loops)$YAPP_HEK_sorbitol_6_2_inter_30.hic

## create aggregated cont counts mcol
mcols(noDroso_loops)$cont_contacts <- mcols(noDroso_loops)$YAPP_HEK_control_1_2_inter_30.hic +
  mcols(noDroso_loops)$YAPP_HEK_control_2_2_inter_30.hic + mcols(noDroso_loops)$YAPP_HEK_control_3_2_inter_30.hic

## create aggregated sorb+cont counts mcol
mcols(noDroso_loops)$agg_contacts <- mcols(noDroso_loops)$sorb_contacts + mcols(noDroso_loops)$cont_contacts

## log transform loop_size to get a normal distribution for matchRanges
mcols(noDroso_loops)$loop_size <- log(mcols(noDroso_loops)$loop_size)
hist(mcols(noDroso_loops)$loop_size) # check for normality

## convert all `0` agg_contact values to NAs and log transform for matchRanges
mcols(noDroso_loops)$agg_contacts[mcols(noDroso_loops)$agg_contacts == 0] <- NA
noDroso_loops <- interactions(noDroso_loops) |> # no longer InteractionMatrix class
  as.data.frame() |>
  na.omit() |>
  mariner::as_ginteractions()

mcols(noDroso_loops)$agg_contacts <- log((mcols(noDroso_loops)$agg_contacts + 1))
hist(mcols(noDroso_loops)$agg_contacts) # check for normality

# Matched set for gained YAPP loops ---------------------------------------
## use matchRanges to select a null set of control sample loops that is matched for size & contact frequency
# noDroso_loops <- readRDS('data/processed/hic/hg38/noDroso_loops.rds')
focal <- noDroso_loops[!noDroso_loops$loop_type %in% c("static", "lost","other")] 
pool <- noDroso_loops[noDroso_loops$loop_type %in% c("static", "lost","other")] 

nullSet <- matchRanges(focal = focal,
                       pool = pool,
                       covar = ~ agg_contacts + loop_size, 
                       method = 'stratified',
                       replace = F)

plotCovariate(nullSet)
plotCovariate(nullSet, covar = "loop_size")
plotPropensity(nullSet, sets = c('f', 'p', 'm'), log = 'x')


gainedLoops <- focal 
ctcfLoops <- pool

gained_bed <- gainedLoops |> 
  as.data.frame() |> 
  dplyr::select(c(1,2,3,6,7,8))

nullSet <- nullSet |> 
  as.data.frame() |> 
  dplyr::select(c(1,2,3,6,7,8))

#------------------------- Plot ATAs of existing vs gained loops in HEKs and T47Ds ---------------------#

aggtad_gain <- aggregateTAD(loops = gained_bed,
                            hic = merged_hicFiles[2],
                            res = 10e3,
                            buffer = 0.5,
                            norm = "SCALE")

aggtad_match <- aggregateTAD(loops = nullSet,
                             hic = merged_hicFiles[1],
                             res = 10e3,
                             buffer = 0.5,
                             norm = "SCALE")

plotAggTAD(aggtad_gain,maxval = aggtad_gain[26,75]*1.2,title="Gained")

plotAggTAD(aggtad_match,maxval = aggtad_match[26,75]*1.2,title="Lost/Existing Loops")

#------------------------- Plot TAD enrichment ---------------------#

## palette of choice
# brewer.pal(n = 9, name = "YlGnBu")
# "#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0" "#225EA8" "#253494" "#081D58"

par(mgp=c(3,.3,0))
plot(1,type="n",xaxt="n",yaxt="n",ann=F,xlim=c(1,76),ylim=c(.9,1.6))
abline(v=c(26,51),lty=2,col="grey90")

data = aggtad_match[cbind(1:76, 25:100)]
bg   = median(data[c(1:20,56:76)])
lines(data/bg,type="l",col="#33A02C")

data = aggtad_gain[cbind(1:76, 25:100)]
bg   = median(data[c(1:20,56:76)])
lines(data/bg,type="l",col="#1D91C0")

axis(side=2,las=2,tcl=0.2,at=seq(1,1.6,.2))
axis(side=1,at = c(26,51),labels=c("left","right"))
legend("topright",bty="n",col=c("#33A02C","#1D91C0"),legend=c("matchedSet","Gained"),lty=1)


#------------------------- PLot APAs of gained loops in HEKs ---------------------#

buffer = 10
res = 10e3
norm="NONE"

apa_gain = aggregateLoops(loops = gained_bed,
                          hic = merged_hicFiles[2],
                          res = res,buffer = buffer,norm=norm)
apa_nullSet = aggregateLoops(loops = nullSet,
                             hic = merged_hicFiles[1],
                             res = res,buffer = buffer,norm=norm)
# plot HEK data
plotAggTAD(apa_gain,maxval = 50,title = "")
plotAggTAD(apa_nullSet,maxval = 50,title = "")

#------------------------- Plot Loop enrichment ---------------------#


plot(1,type="n",xaxt="n",yaxt="n",ann=F,xlim=c(1,50),ylim=c(.9,3.55))
abline(v=26,lty=2,col="grey90")

data = aggtad_match[cbind(1:51, 50:100)]
bg   = median(data[c(1:10,40:50)])
lines(data/bg,type="l",col="#33A02C")

data = aggtad_gain[cbind(1:50, 50:99)]
bg   = median(data[c(1:10,40:50)])
lines(data/bg,type="l",col="#1D91C0")

axis(side=2,las=2,tcl=0.2,at=seq(1,3.5,.5))
axis(side=1,at = 26,labels="loop pixel")
legend("topright",bty="n",col=c("#33A02C","#1D91C0"),legend=c("matchedSet", "Gained"),lty=1)
