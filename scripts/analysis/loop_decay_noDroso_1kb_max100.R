## Calculate (and plot) loop decay for gained PS loops, & nullSet loops

library(InteractionSet)
library(strawr)
library(tidyverse)
library(glue)
library(nullranges)
# if (!requireNamespace("remotes", quietly = TRUE))
#   install.packages("remotes")
# remotes::install_github("EricSDavis/mariner@dev", force = T)
library(mariner)
library(data.table)

## Load in each HiC replicate 
cond <- c("cont", "sorb")
hicFiles <- list.files(glue("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220722_dietJuicerCore/{cond}"),
                       full.names = T)

## Load in merged HiC files for human and d. mel
merged_hicFiles <- list.files(c("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb",
                                "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont"),
                              full.names = T)

## Load in non-d. mel loops and convert to 1kb resolution
noDroso_loops <- readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |>
  interactions() |>
  pullHicPixels(binSize = 10e3,
                files = hicFiles,
                half = "both",
                norm = "NONE",
                matrix = "observed")

# ## replace 10kb extracted counts from original .rds with 1kb extracted counts
# mcols(noDroso_loops) <- cbind(as.matrix(assay(noDroso_loops)), mcols(noDroso_loops)[-c(1:6)])

## create an mcol for loop size
mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

## create an mcol for loop type
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 ~ "gained",
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

# Calculate loop decay for all loop sets (S/O Manjari!)
source("scripts/utils/mh_index.R")

gainedLoops <- focal
ctcfLoops <- pool

# loop decay for gained loops within sorbitol .hic file ---------------------------------------

# calculate APA matrices for gained loops in sorb .hic file
normalized_gained = data.frame()

for(j in seq_along(gainedLoops)){
  l <- 
    gainedLoops[j] |>
    pixelsToMatrices(buffer=10) |>
    pullHicMatrices(binSize = 10e3,
                    files = merged_hicFiles[2],
                    half = "upper",
                    norm = "NONE",
                    matrix = "observed",
                    onDisk = F) |> 
    aggHicMatrices(FUN=sum) |> 
    as.matrix()
  for(i in 0:10){
    normalized_gained[j,i+1]<- median(mh_index(buffer = 10, loop = l, inner = i))
  }
}
head(normalized_gained)
colnames(normalized_gained) <- as.character(c(0:10))

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero_gained=data.frame()
for(i in seq_along(gainedLoops)){
  for(j in 1:11){
    a <- normalized_gained[i,1]-normalized_gained[i,11]
    b <- normalized_gained[i,j]-normalized_gained[i,11]
    zero_gained[i,j] <- round((b/a),2)
  }
  
}
head(zero_gained)
colnames(zero_gained) <- as.character(c(0:10))


Group_gained <- factor(0:10)
normalized_gained <- na.omit(zero_gained)
saveRDS(normalized_gained, "data/processed/gainedLoops_noDroso_1kb_max100_mh_index.rds")

# loop decay for nullSet loops within control .hic file ---------------------------------------
## MatchedGInteractions to GInteractions
nullSet <- nullSet |>
  data.frame() |>
  mariner::as_ginteractions()

# calculate APA matrices for nullSet within control .hic file
normalized_nullSet = data.frame()
for(j in seq_along(nullSet)){
  l <- 
    nullSet[j] |>
    pixelsToMatrices(buffer=10) |>
    pullHicMatrices(binSize = 10e3,
                    files = merged_hicFiles[1],
                    half = "upper",
                    norm = "NONE",
                    matrix = "observed",
                    onDisk = F) |> 
    aggHicMatrices(FUN=sum) |> 
    as.matrix()
  for(i in 0:10){
    normalized_nullSet[j,i+1]<- median(mh_index(buffer = 10, loop = l, inner = i))
  }
}

head(normalized_nullSet)
colnames(normalized_nullSet) <- as.character(c(0:10))

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero_nullSet=data.frame()
for(i in seq_along(nullSet)){
  for(j in 1:11){
    a <- normalized_nullSet[i,1]-normalized_nullSet[i,11]
    b <- normalized_nullSet[i,j]-normalized_nullSet[i,11]
    zero_nullSet[i,j] <- round((b/a),2)
  }
}
head(zero_nullSet)
colnames(zero_nullSet) <- as.character(c(0:10))


Group_nullSet<- factor(0:10)
normalized_nullSet <- na.omit(zero_nullSet)
saveRDS(normalized_nullSet, "data/processed/nullSet_noDroso_1kb_max100_mh_index.rds")

# loop decay for ctcfLoops within control .hic file -----------------------

# calculate APA matrices for ctcf loops in control .hic file
normalized_ctcf = data.frame()

for(j in seq_along(ctcfLoops)){
  l <-
    ctcfLoops[j] |>
    pixelsToMatrices(buffer=10) |>
    pullHicMatrices(binSize = 10e3,
                    files = merged_hicFiles[1],
                    half = "upper",
                    norm = "NONE",
                    matrix = "observed",
                    onDisk = F) |>
    aggHicMatrices(FUN=sum) |>
    as.matrix()
  for(i in 0:10){
    normalized_ctcf[j,i+1]<- median(mh_index(buffer = 10, loop = l, inner = i))
  }
}
head(normalized_ctcf)
colnames(normalized_ctcf) <- as.character(c(0:10))

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero_ctcf=data.frame()
for(i in seq_along(ctcfLoops)){
  for(j in 1:11){
    a <- normalized_ctcf[i,1]-normalized_ctcf[i,11]
    b <- normalized_ctcf[i,j]-normalized_ctcf[i,11]
    zero_ctcf[i,j] <- round((b/a),2)
  }
  
}
head(zero_ctcf)
colnames(zero_ctcf) <- as.character(c(0:10))


Group_ctcf <- factor(0:10)
normalized_ctcf <- na.omit(zero_ctcf)
saveRDS(normalized_ctcf, "data/processed/ctcfLoops_noDroso_1kb_max100_mh_index.rds")


# load processed data and convert to dataframes ---------------------------

normalized_gained <- readRDS("data/processed/gainedLoops_noDroso_1kb_max100_mh_index.rds")
## convert to dataframe
Group_gained <- factor(0:10)
Mean_gained <- as.vector(colMeans(normalized_gained))
df_gained <- data.frame(Group_gained,Mean_gained)
df_gained

normalized_nullSet <- readRDS("data/processed/nullSet_noDroso_1kb_max100_mh_index.rds")
## convert to dataframe
Group_nullSet <- factor(0:10)
Mean_nullSet <- as.vector(colMeans(normalized_nullSet))
df_nullSet <- data.frame(Group_nullSet,Mean_nullSet)
df_nullSet

normalized_ctcf <- readRDS("data/processed/ctcfLoops_noDroso_1kb_max100_mh_index.rds")
## convert to dataframe
Group_ctcf <- factor(0:10)
Mean_ctcf <- as.vector(colMeans(normalized_ctcf))
df_ctcf <- data.frame(Group_ctcf,Mean_ctcf)
df_ctcf

# final dataframe for plotting --------------------------------------------

df <- cbind(df_gained, df_nullSet) |>
  dplyr::select(Mean_nullSet, Mean_gained)

df <- cbind(df_gained, df_nullSet) |>
  dplyr::select(Mean_nullSet, Mean_gained)

df$Group <- factor(seq(0,10))


# visualization -----------------------------------------------------------
library(ggplot2)

## might want to add 'n' of loops per point

(loop_decay_plot <- ggplot(df,
                           aes(x = Group, group = 1)) +
    # geom_point(aes(y = Mean_gained),
    #            color="black") +
    geom_line(aes(y = Mean_gained),
              color="#0076BA") +
    # geom_point(aes(y = Mean_nullSet),
    #            color= "blue") +
    geom_line(aes(y = Mean_nullSet),
              color= "#4F7942") +
    # geom_point(aes(y = Mean_ctcf),
    #            color= "#4F7942") +
    # geom_line(aes(y = Mean_ctcf),
    #           color= "#4F7942") +
    xlab("Distance from the center (kb)") + 
    ylab("Signal relative to the center") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic() +
    theme(
      axis.text.y=element_text(size=8),
      axis.title.x = element_text(color="steelblue", size=12, face="bold"),
      axis.title.y = element_text(color="steelblue", size=12, face="bold")))

ggsave("plots/loop_decay_noDroso_1kb_max100.pdf", loop_decay_plot, device = "pdf")
