library(InteractionSet)
library(strawr)
library(tidyverse)
library(nullranges)
library(hictoolsr, exclude = c("as_ginteractions"))
library(mariner)
library(data.table)
library(RColorBrewer)

source("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/scripts/utils/calcLoopEnrichment.R")

## Load in non-D. mel data
noDroso_loops <- readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  mariner::as_ginteractions()

mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.05 & mcols(noDroso_loops)$log2FoldChange > 1 ~ "truegained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

hicFiles <- list.files(c("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb",
                         "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont",
                         "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont",
                         "/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/YAPP/external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl"),
                       full.names = T)

ctcfLoops <- noDroso_loops |> 
  subset(!loop_type %in% c("truegained", "gained"))

gainedLoops <- noDroso_loops |> 
  subset(loop_type %in% c("truegained", "gained"))


# pull HiC matrices to pull pixels of interest from -----------------------


ctcf_cont_mats <- calcLoopEnrichment(ctcfLoops, hicFiles[1])
ctcf_sorb_mats <- calcLoopEnrichment(ctcfLoops, hicFiles[2])

amat_ctcf_cont_mats <- calcLoopEnrichment(ctcfLoops, hicFiles[3])
amat_ctcf_sorb_mats <- calcLoopEnrichment(ctcfLoops, hicFiles[4])

gained_cont_mats <- calcLoopEnrichment(gainedLoops, hicFiles[1])
gained_sorb_mats <- calcLoopEnrichment(gainedLoops, hicFiles[2])

amat_gained_cont_mats <- calcLoopEnrichment(gainedLoops, hicFiles[3])
amat_gained_sorb_mats <- calcLoopEnrichment(gainedLoops, hicFiles[4])


# create datasets ---------------------------------------------------------
flores_ctcf <- tibble(lfc = log2((ctcf_sorb_mats/ctcf_cont_mats)))
flores_ctcf$author <- "Flores et al"

amat_ctcf <- tibble(lfc= log2((amat_ctcf_sorb_mats/amat_ctcf_cont_mats)))
amat_ctcf$author <- "Amat et al"

ctcf_df <- rbind(flores_ctcf, amat_ctcf)

ctcf_df |> 
  group_by(author) %>%
  summarise(
    count = n(),
    median = median(lfc, na.rm = TRUE),
    IQR = IQR(lfc, na.rm = TRUE)
  )

flores_gained <- tibble(lfc = log2((gained_sorb_mats/gained_cont_mats)))
flores_gained$author <- "Flores et al"

amat_gained <- tibble(lfc = log2((amat_gained_sorb_mats/amat_gained_cont_mats)))
amat_gained$author <- "Amat et al"

gained_df <- rbind(flores_gained, amat_gained)

gained_df |> 
  group_by(author) %>%
  summarise(
    count = n(),
    median = median(lfc, na.rm = TRUE),
    IQR = IQR(lfc, na.rm = TRUE)
  )

# data visualization ------------------------------------------------------
pdf("plots/ctcf_log2FC.pdf")

ctcf_df |> 
  ggplot(aes(x = fct_rev(author), y = lfc, fill = fct_rev(author))) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(x = "",
       fill = "") +
  theme_minimal() +
  ylim(c(-3,3)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  theme(legend.position="right",
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(y = "log2FC(treated/untreated)") +
  scale_fill_manual(values = c("#1D91C0", "#C7E9B4"))

## palette of choice
# brewer.pal(n = 9, name = "YlGnBu")
# "#FFFFD9" "#EDF8B1" "#C7E9B4" "#7FCDBB" "#41B6C4" "#1D91C0" "#225EA8" "#253494" "#081D58" "#33A02C"

dev.off()

pdf("plots/gained_log2FC.pdf")

gained_df |> 
  ggplot(aes(x = fct_rev(author), y = lfc, fill = fct_rev(author))) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  labs(x = "",
       fill = "")  + 
  theme_minimal() +
  ylim(c(-3,3)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  theme(legend.position="right",
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(y = "log2FC(treated/untreated)") +
  scale_fill_manual(values = c("#1D91C0", "#C7E9B4"))


dev.off()
