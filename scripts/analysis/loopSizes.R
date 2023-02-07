## Create plots comparing loop sizes of gained loops vs. existing loops

# library
library(ggplot2)
library(mariner)
library(dplyr)
library(hrbrthemes)

## load in YAPP Hi-C differential loops
noDroso_loops <- readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  mariner::as_ginteractions()

## create new metadata columns for loop type & loop size
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "ctcfLoop",
  mcols(noDroso_loops)$padj > 0.1 ~ "ctcfLoop",
  is.character("NA") ~ "ctcfLoop")

mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

noDroso_loops_df <- noDroso_loops |> 
  as.data.frame() |> 
  dplyr::select(loop_type, loop_size)

gainedLoops <- noDroso_loops_df |> 
  filter(loop_type == "gained")

set.seed(420)
ctcfLoops <- noDroso_loops_df |> 
  filter(loop_type != "gained") |> 
  sample_n(nrow(gainedLoops))

noDroso_loops <- rbind(gainedLoops, ctcfLoops)

# Represent it
(p <- noDroso_loops |> 
  ggplot(aes(x=loop_size, fill=loop_type)) +
  geom_density(position = 'identity', alpha = 0.4) +
  scale_fill_manual(values=c("#1D91C0", "#33A02C")) +
  theme_minimal() +
  theme(legend.position="right",
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank()) +
  labs(fill=""))

ggsave("loopSizes.pdf", plot = p, device = "pdf", path = "plots/")
