## Obj: Look at which genes have promoters at the ends of gained loop anchors

library(InteractionSet)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RMariaDB)
library(AnnotationHub)
library(GO.db)
library(scales)
library(forcats)

## load in YAPP Hi-C differential loops
noDroso_loops <- readRDS("/Users/jpflores/Phanstiel Lab Dropbox/JP Flores/projects/YAPP/MYAP/data/processed/hic/hg38/diffLoops/noDroso/diffLoops_noDroso_10kb.rds") |> 
  mariner::as_ginteractions()

## create new metadata columns for loop type & loop size
mcols(noDroso_loops)$loop_type <- case_when(
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange > 0 ~ "gained",
  mcols(noDroso_loops)$padj < 0.1 & mcols(noDroso_loops)$log2FoldChange < 0 ~ "lost",
  mcols(noDroso_loops)$padj > 0.1 ~ "static",
  is.character("NA") ~ "other")

mcols(noDroso_loops)$loop_size <- pairdist(noDroso_loops)

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg38.knownGene |>
  keepStandardChromosomes()

seqlevels(noDroso_loops) <- seqlevels(txdb)
seqinfo(noDroso_loops) <- seqinfo(txdb)

## subset gained loops
gainedLoops <- noDroso_loops |> 
  subset(loop_type == "gained")

ctcfLoops <- noDroso_loops |> 
  subset(!loop_type == "gained")

## genes
txdb_UCSC <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene",
                              transcript_ids=NULL,
                              circ_seqs=NULL,
                              url="http://genome.ucsc.edu/cgi-bin/",
                              goldenPath.url=getOption("UCSC.goldenPath.url"),
                              taxonomyId=NA,
                              miRBaseBuild=NA)
genes <- genes(txdb_UCSC)
genes <- data.frame(genes)
genes <- GRanges(seqnames = Rle(genes$seqnames), 
                 ranges = IRanges(start = genes$start, end = genes$end), 
                 strand = genes$strand,
                 UCSC_gene_id = genes$gene_id)

##  overlap promoter regions & gained loop anchors
promoters_gained <- subsetByOverlaps(gainedLoops, promoters(genes))
promoters_ctcf <- subsetByOverlaps(ctcfLoops, promoters(genes))

## get % of gained loop anchors @ promoters
library(ggplot2)

# create a dataset
loop_type <- c("Gained", "Gained", "Lost/Existing", "Lost/Existing")
condition <- c("Gained", "notGained", "Lost/Existing", "notLE")
value <- c((length(promoters_gained)/length(gainedLoops)),
           (length(promoters_ctcf)/length(ctcfLoops)),
           (1 - (length(promoters_gained)/length(gainedLoops))),
           (1 - (length(promoters_ctcf)/length(ctcfLoops))))
data <- data.frame(loop_type, condition, value)

# Stacked + percent
ggplot(data, aes(fill=condition, y=value, x = loop_type)) + 
  geom_bar(position="fill", stat="identity", width = 0.3) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 1, hjust = 1),
        legend.position = "none") +
  labs(y = "",
       x = "") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#1D91C0", "#33A02C", "lightgrey", "lightgrey"))


## get gene symbols
ah <- AnnotationHub()
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
genes_gainedLoops <- select(orgdb, keys = promoters_gained$UCSC_gene_id, 
                            columns = c("SYMBOL", "GENENAME")) |> 
  as.data.frame()

genes_ctcfLoops <- select(orgdb, keys = promoters_ctcf$UCSC_gene_id, 
                            columns = c("SYMBOL", "GENENAME")) |> 
  as.data.frame()


nrow(genes_gainedLoops)

## get GO terms
go_gained <- select(orgdb, keys = genes_gainedLoops$SYMBOL, "GO", "SYMBOL")
go2_gained <- select(GO.db, genes_gainedLoops$GO, c("TERM","DEFINITION"), "GOID") 
go2_gained$TERM

go_ctcf <- select(orgdb, keys = genes_ctcfLoops$SYMBOL, "GO", "SYMBOL")
go2_ctcf <- select(GO.db, genes_ctcfLoops$GO, c("TERM","DEFINITION"), "GOID") 
go2_ctcf$TERM

write_csv(genes_ctcfLoops, "tables/genes_atGainedLoopsAnchors.csv")
write_csv(genes_gainedLoops, "tables/genes_atGainedLoopsAnchors.csv")

# ## do any of these genes have multiple gained loops at their promoters?
# loi <- promoters |> 
#   subset(UCSC_gene_id %in% c("7004", "83937", "8463", "163732", "79142"))
# 
# countOverlaps(loi, gainedLoops)
