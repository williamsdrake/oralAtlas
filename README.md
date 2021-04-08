Oral Atlas
================

# Introduction

The following document includes the code used to analyze datasets and
generate figures related to single-cell RNA sequencing from Williams et
al., 2021 (doi: [doi number here](actual.address.here.org))

## Setup

> This section covers dataset setup in R and **Figure S1A-C**

### Data aquisition

The following datasets were used in the manuscript:

  - Oral (Williams et al., 2021):
    [GSE164241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164241)  
  - Skin (He et al., 2020):
    [GSE147424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147424)  
  - Lung (Reyfman et al., 2019):
    [GSE122960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122960)  
  - Ileum (Martin et al., 2019):
    [GSE134809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134809)  
  - Salivary Gland (Huang et al., 2021):
    [covid19cellatlas](https://www.covid19cellatlas.org/index.healthy.html)

### Libraries

``` r
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)
library(reshape2)
library(scales)
library(gdata)
library(ggraph)
library(scater)
library(cowplot)
library(SingleR)
library(data.table)
library(paletteer)
library(R.utils)
library(ComplexHeatmap)
library(clustree)
library(viridis)
library(pheatmap)
library(devtools)
library(Scillus)
library(scFunctions)
library(ggcharts)
library(schex)
library(Cairo)
library(cartography)
library(kableExtra)
library(knitr)
library(patchwork)
library(gsfisher)
library(BiocParallel)
library(nichenetr)
library(easyalluvial)
library(purrr)
```

### Data import

Deposited data was downloaded and converted to Seurat objects (one
example per dataset shown for brevity):

1)  Oral/Ileum:

<!-- end list -->

``` r
sample1.data <- Read10X(data.dir = "/path/to/data/folder/")
sample1 <- CreateSeuratObject(counts = sample1.data, project = "sample1")
```

2)  Skin:

<!-- end list -->

``` r
sample1count <- read.table(file="/path/to/data/folder/sample1.txt",sep=",",header = TRUE, row.names = NULL)
names <- make.unique(as.character(sample1count$X))
rownames(sample1count) <- names
sample1count <- sample1count[,-1] # get rid of old names
sample1 <- CreateSeuratObject(counts = sample1count, min.cells = 3, min.features = 200)
```

3)  Lung:

<!-- end list -->

``` r
sample1.data <- Read10X_h5("/path/to/data/folder/sample1.h5")
sample1 <- CreateSeuratObject(counts = sample1.data, project = "sample1")
```

4)  Salivary gland:

<!-- end list -->

``` r
library(SeuratDisk)
library("reticulate")

anndata<-import ("anndata")
adata = anndata$read("warner20_salivary_gland.processed.h5ad")
adata$write("warner20_salivary_gland.processed.gzip.h5ad", compression="gzip")
Convert("warner20_salivary_gland.processed.gzip.h5ad", dest = "h5seurat", overwrite = TRUE)
SG <- LoadH5Seurat("warner20_salivary_gland.processed.gzip.h5seurat")
```

Individual Seurat objects were merged according to tissue type/disease
status (one example is shown for all tissue types). *The data deposited
at covid19atlas.org was pre-processed and pre-merged, and therefore did
not require this step.*

``` r
GM <- merge(GM136, 
            y = c(GM143, GM144, GM147, GM148, GM183, GM184a, GM238, GM241, GM242, GM169, GM283, GM289), 
            add.cell.ids = c("GM136", "GM143", "GM144", "GM147", "GM148", "GM183", "GM184a", 
                             "GM238", "GM241", "GM242", "GM169", "GM283", "GM289"), 
            project = "GM")
```

In subsequent sections, the following variable names were used for each
merged dataset:

  - `GM` - Healthy gingival mucosa
  - `BM` - Healthy buccal mucosa
  - `PD` - Disease (periodontitis) gingival mucosa
  - `HO` - Healthy buccal and gingival mucosa
  - `HL` - Healthy lung
  - `HI` - Healthy ileum
  - `HS` - Healthy skin
  - `SG` - Healthy salivary gland

### Data quality control

Datasets underwent quality control as follows (only oral datasets shown
for brevity):

``` r
# save pre-qc data
preQC_GM <- as.data.frame(table(GM$orig.ident))
preQC_BM <- as.data.frame(table(BM$orig.ident))
preQC_PD <- as.data.frame(table(PD$orig.ident))

# generate mitochrondial gene metadata 
GM[["percent.mt"]] <- PercentageFeatureSet(GM, 
                                           pattern = "^MT-")
BM[["percent.mt"]] <- PercentageFeatureSet(BM, 
                                           pattern = "^MT-")
PD[["percent.mt"]] <- PercentageFeatureSet(PD, 
                                           pattern = "^MT-")

# generate pre-qc violin plots
preQC.Vln_GM <- VlnPlot(GM, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        ncol = 3, 
                        group.by = "orig.ident", 
                        pt.size = 0)
preQC.Vln_BM <- VlnPlot(BM, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        ncol = 3, 
                        group.by = "orig.ident", 
                        pt.size = 0)
preQC.Vln_PD <- VlnPlot(PD, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        ncol = 3,
                        group.by = "orig.ident", 
                        pt.size = 0)

# filter low quality and potential doublet cells
GM <- subset(GM, 
             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
BM <- subset(BM, 
             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
PD <- subset(PD, 
             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

# save post-qc data
postQC_GM <- as.data.frame(table(GM$orig.ident))
postQC_BM <- as.data.frame(table(BM$orig.ident))
postQC_PD <- as.data.frame(table(PD$orig.ident))

# generate post-qc violin plots
postQC.Vln_GM <- VlnPlot(GM, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                         ncol = 3, 
                         group.by = "orig.ident", 
                         pt.size = 0)
postQC.Vln_BM <- VlnPlot(BM, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                         ncol = 3, 
                         group.by = "orig.ident", 
                         pt.size = 0)
postQC.Vln_PD <- VlnPlot(PD, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                         ncol = 3, 
                         group.by = "orig.ident", 
                         pt.size = 0)

# format pre/post-qc data for graphing
qcSum_GM <- merge(preQC_GM, 
                  postQC_GM, 
                  by = "Var1")
names(qcSum_GM)[1] <- "Patient"
names(qcSum_GM)[2] <- "preQC"
names(qcSum_GM)[3] <- "postQC"
qcSum_GM <- reshape2::melt(qcSum_GM[,c('Patient','preQC','postQC')],
                           id.vars = 1)
names(qcSum_GM)[3] <- "Cells"
qcSum_BM <- merge(preQC_BM, 
                  postQC_BM, 
                  by = "Var1")
names(qcSum_BM)[1] <- "Patient"
names(qcSum_BM)[2] <- "preQC"
names(qcSum_BM)[3] <- "postQC"
qcSum_BM <- reshape2::melt(qcSum_BM[,c('Patient','preQC','postQC')],
                           id.vars = 1)
names(qcSum_BM)[3] <- "Cells"
qcSum_PD <- merge(preQC_PD, 
                  postQC_PD, 
                  by = "Var1")
names(qcSum_PD)[1] <- "Patient"
names(qcSum_PD)[2] <- "preQC"
names(qcSum_PD)[3] <- "postQC"
qcSum_PD <- reshape2::melt(qcSum_PD[,c('Patient','preQC','postQC')],
                           id.vars = 1)
names(qcSum_PD)[3] <- "Cells"

# plot cell numbers pre/post-qc
cellNumSummary_GM <- ggplot(qcSum_GM, 
                            aes(x = Patient, 
                                y = Cells)) + 
  geom_bar(aes(fill = variable), 
           stat = "identity", 
           position = "dodge") +
  theme_cowplot() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  scale_fill_manual(values = c("grey", "black"))

cellNumSummary_BM <- ggplot(qcSum_BM, 
                            aes(x = Patient, 
                                y = Cells)) + 
  geom_bar(aes(fill = variable), 
           stat = "identity", 
           position = "dodge") +
  theme_cowplot() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  scale_fill_manual(values = c("grey", "black"))

cellNumSummary_PD <- ggplot(qcSum_PD, 
                            aes(x = Patient, 
                                y = Cells)) + 
  geom_bar(aes(fill = variable), 
           stat = "identity", 
           position = "dodge") +
  theme_cowplot() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  scale_fill_manual(values = c("grey", "black"))
```

#### Figure S1A-C

``` r
library(patchwork)

FigS1A <- (preQC.Vln_GM / preQC.Vln_BM / preQC.Vln_PD) 
FigS1B <- (postQC.Vln_GM / postQC.Vln_BM / postQC.Vln_PD) 
FigS1C <- cellNumSummary_GM / (cellNumSummary_BM | cellNumSummary_PD) + 
  plot_layout(guides = "collect") 
```

## Gingival vs. Buccal

> This section covers **Figure 1B,D, Figure S1D, Figure 2, Figure S2,
> Figure 3A-D, Figure 7B Table S5**

### Data preparation

#### Transformation

``` r
BM <- NormalizeData(BM)
BM <- FindVariableFeatures(BM, nfeatures = 4000)
GM <- NormalizeData(GM)
GM <- FindVariableFeatures(GM, nfeatures = 4000)
```

#### Integration and batch correction

``` r
sampleList <- list(GM, BM)
oral.anchors <- FindIntegrationAnchors(object.list = sampleList, dims = 1:50)
oralIntegrated <- IntegrateData(anchorset = oral.anchors, dims = 1:50)
```

#### Dimensionality reduction and cell clustering

``` r
# cell-cycle scoring and regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
oralIntegrated <- CellCycleScoring(oralIntegrated, 
                                   s.features = s.genes, 
                                   g2m.features = g2m.genes,
                                   set.ident = TRUE)
# standard workflow for clustering
oralIntegrated <- ScaleData(oralIntegrated, 
                            vars.to.regress = c("S.Score", "G2M.Score"), 
                            verbose = FALSE)
oralIntegrated <- RunPCA(oralIntegrated, 
                         npcs = 50, 
                         verbose = FALSE)
oralIntegrated <- FindNeighbors(oralIntegrated, 
                                reduction = "pca", 
                                dims = 1:50)
oralIntegrated <- FindClusters(oralIntegrated, 
                               resolution = seq(from = 0.1, 
                                                to = 1.0, 
                                                by = 0.1))
oralIntegrated <- RunUMAP(oralIntegrated, 
                          reduction = "pca", 
                          dims = 1:50)
# assess cluster tree
clustree(oralIntegrated)
```

### High-level summary

#### Cell classification

``` r
# set correct assay and idents for GEX
DefaultAssay(oralIntegrated) <- "RNA"
Idents(oralIntegrated) <- "integrated_snn_res.1"
all.markers.1 <- FindAllMarkers(oralIntegrated, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.75)
significant.markers.1 <- all.markers.1[all.markers.1$p_val_adj < 0.2, ]

# remove RBC contaminated clusters and clusters with multiple genotypes (likely doublets)
oralIntegrated <- subset(oralIntegrated, 
                         idents = c("21", "22", "26", "31", "33"), 
                         invert = TRUE)
oralIntegrated <- RenameIdents(oralIntegrated, 
                               "0"="Endothelial", 
                               "1"="Immune", 
                               "2"="Endothelial", 
                               "3"="Immune", 
                               "4"="Immune", 
                               "5"="Endothelial", 
                               "6"="Endothelial", 
                               "7"="Fibroblast", 
                               "8"="Fibroblast", 
                               "9"="Endothelial", 
                               "10"="Fibroblast", 
                               "11"="Fibroblast", 
                               "12"="Endothelial", 
                               "13"="Epithelial", 
                               "14"="Fibroblast", 
                               "15"="Epithelial", 
                               "16"="Endothelial", 
                               "17"="Endothelial", 
                               "18"="Immune", 
                               "19"="Immune", 
                               "20"="Immune", 
                               "23"="Immune", 
                               "24"="Epithelial", 
                               "25"="Epithelial", 
                               "27"="Endothelial", 
                               "28"="Endothelial", 
                               "29"="Other", 
                               "30"="Immune", 
                               "32" = "Immune")

# save general cell type names to metadata and relevel based on order in figures
levels(oralIntegrated) <- c("Endothelial", 
                            "Fibroblast", 
                            "Immune", 
                            "Epithelial", 
                            "Other")
oralIntegrated$generalCellTypes <- oralIntegrated@active.ident
```

##### Figure 1B

``` r
pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5)]
Idents(oralIntegrated) <- 'generalCellTypes'
Fig1B_dimPlot <- DimPlot(oralIntegrated, 
                   reduction = "umap", 
                   cols = pal) + 
             theme_void()
Fig1B_prop <- plot_stat(oralIntegrated, 
                        plot_type = "prop_fill", 
                        group_by = "project") + 
  scale_fill_manual(values = pal) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 4)) +
  guides(fill = FALSE)

# proportion of major cell types by sample for stats
Fig1B_propBySample <- plot_stat(oralIntegrated, 
                                       plot_type = "prop_fill", 
                                       group_by = "orig.ident") + 
  scale_fill_manual(values = pal) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 4), 
        axis.text.x = element_text(angle = 90)) +
  guides(fill = FALSE)
Fig1BStats <- as.data.frame(Fig1B_propBySample[["data"]])
Fig1BStats <- Fig1BStats[order(Fig1BStats$cluster),]
write.csv(Fig1BStats, "Fig1B_forStatistics.csv")
```

##### Figure S1D

``` r
# find cell-type defining markers
generalCellTypeMarkers <- FindAllMarkers(oralIntegrated, 
                                         only.pos = TRUE, 
                                         min.pct = 0.25, 
                                         logfc.threshold = 0.5)
generalSignificant <- generalCellTypeMarkers[generalCellTypeMarkers$p_val_adj < 0.2, ]
generalTop20 <- generalSignificant %>% 
  group_by(cluster) %>% 
  top_n(n=20, wt=avg_logFC)
FigS1D <- plot_heatmap(dataset = oralIntegrated,
                       markers = generalTop20,
                       sort_var = c("generalCellTypes", "orig.ident"), 
                       anno_var = c("generalCellTypes"), 
                       anno_colors = c("Reds"), 
                       n = 20)
```

#### Gene expression

``` r
# isolate different cell types
Idents(oralIntegrated) <- "generalCellTypes"
endo <- subset(oralIntegrated, 
               idents = c("Endothelial"))
epi <- subset(oralIntegrated, 
              idents = c("Epithelial"))
fib <- subset(oralIntegrated, 
              idents = c("Fibroblast"))
im <- subset(oralIntegrated, 
             idents = c("Immune"))

# selected genes to display
curatedGenesEndo <- c("ACKR1", 
                      "RAMP2", 
                      "SELE", 
                      "VWF", 
                      "PECAM1")
curatedGenesFib <- c("LUM", 
                     "COL3A1", 
                     "DCN", 
                     "COL1A1", 
                     "CFD")
curatedGenesIm <- c("CD69", 
                    "CD52", 
                    "CXCR4", 
                    "PTPRC", 
                    "HCST")
curatedGenesEpi <- c("KRT14", 
                     "KRT5", 
                     "S100A2", 
                     "CSTA", 
                     "SPRR1B")
```

##### Figure 1D

``` r
# set ident (BM, GM)
Idents(endo) <- "project"
Idents(im) <- "project"
Idents(fib) <- "project"
Idents(epi) <- "project"

# violin plots
bigPlot1 <- VlnPlot(endo, 
                    features = curatedGenesEndo, 
                    ncol = 5, 
                    pt.size = 0, 
                    assay = "RNA", 
                    combine = FALSE, 
                    cols = c("black", "grey"))
bigPlot2 <- VlnPlot(fib, 
                    features = curatedGenesFib, 
                    ncol = 5, 
                    pt.size = 0, 
                    assay = "RNA", 
                    combine = FALSE, 
                    cols = c("black", "grey"))
bigPlot3 <- VlnPlot(im, 
                    features = curatedGenesIm, 
                    ncol = 5, 
                    pt.size = 0,
                    assay = "RNA", 
                    combine = FALSE, 
                    cols = c("black", "grey"))
bigPlot4 <- VlnPlot(epi, 
                    features = curatedGenesEpi, 
                    ncol = 5, 
                    pt.size = 0, 
                    assay = "RNA", 
                    combine = FALSE, 
                    cols = c("black", "grey"))
bigPlot <- c(bigPlot1, bigPlot2, bigPlot3, bigPlot4)

# format graphs
for(i in 1:length(bigPlot)) {
  bigPlot[[i]] <- bigPlot[[i]] + 
    theme(axis.text.x = element_text(size = 7, angle = 0), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=7, angle = 0)) + 
    guides(fill = FALSE) +
    coord_flip() 
}

Fig1D <- cowplot::plot_grid(plotlist = bigPlot, ncol = 5)
```

### Endothelial subset

#### Cell classification

``` r
DefaultAssay(endo) <- "RNA"
Idents(endo) <- "integrated_snn_res.1"

# endo GEX
endoMarkers <- FindAllMarkers(endo, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.75)
endoSig <- endoMarkers[endoMarkers$p_val_adj < 0.2, ]

# rename
endo <- RenameIdents(endo, 
                     "0" = "H.VEC 1.1", 
                     "2" = "H.VEC 1.2", 
                     "5" = "H.VEC 1.3", 
                     "6" = "H.VEC 1.4",
                    "9" = "H.SMC", 
                    "12" = "H.VEC 1.1", 
                    "16" = "H.LEC", 
                    "17" = "H.SMC", 
                    "27" = "H.SMC",
                    "28" = "H.VEC 1.1")

# relevel
levels(endo) <- c("H.VEC 1.1", 
                  "H.VEC 1.2", 
                  "H.VEC 1.3", 
                  "H.VEC 1.4",
                  "H.SMC", 
                  "H.LEC")

# endo metadata
endo$endoTypes <- Idents(endo)
```

##### Figure 2A

``` r
endo <- RunUMAP(endo, dims = 1:50)
endoColors <- c(carto.pal(pal1 = "red.pal", 
                          n1 = 20)[c(20,8,2,14)],
                carto.pal(pal1 = "wine.pal", 
                          n1 = 10)[c(3,7)])
Idents(endo) <- "endoTypes"
Fig2A_dimPlot <- DimPlot(object = endo, 
                    label = FALSE, 
                    reduction = "umap", 
                    cols = endoColors) + 
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
Fig2A_prop <- plot_stat(endo, 
                        plot_type = "prop_fill", 
                        group_by = "project") +
  scale_fill_manual(values = endoColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 4)) +
  guides(fill = FALSE)

# stats for 2A 
Fig2A_propBySample <- plot_stat(endo, 
                                plot_type = "prop_fill", 
                                group_by = "orig.ident") 
Fig2AStats <- as.data.frame(Fig2A_propBySample[["data"]])
Fig2AStats <- Fig2AStats[order(Fig2AStats$cluster),]
write.csv(Fig2AStats, "Fig2A_forStats.csv")

# prop plots for cartoon
endoPropBySubset <- Fig2A_prop[["data"]]
Fig2B_cartoonProp <- ggplot(endoPropBySubset, 
                           aes(fill=group, 
                           y = freq,
                           x = cluster)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
```

#### Gene expression

``` r
# endo GEX based on new cluster names
Idents(endo) <- "endoTypes"
endoSig <- endoMarkers[endoMarkers$p_val_adj < 0.2, ]
endoTop <- endoSig %>% 
  group_by(cluster) %>% 
  top_n(n=3, wt=avg_logFC)
endoTop <- endoTop[!duplicated(endoTop$gene),] %>%
  map_df(rev)
write.csv(endoSig, "endo_clusterMarkers_significant.csv")
```

##### Figure 2B

``` r
Fig2B <- DotPlot(endo, 
                    features = endoTop$gene, 
                    col.min = 0, 
                    scale.by = "size", 
                    dot.min = 0.01, 
                    dot.scale = 3) +
  coord_flip() + 
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5)) +
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp", 
                             order = 2))
```

#### Pathway analysis

``` r
# helper function
endo.expressed <- getExpressedGenesFromSeuratObject(endo, 
                                                    unique(endo@active.ident), 
                                                    min.pct = 0.25)
annotation <- fetchAnnotation(species = "hs")

endoMarkers.gsf <- FindAllMarkers(endo, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)
endoMarkers.gsf$entrezID <- as.character(annotation$entrez_id[match(endoMarkers.gsf$gene, 
                                                                    annotation$gene_name)])
endoMarkers.gsf <- endoMarkers.gsf[!is.na(endoMarkers.gsf$entrezID),]
background_entrez <- as.character(annotation$entrez_id[match(endo.expressed, 
                                                             annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]

endoMarkers.gsf.filtered <- endoMarkers.gsf[endoMarkers.gsf$p_val_adj < 0.05,]
endoGO <- runGO.all(results=endoMarkers.gsf.filtered,
                    species = "hs",
                    background_ids = background_entrez,
                    gene_id_col="entrezID",
                    gene_id_type="entrez",
                    sample_col="cluster",
                    p_col="p_val_adj",
                    p_threshold=0.05)

endoGO.filtered <- endoGO[endoGO$ontology=="BP",]
endoGO.filtered <- filterGenesets(endoGO.filtered,
                                  min_foreground_genes = 2,
                                  max_genes_geneset = 500,
                                  min_odds_ratio = 2,
                                  p_col = "p.val",
                                  padjust_method = "BH",
                                  use_adjusted_pvalues = FALSE,
                                  pvalue_threshold = 0.05)
# pathways to highlight
endopathways <- c("antigen processing and presentation",
                  "endothelium development", 
                  "response to lipopolysaccharide",
                  "cellular response to type I interferon", 
                  "muscle contraction",
                  "lymphatic endothelial cell fate commitment")

# select genes associated w/ pathways of interest
endoGO.genes <- endoGO.filtered[endoGO.filtered$description %in% endopathways,11]
endoGO.genes2 <- unlist(strsplit(endoGO.genes, ","))
endoGO.genes2 <- endoGO.genes2[!duplicated(endoGO.genes2)]

# new object w/ average expression for all the genes
endoAvgExpression <- AverageExpression(object = endo, 
                                       return.seurat = TRUE)
```

##### Figure S2A

``` r
FigureS2A <- DoHeatmap(endoAvgExpression, 
                       features = endoGO.genes2, 
                       disp.min = -2,
                       group.bar = FALSE, 
                       label = FALSE, 
                       raster = FALSE) + 
  scale_fill_gradient2(low = "steelblue4", 
                       mid = "black", 
                       high = "cornsilk",
                       breaks=c(-2.5,0,2.5), 
                       labels=c(-2.5,0,2.5), 
                       limits = c(-2.5,2.5)) + 
  theme(axis.text.y = element_text(size = 6))
```

### Epithelial subset

#### Cell classification

``` r
Idents(epi) <- "integrated_snn_res.1"

# find known clusters (cornified) are combined into one based on the global clustering @ resolution 1
# recluster to see if they can be resolved
epi2 <- epi
DefaultAssay(epi2) <- "integrated"
epi2 <- FindNeighbors(epi2, 
                      dims = 1:50)
epi2 <- FindClusters(epi2, 
                     resolution = seq(from = 0.1, 
                                      to = 1.0, 
                                      by = 0.1))
Idents(epi2) <- "integrated_snn_res.0.1"
epi2Markers <- FindAllMarkers(epi2, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.75)
epi2Sig <- epi2Markers[epi2Markers$p_val_adj < 0.2, ]
epi2Top <- epi2Sig %>%
  group_by(cluster) %>%
  top_n(n=5, wt=avg_logFC)
print(epi2Top, n=40)
epi <- epi2
epi2 <- NULL

# rename
epi <- RenameIdents(epi, "0"="H.Epi 1.1", "1"="H.Epi 1.1", "2"="H.Epi 1.2", "3"="H.Epi 1.3", "4"="H.Mel")

# relevel
levels(epi) <- c("H.Epi 1.1", "H.Epi 1.2", "H.Epi 1.3", "H.Epi 1.4", "H.Mel")

# save metadata
epi$epiTypes <- Idents(epi)
```

##### Figure 2C

``` r
epi <- RunUMAP(epi, dims = 1:50)
epiColors <- carto.pal(pal1 = "brown.pal", 
                       n1 = 20)[c(20, 15, 10, 5)]

Idents(epi) <- "epiTypes"
Fig2C <- DimPlot(object = epi, 
                    label = FALSE, 
                    reduction = "umap", 
                    cols = epiColors) + 
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
Fig2C_prop <- plot_stat(epi, 
                               plot_type = "prop_fill", 
                               group_by = "project") +
  scale_fill_manual(values = epiColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 4)) +
  guides(fill = FALSE)

# stats for 2C 
Fig2C_propBySample <- plot_stat(epi, 
                                plot_type = "prop_fill", 
                                group_by = "orig.ident")
Fig2CStats <- as.data.frame(Fig2C_propBySample[["data"]])
Fig2CStats <- Fig2CStats[order(Fig2CStats$cluster),]
write.csv(Fig2CStats, "Fig2C_forStats.csv")

# prop plots for cartoon
endoPropBySubset <- Fig2C_prop[["data"]]
Fig2D_cartoonProp <- ggplot(epiPropBySubset, 
       aes(fill=group, 
           y = freq, 
           x = cluster)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
```

#### Gene expression

``` r
# epi GEX based on new cluster names
epiMarkers <- FindAllMarkers(epi, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.75)
epiSig <- epiMarkers[epiMarkers$p_val_adj < 0.2, ]
write.csv(epiSig, "epi_clusterMarkers_significant_R_1.csv")
epiTop <- epiSig %>%
  group_by(cluster) %>%
  top_n(n=5, wt=avg_logFC)
epiTop <- epiTop[!duplicated(epiTop$gene),] 

customEpi <- rev(c("KRT14", "KRT5", "SFN", "KRT15", "DSP",
                   "FDCSP", "IL36G", "MMP12", "ODAM", "SLPI",
                   "CRNN", "CNFN", "TGM3", "SPRR3", "SBSN",
                   "MITF", "DCT", "PMEL", "TYR", "MLANA"))
```

##### Figure 2D

``` r
Fig2D <- DotPlot(epi,
                 features = customEpi, 
                 col.min = 0, 
                 scale.by = "size", 
                 dot.min = 0.01, 
                 dot.scale = 3) +
  coord_flip() + 
  scale_color_gradient2(low = "grey90", high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp", 
                             order = 2))
```

#### Pathway analysis

``` r
# helper function
Idents(epi) <- "epiTypes"
epi.expressed <- getExpressedGenesFromSeuratObject(epi, 
                                                   unique(epi@active.ident), 
                                                   min.pct = 0.25)
annotation <- fetchAnnotation(species = "hs")

epiMarkers.gsf <- FindAllMarkers(epi, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)
epiMarkers.gsf$entrezID <- as.character(annotation$entrez_id[match(epiMarkers.gsf$gene, 
                                                                   annotation$gene_name)])
epiMarkers.gsf <- epiMarkers.gsf[!is.na(epiMarkers.gsf$entrezID),]
background_entrez <- as.character(annotation$entrez_id[match(epi.expressed, 
                                                             annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]

epiMarkers.gsf.filtered <- epiMarkers.gsf[epiMarkers.gsf$p_val_adj < 0.05,]
epiGO <- runGO.all(results=epiMarkers.gsf.filtered,
                   species = "hs",
                   background_ids = background_entrez,
                   gene_id_col="entrezID",
                   gene_id_type="entrez",
                   sample_col="cluster",
                   p_col="p_val_adj",
                   p_threshold=0.05)

epiGO.filtered <- epiGO[epiGO$ontology=="BP",]
epiGO.filtered <- filterGenesets(epiGO.filtered,
                                 min_foreground_genes = 3,
                                 max_genes_geneset = 500,
                                 min_odds_ratio = 2,
                                 p_col = "p.val",
                                 padjust_method = "BH",
                                 use_adjusted_pvalues = TRUE,
                                 pvalue_threshold = 0.05)
# pathways to highlight
epipathways <- c("keratinocyte differentiation",
                 "leukocyte chemotaxis",
                 "cornification",
                 "pigmentation")

# select genes associated w/ pathways of interest
epiGO.genes <- epiGO.filtered[epiGO.filtered$description %in% epipathways,11]
epiGO.genes2 <- unlist(strsplit(epiGO.genes, ","))
epiGO.genes2 <- epiGO.genes2[!duplicated(epiGO.genes2)]

# new object w/ average expression for all the genes
epiAvgExpression <- AverageExpression(object = epi, return.seurat = TRUE)
```

##### Figure S2B

``` r
FigureS2B <- DoHeatmap(epiAvgExpression, 
                       features = epiGO.genes2, 
                       disp.min = -2,
                       group.bar = FALSE, 
                       label = FALSE, 
                       raster = FALSE) + 
  scale_fill_gradient2(low = "steelblue4", 
                       mid = "black", 
                       high = "cornsilk",
                       breaks=c(-2.5,0,2.5), 
                       labels=c(-2.5,0,2.5), 
                       limits = c(-2.5,2.5)) + 
  theme(axis.text.y = element_text(size = 6))
```

### Fibroblast subset

#### Cell classification

``` r
DefaultAssay(fib) <- "RNA"
Idents(fib) <- "integrated_snn_res.1"

# fib GEX
fibMarkers <- FindAllMarkers(fib, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.75)
fibSig <- fibMarkers[fibMarkers$p_val_adj < 0.2, ]

# rename
fib <- RenameIdents(fib, 
                    "7" = "H.Fib 1.1", 
                    "8" = "H.Fib 1.2", 
                    "10" = "H.Fib 1.3", 
                    "11" = "H.Fib 1.4",
                    "14" = "H.Fib 1.5")

# relevel
levels(fib) <- c("H.Fib 1.1", 
                 "H.Fib 1.2", 
                 "H.Fib 1.3", 
                 "H.Fib 1.4", 
                 "H.Fib 1.5")

# fib metadata
fib$fibTypes <- Idents(fib)
```

##### Figure 2E

``` r
fib <- RunUMAP(fib, dims = 1:50)
fibColors <- c(paletteer_d("ggsci::green_material")[c(10, 7, 4, 2)], 
               paletteer_d("ggsci::lime_material")[c(6)])
Idents(fib) <- "fibTypes"
Fig2E_dimPlot <- DimPlot(object = fib, 
                         label = FALSE, 
                         reduction = "umap", 
                         cols = fibColors) + 
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
Fig2E_prop <- plot_stat(fib, 
                        plot_type = "prop_fill", 
                        group_by = "project") +
  scale_fill_manual(values = fibColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 4)) +
  guides(fill = FALSE)

# stats for 2E 
Fig2E_propBySample <- plot_stat(fib, 
                                plot_type = "prop_fill", 
                                group_by = "orig.ident") 
Fig2EStats <- as.data.frame(Fig2E_propBySample[["data"]])
Fig2EStats <- Fig2EStats[order(Fig2EStats$cluster),]
write.csv(Fig2EStats, "Fig2E_forStats.csv")

# prop plots for cartoon
fibPropBySubset <- Fig2E_prop[["data"]]
Fig2F_cartoonProp <- ggplot(endoPropBySubset, 
                           aes(fill=group, 
                           y = freq,
                           x = cluster)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
```

#### Gene expression

``` r
# fib GEX based on new cluster names
fibMarkers <- FindAllMarkers(fib, 
                             only.pos = TRUE, 
                             min.pct = 0.5, 
                             logfc.threshold = 0.25)
fibSig <- fibMarkers[fibMarkers$p_val_adj < 0.2, ]
write.csv(fibSig, "fib_clusterMarkers_significant.csv")
fibTop <- fibSig %>%
  group_by(cluster) %>%
  top_n(n=3, wt=avg_logFC)
fibTop <- fibTop[!duplicated(fibTop$gene),] %>%
  map_df(rev)
```

##### Figure 2F

``` r
Fig2F <- DotPlot(fib, 
                 features = fibTopGenes,
                 col.min = 0, 
                 scale.by = "size", 
                 dot.min = 0.01, 
                 dot.scale = 3) +
  coord_flip() + 
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp", 
                             order = 2))
```

#### Pathway analysis

``` r
# helper function
Idents(fib) <- "fibTypes"
fib.expressed <- getExpressedGenesFromSeuratObject(fib, 
                                                   unique(fib@active.ident), 
                                                   min.pct = 0.25)
annotation <- fetchAnnotation(species = "hs")
fibMarkers.gsf <- FindAllMarkers(fib, 
                                 only.pos = TRUE, 
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)
fibMarkers.gsf$entrezID <- as.character(annotation$entrez_id[match(fibMarkers.gsf$gene, 
                                                                   annotation$gene_name)])
fibMarkers.gsf <- fibMarkers.gsf[!is.na(fibMarkers.gsf$entrezID),]
background_entrez <- as.character(annotation$entrez_id[match(fib.expressed, 
                                                             annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]

fibMarkers.gsf.filtered <- fibMarkers.gsf[fibMarkers.gsf$p_val_adj < 0.05,]
fibGO <- runGO.all(results=fibMarkers.gsf.filtered,
                   species = "hs",
                   background_ids = background_entrez,
                   gene_id_col="entrezID",
                   gene_id_type="entrez",
                   sample_col="cluster",
                   p_col="p_val_adj",
                   p_threshold=0.05)

fibGO.filtered <- fibGO[fibGO$ontology=="BP",]
fibGO.filtered <- filterGenesets(fibGO.filtered,
                                 min_foreground_genes = 3,
                                 max_genes_geneset = 500,
                                 min_odds_ratio = 2,
                                 p_col = "p.val",
                                 padjust_method = "BH",
                                 use_adjusted_pvalues = TRUE,
                                 pvalue_threshold = 0.05)

# pathways to highlight
fibpathways <- c("translational initiation",
                 "tissue remodeling", 
                 "regulation of leukocyte proliferation",
                 "granulocyte migration", 
                 "complement activation")

# select genes associated w/ pathways of interest
fibGO.genes <- fibGO.filtered[fibGO.filtered$description %in% fibpathways,11]
fibGO.genes2 <- unlist(strsplit(fibGO.genes, ","))
fibGO.genes2 <- fibGO.genes2[!duplicated(fibGO.genes2)]

# new object w/ average expression for all the genes
fibAvgExpression <- AverageExpression(object = fib, return.seurat = TRUE)
```

##### Figure S2C

``` r
FigureS2C <- DoHeatmap(fibAvgExpression, 
                       features = fibGO.genes2, 
                       disp.min = -2,
                       group.bar = FALSE, 
                       label = FALSE, 
                       raster = FALSE) + 
  scale_fill_gradient2(low = "steelblue4", 
                       mid = "black", 
                       high = "cornsilk",
                       breaks=c(-2.5,0,2.5), 
                       labels=c(-2.5,0,2.5), 
                       limits = c(-2.5,2.5)) + 
  theme(axis.text.y = element_text(size = 6))
```

##### Table S4

``` r
endoMarkers.top <- endoMarkers.gsf %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_logFC)
epiMarkers.top <- epiMarkers.gsf %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_logFC)
fibMarkers.top <- fibMarkers.gsf %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_logFC)

allMarkers.top <- rbind(endoMarkers.top, epiMarkers.top, fibMarkers.top)
allMarkers.top$cluster <- paste("H", allMarkers.top$cluster, sep=".")

write.csv(allMarkers.top, "TableS4.csv")
```

### Immune subset

#### Cell classification

``` r
# use SingleR
# seurat object -> SingleCellExperiment
immuneSCE <- as.SingleCellExperiment(im)

# reference database
mon.se <- MonacoImmuneData()
blueprint <- BlueprintEncodeData()
commonBlue <- intersect(rownames(im), rownames(blueprint))
blueprint <- blueprint[commonBlue,]
immuneSCE <- immuneSCE[commonBlue,]
immuneSCE <- logNormCounts(immuneSCE)
pred.mon.fine <- SingleR(test=immuneSCE, 
                         ref=mon.se, 
                         labels=mon.se$label.fine,
                         de.method="classic", 
                         fine.tune = TRUE, 
                         assay.type.test = 1, 
                         BPPARAM=MulticoreParam(8))
pred.mon.course <- SingleR(test=immuneSCE, 
                           ref=mon.se, 
                           labels=mon.se$label.main,
                           de.method="classic", 
                           fine.tune = TRUE, 
                           assay.type.test = 1, 
                           BPPARAM=MulticoreParam(8))

# label transfer
im[["mon.fine"]] <- pred.mon.fine$labels
im[["mon.course"]] <- pred.mon.course$labels
im[["course_project"]] <- paste(im$mon.course, im$project, sep = "_")
im[["sample"]] <- im$course_project
Idents(im) <- "mon.fine"
im <- RenameIdents(im, 
                   "Intermediate monocytes"="Monocyte", 
                   "Classical monocytes"="Monocyte", 
                   "Exhausted B cells"="B", 
                   "Effector memory CD8 T cells"="abT (CD8)", 
                   "Myeloid dendritic cells"="mDC", 
                   "T regulatory cells"="Treg", 
                   "Th1 cells"="abT (CD4)", 
                   "Th2 cells"="abT (CD4)", 
                   "Th1/Th17 cells"="abT (CD4)", 
                   "Switched memory B cells"="B", 
                   "Low-density basophils"="Mast",
                   "Terminal effector CD4 T cells"="abT (CD4)", 
                   "Progenitor cells"="Mast", 
                   "Natural killer cells"="NK", 
                   "Non-switched memory B cells"="B", 
                   "Central memory CD8 T cells"="abT (CD8)", 
                   "Terminal effector CD8 T cells"="abT (CD8)", 
                   "Non classical monocytes"="Monocyte", 
                   "Plasmacytoid dendritic cells"="pDC", 
                   "Naive CD4 T cells"="abT (CD4)", 
                   "Vd2 gd T cells"="gd T", 
                   "Th17 cells"="Th17",
                   "Follicular helper T cells"="abT (CD4)", 
                   "Plasmablasts"="Plasma", 
                   "MAIT cells"="MAIT", 
                   "Naive CD8 T cells"="abT (CD8)", 
                   "Low-density neutrophils"="Neutrophil", 
                   "Non-Vd2 gd T cells"="gd T", 
                   "Naive B cells"="B")
levels(im) <- c("abT (CD4)", "Th17", "MAIT", "abT (CD8)", "gd T", "Treg", "NK",
                "pDC", "B", "Plasma",
                "Neutrophil", "Monocyte", "mDC", "Mast")
im$mon.fine <- im@active.ident
DefaultAssay(im) <- "RNA"

im.markers <- FindAllMarkers(im, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# RP and HB genes are not useful for cell ID
im.mkr.sub <- im.markers[-grep("^RP", im.markers$gene),]
im.mkr.sub <- im.mkr.sub[-grep("^HB", im.mkr.sub$gene),] %>%
  group_by(cluster) 
im.mkr.top <- im.mkr.sub %>%
  top_n(n=5, wt=avg_logFC)
```

##### Figure 3A

``` r
Idents(im) <- "mon.fine"
im <- RunUMAP(im, dims = 1:50)
imColors <- c(carto.pal(pal1 = "blue.pal", n1 = 20)[c(20, 17, 14, 11, 8, 5, 2)], 
              carto.pal(pal1 = "purple.pal", n1 = 20)[c(17, 11, 5)], 
              carto.pal(pal1 = "pink.pal", n1 = 20)[c(18, 14, 10, 6)])
Fig3A_dimPlot <- DimPlot(object = im, 
                         label = FALSE, 
                         reduction = "umap", 
                         cols = imColors) +
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
Fig3A_prop <- plot_stat(im, 
                        plot_type = "prop_fill", 
                        group_by = "project") +
  scale_fill_manual(values = imColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 4)) +
  guides(fill = FALSE)

# stats for 2A 
Fig3A_propBySample <- plot_stat(im, 
                                plot_type = "prop_fill", 
                                group_by = "orig.ident") 
Fig3AStats <- as.data.frame(Fig3A_propBySample[["data"]])
Fig3AStats <- Fig3AStats[order(Fig3AStats$cluster),]
write.csv(Fig3AStats, "Fig3A_forStats.csv")
```

##### Figure 3B

``` r
TCgenes <- c(im.mkr.top$gene[1:5], 
             "IL17A", 
             "IL17F", 
             "CD4", 
             "IL2",
             im.mkr.top$gene[6:15], 
             "CD8A", 
             "CD8B",
             im.mkr.top$gene[16:25], 
             "TRDC", 
             "FOXP3", 
             "ICOS",
             im.mkr.top$gene[26:35])
TCgenes <- rev(TCgenes[!duplicated(TCgenes)])
remove <- c("NKG7", 
            "CTSW", 
            "TRBC1", 
            "GZMK")
TCgenes <- setdiff(TCgenes, 
                   remove)
immuneT <- subset(im, 
                  idents = c("abT (CD4)", 
                             "Th17", 
                             "MAIT", 
                             "abT (CD8)", 
                             "gd T", 
                             "Treg", 
                             "NK"))

Figure3B <- DotPlot(immuneT, 
                    features = TCgenes, 
                    col.min = 0, 
                    scale.by = "size", 
                    dot.min = .01) + 
  coord_flip() + 
  scale_color_gradient2(low = "grey90", 
                        high = "black")  +
  theme(axis.title = element_blank(),
        legend.position = "right",
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
```

##### Figure 3C

``` r
Bgenes <- im.mkr.top$gene[c(36:50)]
Bgenes <- Bgenes[!duplicated(Bgenes)]

immuneB <- subset(im, 
                  idents = c("pDC", 
                             "B", 
                             "Plasma"))
Figure3C <- DotPlot(immuneB, 
                    features = rev(c(Bgenes[1:4], 
                                     "IL3RA", 
                                     "CLEC4C", 
                                     "NRP1", 
                                     "AXL", 
                                     Bgenes[6:15], 
                                     Bgenes[5])), 
                    col.min = 0, 
                    scale.by = "size", 
                    dot.min = .01) + 
  coord_flip() + 
  scale_color_gradient2(low = "grey90", 
                        high = "black")  +
  theme(axis.title = element_blank(),
        legend.position = "right",
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
```

##### Figure 3D

``` r
Mgenes <- c(im.mkr.top$gene[c(51:70)])
Mgenes <- Mgenes[!duplicated(Mgenes)]
immuneM <- subset(im, 
                  idents = c("Neutrophil", 
                             "Monocyte", 
                             "mDC", 
                             "Mast"))
Figure3D <- DotPlot(immuneM, 
                    features = rev(c(Mgenes[1:10], 
                                     "CD14", 
                                     "CD163", 
                                     "MARCO", 
                                     "LYVE1",
                                     Mgenes[11:15], 
                                     "CD1C", 
                                     "CLEC9A", 
                                     Mgenes[16:20])), 
                    col.min = 0, 
                    scale.by = "size", 
                    dot.min = .03) + 
  coord_flip() + 
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        legend.position = "right",
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
```

##### Table S5

``` r
Idents(im) <- "mon.fine"
im <- RenameIdents(im, 
                   "abT (CD4)" = "T/NK", 
                   "Th17" = "T/NK", 
                   "MAIT" = "T/NK", 
                   "abT (CD8)" = "T/NK", 
                   "gd T" = "T/NK", 
                   "Treg" = "T/NK", 
                   "NK" = "T/NK",
                   "B" = "B/Plasma", 
                   "Plasma" = "B/Plasma", 
                   "pDC" = "Granulocyte/Myeloid", 
                   "Neutrophil" = "Granulocyte/Myeloid", 
                   "Monocyte" = "Granulocyte/Myeloid", 
                   "mDC" = "Granulocyte/Myeloid", 
                   "Mast" = "Granulocyte/Myeloid")

levels(im) <- c("T/NK", 
                "B/Plasma", 
                "Granulocyte/Myeloid")

im$generalGroups <- im@active.ident

TBMyeloid.markers <- FindAllMarkers(im, 
                                    only.pos = TRUE, 
                                    logfc.threshold = 0.25)
TBMyeloid.top <- TBMyeloid.markers %>% 
  group_by(cluster) %>% 
  top_n(n=20, wt=avg_logFC)
write.csv(TBMyeloid.top, "TableS5.csv")
```

### Microbe Recognition/Entry

``` r
tlrs <- c("TLR1", 
          "TLR2", 
          "TLR3", 
          "TLR4", 
          "TLR5", 
          "TLR6", 
          "TLR7", 
          "TLR8", 
          "TLR9", 
          "TLR10")
nlrs <- c("CIITA", 
          "NAIP", 
          "NOD1", 
          "NLRC4", 
          "NOD2", 
          "NLRC3", 
          "NLRC5", 
          "NLRX1",
          "NLRP1", 
          "NLRP2", 
          "NLRP3", 
          "NLRP4",
          "NLRP5", 
          "NLRP6", 
          "NLRP7",
          "NLRP8", 
          "NLRP9", 
          "NLRP10",
          "NLRP11", 
          "NLRP12", 
          "NLRP13",
          "NLRP14")  
clrs <- c("CLEC7A", 
          "OLR1", 
          "CLEC9A", 
          "CLEC2A", 
          "CLEC2B", 
          "CD69", 
          "CLEC6A",
          "CLEC4D",
          "CLEC4E", 
          "CLEC4C", 
          "CLEC4A",
          "CLEC12A",
          "CLEC1A", 
          "IFI16", 
          "LRRFIP1", 
          "MRC1")
rigs <- c("DDX58", 
          "IFIH1", 
          "DHX58")
prrs <- c(tlrs,
          nlrs, 
          clrs,
          rigs)
damps <- c("HMGB1", 
           "S100A8", 
           "S100A9", 
           "SAA1", 
           "SAA2", 
           "HMGN1", 
           "IL33", 
           "SAP130", 
           "HSPD1", 
           "HSP90B1")
dampR <- c("AGER", 
           "IRAK4", 
           "TREM1", 
           "P2RX1", 
           "CD24")
covidFeats <- c("ACE2", 
                "TMPRSS2", 
                "TMPRSS4",
                "TMPRSS11D", 
                "CTSB", 
                "CTSL", 
                "HNRNPA1", 
                "BSG", 
                "ZCRB1",
                "TOP3B", 
                "ANPEP",
                "CLEC4M",
                "DPP4", 
                "CD209", 
                "FURIN")
oralVirus <- c("ITGA2", "AXL", "TYRO3", "CLEC4G", "F11R", "HSPA1B", "NCAM1", "DAG1", #HSV1&2
               "SDC1", "SDC2", "SDC4", "GPC1", "LN5", "ITGA6", #HPV
               "CD55", #varicella zoster
               "SCARB1", "CLDN9", #epstein barr (mono) 
               "PVR", "TFRC", "LAMP1", #cytomegalo
               "IDE", "EPS15", "HSP1A", "PTX3", #coxsackie
               "CD81", "LDLR", "EGFR", "EPHA2", #Enterovirus
               "CD46", "SLAMF1", "NECTIN4", #measles
               "MOG", #rubella
               "CAV1", "RPSA", "GAS6", "CLDN1", "HAVCR1", "GRK2") #HIV

oralIntegrated$general.site <- paste(oralIntegrated$generalCellTypes, 
                                     oralIntegrated$project, sep = "_")
Idents(oralIntegrated) <- "general.site"
oralIntegrated_covid <- subset(oralIntegrated, 
                               idents = c("Other_GM", "Other_BM"), 
                               invert = TRUE)
levels(oralIntegrated_covid) <- rev(c("Endothelial_GM", "Endothelial_BM",
                                      "Fibroblast_GM", "Fibroblast_BM",
                                      "Immune_GM", "Immune_BM",
                                      "Epithelial_GM", "Epithelial_BM"))

generalViR_Dot_GvB <- DotPlot(oralIntegrated_covid, 
                              features = c(prrs, 
                                           damps,
                                           dampR, 
                                           covidFeats, 
                                           oralVirus), 
                              col.min = 0, 
                              scale.by = "size", 
                          dot.min = 0.00001,
                          dot.scale = 10) + 
  scale_color_gradient2(low = "grey90", high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp", 
                             order = 2))


# find out which factors are differentially regulated in endo cells by oral site
Idents(endo) <- "project"
covidMkrEndo <- FindAllMarkers(endo, 
                               logfc.threshold = 0, 
                               min.pct = 0, 
                               only.pos = TRUE)
significantCOV_endo <- as.data.frame(covidMkrEndo[covidMkrEndo$gene %in% 
                                                    c(covidFeats),])
significantOralV_endo <- as.data.frame(covidMkrEndo[covidMkrEndo$gene %in% 
                                                      c(oralVirus),])
significantPrr_endo <- as.data.frame(covidMkrEndo[covidMkrEndo$gene %in% 
                                                    prrs,])
significantDamp_endo <- as.data.frame(covidMkrEndo[covidMkrEndo$gene %in% 
                                                     c(damps,dampR),])
significantViR_endo <- rbind(significantCOV_endo, 
                             significantOralV_endo, 
                             significantPrr_endo, 
                             significantDamp_endo)
significantViR_endo_Dot <- DotPlot(endo, 
                                   features = rownames(significantViR_endo), 
                                   col.min = 0, 
                                   scale.by = "size", 
                                   dot.min = 0.00001, 
                                   dot.scale = 10) + 
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp",
                             order = 2))
# extract proportions
endoDot_results <- as.data.frame(significantViR_endo_Dot[["data"]])
endoDot_results <- endoDot_results[order(endoDot_results$features.plot),]
endoDot_results <- endoDot_results[endoDot_results$features.plot %in% 
                                     c("IL33", 
                                       "HSPD1", 
                                       "HSP90B1", 
                                       "S100A8", 
                                       "SAA1", 
                                       "IFI16", 
                                       "CLEC2B", 
                                       "CLEC7A", 
                                       "TLR4", 
                                       "NAIP"),]
endoDot_results[is.na(endoDot_results)] = 0


# find out which factors are differentially regulated in epi cells by oral site
Idents(epi) <- "project"
covidMkrepi <- FindAllMarkers(epi, 
                              logfc.threshold = 0, 
                              min.pct = 0, 
                              only.pos = TRUE)
significantCOV_epi <- as.data.frame(covidMkrepi[covidMkrepi$gene %in% 
                                                  c(covidFeats),])
significantOralV_epi <- as.data.frame(covidMkrepi[covidMkrepi$gene %in% 
                                                    c(oralVirus),])
significantPrr_epi <- as.data.frame(covidMkrepi[covidMkrepi$gene %in% 
                                                  prrs,])
significantDamp_epi <- as.data.frame(covidMkrepi[covidMkrepi$gene %in% 
                                                   c(damps,dampR),])
significantViR_epi <- rbind(significantCOV_epi, 
                            significantOralV_epi, 
                            significantPrr_epi, 
                            significantDamp_epi)

significantViR_epi_Dot <- DotPlot(epi, 
                                  features = c(rownames(significantViR_epi)), 
                                  col.min = 0, 
                                  scale.by = "size", 
                                  dot.min = 0.002, 
                                  dot.scale = 10) + 
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", order = 1), 
         size = guide_legend(title = "%cell.exp", order = 2))
# Extract proportion data from dot plot and make a bar graph
epiDot_results <- as.data.frame(significantViR_epi_Dot[["data"]])
epiDot_results <- epiDot_results[order(epiDot_results$features.plot),]
epiDot_results <- epiDot_results[epiDot_results$features.plot %in% 
                                   c("S100A8", 
                                     "S100A9", 
                                     "SAA1", 
                                     "SAA2", 
                                     "HMGN1", 
                                     "CLEC2B", 
                                     "CLEC7A", 
                                     "IFI16", 
                                     "TLR4", 
                                     "NLRP1"),]
epiDot_results[is.na(epiDot_results)] = 0

# find out which factors are differentially regulated in fib cells by oral site
Idents(fib) <- "project"
covidMkrFib <- FindAllMarkers(fib, 
                              logfc.threshold = 0, 
                              min.pct = 0, 
                              only.pos = TRUE)
significantCOV_fib <- as.data.frame(covidMkrFib[covidMkrFib$gene %in% 
                                                  c(covidFeats),])
significantOralV_fib <- as.data.frame(covidMkrFib[covidMkrFib$gene %in% 
                                                    c(oralVirus),])
significantPrr_fib <- as.data.frame(covidMkrFib[covidMkrFib$gene %in% 
                                                  prrs,])
significantDamp_fib <- as.data.frame(covidMkrFib[covidMkrFib$gene %in% 
                                                   c(damps,dampR),])
significantViR_fib <- rbind(significantCOV_fib, 
                            significantOralV_fib, 
                            significantPrr_fib, 
                            significantDamp_fib)

significantViR_fib_Dot <- DotPlot(fib, 
                                  features = c(rownames(significantViR_fib)), 
                                  col.min = 0,
                                  scale.by = "size", 
                                  dot.min = 0.002, 
                                  dot.scale = 10) + 
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp", 
                             order = 2))

# Extract proportion data from dot plot and make a bar graph
fibDot_results <- as.data.frame(significantViR_fib_Dot[["data"]])
fibDot_results <- fibDot_results[order(fibDot_results$features.plot),]
fibDot_results <- fibDot_results[fibDot_results$features.plot %in% 
                                   c("SAP130", 
                                     "HSPD1", 
                                     "HSP90B1", 
                                     "CD24", 
                                     "IL33", 
                                     "IFI16", 
                                     "NLRP1", 
                                     "TLR2", 
                                     "TLR4", 
                                     "TLR5"),]
fibDot_results[is.na(fibDot_results)] = 0

# find out which factors are differentially regulated in immune cells by oral site
Idents(im) <- "project"
covidMkrim <- FindAllMarkers(im, logfc.threshold = 0, min.pct = 0, only.pos = TRUE)
significantCOV_im <- as.data.frame(covidMkrim[covidMkrim$gene %in% c(covidFeats),])
significantOralV_im <- as.data.frame(covidMkrim[covidMkrim$gene %in% c(oralVirus),])
significantPrr_im <- as.data.frame(covidMkrim[covidMkrim$gene %in% prrs,])
significantDamp_im <- as.data.frame(covidMkrim[covidMkrim$gene %in% c(damps,dampR),])
significantViR_im <- rbind(significantCOV_im, significantOralV_im, significantPrr_im, significantDamp_im)

significantViR_im_Dot <- DotPlot(im, 
                                 features = rownames(significantViR_im), 
                                 col.min = 0, 
                                 scale.by = "size", 
                                 dot.min = 0.002, 
                                 dot.scale = 10) + 
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp", 
                             order = 2))

# Extract proportion data from dot plot and make a bar graph
imDot_results <- as.data.frame(significantViR_im_Dot[["data"]])
imDot_results <- imDot_results[order(imDot_results$features.plot),]
imDot_results <- imDot_results[imDot_results$features.plot %in% 
                                 c("SAP130", 
                                   "AGER", 
                                   "IRAK4", 
                                   "TREM1", 
                                   "P2RX1", 
                                   "CD69", 
                                   "CLEC4E", 
                                   "MRC1", 
                                   "CIITA", 
                                   "TLR4"),]
imDot_results[is.na(imDot_results)] = 0
```

##### Figure 7B

``` r
Fig7B_endoGvB <- ggplot(endoDot_results,
                        aes(fill=id,
                            y = pct.exp, 
                            x = features.plot)) +
  geom_bar(position = "fill", 
           stat = "identity") + 
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
Fig7B_epiGvB <- ggplot(epiDot_results, 
                       aes(fill=id, 
                           y = pct.exp, 
                           x = features.plot)) +
  geom_bar(position = "fill", 
           stat = "identity") + 
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
Fig7B_fibGvB <- ggplot(fibDot_results,
                       aes(fill=id, 
                           y = pct.exp, 
                           x = features.plot)) +
  geom_bar(position = "fill",
           stat = "identity") + 
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
Fig7B_immuneGvB <- ggplot(imDot_results, 
                          aes(fill=id, 
                              y = pct.exp,
                              x = features.plot)) +
  geom_bar(position = "fill", 
           stat = "identity") + 
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
```

##### Figure S7A

``` r
FigS7B <- generalViR_Dot_GvB
```

## Oral, Skin, Lung, Ileum

> This section covers **Figure 1E, Figure S3B**

### Dataset preparation

#### Transformation

``` r
HO <- merge(GM, y=BM)
HS <- NormalizeData(HS)
HS <- FindVariableFeatures(HS, nfeatures = 4000)
HL <- NormalizeData(HL)
HL <- FindVariableFeatures(HL, nfeatures = 4000)
HI <- NormalizeData(HI)
HI <- FindVariableFeatures(HI, nfeatures = 4000)
```

#### Downsampling

``` r
# use fewest cells for downsampling
ncol(HS)
ncol(HO)
ncol(HL)
ncol(HI)

# Downsample all datasets
HO.sub <- HO[, sample(colnames(HO), size = ncol(HS), replace = F)]
HL.sub <- HL[, sample(colnames(HL), size = ncol(HS), replace = F)]
HI.sub <- HI[, sample(colnames(HI), size = ncol(HS), replace = F)]
```

#### Integration and batch correction

``` r
sampleList <- list(HS, HO.sub, HI.sub, HL.sub)
oral.anchors <- FindIntegrationAnchors(object.list = sampleList, dims = 1:50)
oralIntegrated_skinLungIleum <- IntegrateData(anchorset = oral.anchors, dims = 1:50)
```

#### Dimensionality reduction and cell clustering

``` r
# cell-cycle scoring and regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
oralIntegrated_skinLungIleum <- CellCycleScoring(oralIntegrated_skinLungIleum, 
                                                 s.features = s.genes, 
                                                 g2m.features = g2m.genes,
                                                 set.ident = TRUE)
# standard workflow for clustering
oralIntegrated_skinLungIleum <- ScaleData(oralIntegrated_skinLungIleum, 
                                          vars.to.regress = c("S.Score", 
                                                              "G2M.Score"), 
                                          verbose = TRUE)
oralIntegrated_skinLungIleum <- RunPCA(oralIntegrated_skinLungIleum, 
                                       npcs = 50, 
                                       verbose = TRUE)
oralIntegrated_skinLungIleum <- FindNeighbors(oralIntegrated_skinLungIleum, 
                                              reduction = "pca", 
                                              dims = 1:50)
oralIntegrated_skinLungIleum <- FindClusters(oralIntegrated_skinLungIleum, 
                                             resolution = seq(from = 0.1, 
                                                              to = 1.0, 
                                                              by = 0.1))
oralIntegrated_skinLungIleum <- RunUMAP(oralIntegrated_skinLungIleum, 
                                        reduction = "pca", 
                                        dims = 1:50)
oralIntegrated_skinLungIleum <- RunTSNE(oralIntegrated_skinLungIleum, 
                                        check_duplicates = FALSE)
# assess cluster tree
clustree(oralIntegrated_skinLungIleum)
```

### High level summary

#### Cell classification

``` r
Idents(oralIntegrated_skinLungIleum) <- "project"
oralIntegrated_skinLungIleum <- RenameIdents(oralIntegrated_skinLungIleum, 
                                             "BM" = "Oral", 
                                             "GM" = "Oral")
oralIntegrated_skinLungIleum$project2 <- oralIntegrated_skinLungIleum@active.ident

# cluster markers
Idents(oralIntegrated_skinLungIleum) <- "integrated_snn_res.0.1"
allmarkers_R01 <- FindAllMarkers(oralIntegrated_skinLungIleum, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.3)

# rename clusters
oralIntegrated_skinLungIleum <- RenameIdents(oralIntegrated_skinLungIleum, "0"="Immune", "1"="Epithelial", "2"="Fibroblast", 
                                             "3"="Endothelial", "4"="Immune", "5"="Immune", "6"="Immune", 
                                             "7"="Endothelial", "8"="Other", "9"="Endothelial", "10"="Epithelial", 
                                             "11"="Endothelial", "12"="Other", "13"="Fibroblast")
levels(oralIntegrated_skinLungIleum) <- c("Endothelial", "Fibroblast", "Immune", "Epithelial", 
                                          "Other")
oralIntegrated_skinLungIleum$extraCourseID <- oralIntegrated_skinLungIleum@active.ident
```

##### Figure 1E

``` r
genCols <- c(paletteer_d("ggsci::nrc_npg")[c(1,3,4,9)], 
             paletteer_d("ggsci::default_jco")[3])


Fig1E_dimPlot <- DimPlot(oralIntegrated_skinLungIleum, 
                 reduction = "tsne", 
                 split.by = "project2", 
                 cols = genCols, 
                 pt.size = 0.01)

Fig1E_prop <- plot_stat(oralIntegrated_skinLungIleum, plot_type = "prop_fill", group_by = "project2")+ 
  scale_fill_manual(values = genCols) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12)) +
  guides(fill = FALSE)
```

### Immune Subset

#### Cell classification

``` r
Idents(oralIntegrated_skinLungIleum) <- "extraCourseID"
immune_OSLI <- subset(oralIntegrated_skinLungIleum, idents = "Immune")

Idents(immune_OSLI) <- "project2"
Idents(oralIntegrated_skinLungIleum) <- "extraCourseID"

# use SingleR
# seurat object -> SingleCellExperiment
immune_SCE <- as.SingleCellExperiment(immune_OSLI)

# reference database
mon.se <- MonacoImmuneData()
blueprint <- BlueprintEncodeData()
commonBlueIm <- intersect(rownames(immune_OSLI), rownames(blueprint))
blueprintIm <- blueprint[commonBlueIm,]
immune_SCE <- immune_SCE[commonBlueIm,]
pred.mon.fine <- SingleR(test=immune_SCE, 
                         ref=mon.se, 
                         labels=mon.se$label.fine,
                         de.method="classic", 
                         fine.tune = TRUE, 
                         assay.type.test = 1, 
                         BPPARAM=MulticoreParam(8))
pred.mon.course <- SingleR(test=immune_SCE, 
                           ref=mon.se, 
                           labels=mon.se$label.main,
                           de.method="classic", 
                           fine.tune = TRUE, 
                           assay.type.test = 1, 
                           BPPARAM=MulticoreParam(8))

# label transfer
immune_OSLI[["mon.fine"]] <- pred.mon.fine$labels
immune_OSLI[["mon.course"]] <- pred.mon.course$labels
Idents(immune_OSLI) <- "mon.fine"
immune_OSLI <- RenameIdents(immune_OSLI, 
                            "Myeloid dendritic cells"="mDC", 
                            "Plasmablasts"="Plasma", 
                            "Classical monocytes"="Monocyte", 
                            "T regulatory cells"="Treg", 
                            "Non-switched memory B cells"="B", 
                            "Th1/Th17 cells"="abT (CD4)", 
                            "Vd2 gd T cells"="gd T", 
                            "Th1 cells"="abT (CD4)", 
                            "Effector memory CD8 T cells"="abT (CD8)", 
                            "Th17 cells"="Th17",
                            "Intermediate monocytes"="Monocyte", 
                            "Switched memory B cells"="B", 
                            "Progenitor cells"="Mast", 
                            "Natural killer cells"="NK", 
                            "MAIT cells"="MAIT", 
                            "Exhausted B cells"="B", 
                            "Central memory CD8 T cells"="abT (CD8)", 
                            "Terminal effector CD8 T cells"="abT (CD8)", 
                            "Th2 cells"="abT (CD4)", 
                            "Non classical monocytes"="Monocyte", 
                            "Low-density basophils"="Mast",
                            "Plasmacytoid dendritic cells"="pDC", 
                            "Naive CD8 T cells"="abT (CD8)", 
                            "Naive B cells"="B",
                            "Follicular helper T cells"="abT (CD4)", 
                            "Naive CD4 T cells"="abT (CD4)", 
                            "Low-density neutrophils"="Neutrophil", 
                            "Non-Vd2 gd T cells"="gd T", 
                            "Terminal effector CD4 T cells"="abT (CD4)")
immune_OSLI$mon.fine.adjusted <- immune_OSLI@active.ident
mylevels <- c("abT (CD4)", "Th17", "MAIT", "abT (CD8)", "gd T", "Treg", "NK", "pDC", "B", "Plasma", "Neutrophil", "Monocyte", "mDC", "Mast")
immune_OSLI$mon.fine.adjusted <- factor(x = immune_OSLI$mon.fine.adjusted, levels = mylevels)
```

##### Figure S3B

``` r
projectCols <- paletteer_d("nationalparkcolors::Redwoods")
Idents(immune_OSLI) <- "project2"
FigS3B <- plot_stat(immune_OSLI, plot_type = "prop_fill", group_by = "mon.fine.adjusted") +
  scale_fill_manual(values = projectCols) + 
  theme_cowplot() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12))
Idents(immune_OSLI) <- "mon.fine.adjusted"
statsPlot <- plot_stat(immune_OSLI, plot_type = "prop_fill", group_by = "orig.ident") +
  scale_fill_manual(values = immuneSubsetColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12)) +
  guides(fill = FALSE)
FigS3ImStats <- as.data.frame(statsPlot[["data"]])
FigS3ImStats <- FigS3ImStats[order(FigS3ImStats$cluster),]
write.csv(FigS3ImStats, "FigS3B_forStats.csv")
```

## Oral vs. Salivary

> This section covers **Figure S3C**

### Dataset preparation

#### Transformation

``` r
SG <- NormalizeData(SG)
SG <- FindVariableFeatures(SG, nfeatures = 4000)
SG$project <- "SG" 
```

#### Downsampling

``` r
# use fewest cells for downsampling
ncol(SG)
ncol(GM)
ncol(BM)

# Downsample all datasets
GM.sub <- GM[, sample(colnames(GM), size = ncol(SG), replace = F)]
BM.sub <- BM[, sample(colnames(BM), size = ncol(SG), replace = F)]
```

#### Integration and batch correction

``` r
sampleList <- list(GM.sub, BM.sub, SG)
oral.anchors <- FindIntegrationAnchors(object.list = sampleList, dims = 1:50)
oralIntegrated <- IntegrateData(anchorset = oral.anchors, dims = 1:50)
```

#### Dimensionality reduction and cell clustering

``` r
# standard workflow for clustering
oralIntegrated_SG <- ScaleData(oralIntegrated_SG, 
                               verbose = TRUE)
oralIntegrated_SG <- RunPCA(oralIntegrated_SG, 
                            npcs = 50, 
                            verbose = TRUE)
oralIntegrated_SG <- FindNeighbors(oralIntegrated_SG, 
                                   reduction = "pca", 
                                   dims = 1:50)
oralIntegrated_SG <- FindClusters(oralIntegrated_SG, 
                                  resolution = seq(from = 0.1, 
                                                   to = 1.0, 
                                                   by = 0.1))
oralIntegrated_SG <- RunUMAP(oralIntegrated_SG, 
                             reduction = "pca", 
                             dims = 1:50)
DefaultAssay(oralIntegrated_SG) <- "RNA"
Idents(oralIntegrated_SG) <- "integrated_snn_res.1"
```

### High level summary

#### Cell classification

``` r
# find cluster defining markers
all.markers.1 <- FindAllMarkers(oralIntegrated_SG, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.75)
significant.markers.1 <- all.markers.1[all.markers.1$p_val_adj < 0.2, ]
top.markers.1 <- significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=20, wt=avg_logFC)
# remove RBC contaminated clusters
oralIntegrated_SG <- subset(oralIntegrated_SG, idents = c("20", "28"), invert = TRUE)
# redo clustering
oralIntegrated_SG <- ScaleData(oralIntegrated_SG, 
                               verbose = TRUE)
oralIntegrated_SG <- RunPCA(oralIntegrated_SG, 
                            npcs = 50, 
                            verbose = TRUE)
oralIntegrated_SG <- FindNeighbors(oralIntegrated_SG, 
                                   reduction = "pca", 
                                   dims = 1:50)
oralIntegrated_SG <- FindClusters(oralIntegrated_SG, 
                                  resolution = seq(from = 0.1, 
                                                   to = 1.0, 
                                                   by = 0.1))
oralIntegrated_SG <- RunUMAP(oralIntegrated_SG, 
                             reduction = "pca", 
                             dims = 1:50)
DefaultAssay(oralIntegrated_SG) <- "RNA"
Idents(oralIntegrated_SG) <- "project"
oralIntegrated_SG <- RenameIdents(oralIntegrated_SG,
                                  "GM"="Mucosa", 
                                  "BM"="Mucosa")
oralIntegrated_SG$project <- oralIntegrated_SG@active.ident
Idents(oralIntegrated_SG) <- "integrated_snn_res.1"
# find cluster defining markers again
all.markers.1 <- FindAllMarkers(oralIntegrated_SG, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.75)
significant.markers.1 <- all.markers.1[all.markers.1$p_val_adj < 0.2, ]
top.markers.1 <- significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=20, wt=avg_logFC)
```

### Immune Subset

#### Cell classification

``` r
# subset immune clusters
immune <- subset(oralIntegrated_SG, idents = c("4", "7", "8", "9", "19", "21", "26"))

# prepare for SingleR
immune_SCE <- as.SingleCellExperiment(immune)
mon.se <- MonacoImmuneData()
blueprint <- BlueprintEncodeData()
commonBlueIm <- intersect(rownames(immune), rownames(blueprint))
blueprintIm <- blueprint[commonBlueIm,]
immune_SCE <- immune_SCE[commonBlueIm,]

# use SingleR to identify cells
pred.mon.fine <- SingleR(test=immune_SCE, ref=mon.se, labels=mon.se$label.fine,
                         de.method="classic", fine.tune = TRUE, assay.type.test = 1, BPPARAM=MulticoreParam(8))
pred.mon.course <- SingleR(test=immune_SCE, ref=mon.se, labels=mon.se$label.main,
                           de.method="classic", fine.tune = TRUE, assay.type.test = 1, BPPARAM=MulticoreParam(8))

# label transfer
immune[["mon.fine"]] <- pred.mon.fine$labels
immune[["mon.course"]] <- pred.mon.course$labels
Idents(immune) <- "mon.fine"
DefaultAssay(immune) <- "integrated"
immune <- RunPCA(immune, npcs = 50)
DimPlot(immune, group.by = "mon.fine", label = TRUE)

immune <- RenameIdents(immune, 
                       "Myeloid dendritic cells"="mDC", 
                       "Plasmablasts"="Plasma", 
                       "Classical monocytes"="Monocyte", 
                       "T regulatory cells"="Treg", 
                       "Non-switched memory B cells"="B", 
                       "Th1/Th17 cells"="abT (CD4)", 
                       "Vd2 gd T cells"="gd T", 
                       "Th1 cells"="abT (CD4)", 
                       "Effector memory CD8 T cells"="abT (CD8)", 
                       "Th17 cells"="Th17",
                       "Intermediate monocytes"="Monocyte", 
                       "Switched memory B cells"="B", 
                       "Progenitor cells"="Mast", 
                       "Natural killer cells"="NK", 
                       "MAIT cells"="MAIT", 
                       "Exhausted B cells"="B", 
                       "Central memory CD8 T cells"="abT (CD8)", 
                       "Terminal effector CD8 T cells"="abT (CD8)", 
                       "Th2 cells"="abT (CD4)", 
                       "Non classical monocytes"="Monocyte", 
                       "Low-density basophils"="Mast",
                       "Plasmacytoid dendritic cells"="pDC", 
                       "Naive CD8 T cells"="abT (CD8)", 
                       "Naive B cells"="B",
                       "Follicular helper T cells"="abT (CD4)", 
                       "Naive CD4 T cells"="abT (CD4)", 
                       "Low-density neutrophils"="Neutrophil", 
                       "Non-Vd2 gd T cells"="gd T", 
                       "Terminal effector CD4 T cells"="abT (CD4)")
immune$mon.fine.adjusted <- immune@active.ident
mylevels <- c("abT (CD4)", 
              "Th17", "MAIT", 
              "abT (CD8)", 
              "gd T", 
              "Treg", 
              "NK", 
              "pDC", 
              "B", 
              "Plasma", 
              "Neutrophil", 
              "Monocyte", 
              "mDC", 
              "Mast")
immune$mon.fine.adjusted <- factor(x = immune$mon.fine.adjusted, 
                                   levels = mylevels)
Idents(immune) <- "mon.fine.adjusted"
immuneSubsetColors <- c(carto.pal(pal1 = "blue.pal", n1 = 20)[c(20, 17, 14, 11, 8, 5, 2)], 
                        carto.pal(pal1 = "purple.pal", n1 = 20)[c(17, 11, 5)], 
                        carto.pal(pal1 = "pink.pal", n1 = 20)[c(18, 14, 10, 6)])
Idents(immune) <- "project"
immune <- RenameIdents(immune, 
                       "GM"="Mucosa", 
                       "BM"="Mucosa")
immune$project <- immune@active.ident
```

##### Figure S3C

``` r
Idents(oralIntegrated_SG) <- "integrated_snn_res.1"
FigS3C_dimPlot <- DimPlot(object = oralIntegrated_SG, 
                          label = TRUE, 
                          reduction = "umap", 
                          split.by = "project") 

Idents(oralIntegrated_SG) <- "project"
FigS3C_prop <- plot_stat(oralIntegrated_SG, 
                         plot_type = "prop_fill", 
                         group_by = "integrated_snn_res.1") + 
  theme_cowplot() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank())
FigS3C_immuneProp <- plot_stat(immune, plot_type = "prop_fill", group_by = "mon.fine.adjusted") +
  theme_cowplot() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12))
```

## Health vs. Perio

> This section covers **Figure 4E, Figure 5A-D, Figure 6, Figure 7,
> Figure S5, Figure S6, Figure S7, Table S6, Table S7**

### Data preparation

#### Transformation

``` r
PD <- NormalizeData(BM)
PD <- FindVariableFeatures(BM, nfeatures = 4000)
```

#### Integration and batch correction

``` r
sampleList <- list(GM, PD)
oral.anchors <- FindIntegrationAnchors(object.list = sampleList, dims = 1:50)
oralIntegrated_HvP <- IntegrateData(anchorset = oral.anchors, dims = 1:50)
```

#### Dimensionality reduction and cell clustering

``` r
# cell-cycle scoring and regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
oralIntegrated <- CellCycleScoring(oralIntegrated, 
                                   s.features = s.genes, 
                                   g2m.features = g2m.genes,
                                   set.ident = TRUE)
# standard workflow for clustering
oralIntegrated_HvP <- CellCycleScoring(oralIntegrated_HvP, 
                                       s.features = s.genes, 
                                       g2m.features = g2m.genes,
                                       set.ident = TRUE)

# Run the standard workflow for clustering
oralIntegrated_HvP <- ScaleData(oralIntegrated_HvP, 
                                vars.to.regress = c("S.Score", 
                                                    "G2M.Score"),
                                verbose = TRUE)
oralIntegrated_HvP <- RunPCA(oralIntegrated_HvP, 
                             npcs = 50, 
                             verbose = FALSE)
oralIntegrated_HvP <- FindNeighbors(oralIntegrated_HvP, 
                                    reduction = "pca", 
                                    dims = 1:50)
oralIntegrated_HvP <- FindClusters(oralIntegrated_HvP, 
                                   resolution = 1)
oralIntegrated_HvP <- RunUMAP(oralIntegrated_HvP, 
                              reduction = "pca", 
                              dims = 1:50)
```

### High-level summary

#### Cell classification

``` r
# set correct assay and idents for GEX
DefaultAssay(oralIntegrated_HvP) <- "RNA"
Idents(oralIntegrated_HvP) <- "integrated_snn_res.1"
all.markers_HvP <- FindAllMarkers(oralIntegrated_HvP, 
                                  only.pos = TRUE)
significant.markers_HvP <- all.markers_HvP[all.markers_HvP$p_val_adj < 0.2, ]
top20.markers_HvP <- significant.markers_HvP %>% 
  group_by(cluster) %>% 
  top_n(n=20, wt=avg_logFC)
write.csv(top20.markers_HvP, "integrated_HvP_top20Markers.csv")
oralIntegrated_HvP <- RenameIdents(oralIntegrated_HvP, 
                                   "0"="Endothelial", 
                                   "1"="Fibroblast", 
                                   "2"="Immune", 
                                   "3"="Endothelial", 
                                   "4"="Immune", 
                                   "5"="Endothelial", 
                                   "6"="Immune", 
                                   "7"="Endothelial", 
                                   "8"="Immune", 
                                   "9"="Epithelial", 
                                   "10"="Fibroblast", 
                                   "11"="Endothelial", 
                                   "12"="Immune", 
                                   "13"="Immune", 
                                   "14"="Fibroblast", 
                                   "15"="Immune", 
                                   "16"="Immune", 
                                   "17"="Fibroblast", 
                                   "18"="RBC", 
                                   "19"="Endothelial", 
                                   "20"="Fibroblast", 
                                   "21"="Other", 
                                   "22"="Immune", 
                                   "23"="Epithelial",
                                   "24"="Immune", 
                                   "25"="Immune", 
                                   "26"="Epithelial", 
                                   "27"="Immune", 
                                   "28"="Endothelial", 
                                   "29"="Immune", 
                                   "30"="Epithelial", 
                                   "31" = "Mixed")

# get rid of RBC and mixed genotype
oralIntegrated_HvP <- subset(oralIntegrated_HvP, 
                             idents = c("RBC", "Mixed"), 
                             invert = TRUE)

levels(oralIntegrated_HvP) <- c("Endothelial", 
                                "Fibroblast", 
                                "Immune", 
                                "Epithelial", 
                                "Other")

# save metadata
oralIntegrated_HvP$celltype <- oralIntegrated_HvP@active.ident
oralIntegrated_HvP$celltype.site <- paste(Idents(oralIntegrated_HvP),
                                          oralIntegrated_HvP$project, 
                                          sep = "_")

# rename the cluster subtypes
Idents(oralIntegrated_HvP) <- "integrated_snn_res.1"
oralIntegrated_HvP <- RenameIdents(oralIntegrated_HvP, "0"="P.VEC 1.1", "1"="P.Fib 1.1", "2"="Immune", "3"="P.VEC 1.2", "4"="Immune", 
                                   "5"="P.SMC", "6"="Immune", "7"="P.VEC 1.3", "8"="Immune", "9"="P.Epi 1", "10"="P.Fib 1.2", 
                                   "11"="P.VEC 1.4", "12"="Immune", "13"="Immune", "14"="P.Fib 1.3", "15"="Immune", 
                                   "16"="Immune", "17"="P.Fib 1.4", "19"="P.LEC", "20"="P.Fib 1.5", "21"="Other", 
                                   "22"="Immune", "23"="P.Epi 2", "24"="Immune", "25"="Immune", "26"="P.Epi 3", "27"="Immune", 
                                   "28"="P.SMC", "29"="Immune", "30"="P.Mel")
# save metadata
oralIntegrated_HvP$fineCellType <- oralIntegrated_HvP@active.ident
```

##### Figure 4E

``` r
genColors <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5)]
Idents(oralIntegrated_HvP) <- "celltype"
Fig4E_dimPlot <- DimPlot(object = oralIntegrated_HvP, 
                    label = FALSE, 
                    reduction = "umap", 
                    cols = genColors) +
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
Fig4E_prop <-  plot_stat(oralIntegrated_HvP_includeOther, 
                            plot_type = "prop_fill", 
                            group_by = "project") +
  scale_fill_manual(values = genColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12)) +
  guides(fill = FALSE)

# proportion of major cell types by sample for stats
Fig4E_propBySample <- plot_stat(oralIntegrated_HvP, 
                                     plot_type = "prop_fill", 
                                     group_by = "orig.ident")
Fig4eStats <- as.data.frame(gmpd_proportionBySample[["data"]])
Fig4eStats <- Fig4eStats[order(Fig4eStats$cluster),]
write.csv(Fig4eStats, "Fig4e_forStatistics.csv")
```

### Immune subset

#### Cell classification

``` r
# use SingleR
# seurat object -> SingleCellExperiment
immune_HvP <- subset(oralIntegrated_HvP, 
                     idents = c("Immune"))
immune_HvPSCE <- as.SingleCellExperiment(immune_HvP)

# reference database
mon.se <- MonacoImmuneData()
blueprint <- BlueprintEncodeData()
commonBlue <- intersect(rownames(immune_HvP), 
                        rownames(blueprint))
blueprint <- blueprint[commonBlue,]
immune_HvPSCE <- immune_HvPSCE[commonBlue,]
immune_HvPSCE <- logNormCounts(immune_HvPSCE)
pred.mon.fine <- SingleR(test=immune_HvPSCE, 
                         ref=mon.se, 
                         labels=mon.se$label.fine,
                         de.method="classic", 
                         fine.tune = TRUE, 
                         assay.type.test = 1, 
                         BPPARAM=MulticoreParam(8))
pred.mon.course <- SingleR(test=immune_HvPSCE, 
                           ref=mon.se, 
                           labels=mon.se$label.main,
                           de.method="classic", 
                           fine.tune = TRUE, 
                           assay.type.test = 1, 
                           BPPARAM=MulticoreParam(8))
immune_HvP[["mon.fine"]] <- pred.mon.fine$labels
immune_HvP[["mon.course"]] <- pred.mon.course$labels

Idents(immune_HvP) <- "mon.fine"
DefaultAssay(immune_HvP) <- "integrated"
immune_HvP <- RunPCA(immune_HvP, 
                     npcs = 50)
immune_HvP <- RunUMAP(immune_HvP,
                      reduction = "pca", 
                      dims = 1:50)
DefaultAssay(immune_HvP) <- "RNA"
immuneMonaco_fine <- FindAllMarkers(immune_HvP)
immuneMonacoFineTop <- immuneMonaco_fine %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC)
immune_HvP <- RenameIdents(immune_HvP, 
                           "Myeloid dendritic cells"="mDC", 
                           "Plasmablasts"="Plasma", 
                           "Classical monocytes"="Monocyte", 
                           "T regulatory cells"="Treg", 
                           "Non-switched memory B cells"="B",
                           "Th1/Th17 cells"="abT (CD4)", 
                           "Vd2 gd T cells"="gd T", 
                           "Th1 cells"="abT (CD4)", 
                           "Effector memory CD8 T cells"="abT (CD8)", 
                           "Th17 cells"="Th17",
                           "Intermediate monocytes"="Monocyte",
                           "Switched memory B cells"="B",
                           "Progenitor cells"="Mast", 
                           "Natural killer cells"="NK", 
                           "MAIT cells"="MAIT",
                           "Exhausted B cells"="B", 
                           "Central memory CD8 T cells"="abT (CD8)", 
                           "Terminal effector CD8 T cells"="abT (CD8)", 
                           "Th2 cells"="abT (CD4)", 
                           "Non classical monocytes"="Monocyte", 
                           "Low-density basophils"="Mast",
                           "Plasmacytoid dendritic cells"="pDC", 
                           "Naive CD8 T cells"="abT (CD8)", 
                           "Naive B cells"="B",
                           "Follicular helper T cells"="abT (CD4)", 
                           "Naive CD4 T cells"="abT (CD4)", 
                           "Low-density neutrophils"="Neutrophil", 
                           "Non-Vd2 gd T cells"="gd T", 
                           "Terminal effector CD4 T cells"="abT (CD4)")
levels(immune_HvP) <- c("abT (CD4)", 
                        "Th17", 
                        "MAIT", 
                        "abT (CD8)", 
                        "gd T", 
                        "Treg", 
                        "NK", 
                        "pDC", 
                        "B", 
                        "Plasma", 
                        "Neutrophil", 
                        "Monocyte", 
                        "mDC", 
                        "Mast")
immune_HvP$mon.fine.adjusted <- immune_HvP@active.ident
immune_HvP$monFine_project <- paste(immune_HvP@active.ident, immune_HvP$project, sep = "_")
Idents(immune_HvP) <- "mon.fine.adjusted"
DefaultAssay(immune_HvP) <- "integrated"
immune_HvP <- RunPCA(immune_HvP, 
                     npcs = 50,
                     verbose = FALSE)
immune_HvP <- RunUMAP(immune_HvP, 
                      reduction = "pca", 
                      dims = 1:50)
DefaultAssay(immune_HvP) <- "RNA"
immuneMonaco_fine <- FindAllMarkers(immune_HvP)
immuneMonacoFineTop <- immuneMonaco_fine %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC)

immuneSubsetColors <- c(carto.pal(pal1 = "blue.pal", n1 = 20)[c(20, 17, 14, 11, 8, 5, 2)], 
                        carto.pal(pal1 = "purple.pal", n1 = 20)[c(17, 11, 5)], 
                        carto.pal(pal1 = "pink.pal", n1 = 20)[c(18, 14, 10, 6)])

immune_HvP <- RenameIdents(immune_HvP, 
                           "abT (CD4)" = "T/NK", 
                           "Th17" = "T/NK", 
                           "MAIT" = "T/NK", 
                           "abT (CD8)" = "T/NK",
                           "gd T" = "T/NK", 
                           "Treg" = "T/NK", 
                           "NK" = "T/NK", 
                           "pDC" = "Myeloid/Granulocyte", 
                           "Neutrophil" = "Myeloid/Granulocyte", 
                           "Monocyte" = "Myeloid/Granulocyte", 
                           "mDC" = "Myeloid/Granulocyte", 
                           "Mast" = "Myeloid/Granulocyte", 
                           "B" = "B/Plasma", 
                           "Plasma" = "B/Plasma")
levels(immune_HvP) <- c("T/NK", 
                        "B/Plasma", 
                        "Myeloid/Granulocyte")
immune_HvP$courseImmune <- immune_HvP@active.ident
immune_HvP$monCourse_project <- paste(immune_HvP$courseImmune, 
                                      immune_HvP$project, sep = "_")
```

##### Figure 5A,B

``` r
Idents(immune_HvP) <- "mon.fine.adjusted"
Fig5A_dimPlot <- DimPlot(immune_HvP, 
                         label = FALSE, 
                         cols = immuneSubsetColors)

Fig5AB_prop <- plot_stat(immune_HvP, plot_type = "prop_fill", group_by = "project") +
  scale_fill_manual(values = immuneSubsetColors) +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12)) +
  guides(fill = FALSE)


################ Figure 5a stats #######################
Idents(immune_HvP) <-"mon.fine.adjusted" 
Fig5AB_propBySample <- plot_stat(immune_HvP, 
                                plot_type = "prop_fill", 
                                group_by = "orig.ident")
Fig5AStats <- as.data.frame(Fig5A_propBySample[["data"]])
Fig5AStats <- Fig5AStats[order(Fig5aStats$cluster),]
write.csv(Fig5AStats, "Fig5A_forStats.csv")
```

##### Figure 5C

``` r
Idents(immune_HvP) <- "mon.fine.adjusted"
immune_HvP_Plasma <- subset(immune_HvP, 
                            idents = "Plasma")
immune_HvP_Plasma <- RunUMAP(immune_HvP_Plasma, 
                             reduction = "pca", 
                             dims = 1:50)

immune_HvP_Plasma <- make_hexbin(immune_HvP_Plasma, 
                                 nbins = 20, 
                                 dimension_reduction = "UMAP")

plasma1 <- plot_hexbin_gene(immune_HvP_Plasma, 
                            type = "data", 
                            gene = "IGHG1", 
                            action = "mean") + 
  scale_fill_viridis_c(option = "A") +
  theme_cowplot()

plasma2 <- plot_hexbin_gene(immune_HvP_Plasma, 
                            type = "data", 
                            gene = "IGHA1", 
                            action = "mean") + 
  scale_fill_viridis_c(option = "A") +
  theme_cowplot()

plasma3 <- plot_hexbin_gene(immune_HvP_Plasma, 
                            type = "data", 
                            gene = "IGKC", 
                            action = "mean") + 
  scale_fill_viridis_c(option = "A") + 
  theme_cowplot()

plasma4 <- plot_hexbin_gene(immune_HvP_Plasma, 
                            type = "data", 
                            gene = "IGLC2", 
                            action = "mean") + 
  scale_fill_viridis_c(option = "A") +
  theme_cowplot()

Fig5C <- (plasma1|plasma2|plasma3|plasma4)
```

##### Figure 5D

``` r
Idents(immune_HvP) <- "mon.fine.adjusted"
immune_HvP_Neut <- subset(immune_HvP, 
                          idents = "Neutrophil")
immune_HvP_Neut <- RunUMAP(immune_HvP_Neut,
                           reduction = "pca", 
                           dims = 1:50)
immune_HvP_Neut <- make_hexbin(immune_HvP_Neut, 
                               nbins = 10, 
                               dimension_reduction = "UMAP")
neut1 <- plot_hexbin_gene(immune_HvP_Neut, 
                          type = "data", 
                          gene = "CXCL8", 
                          action = "mean") + 
  scale_fill_viridis_c(option = "A") +
  theme_cowplot()

neut2 <- plot_hexbin_gene(immune_HvP_Neut, 
                          type = "data", 
                          gene = "SOD2", 
                          action = "mean") + 
  scale_fill_viridis_c(option = "A") +
  theme_cowplot()

neut3 <- plot_hexbin_gene(immune_HvP_Neut,
                          type = "data", 
                          gene = "FCGR3B", 
                          action = "mean") + 
  scale_fill_viridis_c(option = "A") +
  theme_cowplot()

neut4 <- plot_hexbin_gene(immune_HvP_Neut, 
                          type = "data", 
                          gene = "CXCR2", 
                          action = "mean") + 
  scale_fill_viridis_c(option = "A") +
  theme_cowplot()
Fig5D <- (neut1|neut2|neut3|neut4)
```

### Endothelial subset

#### Cell classification

``` r
Idents(oralIntegrated_HvP) <- "fineCellType"
# subset endo
endo_HvP <- subset(oralIntegrated_HvP, 
                   idents = c("P.VEC 1.1", 
                              "P.VEC 1.2", 
                              "P.VEC 1.3", 
                              "P.VEC 1.4",
                              "P.SMC", 
                              "P.LEC"))
levels(endo_HvP) <- c("P.VEC 1.1", 
                      "P.VEC 1.2", 
                      "P.VEC 1.3", 
                      "P.VEC 1.4",
                      "P.SMC", 
                      "P.LEC")
# save metadata
endo_HvP$endoTypes <- Idents(endo_HvP)

# UMAP plotting
endo_HvP <- RunUMAP(endo_HvP, 
                    dims = 1:50)
endo_HvPColors <- c(carto.pal(pal1 = "red.pal", n1 = 20)[c(20,8,2,14)],
                    carto.pal(pal1 = "pink.pal", n1 = 10)[c(2,9)])
```

##### Figure 6A

``` r
Idents(endo_HvP) <- "endoTypes"
Fig6A_dimPlot <- DimPlot(object = endo_HvP, 
                 label = FALSE, 
                 reduction = "umap", 
                 cols = endo_HvPColors) + 
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
Fig6A_prop <- plot_stat(endo_HvP, 
                        plot_type = "prop_fill",
                        group_by = "project") +
  scale_fill_manual(values = endo_HvPColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12)) +
  guides(fill = FALSE)

# stats
Idents(endo_HvP) <- "endoTypes"
Fig6A_propBySample <- plot_stat(endo_HvP, 
                                plot_type = "prop_fill", 
                                group_by = "orig.ident") 
Fig6AStats <- as.data.frame(Fig6A_propBySample[["data"]])
Fig6AStats <- Fig6AStats[order(Fig6AStats$cluster),]
write.csv(Fig6AStats, "Fig6A_forStats.csv")
```

#### Gene expression

``` r
endo_HvPMarkers <- FindAllMarkers(endo_HvP, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.75)
endo_HvPSig <- endo_HvPMarkers[endo_HvPMarkers$p_val_adj < 0.2, ]
write.csv(endo_HvPSig, "endo_HvP_clusterMarkers_significant.csv")
endo_HvPTop <- endo_HvPSig %>%
  group_by(cluster) %>%
  top_n(n=4, wt=avg_logFC)
endo_HvPTop <- endo_HvPTop[!duplicated(endo_HvPTop$gene),] %>%
  map_df(rev)
```

##### Figure S5A,B

``` r
Idents(endo_HvP) <- "endoTypes"
FigureS5A <- DimPlot(object = endo_HvP, label = FALSE, reduction = "umap", cols = endo_HvPColors, split.by = "project") + 
  theme_cowplot() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
FigureS5B <- DotPlot(endo_HvP, 
                     features = rev(c("CLU", "ACKR1", "CD74", "HLA-DRA",
                                      "SELE", "SELP", "CSF2RB",  
                                      "RGCC", "PLVAP", "INSR", "FLT1",
                                      "IGFBP3", "SAT1", "CXCL12", "CLDN5", 
                                      "TAGLN", "ACTA2", "MYL9", "RGS5",
                                      "CCL21", "TFF3", "FABP4", "EFEMP1")), 
                     col.min = 0, 
                     scale.by = "size", 
                     dot.min = 0.01) +
  coord_flip() + scale_color_gradient2(low = "grey90", high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", order = 1), 
         size = guide_legend(title = "%cell.exp", order = 2))
```

#### Pathway analysis

``` r
Idents(endo_HvP) <- "project"
# pathway analysis
endo_HvP.expressed <- getExpressedGenesFromSeuratObject(endo_HvP, 
                                                        unique(endo_HvP@active.ident), 
                                                        min.pct = 0.25)
annotation <- fetchAnnotation(species = "hs")
endo_HvPMarkers.gsf <- FindAllMarkers(endo_HvP, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)
endo_HvPMarkers.gsf$entrezID <- as.character(annotation$entrez_id[match(endo_HvPMarkers.gsf$gene, 
                                                                        annotation$gene_name)])
endo_HvPMarkers.gsf <- endo_HvPMarkers.gsf[!is.na(endo_HvPMarkers.gsf$entrezID),]
background_entrez <- as.character(annotation$entrez_id[match(endo_HvP.expressed, 
                                                             annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]
endo_HvPMarkers.gsf.filtered <- endo_HvPMarkers.gsf[endo_HvPMarkers.gsf$p_val_adj < 0.05,]
endo_HvPGO <- runGO.all(results=endo_HvPMarkers.gsf.filtered,
                        species = "hs",
                        background_ids = background_entrez,
                        gene_id_col="entrezID",
                        gene_id_type="entrez",
                        sample_col="cluster",
                        p_col="p_val_adj",
                        p_threshold=0.01)
endo_HvPGO.filtered <- endo_HvPGO[endo_HvPGO$ontology=="BP",]
endo_HvPGO.filtered <- filterGenesets(endo_HvPGO.filtered,
                                      min_foreground_genes = 2,
                                      max_genes_geneset = 500,
                                      min_odds_ratio = 2,
                                      p_col = "p.val",
                                      padjust_method = "BH",
                                      use_adjusted_pvalues = TRUE,
                                      pvalue_threshold = 0.05)
endo_HvPGO.filtered.PD <- endo_HvPGO.filtered[endo_HvPGO.filtered$cluster == "PD",]

endo_HvPpathways <- c("cell adhesion mediated by integrin",
                      "leukocyte cell-cell adhesion", 
                      "chemokine-mediated signaling pathway")

endo_HvPGO.filtered.PD <- endo_HvPGO.filtered.PD[endo_HvPGO.filtered.PD$description %in% 
                                                   endo_HvPpathways,]
```

##### Figure 6D, Table S6

``` r
Fig6D <- ggplot(endo_HvPGO.filtered.PD, 
                aes(x=-log(p.adj), 
                    y=description)) + 
  geom_point(data=endo_HvPGO.filtered.PD, 
             aes(size=n_fg), 
             shape = 20) +
  theme_cowplot() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 6))
write.csv(endo_HvPGO.filtered.PD, "TableS6_endo.csv")
```

### Epithelial subset

#### Cell classification

``` r
Idents(oralIntegrated_HvP) <- "fineCellType"
# subset epi
epi_HvP <- subset(oralIntegrated_HvP, 
                  idents = c("P.Epi 1", 
                             "P.Epi 2", 
                             "P.Epi 3", 
                             "P.Mel"))
levels(epi_HvP) <- c("P.Epi 1", 
                     "P.Epi 2", 
                     "P.Epi 3", 
                     "P.Mel")

# save metadata
epi_HvP$epiTypes <- Idents(epi_HvP)

# UMAP plotting
epi_HvP <- RunUMAP(epi_HvP, 
                   dims = 1:50)
epi_HvPColors <- carto.pal(pal1 = "brown.pal", n1 = 20)[c(20, 15, 10, 5)]
```

##### Figure 6B

``` r
Idents(epi_HvP) <- "epiTypes"
Fig6B_dimPlot <- DimPlot(object = epi_HvP, 
                         label = FALSE, 
                         reduction = "umap", 
                         cols = epi_HvPColors) + 
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
Fig6B_prop <- plot_stat(epi_HvP, 
                        plot_type = "prop_fill", 
                        group_by = "project") +
  scale_fill_manual(values = epi_HvPColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12)) +
  guides(fill = FALSE)

# stats
Idents(epi_HvP) <- "epiTypes"
Fig6B_propBySample <- plot_stat(epi_HvP, 
                                plot_type = "prop_fill", 
                                group_by = "orig.ident") 
Fig6BStats <- as.data.frame(Fig6B_propBySample[["data"]])
Fig6BStats <- Fig6BStats[order(Fig6BStats$cluster),]
write.csv(Fig6BStats, "Fig6B_forStats.csv")
```

#### Gene expression

``` r
epi_HvPMarkers <- FindAllMarkers(epi_HvP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75)
epi_HvPSig <- epi_HvPMarkers[epi_HvPMarkers$p_val_adj < 0.2, ]
write.csv(epi_HvPSig, "epi_HvP_clusterMarkers_significant_R_1.csv")
epi_HvPTop <- epi_HvPSig %>%
  group_by(cluster) %>%
  top_n(n=5, wt=avg_logFC)
epi_HvPTop <- epi_HvPTop[!duplicated(epi_HvPTop$gene),]
epi_HvPTop <- rev(epi_HvPTop$gene)
```

##### Figure S6C,D

``` r
Idents(epi_HvP) <- "epiTypes"

FigS6C <- DimPlot(object = epi_HvP, 
                  label = FALSE,
                  reduction = "umap", 
                  cols = epi_HvPColors, 
                  split.by = "project") + 
  theme_cowplot() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
FigS6D <- DotPlot(epi_HvP, 
                  features = rev(c("KRT14", "KRT5", "SFN", "KRT15", "DSP",
                                   "IL36G", "IL1A", "TGM1", "ANXA1", "SLPI",
                                   "CRNN", "CNFN", "TGM3", "SPRR3", "SBSN",
                                   "MITF", "DCT", "PMEL", "TYR", "MLANA")), 
                  col.min = 0, 
                  scale.by = "size", 
                  dot.min = 0.01) +
  coord_flip() + scale_color_gradient2(low = "grey90", high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", order = 1), 
         size = guide_legend(title = "%cell.exp", order = 2))
```

#### Pathway analysis

``` r
Idents(epi_HvP) <- "project"

epi_HvP.expressed <- getExpressedGenesFromSeuratObject(epi_HvP, 
                                                       unique(epi_HvP@active.ident), 
                                                       min.pct = 0.25)
annotation <- fetchAnnotation(species = "hs")
epi_HvPMarkers.gsf <- FindAllMarkers(epi_HvP,
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25)
epi_HvPMarkers.gsf$entrezID <- as.character(annotation$entrez_id[match(epi_HvPMarkers.gsf$gene, 
                                                                       annotation$gene_name)])
epi_HvPMarkers.gsf <- epi_HvPMarkers.gsf[!is.na(epi_HvPMarkers.gsf$entrezID),]
background_entrez <- as.character(annotation$entrez_id[match(epi_HvP.expressed,
                                                             annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]
epi_HvPMarkers.gsf.filtered <- epi_HvPMarkers.gsf[epi_HvPMarkers.gsf$p_val_adj < 0.05,]
epi_HvPGO <- runGO.all(results=epi_HvPMarkers.gsf.filtered,
                       species = "hs",
                       background_ids = background_entrez,
                       gene_id_col="entrezID",
                       gene_id_type="entrez",
                       sample_col="cluster",
                       p_col="p_val_adj",
                       p_threshold=0.01)

epi_HvPGO.filtered <- epi_HvPGO[epi_HvPGO$ontology=="BP",]
epi_HvPGO.filtered <- filterGenesets(epi_HvPGO.filtered,
                                     min_foreground_genes = 3,
                                     max_genes_geneset = 500,
                                     min_odds_ratio = 2,
                                     p_col = "p.val",
                                     padjust_method = "BH",
                                     use_adjusted_pvalues = TRUE,
                                     pvalue_threshold = 0.05)
epi_HvPGO.filtered.PD <- epi_HvPGO.filtered[epi_HvPGO.filtered$cluster == "PD",]

epi_HvPpathways <- c("antimicrobial humoral immune response mediated by antimicrobial peptide",
                     "response to lipopolysaccharide",
                     "response to molecule of bacterial origin")

epi_HvPGO.filtered.PD <- epi_HvPGO.filtered.PD[epi_HvPGO.filtered.PD$description %in% 
                                                 epi_HvPpathways,]
```

##### Figure 6E, Table S6

``` r
Fig6E <- ggplot(epi_HvPGO.filtered.PD, 
                aes(x=-log(p.adj), 
                    y=description)) + 
  geom_point(data=epi_HvPGO.filtered.PD, 
             aes(size=n_fg), 
             shape = 20) +
  theme_cowplot() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 4))
write.csv(epi_HvPGO.filtered.PD, "TableS6_epi.csv")
```

### Fibroblast subset

#### Cell classification

``` r
# subset fib
# fib 1.5 was contaminated with epi cells and was not included
fib_HvP <- subset(oralIntegrated_HvP, idents = c("P.Fib 1.1", 
                                                 "P.Fib 1.2", 
                                                 "P.Fib 1.3", 
                                                 "P.Fib 1.4"))
levels(fib_HvP) <- c("P.Fib 1.1", 
                     "P.Fib 1.2",
                     "P.Fib 1.3", 
                     "P.Fib 1.4")

# save metadata
fib_HvP$fibTypes <- fib_HvP@active.ident

# UMAP plotting
fib_HvP <- RunUMAP(fib_HvP, dims = 1:50)
fib_HvPColors <- c(paletteer_d("ggsci::green_material")[c(10, 7, 4)], 
                   paletteer_d("ggsci::lime_material")[c(4,8)])
```

##### Figure 6C

``` r
Idents(fib_HvP) <- "fibTypes"
Fig6C_dimPlot <- DimPlot(object = fib_HvP, 
                         label = FALSE, 
                         reduction = "umap", 
                         cols = fib_HvPColors) + 
  theme_void() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
Fig6C_prop <- plot_stat(fib_HvP, 
                        plot_type = "prop_fill", 
                        group_by = "project") +
  scale_fill_manual(values = fib_HvPColors) + 
  theme_void() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12)) +
  guides(fill = FALSE)

# stats
Idents(fib_HvP) <- "fibTypes"
Fig6C_propBySample <- plot_stat(fib_HvP, 
                                plot_type = "prop_fill", 
                                group_by = "orig.ident") 
Fig6CStats <- as.data.frame(Fig6C_propBySample[["data"]])
Fig6CStats <- Fig6CStats[order(Fig6CStats$cluster),]
write.csv(Fig6CStats, "Fig6C_forStats.csv")
```

#### Gene expression

``` r
fib_HvPMarkers <- FindAllMarkers(fib_HvP, 
                                 only.pos = TRUE, 
                                 min.pct = 0.5, 
                                 logfc.threshold = 0.25)
fib_HvPSig <- fib_HvPMarkers[fib_HvPMarkers$p_val_adj < 0.2, ]
fib_HvPTop <- fib_HvPSig %>%
  group_by(cluster) %>%
  top_n(n=4, wt=avg_logFC)
fib_HvPTop <- fib_HvPTop[!duplicated(fib_HvPTop$gene),] %>%
  map_df(rev)
```

##### Figure S6E,F

``` r
Idents(fib_HvP) <- "fibTypes"
FigS6E <- DimPlot(object = fib_HvP, 
                  label = FALSE, 
                  reduction = "umap", 
                  cols = fib_HvPColors, 
                  split.by = "project") + 
  theme_cowplot() +
  theme(legend.text = element_text(size = 11), 
        legend.justification = c(1,0))
FigS6F <- DotPlot(fib_HvP, 
                  features = rev(c("CXCL2", "CXCL13", "CXCL1", "PHLDA1", 
                                   "APCDD1", "IGFBP2", "MRPS6", "RORB",
                                   "CFD", "APOD", "GSN", "ABCA8",
                                   "SRFP4", "COL11A1", "TIMP3", "ASPN")),
                  col.min = 0, 
                  scale.by = "size", 
                  dot.min = 0.01) +
  coord_flip() + 
  scale_color_gradient2(low = "grey90", high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", order = 1), 
         size = guide_legend(title = "%cell.exp", order = 2))
```

#### Pathway analysis

``` r
Idents(fib_HvP) <- "project"

fib_HvP.expressed <- getExpressedGenesFromSeuratObject(fib_HvP, 
                                                       unique(fib_HvP@active.ident), 
                                                       min.pct = 0.25)
annotation <- fetchAnnotation(species = "hs")
fib_HvPMarkers.gsf <- FindAllMarkers(fib_HvP, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25)
fib_HvPMarkers.gsf$entrezID <- as.character(annotation$entrez_id[match(fib_HvPMarkers.gsf$gene, 
                                                                       annotation$gene_name)])
fib_HvPMarkers.gsf <- fib_HvPMarkers.gsf[!is.na(fib_HvPMarkers.gsf$entrezID),]
background_entrez <- as.character(annotation$entrez_id[match(fib_HvP.expressed, 
                                                             annotation$gene_name)])
background_entrez <- background_entrez[!is.na(background_entrez)]

fib_HvPMarkers.gsf.filtered <- fib_HvPMarkers.gsf[fib_HvPMarkers.gsf$p_val_adj < 0.05,]
fib_HvPGO <- runGO.all(results=fib_HvPMarkers.gsf.filtered,
                       species = "hs",
                       background_ids = background_entrez,
                       gene_id_col="entrezID",
                       gene_id_type="entrez",
                       sample_col="cluster",
                       p_col="p_val_adj",
                       p_threshold=0.01)

fib_HvPGO.filtered <- fib_HvPGO[fib_HvPGO$ontology=="BP",]
fib_HvPGO.filtered <- filterGenesets(fib_HvPGO.filtered,
                                     min_foreground_genes = 3,
                                     max_genes_geneset = 500,
                                     min_odds_ratio = 2,
                                     p_col = "p.val",
                                     padjust_method = "BH",
                                     use_adjusted_pvalues = TRUE,
                                     pvalue_threshold = 0.01)

fib_HvPGO.filtered.PD <- fib_HvPGO.filtered[fib_HvPGO.filtered$cluster == "PD",]

fib_HvPpathways <- c("cellular response to molecule of bacterial origin",
                     "response to lipopolysaccharide", 
                     "cytokine biosynthetic process")
fib_HvPGO.filtered.PD <- fib_HvPGO.filtered.PD[fib_HvPGO.filtered.PD$description %in% 
                                                 fib_HvPpathways,]
```

##### Figure 6F, Table S6

``` r
Fig6F <- ggplot(fib_HvPGO.filtered.PD, 
                aes(x=-log(p.adj), 
                    y=description)) + 
  geom_point(data=fib_HvPGO.filtered.PD, 
             aes(size=n_fg), 
             shape = 20) +
  theme_cowplot() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 6))
write.csv(fib_HvPGO.filtered.PD, "TableS6_fib.csv")
```

### Receptor-Ligand Interactions

``` r
oralIntegrated_HvP_withoutOther <- oralIntegrated_HvP
Idents(oralIntegrated_HvP) <- "celltype"
oralIntegrated_HvP_withoutOther <- subset(oralIntegrated_HvP, idents = "Other", invert = TRUE)

# ligand-target prior model
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

# ligand-receptor network
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

# weighted networks
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

# seurat wrapper from NicheNet
nichenet_output_all = nichenet_seuratobj_aggregate(
  seurat_obj = oralIntegrated_HvP_withoutOther,
  receiver = "Immune",
  condition_colname = "project", condition_oi = "PD", condition_reference = "GM",
  sender = c("Endothelial", "Fibroblast", "Epithelial", "Immune"),
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
  weighted_networks = weighted_networks, organism = "human")

# Set the order of identities to be graphed
Idents(oralIntegrated_HvP_withoutOther) <- "celltype.site"
levels(oralIntegrated_HvP_withoutOther) <- rev(c("Endothelial_GM", "Endothelial_PD", "Epithelial_GM", "Epithelial_PD",
                                    "Fibroblast_GM", "Fibroblast_PD", "Immune_GM", "Immune_PD"))
Idents(immune_HvP) <- "mon.fine.adjusted"

ligands <- DotPlot(oralIntegrated_HvP_withoutOther, 
                   features = nichenet_output_all$top_ligands %>% 
                     rev(), 
                   col.min = 0, 
                   scale.by = "size", 
                   dot.min = 0.01) +
  RotatedAxis() + 
  scale_color_gradient2(low = "grey90", 
                        high = "black")
receptors <- DotPlot(immune_HvP, 
                     features = nichenet_output_all$top_receptors %>% 
                       rev(), 
                     col.min = 0, 
                     scale.by = "size", 
                     dot.min = 0.01) +
  RotatedAxis() + 
  scale_color_gradient2(low = "grey90", 
                        high = "black")
```

##### Figure S6G-I

``` r
FigS6G <- ligands
FigS6H <- receptors
FigS6I <- nichenet_output_all$ligand_receptor_heatmap
```

##### Figure 6E

``` r
# Alluvial plot with truncated list to include only relevant ones
colsAlluvial <- paletteer_d("ggsci::nrc_npg")[c(4,1,3,9)]
# user generated based on L/R prior interaction potential (FigS6I)
LR.alluvial.truncated <- read.csv("LR_alluvial_truncated.csv", sep = ",",header = TRUE)

Fig6E <- alluvial_wide(data = LR.alluvial.truncated,
                       max_variables = 3,
                       fill_by = 'first_variable', 
                       order_levels = c("ITGA1", 
                                        "ITGA6", 
                                        "ITGB4", 
                                        "ITGB2", 
                                        "ITGB5", 
                                        "SDC1", 
                                        "SDC4", 
                                        "CCR6", 
                                        "CCR7", 
                                        "C3AR1", 
                                        "LRP1", 
                                        "FPR1"),
                       col_vector_flow = colsAlluvial[c(2,4,3,1)]) + 
  theme_cowplot()
```

### Cyto-/Chemokine Interactions

``` r
# gene lists
chemoList <- c("CCL4", 
               "CCL8", 
               "CCL17", 
               "CCL26", 
               "CXCL1", 
               "CXCL2", 
               "CXCL3", 
               "CXCL5", 
               "CXCL8", 
               "CCL5", 
               "CCL23", 
               "CXCL16", 
               "CCL28", 
               "CXCL12",
               "CXCL13", 
               "CCL2", 
               "CX3CL1", 
               "CXCL9", 
               "CXCL10", 
               "CXCL11",
               "CCL3", 
               "CCL11", 
               "CCL13", 
               "CCL14", 
               "CCL19", 
               "CCL20", 
               "CCL21")
chemoList2 <-c("CXCL1", 
               "CXCL2", 
               "CXCL3", 
               "CXCL5",
               "CXCL6", 
               "CXCL8",
              "CXCL9", 
              "CXCL10", 
              "CXCL11", 
              "CXCL12", 
              "CXCL13",
              "CXCL16", 
              "XCL1", 
              "XCL2", 
              "CCL1", 
              "CCL2", 
              "CCL3",
              "CCL4", 
              "CCL5", 
              "CCL8", 
              "CCL11",
              "CCL13", 
              "CCL14",
              "CCL15", 
              "CCL16", 
              "CCL17", 
              "CCL19",
              "CCL20", 
              "CCL21", 
              "CCL23", 
              "CCL24", 
              "CCL25", 
              "CCL26", 
              "CCL27", 
              "CCL28",
              "CX3CL1")
chemoColors <- c(rep(c("black"), 7), 
                 rep(c("#003399"), 2),
                 rep(c("black"), 4), 
                 rep(c("#663300"), 2), 
                 rep(c("black"), 2), 
                 rep(c("#cc3300"), 4), 
                 rep(c("black"), 1), 
                 rep(c("#cc3300"), 4), 
                 rep(c("black"), 1),
                 rep(c("chartreuse4"), 3), 
                 rep(c("#cc3300"), 1), 
                 rep(c("black"), 4),
                 rep(c("#CC9900"), 1), 
                 rep(c("#660099"), 1))
chemoRList <- c("CXCR1", 
                "CXCR2", 
                "CXCR3", 
                "CXCR4",
                "CXCR5", 
                "CXCR6", 
                "XCR1", 
                "CCR1", 
                "CCR2",
                "CCR3", 
                "CCR4",
                "CCR5",
                "CCR6", 
                "CCR7", 
                "CCR8", 
                "CCR9",
                "CCR10", 
                "CX3CR1")

chemoRColors <- rev(c(paletteer_d("ggsci::default_jco"), 
                      paletteer_d("ggsci::default_jama")))

Idents(immune_HvP) <- "monFine_project"
levels(immune_HvP) <- rev(c("abT (CD4)_GM", "abT (CD4)_PD", 
                            "Th17_GM", "Th17_PD", 
                            "MAIT_GM", "MAIT_PD", 
                            "abT (CD8)_GM", "abT (CD8)_PD", 
                            "gd T_GM", "gd T_PD", 
                            "Treg_GM", "Treg_PD", 
                            "NK_GM", "NK_PD",
                            "pDC_GM", "pDC_PD", 
                            "B_GM", "B_PD", 
                            "Plasma_GM", "Plasma_PD",
                            "Neutrophil_GM", "Neutrophil_PD", 
                            "Monocyte_GM", "Monocyte_PD", 
                            "mDC_GM", "mDC_PD", 
                            "Mast_GM", "Mast_PD"))

immune_HvP$monFine_project <- immune_HvP@active.ident
Idents(oralIntegrated_HvP) <- "celltype.site"
oralIntegrated_sub3 <- subset(oralIntegrated_HvP, 
                              idents = c("Endothelial_PD",
                                         "Fibroblast_PD", 
                                         "Epithelial_PD",
                                         "Immune_PD"))
levels(oralIntegrated_sub3) <- rev(c("Endothelial_PD", 
                                     "Fibroblast_PD", 
                                     "Epithelial_PD", 
                                     "Immune_PD"))
```

##### Figure S6J

``` r
FigS6J_1 <- DotPlot(immune_HvP, features = rev(chemoRList), col.min = 0,
                      scale.by = "size", 
                     cols = c("white", "black"), 
                     dot.min = 0.01, dot.scale = 3) + 
  theme(axis.title = element_blank(),
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5,
                                   color = rev(chemoRColors))) + 
  coord_flip()

levels(oralIntegrated_HvP) <- c("Endothelial_GM", "Endothelial_PD", "Fibroblast_GM", "Fibroblast_PD", 
                                "Epithelial_GM", "Epithelial_PD", "Immune_GM", "Immune_PD")
FigS6J_2 <- DotPlot(oralIntegrated_HvP, 
                     features = rev(chemoList2), 
                     col.min = 0,
                    scale.by = "size", 
                    cols = c("white", "black"), 
                    dot.min = 0.01) + 
  theme(axis.title = element_blank(),
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90))+
  coord_flip()
```

##### Figure 6F

``` r
Fig6F <- DotPlot(oralIntegrated_sub3, features = chemoList, col.min = 0,
                    scale.by = "size", cols = c("white", "black"), dot.min = 0.01, dot.scale = 10) + 
  scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 90))
```

### Monogenic/GWAS

``` r
perioMonogenic <- c("C1R", 
                    "C1S", 
                    "HAX1", 
                    "LYST", 
                    "CXCR4", 
                    "ITGB2", 
                    "CTSC", 
                    "FPR1")
perioGWAS <- c("SIGLEC5", 
               "AIM2", 
               "NIN", 
               "TENM2",
               "FBXO38", 
               "TSNAX", 
               "DISC1")
oralIntegrated_sub <- subset(oralIntegrated_HvP, 
                             idents = c("Immune_GM", "Immune_PD"), 
                             invert = TRUE)
oralIntegrated_subIm <- merge(oralIntegrated_sub, 
                              y=immune_HvP)
levels(oralIntegrated_subIm) <- rev(c("Endothelial_GM", "Endothelial_PD", 
                                      "Fibroblast_GM", "Fibroblast_PD", 
                                      "Epithelial_GM", "Epithelial_PD", 
                                      "abT (CD4)_GM", "abT (CD4)_PD", 
                                      "Th17_GM", "Th17_PD", 
                                      "MAIT_GM", "MAIT_PD", 
                                      "abT (CD8)_GM", "abT (CD8)_PD", 
                                      "gd T_GM", "gd T_PD", 
                                      "Treg_GM", "Treg_PD", 
                                      "NK_GM", "NK_PD",
                                      "pDC_GM", "pDC_PD", 
                                      "B_GM", "B_PD", 
                                      "Plasma_GM", "Plasma_PD",
                                      "Neutrophil_GM", "Neutrophil_PD", 
                                      "Monocyte_GM", "Monocyte_PD", 
                                      "mDC_GM", "mDC_PD", 
                                      "Mast_GM", "Mast_PD"))
```

##### Figure 7A

``` r
Fig7A <- DotPlot(oralIntegrated_subIm, 
                 features = c(perioMonogenic, 
                              perioGWAS), 
                 col.min = 0,
                 scale.by = "size",
                 dot.min = 0.01, 
                 dot.scale = 3) + 
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        legend.position = "right",
        text = element_text(size = 9),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 11)
  )
```

### Microbe Recognition/Entry

``` r
tlrs <- c("TLR1", 
          "TLR2", 
          "TLR3", 
          "TLR4", 
          "TLR5", 
          "TLR6", 
          "TLR7", 
          "TLR8", 
          "TLR9", 
          "TLR10")
nlrs <- c("CIITA", 
          "NAIP", 
          "NOD1", 
          "NLRC4", 
          "NOD2",
          "NLRC3",
          "NLRC5", 
          "NLRX1",
          "NLRP1", 
          "NLRP2", 
          "NLRP3", 
          "NLRP4", 
          "NLRP5", 
          "NLRP6", 
          "NLRP7", 
          "NLRP8", 
          "NLRP9", 
          "NLRP10", 
          "NLRP11", 
          "NLRP12", 
          "NLRP13", 
          "NLRP14")  
clrs <- c("CLEC7A", 
          "OLR1", 
          "CLEC9A", 
          "CLEC2A", 
          "CLEC2B", 
          "CD69", 
          "CLEC6A", 
          "CLEC4D", 
          "CLEC4E", 
          "CLEC4C", 
          "CLEC4A", 
          "CLEC12A", 
          "CLEC1A", 
          "IFI16", 
          "LRRFIP1", 
          "MRC1")
rigs <- c("DDX58", 
          "IFIH1", 
          "DHX58")
prrs <- c(tlrs, nlrs, clrs, rigs)
damps <- c("HMGB1", 
           "S100A8", 
           "S100A9", 
           "SAA1", 
           "SAA2", 
           "HMGN1", 
           "IL33", 
           "SAP130", 
           "HSPD1", 
           "HSP90B1")
dampR <- c("AGER", 
           "IRAK4", 
           "TREM1", 
           "P2RX1", 
           "CD24")
oralVirus <- c("ITGA2", "AXL", "TYRO3", "CLEC4G", "F11R", "HSPA1B", "NCAM1", "DAG1", #HSV1&2
               "SDC1", "SDC2", "SDC4", "GPC1", "LN5", "ITGA6", #HPV
               "CD55", #varicella zoster
               "SCARB1", "CLDN9", #epstein barr (mono) 
               "PVR", "TFRC", "LAMP1", #cytomegalo
               "IDE", "EPS15", "HSP1A", "PTX3", #coxsackie
               "CD81", "LDLR", "EGFR", "EPHA2", #Enterovirus
               "CD46", "SLAMF1", "NECTIN4", #measles
               "MOG", #rubella
               "CAV1", "RPSA", "GAS6", "CLDN1", "HAVCR1", "GRK2") #HIV

oralIntegrated_HvP_withoutOther$general.site <- paste(oralIntegrated_HvP_withoutOther$celltype, 
                                         oralIntegrated_HvP_withoutOther$project, sep = "_")
Idents(oralIntegrated_HvP_withoutOther) <- "general.site"
oralIntegrated_HvP_covid <- oralIntegrated_HvP_withoutOther
levels(oralIntegrated_HvP_covid) <- rev(c("Endothelial_GM", "Endothelial_PD",
                                          "Fibroblast_GM", "Fibroblast_PD",
                                          "Epithelial_GM", "Epithelial_PD", 
                                          "Immune_GM", "Immune_PD"))
covidFeats <- c("ACE2", 
                "TMPRSS2", 
                "TMPRSS4",
                "TMPRSS11D", 
                "CTSB", 
                "CTSL", 
                "HNRNPA1", 
                "BSG", 
                "ZCRB1", 
                "TOP3B", 
                "ANPEP", 
                "CLEC4M", 
                "DPP4", 
                "CD209", 
                "FURIN")

generalViR_Dot_HvP <- DotPlot(oralIntegrated_HvP_covid, 
                              features = c(prrs, 
                                           damps, 
                                           dampR, 
                                           covidFeats, 
                                           oralVirus), 
                              col.min = 0, 
                              scale.by = "size", 
                              dot.min = 0.002, 
                              dot.scale = 10) + 
  scale_color_gradient2(low = "grey90", high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1),
         size = guide_legend(title = "%cell.exp", 
                             order = 2))

# find out which factors are differentially regulated in endo cells in disease
Idents(endo_HvP) <- "project"
covidMkrEndo <- FindAllMarkers(endo_HvP, 
                               logfc.threshold = 0, 
                               min.pct = 0, 
                               only.pos = TRUE)
significantCOV_endo <- as.data.frame(covidMkrEndo[covidMkrEndo$gene %in% 
                                                    c(covidFeats),])
significantOralV_endo <- as.data.frame(covidMkrEndo[covidMkrEndo$gene %in% 
                                                      c(oralVirus),])
significantPrr_endo <- as.data.frame(covidMkrEndo[covidMkrEndo$gene %in% 
                                                    prrs,])
significantDamp_endo <- as.data.frame(covidMkrEndo[covidMkrEndo$gene %in% 
                                                     c(damps,dampR),])
significantViR_endo <- rbind(significantCOV_endo, 
                             significantOralV_endo, 
                             significantPrr_endo, 
                             significantDamp_endo)

significantViR_endo_Dot <- DotPlot(endo_HvP, 
                                   features = rownames(significantViR_endo), 
                                   col.min = 0, 
                                   scale.by = "size",
                                   dot.min = 0.002, 
                                   dot.scale = 10) +
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp", 
                             order = 2))

# extract proportions
endoDot_results <- as.data.frame(significantViR_endo_Dot[["data"]])
endoDot_results <- endoDot_results[order(endoDot_results$features.plot),]
endoDot_results <- endoDot_results[endoDot_results$features.plot %in% 
                                     c("TLR4", 
                                       "NLRX1", 
                                       "NAIP", 
                                       "NLRP3", 
                                       "NOD1",
                                       "IL33", 
                                       "SAP130", 
                                       "HSPD1", 
                                       "HSP90B1", 
                                       "AGER"),]
endoDot_results[is.na(endoDot_results)] = 0


# find out which factors are differentially regulated in epi cells in disease
Idents(epi_HvP) <- "project"
covidMkrepi <- FindAllMarkers(epi_HvP, 
                              logfc.threshold = 0,
                              min.pct = 0, 
                              only.pos = TRUE)
significantCOV_epi <- as.data.frame(covidMkrepi[covidMkrepi$gene %in%
                                                  c(covidFeats),])
significantOralV_epi <- as.data.frame(covidMkrepi[covidMkrepi$gene %in% 
                                                    c(oralVirus),])
significantPrr_epi <- as.data.frame(covidMkrepi[covidMkrepi$gene %in%
                                                  prrs,])
significantDamp_epi <- as.data.frame(covidMkrepi[covidMkrepi$gene %in%
                                                   c(damps,dampR),])
significantViR_epi <- rbind(significantCOV_epi, 
                            significantOralV_epi, 
                            significantPrr_epi, 
                            significantDamp_epi)

significantViR_epi_Dot <- DotPlot(epi_HvP, 
                                  features = c(rownames(significantViR_epi)), 
                                  col.min = 0, 
                                  scale.by = "size",
                                  dot.min = 0.002, 
                                  dot.scale = 10) +
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp",
                                 order = 1), 
         size = guide_legend(title = "%cell.exp", 
                             order = 2))
# extract proportion
epiDot_results <- as.data.frame(significantViR_epi_Dot[["data"]])
epiDot_results <- epiDot_results[order(epiDot_results$features.plot),]
epiDot_results <- epiDot_results[epiDot_results$features.plot %in% 
                                   c("S100A9", 
                                     "S100A8", 
                                     "SAA1", 
                                     "SAA2", 
                                     "HSP90B1",
                                     "CD69", 
                                     "CLEC1A", 
                                     "CLEC2B", 
                                     "OLR1", 
                                     "CLEC4D"),]
epiDot_results[is.na(epiDot_results)] = 0

# find out which factors are differentially regulated in fib cells in disease
Idents(fib_HvP) <- "project"
covidMkrFib <- FindAllMarkers(fib_HvP, 
                              logfc.threshold = 0, 
                              min.pct = 0, 
                              only.pos = TRUE)
significantCOV_fib <- as.data.frame(covidMkrFib[covidMkrFib$gene %in% 
                                                  c(covidFeats),])
significantOralV_fib <- as.data.frame(covidMkrFib[covidMkrFib$gene %in% 
                                                    c(oralVirus),])
significantPrr_fib <- as.data.frame(covidMkrFib[covidMkrFib$gene %in%
                                                  prrs,])
significantDamp_fib <- as.data.frame(covidMkrFib[covidMkrFib$gene %in% 
                                                   c(damps,dampR),])
significantViR_fib <- rbind(significantCOV_fib, 
                            significantOralV_fib,
                            significantPrr_fib, 
                            significantDamp_fib)

significantViR_fib_Dot <- DotPlot(fib_HvP, 
                                  features = rownames(significantViR_fib),
                                  col.min = 0, 
                                  scale.by = "size",
                                  dot.min = 0.002, 
                                  dot.scale = 10) +
  scale_color_gradient2(low = "grey90", 
                        high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1), 
         size = guide_legend(title = "%cell.exp",
                             order = 2))

# extract proportions 
fibDot_results <- as.data.frame(significantViR_fib_Dot[["data"]])
fibDot_results <- fibDot_results[order(fibDot_results$features.plot),]
fibDot_results <- fibDot_results[fibDot_results$features.plot %in% 
                                   c("SAA1", 
                                     "SAA2", 
                                     "HSP90B1", 
                                     "HSPD1", 
                                     "SAP130",
                                     "NLRC5", 
                                     "NOD2", 
                                     "NLRP1", 
                                     "TLR2",
                                     "TLR4"),]
fibDot_results[is.na(fibDot_results)] = 0

# find out which factors are differentially regulated in immune cells in disease
Idents(immune_HvP) <- "project"
covidMkrImmune <- FindAllMarkers(immune_HvP, 
                                 logfc.threshold = 0, 
                                 min.pct = 0, 
                                 only.pos = TRUE)
significantCOV_immune <- as.data.frame(covidMkrImmune[covidMkrImmune$gene %in% 
                                                        c(covidFeats),])
significantOralV_immune <- as.data.frame(covidMkrImmune[covidMkrImmune$gene %in% 
                                                          c(oralVirus),])
significantPrr_immune <- as.data.frame(covidMkrImmune[covidMkrImmune$gene %in% 
                                                        prrs,])
significantDamp_immune <- as.data.frame(covidMkrImmune[covidMkrImmune$gene %in% 
                                                         c(damps,dampR),])
significantViR_immune <- rbind(significantCOV_immune, 
                               significantOralV_immune, 
                               significantPrr_immune, 
                               significantDamp_immune)

significantViR_immune_Dot <- DotPlot(immune_HvP, 
                                     features = c(rownames(significantViR_immune)), 
                                     col.min = 0, 
                                     scale.by = "size",
                                     dot.min = 0.002, 
                                     dot.scale = 10) +
  scale_color_gradient2(low = "grey90", high = "black") +
  theme(axis.title = element_blank(),
        text = element_text(size = 9.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5))+
  guides(colour = guide_colorbar(title = "avg.exp", 
                                 order = 1),
         size = guide_legend(title = "%cell.exp", 
                             order = 2))

# extract proportions
immuneDot_results <- as.data.frame(significantViR_immune_Dot[["data"]])
immuneDot_results <- immuneDot_results[order(immuneDot_results$features.plot),]
immuneDot_results <- immuneDot_results[immuneDot_results$features.plot %in% 
                                         c("HSP90B1", 
                                           "HSPD1",
                                           "P2RX1", 
                                           "CD24", 
                                           "HMGN1",
                                           "CLEC2B", 
                                           "IFI16", 
                                           "TLR6", 
                                           "TLR4", 
                                           "IFIH1"),]
immuneDot_results[is.na(immuneDot_results)] = 0
```

##### Figure 7B

``` r
Fig7B_endoHvP <- ggplot(endoDot_results, 
                          aes(fill=id, 
                              y = pct.exp, 
                              x = features.plot)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
Fig7B_epiHvP <- ggplot(epiDot_results, 
                         aes(fill=id, 
                             y = pct.exp, 
                             x = features.plot)) +
  geom_bar(position = "fill",
           stat = "identity") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
Fig7B_fibHvP <- ggplot(fibDot_results, 
                         aes(fill=id, 
                             y = pct.exp, 
                             x = features.plot)) +
  geom_bar(position = "fill",
           stat = "identity") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
Fig7B_immuneHvP <- ggplot(immuneDot_results, 
                            aes(fill=id, 
                                y = pct.exp, 
                                x = features.plot)) +
  geom_bar(position = "fill", 
           stat = "identity") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
```

##### Figure S7B

``` r
FigS7B <- generalViR_Dot_HvP
```

## SessionInfo

``` r
sessionInfo()

R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8
[12] LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] SoupX_1.4.8                 purrr_0.3.4                 easyalluvial_0.2.3          nichenetr_1.0.0             gsfisher_0.2                knitr_1.30                  kableExtra_1.3.1            cartography_2.4.2           Cairo_1.5-12.2              schex_1.0.55                ggcharts_0.2.1             
[12] Scillus_0.3.0               devtools_2.3.2              usethis_1.6.3               pheatmap_1.0.12             viridis_0.5.1               viridisLite_0.3.0           clustree_0.4.3              ComplexHeatmap_2.2.0        R.utils_2.10.1              R.oo_1.24.0                 R.methodsS3_1.8.1          
[23] paletteer_1.2.0             data.table_1.13.4           SingleR_1.0.6               SingleCellSignalR_0.0.1.8   cowplot_1.1.0               scater_1.14.6               ggraph_2.0.4                gdata_2.18.0                scales_1.1.1                reshape2_1.4.4              tidyr_1.1.2                
[34] patchwork_1.1.0             ggplot2_3.3.3.9000          MAST_1.12.0                 SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1 DelayedArray_0.12.3         BiocParallel_1.20.1         matrixStats_0.57.0          Biobase_2.46.0              GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
[45] IRanges_2.20.2              S4Vectors_0.24.4            BiocGenerics_0.32.0         Seurat_3.2.2                dplyr_1.0.4                

loaded via a namespace (and not attached):
  [1] rsvd_1.0.3                    Hmisc_4.4-2                   ica_1.0-2                     class_7.3-17                  ps_1.5.0                      foreach_1.5.1                 lmtest_0.9-38                 rprojroot_2.0.2               crayon_1.3.4                  MASS_7.3-53                  
 [11] backports_1.2.0               nlme_3.1-150                  rlang_0.4.10                  XVector_0.26.0                caret_6.0-86                  ROCR_1.0-11                   irlba_2.3.3                   callr_3.5.1                   limma_3.42.2                  rjson_0.2.20                 
 [21] bit64_4.0.5                   glue_1.4.2                    sctransform_0.3.2             processx_3.4.5                vipor_0.4.5                   AnnotationDbi_1.48.0          tidyselect_1.1.0              fitdistrplus_1.1-3            XML_3.99-0.3                  zoo_1.8-8                    
 [31] org.Mm.eg.db_3.10.0           xtable_1.8-4                  formattable_0.2.0.1           magrittr_2.0.1                evaluate_0.14                 cli_2.2.0                     zlibbioc_1.32.0               rstudioapi_0.13               miniUI_0.1.1.1                sp_1.4-4                     
 [41] rpart_4.1-15                  tinytex_0.27                  shiny_1.5.0                   BiocSingular_1.2.2            xfun_0.19                     askpass_1.1                   clue_0.3-58                   pkgbuild_1.2.0                multtest_2.42.0               cluster_2.1.0                
 [51] caTools_1.18.0                tidygraph_1.2.0               tibble_3.0.5                  interactiveDisplayBase_1.24.0 ggrepel_0.8.2                 listenv_0.8.0                 png_0.1-7                     future_1.20.1                 ipred_0.9-9                   withr_2.3.0                  
 [61] bitops_1.0-6                  ggforce_0.3.2                 plyr_1.8.6                    e1071_1.7-4                   pracma_2.2.9                  pROC_1.16.2                   pillar_1.4.7                  gplots_3.1.1                  GlobalOptions_0.1.2           fs_1.5.0                     
 [71] GetoptLong_1.0.4              DelayedMatrixStats_1.8.0      vctrs_0.3.6                   ellipsis_0.3.1                generics_0.1.0                lava_1.6.8.1                  tools_3.6.3                   foreign_0.8-76                beeswarm_0.2.3                entropy_1.2.1                
 [81] munsell_0.5.0                 tweenr_1.0.1                  fastmap_1.0.1                 compiler_3.6.3                pkgload_1.1.0                 abind_1.4-5                   httpuv_1.5.4                  ExperimentHub_1.12.0          sessioninfo_1.1.1             plotly_4.9.2.9000            
 [91] rgeos_0.5-5                   GenomeInfoDbData_1.2.2        prodlim_2019.11.13            gridExtra_2.3                 edgeR_3.28.1                  lattice_0.20-41               deldir_0.2-3                  visNetwork_2.0.9              later_1.1.0.1                 BiocFileCache_1.10.2         
[101] recipes_0.1.15                jsonlite_1.7.2                pbapply_1.4-3                 lazyeval_0.2.2                promises_1.1.1                spatstat_1.64-1               latticeExtra_0.6-29           goftest_1.2-2                 checkmate_2.0.0               spatstat.utils_1.17-0        
[111] reticulate_1.18               rmarkdown_2.5                 webshot_0.5.2                 Rtsne_0.15                    forcats_0.5.0                 uwot_0.1.9                    igraph_1.2.6                  survival_3.2-7                yaml_2.2.1                    htmltools_0.5.0              
[121] memoise_1.1.0                 locfit_1.5-9.4                graphlayouts_0.7.1            digest_0.6.27                 assertthat_0.2.1              mime_0.9                      rappdirs_0.3.1                SIMLR_1.12.0                  RSQLite_2.2.1                 future.apply_1.6.0           
[131] remotes_2.2.0                 blob_1.2.1                    DiagrammeR_1.0.6.1            splines_3.6.3                 Formula_1.2-4                 rematch2_2.1.2                AnnotationHub_2.18.0          RCurl_1.98-1.2                hms_0.5.3                     colorspace_2.0-0             
[141] base64enc_0.1-3               BiocManager_1.30.10           ggbeeswarm_0.6.0              shape_1.4.5                   nnet_7.3-14                   Rcpp_1.0.6                    RANN_2.6.1                    circlize_0.4.11               fansi_0.4.2                   parallelly_1.21.0            
[151] ModelMetrics_1.2.2.2          R6_2.5.0                      ggridges_0.5.2                lifecycle_0.2.0               curl_4.3                      leiden_0.3.6                  testthat_3.0.1                Matrix_1.2-18                 desc_1.2.0                    RcppAnnoy_0.0.17             
[161] org.Hs.eg.db_3.10.0           RColorBrewer_1.1-2            iterators_1.0.13              stringr_1.4.0                 gower_0.2.2                   htmlwidgets_1.5.2.9000        polyclip_1.10-0               biomaRt_2.42.1                rvest_0.3.6                   mgcv_1.8-33                  
[171] globals_0.14.0                openssl_1.4.3                 htmlTable_2.1.0               codetools_0.2-18              lubridate_1.7.9.2             randomForest_4.6-14           gtools_3.8.2                  prettyunits_1.1.1             dbplyr_1.3.0                  RSpectra_0.16-0              
[181] gtable_0.3.0                  DBI_1.1.1                     ggalluvial_0.12.2             tensor_1.5                    httr_1.4.2                    KernSmooth_2.23-18            stringi_1.5.3                 progress_1.2.2                farver_2.0.3                  fdrtool_1.2.15               
[191] hexbin_1.28.1                 timeDate_3043.102             xml2_1.3.2                    BiocNeighbors_1.4.2           readr_1.4.0                   BiocVersion_3.10.1            bit_4.0.4                     jpeg_0.1-8.1                  spatstat.data_1.5-2           pkgconfig_2.0.3     
```
