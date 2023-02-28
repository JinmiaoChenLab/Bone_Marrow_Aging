setwd("/disco_500t/mengwei/backup/BM_aging_revision/github/")
library(Seurat)
library(pbmcapply)
library(scDblFinder)
library(BiocParallel)
library(stringr)
library(FastIntegration)


options(bitmapType = "cairo")
`%!in%` <- Negate('%in%')
meta = FindSampleByMetadata()

meta = meta[which(meta$tissue == "bone marrow"),]
meta = meta[grep("10", meta$platform),]
meta = meta[which(meta$sampleType == "Normal"),]
DownloadDiscoData(meta, dir = "./disco") 


rna.list = pbmclapply(
  1:nrow(meta), function(i) {
    rna = readRDS(paste0("disco/", meta$sampleId[i], ".rds"))
    rna = RenameCells(rna, new.names = paste0(Cells(rna), "--", meta$sampleId[i]))
    rna@assays$RNA@counts = sweep(expm1(rna@assays$RNA@data), 2, rna$nCount_RNA, `*`) /10000
    rna[["percent.mt"]] = PercentageFeatureSet(rna, pattern = "^MT-")
    rna[["percent.rp"]] = PercentageFeatureSet(rna, pattern = "^RP[S|L]")
    sce = scDblFinder(as.SingleCellExperiment(rna))
    rna = subset(rna, cells = Cells(rna)[which(sce@colData@listData[["scDblFinder.class"]] == "singlet")])
    rna = NormalizeData(rna)
    rna = FindVariableFeatures(rna)
    return(rna)
  }, mc.cores = 10
)


rna.list = rna.list[which(unlist(lapply(rna.list, ncol)) > 500)]
rna.list = rna.list[which(unlist(lapply(rna.list, nrow)) == 33538)]

saveRDS(rna.list, "rna_list.rds", compress = F)



BuildIntegrationFile(rna.list = rna.list, tmp.dir = "./", nCores = 30, nfeatures = 3000)
FastFindAnchors(tmp.dir = "./", nCores = 100)

genes = readRDS("FastIntegrationTmp/raw/1.rds")
genes = rownames(genes)
idx = split(1:length(genes), cut(1:length(genes), 20, labels = FALSE))

pbmclapply(
  1:20, function(i) {
    rna.integrated = FastIntegration(tmp.dir = "./", npcs = 1:30, slot = "data",cut.low = -0.1,
                                     features.to.integrate = genes[idx[[i]]])
    saveRDS(rna.integrated, paste0("./FastIntegrationTmp/inte/inte_", i, ".rds"), 
            compress = F)
  }, mc.cores = 20
)

