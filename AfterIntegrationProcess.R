library(Seurat)
library(pbmcapply)
library(stringr)
library(Matrix)
library(FastIntegration)
options(bitmapType = "cairao")
`%!in%` <- Negate('%in%')

features = readRDS("FastIntegrationTmp/others/features.rds")
rna.data = pbmclapply(
  1:20, function(i) {
    rna = readRDS(paste0("./FastIntegrationTmp/inte/inte_", i, ".rds"))
    rna = rna[intersect(rownames(rna), features),,drop=F]
    return(rna)
  }, mc.cores = 20
)
rna.data = do.call(rbind, rna.data)
rna.data = CreateSeuratObject(rna.data)
rna.data = ScaleData(rna.data, features = features)
rna.data = RunPCA(rna.data, npcs = 30, features = features)
rna.data = RunUMAP(rna.data, dims = 1:30)
rna.data = FindNeighbors(rna.data, dims = 1:30)
rna.data = FindClusters(rna.data, graph.name = "RNA_snn", algorithm = 2)

rna.list = readRDS("rna_list.rds")
idx = split(sample(1:length(rna.list), size = length(rna.list)), 
            cut(1:length(rna.list), round(length(rna.list)/5), labels = FALSE))

rna.bind = pbmcapply::pbmclapply(
  idx, function(i) {
    rna.bind = rna.list[[i[1]]]@assays$RNA@data
    for (j in i[2:length(i)]) {
      rna.bind = cbind(rna.bind, rna.list[[j]]@assays$RNA@data)
    }
    return(rna.bind)
  }, mc.cores = length(idx)
)
rna.bind = do.call(cbind, rna.bind)
rna.data[["raw"]] = CreateAssayObject(rna.bind[,Cells(rna.data)])


rna.data$sample = str_match(Cells(rna.data), "--(.*$)")[,2]
meta = FindSampleByMetadata()
rownames(meta) = meta$sampleId
meta = meta[rna.data$sample,]
rownames(meta) = Cells(rna.data)
rna.data = AddMetaData(rna.data, meta)

rna.data@assays$RNA@counts = matrix(0)
rna.data@assays$RNA@scale.data = matrix(0)
rna.data@assays$raw@counts = matrix(0)
saveRDS(rna.data, "rna_processed.rds", compress = F)
