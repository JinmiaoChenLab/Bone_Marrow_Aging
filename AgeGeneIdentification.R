library(Seurat)
library(ktools)
library(pbmcapply)
library(ComplexHeatmap)
library(randomForest)
library(stringr)
library(plyr)
library(ggpubr)
library(glmnet)
library(caret)
library(enrichR)
library(e1071)
library(dplyr)
library(data.table)
options(bitmapType = "cairo")

rna.data = readRDS("rna_data_annotated.rds")
meta = read.csv("metadata.txt", sep = "\t")
rownames(meta) = meta$sample_id
meta = meta[rna.data$sample,]
rownames(meta) = Cells(rna.data)
rna.data = AddMetaData(rna.data, meta)
rna.data = subset(rna.data, age > 2)

DefaultAssay(rna.data) = "RNA"

a = data.frame(table(rna.data$sample, rna.data$ct_rename))
a = a[which(a$Freq > 20),]
a = data.frame(table(a$Var2))
a = a[which(a$Freq > 20),]

meta = read.csv("metadata.txt", sep = "\t")
rownames(meta) = meta$sample_id
a$Var1 = as.character(a$Var1)

select.feature.integrated = pbmclapply(
  a$Var1, function(jj) {
    rna.sub = subset(rna.data, ct_rename == jj)
    DefaultAssay(rna.sub) = "raw"
    
    nFeature = rowSums(x = GetAssayData(object = rna.sub, slot = "data") > 0) 
    features = rownames(rna.sub)[which(nFeature > ncol(rna.sub)/10)]
    
    rna.bind = pbmclapply(
      1:50, function(i) {
        rna = readRDS(paste0("./FastIntegrationTmp/inte/inte_", i, ".rds"))
        rna = rna[,Cells(rna.sub)]
        rna = rna[intersect(rownames(rna), features),]
        return(rna)
      }, mc.cores = 20
    )
    rna.bind = do.call(rbind, rna.bind)
    rna.sub[["RNA"]] = CreateAssayObject(rna.bind)
    
    rna.avergae = AverageExpression(rna.sub, group.by = "sample")
    rna.avergae = rna.avergae$RNA
    
    age = meta[colnames(rna.avergae), "age"]
    
    cor.res = apply(rna.avergae, 1, function(i){
      a = cor.test(i, age, method = "spearman")
      return(c(a$estimate, a$p.value))
    })
    res = data.frame(gene = rownames(rna.avergae), cor = cor.res[1,], cell.type = jj, p = cor.res[2,])
    res$q = p.adjust(res$p)
    return(res)
  }, mc.cores = 5
)
select.feature.integrated = do.call(rbind, select.feature.integrated)
select.feature.integrated$cell.type = as.character(select.feature.integrated$cell.type)

select.feature.integrated = select.feature.integrated[which(select.feature.integrated$p < 0.01),]

saveRDS(select.feature.integrated, "rds/age_cor_integrated.rds", compress = F)

select.feature = select.feature.integrated[which(select.feature.integrated$p < 0.01),]
select.feature = readRDS("rds/age_cor_integrated.rds")
select.feature %>%
  group_by(cell.type) %>%
  top_n(n = 50, wt = -p) -> select.feature


plot.list = pbmclapply(
  unique(select.feature$cell.type), function(ct) {
    feature = select.feature[which(select.feature$cell.type == ct),]
    feature = feature$gene
    dbs <- listEnrichrDbs()
    dbs <- c("Reactome_2022")
    enriched = enrichr(feature, dbs)
    enriched = do.call(rbind, enriched)
    enriched = enriched[which(enriched$Adjusted.P.value < 0.05),]
    enriched$Type = str_match(rownames(enriched), "^(.*?)_\\d+")[,2]
    
    enriched$Term = sub("\\(GO:\\d+\\)", "", enriched$Term)
    enriched$Term = sub("R\\-HSA\\-\\d+", "", enriched$Term)
    enriched$Term = sub("\\(.*\\)", "", enriched$Term)
    
    if(nrow(enriched) == 0) {
      return(NULL)
    }
    enriched = enriched[order(enriched$P.value),]
    enriched = enriched[1:min(nrow(enriched), 10),]
    
    plot.data = data.frame(term = enriched$Term, logp = -log(enriched$P.value), type = enriched$Type)
    plot.data$ct = ct
    return(plot.data)
  }, mc.cores = 10
)

plot.list = do.call(rbind, plot.list)


plot.data = dcast(plot.list, term ~ ct, value.var = "logp", fill = 0)
rownames(plot.data) = plot.data$term
plot.data = plot.data[,-1]


pdf("plot/Figure_3_pathway_age_gene_Reactome.pdf", width = 16, height = 8)
Heatmap(t(plot.data), name = "Error", column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 7), col = c("#dee3ed", "red"))
dev.off()

library(UpSetR)

select.feature = readRDS("rds/age_cor_integrated.rds")
select.feature = select.feature[which(select.feature$p < 1e-2),]


gene = unique(select.feature$gene)
gene = data.frame(gene = gene)
rownames(gene) = gene$gene

ct = unique(select.feature$cell.type)

for (i in 1:length(ct)) {
  print(i)
  a = select.feature[which(select.feature$cell.type == ct[i]),]
  gene[,i+1] = 0
  gene[a$gene,i+1] = 1
}

colnames(gene)[2:(length(ct) + 1)] = ct

pdf("plot/Figure_3A_upset.pdf", width = 10, height = 6)
upset(gene, sets = colnames(gene)[2:(length(ct) + 1)], mb.ratio = c(0.3, 0.7), show.numbers = F,  
      order.by = "freq")
dev.off()
