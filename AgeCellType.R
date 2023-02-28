library(Seurat)
library(pbmcapply)
library(stringr)
library(Matrix)
library(ggplot2)
library(MASS)
library(dplyr)
library(reshape2)
options(bitmapType = "cairo")
`%!in%` <- Negate('%in%')


rna.data = readRDS("rna_data_annotated.rds")
rna.sub = subset(rna.data, project_id %in% c("GSE185381", "GSE120221", "HCA_HematopoieticImmuneCellAtlas"))

plot.data = data.frame(table(rna.sub$sample, rna.sub$ct_rename))
plot.data = dcast(plot.data, Var1 ~ Var2, value.var = "Freq", fill = 0)

rownames(plot.data) = plot.data$Var1
plot.data = plot.data[,-1]
plot.data.pct = sweep(plot.data, 1, rowSums(plot.data),`/`)

res = pbmclapply(
  colnames(plot.data), function(i){
    model.data = data.frame(pct = as.numeric(plot.data.pct[,i]),
                            age = meta[rownames(plot.data),"age"],
                            count = as.numeric(plot.data[,i]), 
                            project = meta[rownames(plot.data),"project_id"])
    # model.data = model.data %>% 
    #   mutate(zscore = (pct - mean(pct))/sd(pct))
    
    # model.data = model.data[which(model.data$zscore < 1.64 & model.data$zscore > -1.64),]
    
    m1 = lm(pct ~ project, data = model.data)
    m2 = lm(pct ~ project + age, data = model.data)
    a = anova(m1, m2)
    return(a$`Pr(>F)`[2])
  }, mc.cores = 54
)
res = data.frame(ct = colnames(plot.data), p = unlist(res))


cor.ct = pbmclapply(
  colnames(plot.data), function(i){
    model.data = data.frame(pct = as.numeric(plot.data.pct[,i]), 
                            count = as.numeric(plot.data[,i]), 
                            age = meta[rownames(plot.data),"age"],
                            project = meta[rownames(plot.data),"project_id"])
    unique(model.data$project)
    
    a1 = cor.test(model.data$pct[which(model.data$project == "GSE185381")], model.data$age[which(model.data$project == "GSE185381")])
    a2 = cor.test(model.data$pct[which(model.data$project == "GSE120221")], model.data$age[which(model.data$project == "GSE120221")])
    a3 = cor.test(model.data$pct[which(model.data$project == "HCA_HematopoieticImmuneCellAtlas")], model.data$age[which(model.data$project == "HCA_HematopoieticImmuneCellAtlas")])
    return(c(a1$estimate, a1$p.value,a2$estimate, a2$p.value,a3$estimate, a3$p.value))
  }, mc.cores = 10
)
cor.ct = do.call(rbind, cor.ct)
res$cor_GSE185381 = cor.ct[,1]
res$p_GSE185381 = cor.ct[,2]
res$cor_GSE120221 = cor.ct[,3]
res$p_GSE120221 = cor.ct[,4]
res$cor_HCA = cor.ct[,5]
res$p_HCA = cor.ct[,6]


write.table(res, "plot/Table_age_cor.txt", sep = "\t", quote = F, col.names = T, row.names = T)


plot.data = sweep(plot.data, 1, rowSums(plot.data),`/`)

model.data = data.frame(pct = as.numeric(plot.data[,"Naive CD8 T cell"]), 
                        age = meta[rownames(plot.data),"age"],
                        project = meta[rownames(plot.data),"project_id"])

pdf("plot/Figure_2_naiveCD8.pdf", width = 6, height = 6)
ggplot(model.data, aes(x=age, y=pct, color=project)) +
  geom_point() + ylab("percentage") + 
  geom_smooth(method=lm) + theme_minimal() + theme(legend.position="bottom") + 
  ggtitle("Naive CD8 T cell")
dev.off()


model.data = data.frame(pct = as.numeric(plot.data[,"Effector/Memory CD4 T cell"]), 
                        age = meta[rownames(plot.data),"age"],
                        project = meta[rownames(plot.data),"project_id"])

pdf("plot/Figure_2_memoryCD4.pdf", width = 6, height = 6)
ggplot(model.data, aes(x=age, y=pct, color=project)) +
  geom_point() + ylab("percentage") + 
  geom_smooth(method=lm) + theme_minimal() + theme(legend.position="bottom") + 
  ggtitle("Effector/Memory CD4 T cell")
dev.off()


model.data = data.frame(pct = as.numeric(plot.data[,"CD16+ NK"]), 
                        age = meta[rownames(plot.data),"age"],
                        project = meta[rownames(plot.data),"project_id"])

pdf("plot/Figure_2_CD16_NK.pdf", width = 6, height = 6)
ggplot(model.data, aes(x=age, y=pct, color=project)) +
  geom_point() + ylab("percentage") + 
  geom_smooth(method=lm) + theme_minimal() + theme(legend.position="bottom") + 
  ggtitle("CD16+ NK")
dev.off()


meta = read.csv("metadata_used.txt", header = T, sep = "\t")
model.data = data.frame(pct = as.numeric(plot.data.pct[,"Naive CD4 T cell"]), 
                        age = meta[rownames(plot.data.pct),"age"],
                        project = meta[rownames(plot.data.pct),"project_id"])

pdf("plot/Naive_CD4.pdf", width = 6, height = 6)
ggplot(model.data, aes(x=age, y=pct, color=project)) +
  geom_point() + ylab("percentage") + 
  geom_smooth(method=lm) + theme_minimal() + theme(legend.position="bottom") + 
  ggtitle("Naive CD4 T cell")
dev.off()

