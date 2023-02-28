library(Seurat)
library(ktools)
library(pbmcapply)
library(randomForest)
library(plyr)
library(ggpubr)
library(glmnet)
library(caret)
library(e1071)
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
a = a[which(a$Freq > 30),]
meta = read.csv("metadata.txt", sep = "\t")
rownames(meta) = meta$sample_id


select.feature = readRDS("rds/age_cor_integrated.rds")
select.feature$cell.type = as.character(select.feature$cell.type)
select.feature$gene = as.character(select.feature$gene)
select.feature = select.feature[which(select.feature$p < 0.001),]
ct.select = intersect(names(which(table(select.feature$cell.type) > 5)),a$Var1)


rna.bind = pbmclapply(
  1:50, function(i) {
    rna = readRDS(paste0("./FastIntegrationTmp/inte/inte_", i, ".rds"))
    rna = rna[intersect(unique(select.feature$gene), rownames(rna)),]
    return(rna)
  }, mc.cores = 20
)


res = pbmclapply(
  ct.select, function(jj) {
    rna.sub = subset(rna.data, ct_rename == jj)
    rna.sub = subset(rna.sub, sample %in% names(which(table(rna.sub$sample) > 20)))
    
    
    DefaultAssay(rna.sub) = "raw"
    
    features = select.feature$gene[which(select.feature$cell.type == jj)]
    
    rna.bind2 = pbmclapply(
      rna.bind, function(rna) {
        rna = rna[,Cells(rna.sub)]
        rna = rna[intersect(rownames(rna), features),,drop= F]
        return(rna)
      }, mc.cores = 20
    )
    
    
    rna.bind2 = do.call(rbind, rna.bind2)
    rna.sub[["RNA"]] = CreateAssayObject(rna.bind2)
    
    rna.avergae = AverageExpression(rna.sub, group.by = "sample")
    rna.avergae = rna.avergae$RNA
    rna.avergae = rna.avergae[features,]
    age = meta[colnames(rna.avergae), "age"]
    
    train = data.frame(t(rna.avergae))
    train$age = age
    
    #####glmnet######
    control <- trainControl(method = "repeatedcv",
                            number = 5,
                            repeats = 5,
                            search = "random",
                            verboseIter = TRUE)
    rf <- train(age ~ ., data = train,
                method = "glmnet",
                preProcess = c("center", "scale"),
                tuneLength = 25,
                trControl = control)
    # a = data.frame(predict(rf, train[,1:(ncol(train)-1)]))
    #####rf######
    # rf = randomForest(age~., data=train, proximity=TRUE)
    # a = data.frame(rf$predicted)
    
    #####svm######
    # svmfit = svm(age ~ ., data = train, kernel = "linear", cost = 10, scale = FALSE)
    # a = data.frame(svmfit$fitted)
    
    # colnames(a) = jj
    return(rf)
  }, mc.cores = 10
)
names(res) = ct.select
saveRDS(res, "rds/age_model.rds", compress = F)

