library(readr)
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/MRI/results/RPE_category/2-5category/all")

file <- list.files(path=".", pattern="Result_PRE_")
file_name <- strsplit(file, ".csv")
file_name <- substr(file_name, start=12, stop=30)

result <- read_csv(file[1])
Attribute <- data.frame(Gene.Name=result$Gene.Name)
Attribute_pval <- data.frame(Gene.Name=result$Gene.Name)
for (i in 1:length(file)){
  result <- read_csv(file[i])
  
  sig_result <- subset(result, result$Pval < 0.05)
  sig_result$NetAttr <- c(1)
  
  result <- merge(result, sig_result[,c("Gene.Name","NetAttr")], by.x="Gene.Name", by.y="Gene.Name", all.x=TRUE)
  result[is.na(result$NetAttr),]$NetAttr <- c(0)
  
  #write.csv(result, file=paste0("NetworkAttribute/",file[i],"_NetAttr.csv"))
  Attribute <- merge(Attribute, result[,c("Gene.Name","NetAttr")], by.x="Gene.Name", by.y="Gene.Name", all.x=TRUE, all.y=TRUE)
  colnames(Attribute)[i+1] <- file_name[i]
  Attribute_pval <- merge(Attribute_pval, result[,c("Gene.Name","Pval")], by.x="Gene.Name", by.y="Gene.Name", all.x=TRUE, all.y=TRUE)
  colnames(Attribute_pval)[i+1] <- file_name[i]
}
write.csv(Attribute, "Attribute_sigGenes.csv")
write.csv(Attribute_pval, "Attribute_sigGenes_pvalue.csv")


###
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/OCT/results/PRE_category/PRE_factor/all")

file <- list.files(path=".", pattern="Result_PRE_")
file_name <- strsplit(file, ".csv")
file_name <- substr(file_name, start=12, stop=30)

result <- read_csv(file[1])
Attribute <- data.frame(Gene.Name=result$Gene.Name)
for (i in 1:length(file)){
  result <- read_csv(file[i])
  
  sig_result <- subset(result, result$Pval1 < 0.05 | result$Pval2 < 0.05 | result$Pval3 < 0.05)
  sig_result$NetAttr <- c(1)
  
  result <- merge(result, sig_result[,c("Gene.Name","NetAttr")], by.x="Gene.Name", by.y="Gene.Name", all.x=TRUE)
  result[is.na(result$NetAttr),]$NetAttr <- c(0)
  
  #write.csv(result, file=paste0("NetworkAttribute/",file[i],"_NetAttr.csv"))
  Attribute <- merge(Attribute, result[,c("Gene.Name","NetAttr")], by.x="Gene.Name", by.y="Gene.Name", all.x=TRUE, all.y=TRUE)
  colnames(Attribute)[i+1] <- file_name[i]
}
write.csv(Attribute, "Attribute_sigGenes(0.05).csv")

###
file <- list.files(path=".", pattern="Result_PRE_")
file_name <- strsplit(file, ".csv")
file_name <- substr(file_name, start=12, stop=30)

result <- read_csv(file[1])
Attribute <- data.frame(Gene.Name=result$Gene.Name)
for (i in 1:length(file)){
  result <- read_csv(file[i])
  
  sig_result <- subset(result, result$Pval1 < 0.01 | result$Pval2 < 0.01 | result$Pval3 < 0.01)
  sig_result$NetAttr <- c(1)
  
  result <- merge(result, sig_result[,c("Gene.Name","NetAttr")], by.x="Gene.Name", by.y="Gene.Name", all.x=TRUE, all.y=TRUE)
  result[is.na(result$NetAttr),]$NetAttr <- c(0)
  
  #write.csv(result, file=paste0("NetworkAttribute/",file[i],"_NetAttr.csv"))
  Attribute <- merge(Attribute, result[,c("Gene.Name","NetAttr")], by.x="Gene.Name", by.y="Gene.Name", all.x=TRUE, all.y=TRUE)
  colnames(Attribute)[i+1] <- file_name[i]
}
write.csv(Attribute, "Attribute_sigGenes(0.01).csv")