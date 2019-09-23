#####  OCT  - logistic regression analysis using slope and genotypes  #####
library(SNPassoc)
library(lmerTest)
#library(lme4)
#library(nlme)
EPIC_OCT_slope_snp

##### check basic numbers (Gender, AgeAtExam, DiseaseCourse, DiseaseDuration) #####
EPIC_OCT_slope_snp
tmp <- merge(EPIC_OCT_slope_snp, EPIC_OCT_metadata, by.x=c("EPICID","Eye"), by.y=c("EPIC_ID","Eye"))
tmp_matrices <- subset(tmp, tmp$Eye == "L")  # for L
tmp_matrices <- subset(tmp, tmp$Eye == "R")  # for R

length(unique(tmp_matrices$EPICID))                # total Left n=354, Right n=357
colnames(OCT_matrices)[1:20]

length(unique(subset(tmp, Gender=="M")$EPICID))     # n=121
length(unique(subset(tmp, Gender=="F")$EPICID))     # n=242
# CIS PP PR RR SP UNC UNK
length(unique(subset(tmp, DiseaseCourse=="CIS")$EPICID))  # n=14
length(unique(subset(tmp, DiseaseCourse=="PP")$EPICID))   # n=10
length(unique(subset(tmp, DiseaseCourse=="RR")$EPICID))   # n=284
length(unique(subset(tmp, DiseaseCourse=="SP")$EPICID))   # n=66
length(unique(subset(tmp, DiseaseCourse=="UNC")$EPICID))  # n=10

mean(tmp$AgeAtExam)                  # 49.92562
sd(tmp$AgeAtExam)                    # +- 9.874638
mean(tmp$DiseaseDuration)            # 16.34389
sd(tmp$DiseaseDuration)              # +- 8.540359



##### check outliers #####
# MAC.VisitSeq, GCC.VisitSeq, GCL.VisitSeq, INL.VisitSeq, IPL.VisitSeq, pRNFL.VisitSeq
hist(OCT_slope_LR$MAC.VisitSeq, breaks=100)
hist(OCT_slope_LR[OCT_slope_LR$MAC.VisitSeq < 50,]$MAC.VisitSeq, breaks=100)        # MAC: select > < 50
hist(OCT_slope_LR$GCC.VisitSeq, breaks=100)
hist(OCT_slope_LR[OCT_slope_LR$GCC.VisitSeq > -0.2,]$GCC.VisitSeq, breaks=100)      # GCC: select > -0.2
hist(OCT_slope_LR$GCL.VisitSeq, breaks=100)
hist(OCT_slope_LR[OCT_slope_LR$GCL.VisitSeq > -0.2,]$GCL.VisitSeq, breaks=100)      # GCL: select > -0.2
hist(OCT_slope_LR$INL.VisitSeq, breaks=100)                                         # INL
hist(OCT_slope_LR$IPL.VisitSeq, breaks=100)
hist(OCT_slope_LR[OCT_slope_LR$IPL.VisitSeq > -0.05,]$IPL.VisitSeq, breaks=100)     # IPL:  select > -0.05
hist(OCT_slope_LR$pRNFL.VisitSeq, breaks=100)                                       # pRNFL


#####  test for first SNP #####
tmp <- subset(EPIC_OCT_slope_snp, MAC.VisitSeq < 50)
i <- 9

do <- dominant(snp(tmp[,i], sep=""))  #dominant model
ad <- additive(snp(tmp[,i], sep=""))  #additive model
re <- recessive(snp(tmp[,i], sep=""))  #recessive model
co <- codominant(snp(tmp[,i], sep=""))  #codominant model = genotype

# 'lmer' function of 'lme4' package
library(lme4)
mod <- lmer(MAC.VisitSeq ~ factor(do) + (1|EPICID), tmp)
coefs <- data.frame(coef(summary(mod)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs                  # calculation p-value
mod <- lmer(MAC.VisitSeq ~ factor(re) + (1|EPICID), tmp)
coefs <- data.frame(coef(summary(mod)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs                  # calculation p-value


# 'lmer' function of 'lmerTest' package
library(lmerTest)
mod_lmer <- summary(lmer(MAC.VisitSeq ~ factor(re) + (1|EPICID), tmp))
dfit <- summary(mod_lmer)
coef(dfit)[,5]

##########

lmData <- subset(EPIC_OCT_slope_snp, MAC.VisitSeq < 50)
vol <- "MAC.VisitSeq"


###################   Longitudinal analysis - logistic regression using linear model  #####################33
##### LogiRegression_noCovar: function for logistic regression without covariates for both eyes
Regression_noCovar <- function(lmData, vol){
  library(lmerTest)
  library(SNPassoc)
  Result <- data.frame("SNP.No"=numeric(), "SNP.Name"=character(), 
                       "Domi"=numeric(), "Additi"=numeric(), "Recess"=numeric(), "Codo1"=numeric(), "Codo2"=numeric(),
                       stringsAsFactors = F)
  slope_value <- lmData[,vol]         # slope value
  j <- 1
  for (i in 9:ncol(lmData)){
    dfit <- try(summary(lmer(slope_value ~ factor(dominant(snp(lmData[,i], sep=""))) + (1|EPICID), data=lmData)), silent=TRUE)
    if (isTRUE(class(dfit)=="try-error")){
      p.dfit <- NA
    } else {p.dfit <- coef(dfit)[,5]}
    
    afit <- try(summary(lmer(slope_value ~ additive(snp(lmData[,i], sep="")) + (1|EPICID), data=lmData)), silent=TRUE)
    if (isTRUE(class(afit)=="try-error")){
      p.afit <- NA
    } else {p.afit <- coef(afit)[,5]}
    
    rfit <- try(summary(lmer(slope_value ~ factor(recessive(snp(lmData[,i], sep=""))) + (1|EPICID), data=lmData)), silent=TRUE)
    if (isTRUE(class(rfit)=="try-error")){
      p.rfit <- NA
    } else {p.rfit <- coef(rfit)[,5]}
    
    cfit <- try(summary(lmer(slope_value ~ factor(codominant(snp(lmData[,i], sep=""))) + (1|EPICID), data=lmData)), silent=TRUE)     ## same as genotype
    if (isTRUE(class(cfit)=="try-error")){
      p.cfit <- NA
    } else {p.cfit <- coef(cfit)[,5]}
    
    Result <- rbind(Result, data.frame("SNP.No"=j, "SNP.Name"=colnames(lmData)[i], 
                               "Domi"=p.dfit[2], "Additi"=p.afit[2], "Recess"=p.rfit[2], "Codo1"=p.cfit[2], "Codo2"=p.cfit[3],
                               stringsAsFactors = F))
    j <- j+1
    rm(p.dfit, p.afit, p.rfit, p.cfit)
  }
  return(Result)
}


##### test for each volume slope (lme slope) #####
# MAC.VisitSeq, GCC.VisitSeq, GCL.VisitSeq, INL.VisitSeq, IPL.VisitSeq, pRNFL.VisitSeq
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/Analysis_results/Results_OCT/194SNP_longitudinal_results")

Sig.genes <- data.frame()
# 1 MAC.VisitSeq
lmData <- EPIC_OCT_slope_snp[EPIC_OCT_slope_snp$MAC.VisitSeq < 50,]
nrow(subset(lmData, Eye == "L")); nrow(subset(lmData, Eye == "R"))      # left n=354; right n=356
Result_table <- Regression_noCovar(lmData=lmData, vol="MAC.VisitSeq")
row.names(Result_table) <- Result_table$SNP.No
Result_table$Do.adj <- p.adjust(Result_table$Domi, method="fdr")
Result_table$Ad.adj <- p.adjust(Result_table$Additi, method="fdr")
Result_table$Re.adj <- p.adjust(Result_table$Recess, method="fdr")
Result_table$Co1.adj <- p.adjust(Result_table$Codo1, method="fdr")
Result_table$Co2.adj <- p.adjust(Result_table$Codo2, method="fdr")
Result_table$region <- c("MAC")
write.csv(Result_table, "Result_MAC_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result_table, 
                                     Domi < 0.05 | Additi < 0.05 | Recess < 0.05 | Codo1 < 0.05 | Codo2 < 0.05))
rm(Result_table, lmData)
# 2 GCC.VisitSeq
lmData <- EPIC_OCT_slope_snp[EPIC_OCT_slope_snp$GCC.VisitSeq > -0.2,]
nrow(subset(lmData, Eye == "L")); nrow(subset(lmData, Eye == "R"))      # left n=318; right n=323
Result_table <- Regression_noCovar(lmData=lmData, vol="GCC.VisitSeq")
row.names(Result_table) <- Result_table$SNP.No
Result_table$Do.adj <- p.adjust(Result_table$Domi, method="fdr")
Result_table$Ad.adj <- p.adjust(Result_table$Additi, method="fdr")
Result_table$Re.adj <- p.adjust(Result_table$Recess, method="fdr")
Result_table$Co1.adj <- p.adjust(Result_table$Codo1, method="fdr")
Result_table$Co2.adj <- p.adjust(Result_table$Codo2, method="fdr")
Result_table$region <- c("GCC")
write.csv(Result_table, "Result_GCC_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result_table, 
                                     Domi < 0.05 | Additi < 0.05 | Recess < 0.05 | Codo1 < 0.05 | Codo2 < 0.05))
rm(Result_table, lmData)

# 3 GCL.VisitSeq
lmData <- EPIC_OCT_slope_snp[EPIC_OCT_slope_snp$GCL.VisitSeq > -0.2,]
nrow(subset(lmData, Eye == "L")); nrow(subset(lmData, Eye == "R"))      # left n=351; right n=355
Result_table <- Regression_noCovar(lmData=lmData, vol="GCL.VisitSeq")
row.names(Result_table) <- Result_table$SNP.No
Result_table$Do.adj <- p.adjust(Result_table$Domi, method="fdr")
Result_table$Ad.adj <- p.adjust(Result_table$Additi, method="fdr")
Result_table$Re.adj <- p.adjust(Result_table$Recess, method="fdr")
Result_table$Co1.adj <- p.adjust(Result_table$Codo1, method="fdr")
Result_table$Co2.adj <- p.adjust(Result_table$Codo2, method="fdr")
Result_table$region <- c("GCL")
write.csv(Result_table, "Result_GCL_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result_table, 
                                     Domi < 0.05 | Additi < 0.05 | Recess < 0.05 | Codo1 < 0.05 | Codo2 < 0.05))
rm(Result_table, lmData)

# 4 INL.VisitSeq
lmData <- EPIC_OCT_slope_snp
nrow(subset(lmData, Eye == "L")); nrow(subset(lmData, Eye == "R"))      # left n=354; right n=357
Result_table <- Regression_noCovar(lmData=lmData, vol="INL.VisitSeq")
row.names(Result_table) <- Result_table$SNP.No
Result_table$Do.adj <- p.adjust(Result_table$Domi, method="fdr")
Result_table$Ad.adj <- p.adjust(Result_table$Additi, method="fdr")
Result_table$Re.adj <- p.adjust(Result_table$Recess, method="fdr")
Result_table$Co1.adj <- p.adjust(Result_table$Codo1, method="fdr")
Result_table$Co2.adj <- p.adjust(Result_table$Codo2, method="fdr")
Result_table$region <- c("INL")
write.csv(Result_table, "Result_INL_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result_table, 
                                     Domi < 0.05 | Additi < 0.05 | Recess < 0.05 | Codo1 < 0.05 | Codo2 < 0.05))
rm(Result_table, lmData)

# 5 IPL.VisitSeq
lmData <- EPIC_OCT_slope_snp[EPIC_OCT_slope_snp$IPL.VisitSeq > -0.05,]
nrow(subset(lmData, Eye == "L")); nrow(subset(lmData, Eye == "R"))      # left n=354; right n=355
Result_table <- Regression_noCovar(lmData=lmData, vol="IPL.VisitSeq")
row.names(Result_table) <- Result_table$SNP.No
Result_table$Do.adj <- p.adjust(Result_table$Domi, method="fdr")
Result_table$Ad.adj <- p.adjust(Result_table$Additi, method="fdr")
Result_table$Re.adj <- p.adjust(Result_table$Recess, method="fdr")
Result_table$Co1.adj <- p.adjust(Result_table$Codo1, method="fdr")
Result_table$Co2.adj <- p.adjust(Result_table$Codo2, method="fdr")
Result_table$region <- c("IPL")
write.csv(Result_table, "Result_IPL_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result_table, 
                                     Domi < 0.05 | Additi < 0.05 | Recess < 0.05 | Codo1 < 0.05 | Codo2 < 0.05))
rm(Result_table, lmData)

# 6 pRNFL.VisitSeq
lmData <- EPIC_OCT_slope_snp
nrow(subset(lmData, Eye == "L")); nrow(subset(lmData, Eye == "R"))      # left n=354; right n=357
Result_table <- Regression_noCovar(lmData=lmData, vol="pRNFL.VisitSeq")
row.names(Result_table) <- Result_table$SNP.No
Result_table$Do.adj <- p.adjust(Result_table$Domi, method="fdr")
Result_table$Ad.adj <- p.adjust(Result_table$Additi, method="fdr")
Result_table$Re.adj <- p.adjust(Result_table$Recess, method="fdr")
Result_table$Co1.adj <- p.adjust(Result_table$Codo1, method="fdr")
Result_table$Co2.adj <- p.adjust(Result_table$Codo2, method="fdr")
Result_table$region <- c("pRNFL")
write.csv(Result_table, "Result_pRNFL_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result_table, 
                                     Domi < 0.05 | Additi < 0.05 | Recess < 0.05 | Codo1 < 0.05 | Codo2 < 0.05))
rm(Result_table, lmData)
#

write.csv(Sig.genes, "SigGenes_slope.csv")

##########

