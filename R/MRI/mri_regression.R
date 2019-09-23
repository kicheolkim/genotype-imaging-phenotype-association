#####  MRI  - logistic regression analysis using slope and genotypes  #####
library(SNPassoc)
EPIC_MRI_snp
save(MRI_slope_lme, EPIC_metadata, EPIC_194snp, EPIC_MRI_dmt, EPIC_MRI_Normalized, file="EPIC_MRI_dataset.RData")
load(file="C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/MRI/RData/EPIC_MRI_dataset.RData")

##### check basic numbers (Gender, AgeAtExam, DiseaseCourse, DiseaseDuration) #####
length(unique(EPIC_MRI_snp$EPICID))
colnames(EPIC_MRI_snp)[1:20]

length(unique(subset(EPIC_MRI_snp, Gender=="M")$EPICID))
length(unique(subset(EPIC_MRI_snp, Gender=="F")$EPICID))
length(unique(subset(EPIC_MRI_snp, DiseaseCourse=="CIS")$EPICID))  # CIS PP PR RR SP UNC UNK
length(unique(subset(EPIC_MRI_snp, DiseaseCourse=="PP")$EPICID))
length(unique(subset(EPIC_MRI_snp, DiseaseCourse=="PR")$EPICID))
length(unique(subset(EPIC_MRI_snp, DiseaseCourse=="RR")$EPICID))
length(unique(subset(EPIC_MRI_snp, DiseaseCourse=="SP")$EPICID))
length(unique(subset(EPIC_MRI_snp, DiseaseCourse=="UNC")$EPICID))
length(unique(subset(EPIC_MRI_snp, DiseaseCourse=="UNK")$EPICID))

mean(EPIC_MRI_snp$AgeAtExam)
sd(EPIC_MRI_snp$AgeAtExam)
mean(EPIC_MRI_snp$DiseaseDuration)
sd(EPIC_MRI_snp$DiseaseDuration)
#
#####
##### check outliers #####
# check outliers for lm slope
hist(EPIC_MRI_slope$NBV_Slop, breaks=100)         # select > -300 | < 250
hist(EPIC_MRI_slope$NGMV_Slop, breaks=100)
hist(EPIC_MRI_slope$NWMV_Slop, breaks=100)        # select > -150 | < 150
hist(EPIC_MRI_slope$NT2Vlog_Slop, breaks=100)     # select > -1.5
hist(EPIC_MRI_slope$NVCSFlog_Slop, breaks=100)
hist(EPIC_MRI_slope$NCGM_Slop, breaks=100)
hist(EPIC_MRI_slope$PBVC_since_base_Slop, breaks=100)        # select < 5

# check outliers for lme slope
hist(EPIC_MRI_snp$NBV.Time, breaks=100)         # select > -100
hist(EPIC_MRI_snp$NGMV.Time, breaks=100)        # select > -60
hist(EPIC_MRI_snp$NWMV.Time, breaks=100)        # select > -55
hist(EPIC_MRI_snp$NT2Vlog.Time, breaks=100)     # select > -1 | < 1
hist(EPIC_MRI_snp$NVCSFlog.Time, breaks=100)
hist(EPIC_MRI_snp$NCGM.Time, breaks=100)        # select > -60

#####  test for first SNP #####
lmData <- subset(EPIC_MRI_snp, NBV.Time > -100)
i <- 30
do <- dominant(snp(EPIC_MRI_snp[,"chr11:118762073"], sep=""))  #dominant model
ad <- additive(snp(EPIC_MRI_snp[,i], sep=""))  #additive model
re <- recessive(snp(EPIC_MRI_snp[,i], sep=""))  #recessive model
co <- codominant(snp(EPIC_MRI_snp[,i], sep=""))  #codominant model = genotype
summary(lm(NBV_Slop ~ ad + DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=EPIC_MRI_snp))   # with covariates
summary(lm(NBV.Time ~ ad, data=EPIC_MRI_snp))   # without covariates

# multiple SNPs
summary(lm(slope_value ~ factor(dominant(snp(lmData[,157], sep="")))+factor(dominant(snp(lmData[,129], sep="")))+factor(dominant(snp(lmData[,46], sep=""))) , data=lmData))
#####


##########
###################   Longitudinal analysis - logistic regression using linear model  #####################33
##### TestRegression_noCovar: function for logistic regression without covariates  #####
TestRegression_noCovar <- function(lmData, vol){
  library(SNPassoc)
  Result <- data.frame("SNP.No"=numeric(), "SNP.Name"=character(), 
                       "Domi"=numeric(), "Additi"=numeric(), "Recess"=numeric(), "Codo1"=numeric(), "Codo2"=numeric(),
                       stringsAsFactors = F)
  slope_value <- lmData[,vol]         # slope value
  j <- 1
  for (i in 28:ncol(lmData)){
    dfit <- try(summary(lm(slope_value ~ factor(dominant(snp(lmData[,i], sep=""))) , data=lmData)), silent=TRUE)
    if (isTRUE(class(dfit)=="try-error")){
      p.dfit <- NA
    } else {p.dfit <- coef(dfit)[,4]}
    
    afit <- try(summary(lm(slope_value ~ additive(snp(lmData[,i], sep="")) , data=lmData)), silent=TRUE)
    if (isTRUE(class(afit)=="try-error")){
      p.afit <- NA
    } else {p.afit <- coef(afit)[,4]}
    
    rfit <- try(summary(lm(slope_value ~ factor(recessive(snp(lmData[,i], sep=""))) , data=lmData)), silent=TRUE)
    if (isTRUE(class(rfit)=="try-error")){
      p.rfit <- NA
    } else {p.rfit <- coef(rfit)[,4]}
    
    cfit <- try(summary(lm(slope_value ~ factor(codominant(snp(lmData[,i], sep=""))) , data=lmData)), silent=TRUE)     ## same as genotype
    if (isTRUE(class(cfit)=="try-error")){
      p.cfit <- NA
    } else {p.cfit <- coef(cfit)[,4]}
    
    Result <- rbind(Result, 
                    data.frame("SNP.No"=j, "SNP.Name"=colnames(lmData)[i], 
                               "Domi"=p.dfit[2], "Additi"=p.afit[2], "Recess"=p.rfit[2], "Codo1"=p.cfit[2], "Codo2"=p.cfit[3],
                               stringsAsFactors = F))
    j <- j+1
    rm(p.dfit, p.afit, p.rfit, p.cfit)
  }
  rm(slope_value, lmData)
  return(Result)
}
#####
##### test for each volume slope (lme slope) #####
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/Analysis_results/Results_MRI/194SNP_longitudinal_results")


Sig.genes <- data.frame()
# 1 NBV_Slop
lmData <- subset(EPIC_MRI_snp, NBV.Time > -100)     # filter outlier: slope of normalized brain volume
Result <- TestRegression_noCovar(lmData=lmData, vol="NBV.Time")
row.names(Result) <- Result$SNP.No
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Co1.adj <- p.adjust(Result$Codo1, method="fdr")
Result$Co2.adj <- p.adjust(Result$Codo2, method="fdr")
Result$region <- c("Total brain")
write.csv(Result, "Result_NBV_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
rm(Result)
# 2 NGMV_Slop
lmData <- subset(EPIC_MRI_snp, NGMV.Time > -60)                                   # filter outlier: slope of normalized gray matter volume
Result <- TestRegression_noCovar(lmData=lmData, vol="NGMV.Time")
row.names(Result) <- Result$SNP.No
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Co1.adj <- p.adjust(Result$Codo1, method="fdr")
Result$Co2.adj <- p.adjust(Result$Codo2, method="fdr")
Result$region <- c("Gray matter")
write.csv(Result, "Result_NGMV_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
rm(Result)
# 3 NWMV_Slop
lmData <- subset(EPIC_MRI_snp, NWMV.Time > -55)    # filter outlier: slope of normalized white matter volume
Result <- TestRegression_noCovar(lmData=lmData, vol="NWMV.Time")
row.names(Result) <- Result$SNP.No
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Co1.adj <- p.adjust(Result$Codo1, method="fdr")
Result$Co2.adj <- p.adjust(Result$Codo2, method="fdr")
Result$region <- c("White matter")
write.csv(Result, "Result_NWMV_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
rm(Result)
# 4 NT2Vlog_Slop
lmData <- subset(EPIC_MRI_snp, NT2Vlog.Time > -1 | NT2Vlog.Time < 1)     # filter outlier: slope of normalized T2 volume
Result <- TestRegression_noCovar(lmData=lmData, vol="NT2Vlog.Time")
row.names(Result) <- Result$SNP.No
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Co1.adj <- p.adjust(Result$Codo1, method="fdr")
Result$Co2.adj <- p.adjust(Result$Codo2, method="fdr")
Result$region <- c("T2 lesion (log)")
write.csv(Result, "Result_NT2Vlog_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
rm(Result)
# 5 NVCSFlog_Slop
lmData <- EPIC_MRI_snp                                 # filter outlier: slope of normalized ventricular csf volume
Result <- TestRegression_noCovar(lmData=lmData, vol="NVCSFlog.Time")
row.names(Result) <- Result$SNP.No
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Co1.adj <- p.adjust(Result$Codo1, method="fdr")
Result$Co2.adj <- p.adjust(Result$Codo2, method="fdr")
Result$region <- c("VCSF (log)")
write.csv(Result, "Result_NVCSFlog_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
rm(Result)
# 6 NCGM_Slop
lmData <- subset(EPIC_MRI_snp, NCGM.Time > -60)     # filter outlier: slope of normalized ventricular csf volume
Result <- TestRegression_noCovar(lmData=lmData, vol="NCGM.Time")
row.names(Result) <- Result$SNP.No
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Co1.adj <- p.adjust(Result$Codo1, method="fdr")
Result$Co2.adj <- p.adjust(Result$Codo2, method="fdr")
Result$region <- c("Cortical gray matter")
write.csv(Result, "Result_NCGM_slope.csv")
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
rm(Result)
#
write.csv(Sig.genes, "SigGenes_slope.csv")
##########
##########################################################################################################333






###################   Cross-sectional analysis - logistic regression using linear model (covariates: hla, age at exam, gender, disease duration)  #####################33
TestRegression <- function(lmData, vol){
  library(SNPassoc)
  Result <- data.frame("SNP.No"=numeric(), "SNP.Name"=character(), 
                       "Domi"=numeric(), "HLA"=numeric(), "Age"=numeric(), "Gendr"=numeric(), "Dis.Dur"=numeric(),
                       "Additi"=numeric(), "HLA2"=numeric(), "Age2"=numeric(), "Gendr2"=numeric(), "Dis.Dur2"=numeric(),
                       "Recess"=numeric(), "HLA3"=numeric(), "Age3"=numeric(), "Gendr3"=numeric(), "Dis.Dur3"=numeric(),
                       stringsAsFactors = F)
  slope_value <- lmData[,vol]         # slope value
  j <- 1
  for (i in 35:ncol(lmData)){
    dfit <- try(summary(lm(slope_value ~ factor(dominant(snp(lmData[,i], sep=""))) + DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(dfit)=="try-error")){
      p.dfit <- NA
    } else {p.dfit <- coef(dfit)[,4]}
    
    afit <- try(summary(lm(slope_value ~ additive(snp(lmData[,i], sep="")) + DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(afit)=="try-error")){
      p.afit <- NA
    } else {p.afit <- coef(afit)[,4]}
    
    rfit <- try(summary(lm(slope_value ~ factor(recessive(snp(lmData[,i], sep=""))) + DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(rfit)=="try-error")){
      p.rfit <- NA
    } else {p.rfit <- coef(rfit)[,4]}
    
    Result <- rbind(Result, 
                    data.frame("SNP.No"=j, "SNP.Name"=colnames(lmData)[i], 
                               "Domi"=p.dfit[2], "HLA"=p.dfit[3], "Age"=p.dfit[4], "Gendr"=p.dfit[5], "Dis.Dur"=p.dfit[6],
                               "Additi"=p.afit[2], "HLA2"=p.afit[3], "Age2"=p.afit[4], "Gendr2"=p.afit[5], "Dis.Dur2"=p.afit[6],
                               "Recess"=p.rfit[2], "HLA3"=p.rfit[3], "Age3"=p.rfit[4], "Gendr3"=p.rfit[5], "Dis.Dur3"=p.rfit[6],
                               stringsAsFactors = F))
    j <- j+1
    rm(p.dfit, p.afit, p.rfit, p.cfit)
  }
  rm(slope_value, lmData)
  return(Result)
}


# check outliers for lm slope
hist(EPIC_MRI_base_meta_snp$NBV, breaks=100) 
hist(EPIC_MRI_base_meta_snp$NGMV, breaks=100)
hist(EPIC_MRI_base_meta_snp$NWMV, breaks=100)
hist(EPIC_MRI_base_meta_snp$NT2Vlog, breaks=100)
hist(EPIC_MRI_base_meta_snp$NVCSFlog, breaks=100)
hist(EPIC_MRI_base_meta_snp$NCGM, breaks=100)


##### test for each volume slope (lme slope) #####
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/Analysis_results/Results_MRI/194SNP_crosssectional_results")
Sig.genes <- data.frame()
lmData <- EPIC_MRI_base_meta_snp

# 1 NBV
Result <- TestRegression(lmData=lmData, vol="NBV")
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Region <- "BV"
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
write.csv(Result, "Result_NBV_baseline.csv")
rm(Result)

# 2 NGMV
Result <- TestRegression(lmData=lmData, vol="NGMV")
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Region <- "GMV"
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
write.csv(Result, "Result_NGMV_baseline.csv")
rm(Result)

# 3 NWMV
Result <- TestRegression(lmData=lmData, vol="NWMV")
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Region <- "WMV"
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
write.csv(Result, "Result_NWMV_baseline.csv")
rm(Result)

# 4 NT2Vlog
Result <- TestRegression(lmData=lmData, vol="NT2Vlog")
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Region <- "logT2V"
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
write.csv(Result, "Result_NT2Vlog_baseline.csv")
rm(Result)

# 5 NVCSFlog
Result <- TestRegression(lmData=lmData, vol="NVCSFlog")
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Region <- "logVCSF"
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
write.csv(Result, "Result_NVCSFlog_baseline.csv")
rm(Result)

# 6 NCGM
Result <- TestRegression(lmData=lmData, vol="NCGM")
Result$Do.adj <- p.adjust(Result$Domi, method="fdr")
Result$Ad.adj <- p.adjust(Result$Additi, method="fdr")
Result$Re.adj <- p.adjust(Result$Recess, method="fdr")
Result$Region <- "CGM"
Sig.genes <- rbind(Sig.genes, subset(Result, Domi < 0.05 | Additi < 0.05 | Recess < 0.05))
write.csv(Result, "Result_NCGM_baseline.csv")
rm(Result)
#

write.csv(Sig.genes, "SigGenes_baseline.csv")

##########
##########################################################################################################333






###################   Longitudinal analysis - logistic regression (covariates: hla, age at exam, gender, disease duration)  #####################33
##### TestRegression: function for logistic regression with covariates #####
TestRegression <- function(lmData, vol){
  library(SNPassoc)
  Result <- data.frame("SNP.No"=numeric(), "SNP.Name"=character(), 
                       "Domi"=numeric(), "HLA"=numeric(), "Age"=numeric(), "Gendr"=numeric(), "Dis.Dur"=numeric(),
                       "Additi"=numeric(), "HLA2"=numeric(), "Age2"=numeric(), "Gendr2"=numeric(), "Dis.Dur2"=numeric(),
                       "Recess"=numeric(), "HLA3"=numeric(), "Age3"=numeric(), "Gendr3"=numeric(), "Dis.Dur3"=numeric(),
                       "Codo1"=numeric(), "Codo2"=numeric(), "HLA4"=numeric(), "Age4"=numeric(), "Gendr4"=numeric(), "Dis.Dur4"=numeric(),
                       stringsAsFactors = F)
  slope_value <- lmData[,vol]         # slope value
  j <- 1
  for (i in 30:ncol(lmData)){
    dfit <- try(summary(lm(slope_value ~ factor(dominant(snp(lmData[,i], sep=""))) + DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(dfit)=="try-error")){
      p.dfit <- NA
    } else {p.dfit <- coef(dfit)[,4]}
    
    afit <- try(summary(lm(slope_value ~ additive(snp(lmData[,i], sep="")) + DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(afit)=="try-error")){
      p.afit <- NA
    } else {p.afit <- coef(afit)[,4]}
    
    rfit <- try(summary(lm(slope_value ~ factor(recessive(snp(lmData[,i], sep=""))) + DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(rfit)=="try-error")){
      p.rfit <- NA
    } else {p.rfit <- coef(rfit)[,4]}
    
    cfit <- try(summary(lm(slope_value ~ factor(codominant(snp(lmData[,i], sep=""))) + DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)     ## same as genotype
    if (isTRUE(class(cfit)=="try-error")){
      p.cfit <- NA
    } else {p.cfit <- coef(cfit)[,4]}
    
    Result <- rbind(Result, 
                    data.frame("SNP.No"=j, "SNP.Name"=colnames(lmData)[i], 
                               "Domi"=p.dfit[2], "HLA"=p.dfit[3], "Age"=p.dfit[4], "Gendr"=p.dfit[5], "Dis.Dur"=p.dfit[6],
                               "Additi"=p.afit[2], "HLA2"=p.afit[3], "Age2"=p.afit[4], "Gendr2"=p.afit[5], "Dis.Dur2"=p.afit[6],
                               "Recess"=p.rfit[2], "HLA3"=p.rfit[3], "Age3"=p.rfit[4], "Gendr3"=p.rfit[5], "Dis.Dur3"=p.rfit[6],
                               "Codo1"=p.cfit[2], "Codo2"=p.cfit[3], "HLA4"=p.cfit[4], "Age4"=p.cfit[5], "Gendr4"=p.cfit[6], "Dis.Dur4"=p.cfit[7],
                               stringsAsFactors = F))
    j <- j+1
    rm(p.dfit, p.afit, p.rfit, p.cfit)
  }
  rm(slope_value, lmData)
  return(Result)
}

##### test for each volume slope (lm slope) #####
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/MRI/results")

# 1 NBV_Slop
lmData <- subset(EPIC_MRI_snp, NBV_Slop > -300 & NBV_Slop < 250)     # filter outlier: slope of normalized brain volume
Result <- TestRegression(lmData=lmData, vol="NBV_Slop")
result_NBV_slop <- Result
write.csv(Result, "Result_NBV_slope.csv")
rm(Result)

# 2 NGMV_Slop
lmData <- EPIC_MRI_snp                                                 # filter outlier: slope of normalized gray matter volume
Result <- TestRegression(lmData=lmData, vol="NGMV_Slop")
result_NGMV_slop <- Result
write.csv(Result, "Result_NGMV_slope.csv")
rm(Result)

# 3 NWMV_Slop
lmData <- subset(EPIC_MRI_snp, NWMV_Slop > -150 & NWMV_Slop < 150)    # filter outlier: slope of normalized white matter volume
Result <- TestRegression(lmData=lmData, vol="NWMV_Slop")
result_NWMV_slop <- Result
write.csv(Result, "Result_NWMV_slope.csv")
rm(Result)

# 4 NT2Vlog_Slop
lmData <- subset(EPIC_MRI_snp, NT2Vlog_Slop > -1.5)     # filter outlier: slope of normalized T2 volume
Result <- TestRegression(lmData=lmData, vol="NT2Vlog_Slop")
result_NT2Vlog_slop <- Result
write.csv(Result, "Result_NT2Vlog_slope.csv")
rm(Result)

# 5 NVCSFlog_Slop
lmData <- EPIC_MRI_snp     # filter outlier: slope of normalized ventricular csf volume
Result <- TestRegression(lmData=lmData, vol="NVCSFlog_Slop")
result_NVCSFlog_slop <- Result
write.csv(Result, "Result_NVCSFlog_slope.csv")
rm(Result)

# 6 NCGM_Slop
lmData <- EPIC_MRI_snp     # filter outlier: slope of normalized ventricular csf volume
Result <- TestRegression(lmData=lmData, vol="NCGM_Slop")
result_NCGM_slop <- Result
write.csv(Result, "Result_NCGM_slope.csv")
rm(Result)

# 7 PBVC_since_base_Slop
lmData <- subset(EPIC_MRI_snp, PBVC_since_base_Slop < 5)     # filter outlier: slope of PBVC_since_baseline
Result <- TestRegression(lmData=lmData, vol="PBVC_since_base_Slop")
result_PBVC_since_base_slop <- Result
write.csv(Result, "Result_PBVC_since_base_slope.csv")
rm(Result)
#

##########


###################   Longitudinal analysis - Interaction of HLA-DRB*15:01     ##################################################################33
### logistic regression (hla as an interaction term, covariates: age at exam, gender, disease duration)
######  function for logistic regression with interaction of HLA-DRB*15:01 with covariates #####
TestRegressionHLA <- function(lmData, vol){
  Result <- data.frame("SNP.No"=numeric(), "SNP.Name"=character(), 
                       "Domi"=numeric(), "HLA"=numeric(), "Age"=numeric(), "Gendr"=numeric(), "Dis.Dur"=numeric(), "Do+HLA"=numeric(),
                       "Additi"=numeric(), "HLA2"=numeric(), "Age2"=numeric(), "Gendr2"=numeric(), "Dis.Dur2"=numeric(), "Ad+HLA"=numeric(),
                       "Recess"=numeric(), "HLA3"=numeric(), "Age3"=numeric(), "Gendr3"=numeric(), "Dis.Dur3"=numeric(), "Re+HLA"=numeric(),
                       "Codo1"=numeric(), "Codo2"=numeric(), "HLA4"=numeric(), "Age4"=numeric(), "Gendr4"=numeric(), "Dis.Dur4"=numeric(), "Co1+HLA"=numeric(), "Co2+HLA"=numeric(),
                       stringsAsFactors = F)
  slope_value <- lmData[,vol]         # slope value
  j <- 1
  for (i in 30:ncol(lmData)){
    dfit <- try(summary(lm(slope_value ~ factor(dominant(snp(lmData[,i], sep=""))) * DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(dfit)=="try-error")){
      p.dfit <- NA
    } else {p.dfit <- coef(dfit)[,4]}
    
    afit <- try(summary(lm(slope_value ~ additive(snp(lmData[,i], sep="")) * DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(afit)=="try-error")){
      p.afit <- NA
    } else {p.afit <- coef(afit)[,4]}
    
    rfit <- try(summary(lm(slope_value ~ factor(recessive(snp(lmData[,i], sep=""))) * DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)
    if (isTRUE(class(rfit)=="try-error")){
      p.rfit <- NA
    } else {p.rfit <- coef(rfit)[,4]}
    
    cfit <- try(summary(lm(slope_value ~ factor(codominant(snp(lmData[,i], sep=""))) * DRBstat1 + AgeAtExam + Gender + DiseaseDuration, data=lmData)), silent=TRUE)     ## same as genotype
    if (isTRUE(class(cfit)=="try-error")){
      p.cfit <- NA
    } else {p.cfit <- coef(cfit)[,4]}
    
    Result <- rbind(Result, 
                    data.frame("SNP.No"=j, "SNP.Name"=colnames(lmData)[i], 
                               "Domi"=p.dfit[2], "HLA"=p.dfit[3], "Age"=p.dfit[4], "Gendr"=p.dfit[5], "Dis.Dur"=p.dfit[6], "Do+HLA"=p.dfit[7],
                               "Additi"=p.afit[2], "HLA2"=p.afit[3], "Age2"=p.afit[4], "Gendr2"=p.afit[5], "Dis.Dur2"=p.afit[6], "Ad+HLA"=p.afit[7],
                               "Recess"=p.rfit[2], "HLA3"=p.rfit[3], "Age3"=p.rfit[4], "Gendr3"=p.rfit[5], "Dis.Dur3"=p.rfit[6], "Re+HLA"=p.rfit[7],
                               "Codo1"=p.cfit[2], "Codo2"=p.cfit[3], "HLA4"=p.cfit[4], "Age4"=p.cfit[5], "Gendr4"=p.cfit[6], "Dis.Dur4"=p.cfit[7], "Co1+HLA"=p.cfit[8], "Co2+HLA"=p.cfit[9],
                               stringsAsFactors = F))
    j <- j+1
    rm(p.dfit, p.afit, p.rfit, p.cfit)
  }
  rm(slope_value, lmData)
  return(Result)
}


######  test for each volume slope  #####
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/MRI/results/withHLA")

# 1 NBV_Slop
lmData <- subset(EPIC_MRI_snp, NBV_Slop > -300 & NBV_Slop < 250)     # filter outlier: slope of normalized brain volume
Result <- TestRegressionHLA(lmData=lmData, vol="NBV_Slop")
result_NBV_slop_hla <- Result
write.csv(Result, "Result_NBV_slope_hla.csv")
rm(Result)

# 2 NGMV_Slop
lmData <- subset(EPIC_MRI_snp, NGMV_Slop > -200)                     # filter outlier: slope of normalized gray matter volume
Result <- TestRegressionHLA(lmData=lmData, vol="NGMV_Slop")
result_NGMV_slop_hla <- Result
write.csv(Result, "Result_NGMV_slope_hla.csv")
rm(Result)

# 3 NWMV_Slop
lmData <- subset(EPIC_MRI_snp, NWMV_Slop > -150 & NWMV_Slop < 150)    # filter outlier: slope of normalized white matter volume
Result <- TestRegressionHLA(lmData=lmData, vol="NWMV_Slop")
result_NWMV_slop_hla <- Result
write.csv(Result, "Result_NWMV_slope_hla.csv")
rm(Result)

# 4 NT2Vlog_Slop
lmData <- subset(EPIC_MRI_snp, NT2Vlog_Slop > -2.5 & NT2Vlog_Slop < 2)     # filter outlier: slope of normalized T2 volume
Result <- TestRegressionHLA(lmData=lmData, vol="NT2Vlog_Slop")
result_NT2Vlog_slop_hla <- Result
write.csv(Result, "Result_NT2Vlog_slope_hla.csv")
rm(Result)

# 5 NVCSFlog_Slop
lmData <- subset(EPIC_MRI_snp, NVCSFlog_Slop < 0.6)     # filter outlier: slope of normalized ventricular csf volume
Result <- TestRegressionHLA(lmData=lmData, vol="NVCSFlog_Slop")
result_NVCSFlog_slop_hla <- Result
write.csv(Result, "Result_NVCSFlog_slope_hla.csv")
rm(Result)

# 6 NCGM_Slop
lmData <- subset(EPIC_MRI_snp, NCGM_Slop > -150 & NCGM_Slop < 120)     # filter outlier: slope of normalized ventricular csf volume
Result <- TestRegressionHLA(lmData=lmData, vol="NCGM_Slop")
result_NCGM_slop_hla <- Result
write.csv(Result, "Result_NCGM_slope_hla.csv")
rm(Result)

# 7 PBVC_since_base_Slop
lmData <- subset(EPIC_MRI_snp, PBVC_since_base_Slop < 5)     # filter outlier: slope of PBVC_since_baseline
Result <- TestRegressionHLA(lmData=lmData, vol="PBVC_since_base_Slop")
result_PBVC_since_base_slop_hla <- Result
write.csv(Result, "Result_PBVC_since_base_slope_hla.csv")
rm(Result)
#

##########
##########  ###
