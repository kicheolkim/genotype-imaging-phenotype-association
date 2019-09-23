EPIC_OCT_new2_snp_slope
EPIC_OCT_7SNPs_slopes

both <- merge(EPIC_OCT_new2_snp_slope[,c(1,2,7,9,11,13,15,17,19,21,22,23,24,25,26,27,33,34)], 
              EPIC_OCT_7SNPs_slopes[,c(1,2,7,9,11,12,13,19,20,21,22,23,24,25,27,29,31,32,35,36,40)], 
              by.x=c("EPIC_ID","Eye"), by.y=c("EPIC_ID","Eye"), all.x=TRUE, all.y=TRUE)
write.csv(both, "EPIC_dataset.csv")

length(unique(EPIC_OCT_new2_snp_slope$EPIC_ID)) # 386 patients

# MAC, GCC, GCL, INL, IPL, pRNFL
hist(EPIC_OCT_new2_snp_slope$MAC_Slope, breaks=100)
hist(EPIC_OCT_7SNPs_slopes$MAC_Slope, breaks=100)
hist(EPIC_OCT_new2_snp_slope$GCC_Slope, breaks=100)
hist(EPIC_OCT_7SNPs_slopes$GCC_Slope, breaks=100)
hist(EPIC_OCT_new2_snp_slope$GCL_Slope, breaks=100)
hist(EPIC_OCT_7SNPs_slopes$GCL_Slope, breaks=100)
hist(EPIC_OCT_new2_snp_slope$INL_Slope, breaks=100)
hist(EPIC_OCT_7SNPs_slopes$INL_Slope, breaks=100)
hist(EPIC_OCT_new2_snp_slope$IPL_Slope, breaks=100)
hist(EPIC_OCT_7SNPs_slopes$IPL_Slope, breaks=100)
hist(EPIC_OCT_new2_snp_slope$pRNFL_Slope, breaks=100)
hist(EPIC_OCT_7SNPs_slopes$pRNFL_Slope, breaks=100)


setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/Collabo_Calabresi/ReAnalysis/MSChip_results_withoutHLA")



### Filtering outliers
# 1. MAC
hist(EPIC_OCT_new2_snp_slope$MAC_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$MAC_Slope < 1.5 & EPIC_OCT_new2_snp_slope$MAC_Slope > -1.5),]
hist(OCT_dataset$MAC_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 383 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$MAC_Slope),]
vol <- as.numeric(as.character(lme_data$MAC_Slope))
filename <- c("EPIC_OCT_new_slope-lme_MAC.csv")


# 2. GCC
hist(EPIC_OCT_new2_snp_slope$GCC_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$GCC_Slope > -5),]
hist(OCT_dataset$GCC_Slope, breaks=100)
length(unique(OCT_dataset$EPIC_ID)) # 348 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$GCC_Slope),]
vol <- as.numeric(as.character(lme_data$GCC_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_GCC.csv")


# 3. GCL
hist(EPIC_OCT_new2_snp_slope$GCL_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$GCL_Slope > -2),]
hist(OCT_dataset$GCL_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 384 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$GCL_Slope),]
vol <- as.numeric(as.character(lme_data$GCL_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_GCL.csv")


# 4. INL
hist(EPIC_OCT_new2_snp_slope$INL_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$INL_Slope > -6e-15 & EPIC_OCT_new2_snp_slope$INL_Slope < 6e-15),]
hist(OCT_dataset$INL_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 370 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$INL_Slope),]
vol <- as.numeric(as.character(lme_data$INL_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_INL.csv")


# 5. IPL
hist(EPIC_OCT_new2_snp_slope$IPL_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$IPL_Slope > -2),]
hist(OCT_dataset$IPL_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 384 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$IPL_Slope),]
vol <- as.numeric(as.character(lme_data$IPL_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_IPL.csv")


# 6. pRNFL
hist(EPIC_OCT_new2_snp_slope$pRNFL_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$pRNFL_Slope > -50 & EPIC_OCT_new2_snp_slope$pRNFL_Slope < 50),]
hist(OCT_dataset$pRNFL_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 361 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$pRNFL_Slope),]
vol <- as.numeric(as.character(lme_data$pRNFL_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_pRNFL.csv")




#########################################################################################################################=====
#####  Linear mixed effect model, adjusting for inter-eye correlation  (with optic neuritis information)  #####
#####  WITHOUT HLA-DRB
Result <- data.frame("SNP.No"=numeric(), "SNP.name"=character(), 
                     "Domi"=numeric(), "ON"=numeric(), "Age"=numeric(), "Gendr"=numeric(), "Dis.Dur"=numeric(), 
                     "Additi"=numeric(), "ON2"=numeric(), "Age2"=numeric(), "Gendr2"=numeric(), "Dis.Dur2"=numeric(), 
                     "Recess"=numeric(), "ON3"=numeric(), "Age3"=numeric(), "Gendr3"=numeric(), "Dis.Dur3"=numeric(), 
                     "Codo1"=numeric(), "Codo2"=numeric(), "ON4"=numeric(), "Age4"=numeric(), "Gendr4"=numeric(), "Dis.Dur4"=numeric(),  
                     stringsAsFactors = F)
j <- 1

for (i in 36:55){
  snp.name <- colnames(lme_data)[i]
  snp <- lme_data[,i]               # genotype
  
  do <- dominant(snp(snp, sep=""))  #dominant model
  ad <- additive(snp(snp, sep=""))  #additive model
  re <- recessive(snp(snp, sep=""))  #recessive model
  co <- codominant(snp(snp, sep=""))  #codominant model = genotype
  
  ## some patients are diagnosed in same year, so including month.   
  dfit <- try(summary(lme(vol ~ do+ON+AgeAtExam+Gender+DiseaseDuration, random=~1|Vision.PIDN, data=lme_data, na.action=na.exclude)), silent=TRUE)
  if (isTRUE(class(dfit)=="try-error")){
    p.dfit <- NA
  } else {p.dfit <- dfit$tTable[,5]}
  
  afit <- try(summary(lme(vol ~ ad+ON+AgeAtExam+Gender+DiseaseDuration, random=~1|Vision.PIDN, data=lme_data, na.action=na.exclude)), silent=TRUE)
  if (isTRUE(class(afit)=="try-error")){
    p.afit <- NA
  } else {p.afit <- afit$tTable[,5]}
  
  rfit <- try(summary(lme(vol ~ re+ON+AgeAtExam+Gender+DiseaseDuration, random=~1|Vision.PIDN, data=lme_data, na.action=na.exclude)), silent=TRUE)
  if (isTRUE(class(rfit)=="try-error")){
    p.rfit <- NA
  } else {p.rfit <- rfit$tTable[,5]}
  
  cfit <- try(summary(lme(vol ~ co+ON+AgeAtExam+Gender+DiseaseDuration, random=~1|Vision.PIDN, data=lme_data, na.action=na.exclude)), silent=TRUE)
  if (isTRUE(class(cfit)=="try-error")){
    p.cfit <- NA
  } else {p.cfit <- cfit$tTable[,5]}
  
  Result <- rbind(Result, data.frame("SNP.No"=j, "SNP.name"=snp.name, 
                                     "Domi"=p.dfit[2], "ON"=p.dfit[3], "Age"=p.dfit[4], "Gendr"=p.dfit[5], "Dis.Dur"=p.dfit[6],
                                     "Additi"=p.afit[2], "ON2"=p.afit[3], "Age2"=p.afit[4], "Gendr2"=p.afit[5], "Dis.Dur2"=p.afit[6],
                                     "Recess"=p.rfit[2], "ON3"=p.rfit[3], "Age3"=p.rfit[4], "Gendr3"=p.rfit[5], "Dis.Dur3"=p.rfit[6],
                                     "Codo1"=p.cfit[2], "Codo2"=p.cfit[3], "ON4"=p.cfit[4], "Age4"=p.cfit[5], "Gendr4"=p.cfit[6], "Dis.Dur4"=p.cfit[7],
                                     stringsAsFactors = F))
  
  j <- j+1
}
write.csv(Result, filename)
#

rm(ad, co, do, re, snp.name, snp, vol, p.afit, p.cfit, p.dfit, p.rfit, i, j, afit, cfit, dfit, rfit, Result, filename, lme_data)
#




#########################################################################################################################=====
#####  Linear mixed effect model, adjusting for inter-eye correlation  (with optic neuritis information)  #####
#####  WITH HLA-DRB
library(nlme)
library(SNPassoc)

setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/Collabo_Calabresi/ReAnalysis/MSChip_results_withHLA")

OCT_dataset <- EPIC_OCT_new2_snp_slope
i <- 36     # genotype (36 ~ 55)
colnames(OCT_dataset)[i]

###
# 1. MAC
hist(EPIC_OCT_new2_snp_slope$MAC_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$MAC_Slope < 1.5 & EPIC_OCT_new2_snp_slope$MAC_Slope > -1.5),]
hist(OCT_dataset$MAC_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 383 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$MAC_Slope),]
vol <- as.numeric(as.character(lme_data$MAC_Slope))
filename <- c("EPIC_OCT_new_slope-lme_MAC.csv")


# 2. GCC
hist(EPIC_OCT_new2_snp_slope$GCC_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$GCC_Slope > -5),]
hist(OCT_dataset$GCC_Slope, breaks=100)
length(unique(OCT_dataset$EPIC_ID)) # 348 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$GCC_Slope),]
vol <- as.numeric(as.character(lme_data$GCC_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_GCC.csv")


# 3. GCL
hist(EPIC_OCT_new2_snp_slope$GCL_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$GCL_Slope > -2),]
hist(OCT_dataset$GCL_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 384 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$GCL_Slope),]
vol <- as.numeric(as.character(lme_data$GCL_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_GCL.csv")


# 4. INL
hist(EPIC_OCT_new2_snp_slope$INL_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$INL_Slope > -6e-15 & EPIC_OCT_new2_snp_slope$INL_Slope < 6e-15),]
hist(OCT_dataset$INL_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 370 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$INL_Slope),]
vol <- as.numeric(as.character(lme_data$INL_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_INL.csv")


# 5. IPL
hist(EPIC_OCT_new2_snp_slope$IPL_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$IPL_Slope > -2),]
hist(OCT_dataset$IPL_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 384 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$IPL_Slope),]
vol <- as.numeric(as.character(lme_data$IPL_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_IPL.csv")


# 6. pRNFL
hist(EPIC_OCT_new2_snp_slope$pRNFL_Slope, breaks=100)
OCT_dataset <- EPIC_OCT_new2_snp_slope[(EPIC_OCT_new2_snp_slope$pRNFL_Slope > -50 & EPIC_OCT_new2_snp_slope$pRNFL_Slope < 50),]
hist(OCT_dataset$pRNFL_Slope, breaks=100)     # 
length(unique(OCT_dataset$EPIC_ID)) # 361 patients

lme_data <- OCT_dataset[!is.na(OCT_dataset$pRNFL_Slope),]
vol <- as.numeric(as.character(lme_data$pRNFL_Slope))             # volume variable                   ##########
filename <- c("EPIC_OCT_new_slope-lme_pRNFL.csv")


###
Result <- data.frame("SNP.No"=numeric(), "SNP.name"=character(), 
                     "Domi"=numeric(), "HLA"=numeric(), "ON"=numeric(), "Age"=numeric(), "Gendr"=numeric(), "Dis.Dur"=numeric(), "Do+HLA"=numeric(),
                     "Additi"=numeric(), "HLA2"=numeric(), "ON2"=numeric(), "Age2"=numeric(), "Gendr2"=numeric(), "Dis.Dur2"=numeric(), "Ad+HLA2"=numeric(),
                     "Recess"=numeric(), "HLA3"=numeric(), "ON3"=numeric(), "Age3"=numeric(), "Gendr3"=numeric(), "Dis.Dur3"=numeric(), "Re+HLA3"=numeric(),
                     "Codo1"=numeric(), "Codo2"=numeric(), "HLA4"=numeric(), "ON4"=numeric(), "Age4"=numeric(), "Gendr4"=numeric(), "Dis.Dur4"=numeric(), "Codo1+HLA4"=numeric(), "Codo2+HLA4"=numeric(), 
                     stringsAsFactors = F)
j <- 1
for (i in 36:55){
  snp.name <- colnames(lme_data)[i]
  snp <- lme_data[,i]               # genotype
  
  do <- dominant(snp(snp, sep=""))  #dominant model
  ad <- additive(snp(snp, sep=""))  #additive model
  re <- recessive(snp(snp, sep=""))  #recessive model
  co <- codominant(snp(snp, sep=""))  #codominant model = genotype
  
  ## some patients are diagnosed in same year, so including month.   
  dfit <- try(summary(lme(vol ~ do*DRBstat+ON+AgeAtExam+Gender+DiseaseDuration, random=~1|Vision.PIDN, data=lme_data, na.action=na.exclude)), silent=TRUE)
  if (isTRUE(class(dfit)=="try-error")){
    p.dfit <- NA
  } else {p.dfit <- dfit$tTable[,5]}
  
  afit <- try(summary(lme(vol ~ ad*DRBstat+ON+AgeAtExam+Gender+DiseaseDuration, random=~1|Vision.PIDN, data=lme_data, na.action=na.exclude)), silent=TRUE)
  if (isTRUE(class(afit)=="try-error")){
    p.afit <- NA
  } else {p.afit <- afit$tTable[,5]}
  
  rfit <- try(summary(lme(vol ~ re*DRBstat+ON+AgeAtExam+Gender+DiseaseDuration, random=~1|Vision.PIDN, data=lme_data, na.action=na.exclude)), silent=TRUE)
  if (isTRUE(class(rfit)=="try-error")){
    p.rfit <- NA
  } else {p.rfit <- rfit$tTable[,5]}
  
  cfit <- try(summary(lme(vol ~ co*DRBstat+ON+AgeAtExam+Gender+DiseaseDuration, random=~1|Vision.PIDN, data=lme_data, na.action=na.exclude)), silent=TRUE)
  if (isTRUE(class(cfit)=="try-error")){
    p.cfit <- NA
  } else {p.cfit <- cfit$tTable[,5]}
  
  Result <- rbind(Result, data.frame("SNP.No"=j, "SNP.name"=snp.name, 
                                     "Domi"=p.dfit[2], "HLA"=p.dfit[3], "ON"=p.dfit[4], "Age"=p.dfit[5], "Gendr"=p.dfit[6], "Dis.Dur"=p.dfit[7], "Do+HLA"=p.dfit[8],
                                     "Additi"=p.afit[2], "HLA2"=p.afit[3], "ON2"=p.afit[4], "Age2"=p.afit[5], "Gendr2"=p.afit[6], "Dis.Dur2"=p.afit[7], "Ad+HLA2"=p.afit[8],
                                     "Recess"=p.rfit[2], "HLA3"=p.rfit[3], "ON3"=p.rfit[4], "Age3"=p.rfit[5], "Gendr3"=p.rfit[6], "Dis.Dur3"=p.rfit[7], "Re+HLA3"=p.rfit[8],
                                     "Codo1"=p.cfit[2], "Codo2"=p.cfit[3], "HLA4"=p.cfit[4], "ON4"=p.cfit[5], "Age4"=p.cfit[6], "Gendr4"=p.cfit[7], "Dis.Dur4"=p.cfit[8], "Codo1+HLA4"=p.cfit[9], "Codo2+HLA4"=p.cfit[10],
                                     stringsAsFactors = F))
  j <- j+1
}
write.csv(Result, filename)
#

rm(ad, co, do, re, snp.name, snp, vol, p.afit, p.cfit, p.dfit, p.rfit, i, j, afit, cfit, dfit, rfit, Result, filename)
#