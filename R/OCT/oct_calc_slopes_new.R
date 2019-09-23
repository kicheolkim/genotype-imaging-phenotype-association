#####  MRI  - calculate slope -> test using logistic regression model  #####
# https://stats.stackexchange.com/questions/122009/extracting-slopes-for-cases-from-a-mixed-effects-model-lme4
## *** Do not use factor variable to calculate slopes


##### check basic numbers (Gender, AgeAtExam, DiseaseCourse, DiseaseDuration) #####
OCT_matrices_longitud
OCT_matrices <- subset(OCT_matrices_longitud, OCT_matrices_longitud$Eye == "L")  # for L
OCT_matrices <- subset(OCT_matrices_longitud, OCT_matrices_longitud$Eye == "R")  # for R

length(unique(OCT_matrices_longitud$EPIC_ID))                # total n=437
colnames(OCT_matrices)[1:20]

length(unique(subset(OCT_matrices_longitud, Gender=="M")$EPIC_ID))     # n=138
length(unique(subset(OCT_matrices_longitud, Gender=="F")$EPIC_ID))     # n=299
# CIS PP PR RR SP UNC UNK
length(unique(subset(OCT_matrices_longitud, DiseaseCourse.x=="CIS")$EPIC_ID))  # n=32
length(unique(subset(OCT_matrices_longitud, DiseaseCourse.x=="PP")$EPIC_ID))   # n=15
length(unique(subset(OCT_matrices_longitud, DiseaseCourse.x=="RR")$EPIC_ID))   # n=336
length(unique(subset(OCT_matrices_longitud, DiseaseCourse.x=="SP")$EPIC_ID))   # n=67
length(unique(subset(OCT_matrices_longitud, DiseaseCourse.x=="UNC")$EPIC_ID))  # n=16

mean(OCT_matrices_longitud$AgeAtExam)                  # 47.96019
sd(OCT_matrices_longitud$AgeAtExam)                    # +- 10.92244
mean(OCT_matrices_longitud$DiseaseDuration)            # 14.44221
sd(OCT_matrices_longitud$DiseaseDuration)              # +- 9.25024



##### calculate slopes using linear mixed effect model adjusted for the effects of age at first MRI, sex #####
# covariates: gender, age at exam, disease duration, disease modifying treatment at exam 
# I did not used HLA-DRB*15:01 because data is not completed
library(lmerTest)
library(readr)
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/OCT/RData")
# MAC, GCC, GCL, INL, IPL, pRNFL
OCT_matrices_longitud

#OCT_matrices <- subset(OCT_matrices_longitud, OCT_matrices_longitud$Eye == "L")  # for L
#OCT_matrices <- subset(OCT_matrices_longitud, OCT_matrices_longitud$Eye == "R")  # for R

hist(OCT_matrices$MAC, breaks=100)
hist(OCT_matrices[OCT_matrices$MAC < 400,]$MAC, breaks=100)       # MAC: remove outlier for Right eye only
hist(OCT_matrices$GCC, breaks=100)                                # GCC: use all  (skewed for left eye, how can I treat it?)
hist(OCT_matrices$GCL, breaks=100)                                # GCL: use all  (skewed... how can I treat it?)
hist(OCT_matrices$INL, breaks=100)                                # INL: use all
hist(OCT_matrices$IPL, breaks=100)                                # IPL: use all
hist(OCT_matrices$pRNFL, breaks=100)         # for L and R (What is pRNFL=0 : it doesn't exist, remove it!!)
hist(OCT_matrices[OCT_matrices$pRNFL > 0,]$pRNFL, breaks=100)
hist(OCT_matrices[OCT_matrices$pRNFL > 0 & OCT_matrices$pRNFL < 160,]$pRNFL, breaks=100)     # pRNFL: remove 0 values and outliers for Left eye
hist(OCT_matrices[OCT_matrices$pRNFL > 0 & OCT_matrices$pRNFL < 150,]$pRNFL, breaks=100)     # pRNFL: remove 0 values and outliers for Right eye


OCT_matrices <- OCT_matrices_longitud
# MAC
flme <- lmer(MAC ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(MAC=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- slopeTable
# GCC
flme <- lmer(GCC ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(GCC=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)
# GCL
flme <- lmer(GCL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(GCL=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)
# INL
flme <- lmer(INL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(INL=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)
# IPL
flme <- lmer(IPL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(IPL=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)
# pRNFL
flme <- lmer(pRNFL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(pRNFL=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)

write.csv(OCT_slope, "OCT_slope_LR.csv")
OCT_slope_LR <- OCT_slope



# remove extra EPICID column
#EPIC_MRI_slope <- read_csv("EPIC_MRI_slope.csv")
#save(EPIC_MRI_slope, file="EPIC_MRI_slopes.RData")

# MAC, GCC, GCL, INL, IPL, pRNFL
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


# select column
OCT_slope_LR_all <- OCT_slope_LR
OCT_slope_LR <- OCT_slope_LR[,c(1,2,4,10,16,22,28,34)]


### merge dataset with metadata and SNPs
EPIC_OCT_slope_snp <- merge(OCT_slope_LR, EPIC_194snp, by.x="EPICID", by.y="EPICID")
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/OCT/RData")
save(EPIC_OCT_slope_snp, file="EPIC_OCT_slopeLME+194genotypes.RData")
write.csv(EPIC_OCT_slope_snp, file="EPIC_OCT_slopeLME+194genotypes.csv")

OCT_slopes_linear_NoAdjust <- read_csv("RData/EPIC_OCT_slopes_linear_noAdjusted.csv")
OCT_slopes_linear_NoAdjust_SNPs <- merge(OCT_slopes_linear_NoAdjust, EPIC_194snp, by.x="EPICID", by.y="EPICID")
write.csv(OCT_slopes_linear_NoAdjust_SNPs, file="EPIC_OCT_slopes_linear_noAdjusted+194genotypes.csv")




##### plot for single patient
library(ggplot2)
OCT_matrices <- subset(OCT_matrices_longitud, OCT_matrices_longitud$Eye == "L")
OCT_matrices <- subset(OCT_matrices_longitud, OCT_matrices_longitud$Eye == "R")
id_list <- unique(OCT_matrices$EPICID)
tmp_volume <- OCT_matrices

indiv <- subset(tmp_volume, tmp_volume$EPIC_ID == "EPIC0002")
indiv <- indiv[!is.na(indiv$MAC),]
coef <- coef(lm(MAC ~ VisitSeq, data=indiv))
ggplot(indiv, aes(x=VisitSeq, y=MAC)) + geom_point() + geom_line() + 
  geom_abline(aes(intercept=coef[1], slope=coef[2]), data=indiv, color="Red") +
  ylab("MAC (macular) thickness") + xlab("Visit sequence") + ggtitle("EPIC0002, left eye")

indiv <- subset(tmp_volume, tmp_volume$EPIC_ID == "EPIC0005")
indiv <- indiv[!is.na(indiv$pRNFL),]
coef <- coef(lm(pRNFL ~ VisitSeq, data=indiv))
ggplot(indiv, aes(x=VisitSeq, y=pRNFL)) + geom_point() + geom_line() + 
  geom_abline(aes(intercept=coef[1], slope=coef[2]), data=indiv, color="Red") +
  ylab("pRNFL thickness") + xlab("Visit sequence") + ggtitle("EPIC0005, left eye")



