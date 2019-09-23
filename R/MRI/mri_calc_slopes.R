#####  MRI  - calculate slope -> test using logistic regression model  #####
## *** Do not use factor variable to calculate slopes

##### calculate slopes using simple linear model #####
EPIC_MRI_Normalized

vols <- c("NBV","NGMV","NWMV","NT2V","NVCSF","NCGM","New_T2_count","PBVC_since_baseline","Gad_Lesions")
id_list <- unique(EPIC_MRI_Normalized$EPICID)

# NBV
tmp_volume <- EPIC_MRI_Normalized
slopeTable <- data.frame()
for (i in 1:length(id_list)){
  indiv <- subset(tmp_volume, tmp_volume$EPICID == id_list[i])
  indiv <- indiv[!is.na(indiv$NBV),]
  if (nrow(indiv) > 1){
    coef <- coef(lm(NBV ~ Time, data=indiv))[2]
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NBV_Slop"=coef, "NBV_n"=nrow(indiv)))
    rm(indiv)
  } else {
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NBV_Slop"=NA, "NBV_n"=nrow(indiv)))
    rm(indiv)
  }
}
MRI_slope <- slopeTable
rm(slopeTable, tmp_volume)

# NGMV
tmp_volume <- EPIC_MRI_Normalized
slopeTable <- data.frame()
for (i in 1:length(id_list)){
  indiv <- subset(tmp_volume, tmp_volume$EPICID == id_list[i])
  indiv <- indiv[!is.na(indiv$NGMV),]
  if (nrow(indiv) > 1){
    coef <- coef(lm(NGMV ~ Time, data=indiv))[2]
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NGMV_Slop"=coef, "NGMV_n"=nrow(indiv)))
    rm(indiv)
  } else {
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NGMV_Slop"=NA, "NGMV_n"=nrow(indiv)))
    rm(indiv)
  }
}
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, tmp_volume)

# NWMV
tmp_volume <- EPIC_MRI_Normalized
slopeTable <- data.frame()
for (i in 1:length(id_list)){
  indiv <- subset(tmp_volume, tmp_volume$EPICID == id_list[i])
  indiv <- indiv[!is.na(indiv$NWMV),]
  if (nrow(indiv) > 1){
    coef <- coef(lm(NWMV ~ Time, data=indiv))[2]
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NWMV_Slop"=coef, "NWMV_n"=nrow(indiv)))
    rm(indiv)
  } else {
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NWMV_Slop"=NA, "NWMV_n"=nrow(indiv)))
    rm(indiv)
  }
}
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, tmp_volume)

# NT2V_log
tmp_volume <- EPIC_MRI_Normalized
slopeTable <- data.frame()
for (i in 1:length(id_list)){
  indiv <- subset(tmp_volume, tmp_volume$EPICID == id_list[i])
  indiv <- indiv[indiv$NT2V_log > 0,]
  indiv <- indiv[!is.na(indiv$NT2V_log),]
  if (nrow(indiv) > 1){
    coef <- coef(lm(NT2V_log ~ Time, data=indiv))[2]
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NT2Vlog_Slop"=coef, "NT2Vlog_n"=nrow(indiv)))
    rm(indiv)
  } else {
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NT2Vlog_Slop"=NA, "NT2Vlog_n"=nrow(indiv)))
    rm(indiv)
  }
}
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, tmp_volume)

# NVCSF_log
tmp_volume <- EPIC_MRI_Normalized
slopeTable <- data.frame()
for (i in 1:length(id_list)){
  indiv <- subset(tmp_volume, tmp_volume$EPICID == id_list[i])
  indiv <- indiv[indiv$NT2V_log > 0,]
  indiv <- indiv[!is.na(indiv$NVCSF_log),]
  if (nrow(indiv) > 1){
    coef <- coef(lm(NVCSF_log ~ Time, data=indiv))[2]
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NVCSFlog_Slop"=coef, "NVCSFlog_n"=nrow(indiv)))
    rm(indiv)
  } else {
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NVCSFlog_Slop"=NA, "NVCSFlog_n"=nrow(indiv)))
    rm(indiv)
  }
}
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, tmp_volume)

# NCGM
tmp_volume <- EPIC_MRI_Normalized
slopeTable <- data.frame()
for (i in 1:length(id_list)){
  indiv <- subset(tmp_volume, tmp_volume$EPICID == id_list[i])
  indiv <- indiv[!is.na(indiv$NCGM),]
  if (nrow(indiv) > 1){
    coef <- coef(lm(NCGM ~ Time, data=indiv))[2]
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NCGM_Slop"=coef, "NCGM_n"=nrow(indiv)))
    rm(indiv)
  } else {
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "NCGM_Slop"=NA, "NCGM_n"=nrow(indiv)))
    rm(indiv)
  }
}
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, tmp_volume)

# PBVC_since_baseline
tmp_volume <- EPIC_MRI_Normalized
slopeTable <- data.frame()
for (i in 1:length(id_list)){
  indiv <- subset(tmp_volume, tmp_volume$EPICID == id_list[i])
  indiv <- indiv[!is.na(indiv$PBVC_since_baseline),]
  if (nrow(indiv) > 1){
    coef <- coef(lm(PBVC_since_baseline ~ Time, data=indiv))[2]
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "PBVC_since_base_Slop"=coef, "PBVC_since_base_n"=nrow(indiv)))
    rm(indiv)
  } else {
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv$EPICID[1], "PBVC_since_base_Slop"=NA, "PBVC_since_base_n"=nrow(indiv)))
    rm(indiv)
  }
}
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, tmp_volume)

# 
MRI_slope$Nsum <- MRI_slope$NBV_n + MRI_slope$NGMV_n + MRI_slope$NWMV_n + MRI_slope$NT2Vlog_n + MRI_slope$NVCSFlog_n + MRI_slope$NCGM_n
MRI_slope <- MRI_slope[!is.na(MRI_slope$EPICID),]
MRI_slope <- MRI_slope[MRI_slope$Nsum > 6,]
EPIC_MRI_slope <- MRI_slope
write.csv(EPIC_MRI_slope, "EPIC_MRI_slope.csv")




################################################################################################################

##### calculate slopes using linear mixed effect model adjusted for the effects of age at first MRI, sex #####
# covariates: gender, presence of HLA-DRB*15:01 allele, age at first MRI, disease duration, disease modifying treatment at exam
library(lme4)
library(readr)
setwd("C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/MRI/RData")

vols <- c("NBV","NGMV","NWMV","NT2V","NVCSF","NCGM","New_T2_count","PBVC_since_baseline","Gad_Lesions")
#EPIC_MRI_Vol <- read_csv("EPIC_MRI_Normalized_withMetadata_removed-NA&baseline.csv", col_types = cols(X1 = col_skip()))
#load(file="EPIC_MRI_volume+metadata+dmt+194genotypes.RData")
EPIC_MRI_Vol <- EPIC_MRI_dmt

# NBV
slopeTable <- data.frame()
flme <- lmer(NBV ~ Time + DRBstat1 + Gender + AgeAtExam + DiseaseDuration + Treatment + (Time|EPICID), EPIC_MRI_Vol)
slopeTable <- data.frame(NBV=coef(flme)$EPICID)
MRI_slope <- slopeTable
rm(slopeTable, flme)

# NGMV
slopeTable <- data.frame()
flme <- lmer(NGMV ~ Time + DRBstat1 + Gender + AgeAtExam + DiseaseDuration + Treatment + (Time|EPICID), EPIC_MRI_Vol)
slopeTable <- data.frame(NGMV=coef(flme)$EPICID)
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, flme)

# NWMV
slopeTable <- data.frame()
flme <- lmer(NWMV ~ Time + DRBstat1 + Gender + AgeAtExam + DiseaseDuration + Treatment + (Time|EPICID), EPIC_MRI_Vol)
slopeTable <- data.frame(NWMV=coef(flme)$EPICID)
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, flme)

# NT2V_log (NT2V+1)
slopeTable <- data.frame()
flme <- lmer(log(NT2V+1) ~ Time + DRBstat1 + Gender + AgeAtExam + DiseaseDuration + Treatment + (Time|EPICID), EPIC_MRI_Vol)
slopeTable <- data.frame(NT2Vlog=coef(flme)$EPICID)
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, flme)

# NVCSF_log (NVCSF+1)
slopeTable <- data.frame()
flme <- lmer(log(NVCSF+1) ~ Time + DRBstat1 + Gender + AgeAtExam + DiseaseDuration + Treatment + (Time|EPICID), EPIC_MRI_Vol)
slopeTable <- data.frame(NVCSFlog=coef(flme)$EPICID)
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, flme)

# NCGM
slopeTable <- data.frame()
flme <- lmer(NCGM ~ Time + DRBstat1 + Gender + AgeAtExam + DiseaseDuration + Treatment + (Time|EPICID), EPIC_MRI_Vol)
slopeTable <- data.frame(NCGM=coef(flme)$EPICID)
MRI_slope <- cbind(MRI_slope, slopeTable)
rm(slopeTable, flme)

#
MRI_slope$EPICID <- row.names(MRI_slope)
MRI_slope_lme <- MRI_slope[,c("EPICID","NBV..Intercept.","NBV.Time","NGMV..Intercept.","NGMV.Time","NWMV..Intercept.","NWMV.Time","NT2Vlog..Intercept.","NT2Vlog.Time","NVCSFlog..Intercept.","NVCSFlog.Time","NCGM..Intercept.","NCGM.Time")]
head(MRI_slope_lme)
write.csv(MRI_slope_lme, "EPIC_MRI_slope_usingLME.csv")

################################################################################################################

# remove extra EPICID column
#EPIC_MRI_slope <- read_csv("EPIC_MRI_slope.csv")
#save(EPIC_MRI_slope, file="EPIC_MRI_slopes.RData")


hist(EPIC_MRI_slope$NBV_Slop, breaks=100)         # select > -300 | < 250
hist(EPIC_MRI_slope$NGMV_Slop, breaks=100)
hist(EPIC_MRI_slope$NWMV_Slop, breaks=100)        # select > -150 | < 150
hist(EPIC_MRI_slope$NT2Vlog_Slop, breaks=100)     # select > -1.5
hist(EPIC_MRI_slope$NVCSFlog_Slop, breaks=100)
hist(EPIC_MRI_slope$NCGM_Slop, breaks=100)


### merge dataset
EPIC_MRI_snp <- merge(MRI_slope_lme, EPIC_metadata, by.x="EPICID", by.y="EPICID")
EPIC_MRI_snp <- merge(EPIC_MRI_snp, EPIC_194snp, by.x="EPID", by.y="EPICID")
save(EPIC_MRI_snp, file="EPIC_MRI_slopeLME+metadata+194genotypes.RData")
write.csv(EPIC_MRI_snp, file="EPIC_MRI_slopeLME+metadata+194genotypes.csv")
#EPIC_MRI <- merge(MRI_slope_lme, EPIC_metadata, by.x="EPICID", by.y="EPICID")
#write.csv(EPIC_MRI, "EPIC_MRI_slopeLME+metadata.csv")

## all NA samples are removed manually
#EPIC_MRI <- merge(MRI_slope_lme, EPIC_194snp, by.x="EPICID", by.y="EPICID")
#save(EPIC_MRI, file="EPIC_MRI_slopeLME+194genotypes.RData")
#write.csv(EPIC_MRI, file="EPIC_MRI_slopeLME+194genotypes.csv")



##### histogram for cross-sectional analyssi
hist(EPIC_MRI_baseline$NBV, breaks=100)         # select > -300 | < 250
hist(EPIC_MRI_baseline$NGMV, breaks=100)
hist(EPIC_MRI_baseline$NWMV, breaks=100)        # select > -150 | < 150
hist(EPIC_MRI_baseline$NT2V, breaks=100)     # select > -1.5
hist(EPIC_MRI_baseline$NVCSF, breaks=100)
hist(EPIC_MRI_baseline$NCGM, breaks=100)
#
EPIC_MRI_baseline$NT2Vlog <- log(EPIC_MRI_baseline$NT2V)
EPIC_MRI_baseline$NVCSFlog <- log(EPIC_MRI_baseline$NVCSF)
hist(EPIC_MRI_baseline$NT2Vlog, breaks=100)
hist(EPIC_MRI_baseline$NVCSFlog, breaks=100)

write.csv(EPIC_MRI_baseline, "EPIC_MRI_baseline.csv")

EPIC_MRI_base_meta <- merge(EPIC_MRI_baseline, EPIC_metadata, by.x="EPICID", by.y="EPICID")
EPIC_MRI_base_meta_snp <- merge(EPIC_MRI_base_meta, EPIC_194snp, by.x="EPID", by.y="EPICID")
save(EPIC_MRI_base_meta_snp, file="EPIC_MRI_baseline+metadata+194genotypes.RData")



##### plot for single patient
library(ggplot2)
id_list <- unique(EPIC_MRI_Normalized$EPICID)
i=1
tmp_volume <- EPIC_MRI_Normalized
indiv <- subset(tmp_volume, tmp_volume$EPICID == "EPIC0001")
indiv <- indiv[!is.na(indiv$NBV),]
coef <- coef(lm(NBV ~ Time, data=indiv))

ggplot(indiv, aes(x=Time, y=NBV)) + geom_point() + geom_line() + 
  geom_abline(aes(intercept=coef[1], slope=coef[2]), data=indiv, color="Red") +
  ylab("Normalized Total Brain volume") + xlab("Baseline - 2 Year - 5 Year")


