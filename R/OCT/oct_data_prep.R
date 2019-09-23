#####  OCT  - calculate slope -> test using logistic regression model  #####

##### OCT data preparation #####
### prepare OCT data with metadata (NAs are removed by manual)
OCT_matrices <- read_csv("C:/Users/kicheol/Desktop/Projects/Phenotype/raw.data/UCSF/OCT/EPIC_OCT_matrices_201801+metadata(filtered).csv")
EPIC_OCT_metadata <- read_csv("C:/Users/kicheol/Desktop/Projects/Phenotype/raw.data/UCSF/OCT/EPIC_OCT_metadata(filtered).csv")

### add DMT information to data
EPIC_OCT_metadata_DMTs <- read_csv("C:/Users/kicheol/Desktop/Projects/Phenotype/raw.data/UCSF/OCT/EPIC_ORIGINS_with_OCT_metadata(20180228)-DMTs.csv")
OCT_matrices_dmt <- merge(OCT_matrices, EPIC_OCT_metadata_DMTs, by.x=c("EPIC_ID","ExamDate"), by.y=c("EPICID","VisitDate2"), all.x=TRUE)
write.csv(OCT_matrices_dmt, "C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/OCT/RData/OCT_matrices_dmt.csv")

  # removed #VALUE and NAs
OCT_matrices_edited <- read_csv("C:/Users/kicheol/Desktop/Projects/Phenotype/R_ImgPheno/OCT/RData/OCT_matrices_edited.csv", col_types = cols(X1 = col_skip()))

  # remove baseline
patient_list <- unique(OCT_matrices_edited$EPIC_ID)
OCT_matrices_longitud <- data.frame()
for (i in 1:length(patient_list)){
  tmp_L <- OCT_matrices_edited[OCT_matrices_edited$EPIC_ID == patient_list[i] & OCT_matrices_edited$Eye == "L",]
  tmp_R <- OCT_matrices_edited[OCT_matrices_edited$EPIC_ID == patient_list[i] & OCT_matrices_edited$Eye == "R",]
  if (nrow(tmp_L) > 1){OCT_matrices_longitud <- rbind(OCT_matrices_longitud, tmp_L)}
  if (nrow(tmp_R) > 1){OCT_matrices_longitud <- rbind(OCT_matrices_longitud, tmp_R)}
}
write.csv(OCT_matrices_longitud, "OCT_matrices_longitudinal.csv")


### check histogram for each OCT values
# MAC, MAClog GCC, GCL, INL, IPL, pRNFL
hist(OCT_matrices_longitud$MAC, breaks=100)
hist(OCT_matrices_longitud$GCC, breaks=100)
hist(OCT_matrices_longitud$GCL, breaks=100)
hist(OCT_matrices_longitud$INL, breaks=100)
hist(OCT_matrices_longitud$IPL, breaks=100)
hist(OCT_matrices_longitud$pRNFL, breaks=100)                 # What is pRNFL=0 : it doesn't exist, remove it!!



### check number of samples
length(unique(OCT_matrices_longitud[!is.na(OCT_matrices_longitud$MAC),]$EPIC_ID))                # MAC = 430
length(unique(OCT_matrices_longitud[!is.na(OCT_matrices_longitud$GCC),]$EPIC_ID))                # GCC = 387
length(unique(OCT_matrices_longitud[!is.na(OCT_matrices_longitud$GCL),]$EPIC_ID))                # GCL = 430
length(unique(OCT_matrices_longitud[!is.na(OCT_matrices_longitud$INL),]$EPIC_ID))                # INL = 429
length(unique(OCT_matrices_longitud[!is.na(OCT_matrices_longitud$IPL),]$EPIC_ID))                # IPL = 430
length(unique(OCT_matrices_longitud[!is.na(OCT_matrices_longitud$pRNFL) & OCT_matrices_longitud$pRNFL > 0,]$EPIC_ID))    # pRNFL = 424


### change data type to factor (ON, Eye, gender, DiseaseCourse, DRBstat1, DRBstat2, Treatment)
OCT_matrices_longitud$ON <- factor(OCT_matrices_longitud$ON)
OCT_matrices_longitud$Eye <- factor(OCT_matrices_longitud$Eye)
OCT_matrices_longitud$Gender <- factor(OCT_matrices_longitud$Gender)
OCT_matrices_longitud$DiseaseCourse.x <- factor(OCT_matrices_longitud$DiseaseCourse.x)
OCT_matrices_longitud$DRBstat <- factor(OCT_matrices_longitud$DRBstat)
OCT_matrices_longitud$DRBstat2 <- factor(OCT_matrices_longitud$DRBstat2)
OCT_matrices_longitud$DMTs_at_Visit <- factor(OCT_matrices_longitud$DMTs_at_Visit)

write.csv(OCT_matrices_longitud, "OCT_matrices_longitudinal.csv")
save(OCT_matrices_longitud, file="OCT_matrices_longitudinal.RData")




### preparing genotype (194 SNPs)
EPIC_194snp <- read_csv("C:/Users/kicheol/Desktop/Projects/Phenotype/raw.data/UCSF/EPIC_genetics_2015(2).csv")
#save(EPIC_194snp, file="EPIC_194SNPs.RData")

## merge dataset with SNPs
EPIC_OCT_snp <- merge(EPIC_OCT, EPIC_194snp, by.x="EPID", by.y="EPICID")

