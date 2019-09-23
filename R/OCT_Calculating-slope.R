########## calculate slopes #####
tmp_slopeCalc <- EPIC_OCT_7SNPs
tmp_slopeCalc <- EPIC_OCT_new2_snp
tmp_slopeCalc$VisitSeq <- as.numeric(as.character(tmp_slopeCalc$VisitSeq))

hist(log(tmp_slopeCalc$MAC))
hist(tmp_slopeCalc$GCC)
hist(log(tmp_slopeCalc$GCC))
hist(tmp_slopeCalc$GCL)
hist(tmp_slopeCalc$INL)
hist(log(tmp_slopeCalc$INL))
hist(tmp_slopeCalc$IPL)
hist(log(tmp_slopeCalc$pRNFL))
hist(tmp_slopeCalc$pRNFL)

id_list <- unique(tmp_slopeCalc$EPIC_ID)

# MAC, GCC, GCL, INL, IPL, pRNFL
slopeTable <- data.frame("EPICID"=character(),"ExamDate"=character(),"Eye"=character(), "pRNFL_Slope"=numeric(), "Count"=numeric())  ########## Volume
for (i in 1:length(id_list)){
  indiv_table <- subset(tmp_slopeCalc, tmp_slopeCalc$EPIC_ID==id_list[i])
  indiv_table <- indiv_table[!is.na(indiv_table$pRNFL),]                                  ########## Volume
  indiv_L <- indiv_table[indiv_table$Eye=="L",]
  indiv_R <- indiv_table[indiv_table$Eye=="R",]
  # calculate slope
  if (nrow(indiv_L) >= 2){
    coef_L <- coef(lm(pRNFL ~ VisitSeq, data=indiv_L))[2]                                 ########## Volume
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv_L$EPIC_ID[1], "Eye"=indiv_L$Eye[1], "pRNFL_Slope"=coef_L, "Count"=nrow(indiv_L)))  ########## Volume
    rm(coef_L, indiv_L)
  }
  if (nrow(indiv_R) >= 2){
    coef_R <- coef(lm(pRNFL ~ VisitSeq, data=indiv_R))[2]                                 ########## Volume
    slopeTable <- rbind(slopeTable, data.frame("EPICID"=indiv_R$EPIC_ID[1],"Eye"=indiv_R$Eye[1], "pRNFL_Slope"=coef_R, "Count"=nrow(indiv_R)))  ########## Volume
    rm(coef_R, indiv_R)
  }
  rm(indiv_table)
}
OCT_Slopes <- merge(OCT_Slopes, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE, all.y=TRUE)
rm(slopeTable)


OCT_Slopes <- slopeTable
write.csv(OCT_Slopes, "EPIC_OCT_slopes_MSchip.csv")



#
i <- 3
indiv_table <- subset(tmp_slopeCalc, tmp_slopeCalc$EPIC_ID==id_list[i])
indiv_table <- indiv_table[!is.na(indiv_table$INL),]                                ########## Volume
indiv_table$INL <- log(indiv_table$INL)
indiv_L <- indiv_table[indiv_table$Eye=="L",]
tmp_fit <- lm(formula = INL ~ VisitSeq, data = indiv_L)
summary(tmp_fit)
plot(GCL ~ VisitSeq, data = indiv_L)


EPIC_OCT_7SNPs_slopes <- merge(EPIC_OCT_7SNPs_baseline, OCT_Slopes, by.x=c("EPIC_ID", "Eye"), by.y=c("EPICID","Eye"))
write.csv(EPIC_OCT_7SNPs_slopes, "EPIC_OCT_7SNPs_slopes.csv")
EPIC_OCT_7SNPs_slopes <- read.csv("C:/Users/kicheol/Desktop/Projects/Phenotype/Collabo_Calabresi/ReAnalysis/EPIC_OCT_7SNPs_slopes.csv")


EPIC_OCT_new2_snp_slope <- merge(EPIC_OCT_new2_snp, OCT_Slopes, by.x=c("EPIC_ID", "Eye"), by.y=c("EPICID","Eye"))
write.csv(EPIC_OCT_new2_snp_slope, "EPIC_OCT_mschip_slopes.csv")


#
rm(ad, co, do, re, snp.name, snp, vol, p.afit, p.cfit, p.dfit, p.rfit, i, j, afit, cfit, dfit, rfit, Result, filename)
sessionInfo()
#
