# ON + Gender + AgeAtExam + DiseaseDuration + 
library(lmerTest)
setwd("C:/Users/kicheol/Desktop/test OCT slopes")
OCT_matrices <- OCT_matrices_longitud

# MAC
flme <- lmer(MAC ~ VisitSeq + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(MAC=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- slopeTable
# GCC
flme <- lmer(GCC ~ VisitSeq + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(GCC=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)
# GCL
flme <- lmer(GCL ~ VisitSeq + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(GCL=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)
# IPL
flme <- lmer(IPL ~ VisitSeq + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(IPL=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)
# pRNFL
flme <- lmer(pRNFL ~ VisitSeq + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(pRNFL=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)

OCT_slope$Cov <- "None"
write.csv(OCT_slope, "OCT_slope_Cov-none.csv")

slopes_covar <- OCT_slope


#########################################################################################################################
library(ggplot2)
library(readr)
slopes_covar <- read_csv("OCT_slope_Cov-AgeExam_ed.csv")
slopes_covar <- rbind(slopes_covar, read_csv("OCT_slope_Cov-DisDur_ed.csv"))
slopes_covar <- rbind(slopes_covar, read_csv("OCT_slope_Cov-Gender_ed.csv"))
slopes_covar <- rbind(slopes_covar, read_csv("OCT_slope_Cov-ON_ed.csv"))
slopes_covar <- rbind(slopes_covar, read_csv("OCT_slope_Cov-none.csv"))

plotData <- slopes_covar[slopes_covar$MAC.VisitSeq < 50,]
ggplot(plotData, aes(MAC.VisitSeq, fill = Cov)) + geom_density(alpha = 0.2)

plotData <- slopes_covar[slopes_covar$GCL.VisitSeq > -0.5,]
ggplot(plotData, aes(GCL.VisitSeq, fill = Cov)) + geom_density(alpha = 0.2)

ggplot(slopes_covar, aes(IPL.VisitSeq, fill = Cov)) + geom_density(alpha = 0.2)

ggplot(slopes_covar, aes(pRNFL.VisitSeq, fill = Cov)) + geom_density(alpha = 0.2)

plotData <- slopes_covar[slopes_covar$GCC.VisitSeq > -1,]
ggplot(plotData, aes(GCC.VisitSeq, fill = Cov)) + geom_density(alpha = 0.2)
#



#########################################################################################################################
# GCC
flme <- lmer(GCC ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(GCC=coef(flme)$EPIC_ID)
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

# INL
flme <- lmer(INL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), OCT_matrices)
slopeTable <- data.frame(INL=coef(flme)$EPIC_ID)
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPICID","Eye")
slopeTable <- cbind(slopeTable, df)
OCT_slope <- merge(OCT_slope, slopeTable, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"), all.x=TRUE)





