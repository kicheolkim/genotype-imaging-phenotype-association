library(ggplot2)

slope_original <- read_csv("C:/Users/kicheol/Desktop/OCT_slope.csv", col_types = cols(X1 = col_skip()))
slopeTable

ggplot(slope_original, aes(x=INL_Slope, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("original slope for OCT")

ggplot(slope_original, aes(x=GCL_Slope, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("original slope for OCT")



ggplot(OCT_slope_LR, aes(x=MAC.VisitSeq, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT")

ggplot(OCT_slope_LR, aes(x=pRNFL.VisitSeq, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT")


OCT_matrices <- OCT_matrices_longitud

ggplot(OCT_matrices[OCT_matrices$GCL < 400,], aes(x=GCL, color=Gender)) + 
  geom_histogram(aes(fill=Gender), alpha=0.3, position="identity", bins = 100) +
  ggtitle("original OCT values")

ggplot(OCT_matrices[OCT_matrices$GCL < 400,], aes(x=GCL, color=Eye, fill=ON)) + 
  geom_histogram(aes(fill=ON), alpha=0.3, position="identity", bins = 100) +
  ggtitle("original OCT values")

ggplot(OCT_matrices[OCT_matrices$GCL < 400,], aes(x=GCL, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("original OCT values")

ggplot(OCT_matrices[OCT_matrices$GCL < 400,], aes(x=GCL, color=ON)) + 
  geom_histogram(aes(fill=ON), alpha=0.3, position="identity", bins = 100) +
  ggtitle("original OCT values")



##################
tmp <- tmp_LR[tmp_LR$EPICID %in% unique(OCT_matrices_longitud[OCT_matrices_longitud$ON == 0,]$EPIC_ID),]

tmp_L <- slopeTable
tmp_L$Eye <- "L"
tmp_R <- slopeTable
tmp_R$Eye <- "R"
tmp_LR <- rbind(tmp_L, tmp_R)
ggplot(tmp, aes(x=MAC.VisitSeq, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT")


ggplot(tmp_LR, aes(x=MAC.VisitSeq)) + 
  geom_histogram(alpha=0.5, position="identity", bins = 100) +
  ggtitle("adjusted slope with both eyes for OCT")


##################
OCT_matrices <- OCT_matrices_longitud
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$MAC) & OCT_matrices$MAC < 400)
flme <- lmer(MAC ~ VisitSeq * Eye + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID), tmp)
slopeTable <- data.frame(MAC=coef(flme)$EPIC_ID)
slopeTable$EPICID <- row.names(slopeTable)
rm(tmp)

ggplot(slopeTable, aes(x=MAC.VisitSeq)) + 
  geom_histogram(alpha=0.5, position="identity", bins = 100) +
  ggtitle("adjusted slope with both eyes for OCT")


# GCL - left eye
OCT_matrices <- subset(OCT_matrices_longitud, Eye == "L")  # for L
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$GCL))
flme <- lmer(GCL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID), tmp)
slopeTable <- data.frame(GCL=coef(flme)$EPIC_ID)
slopeTable$EPICID <- row.names(slopeTable); rm(tmp)
tmp_L <- slopeTable
tmp_L$Eye <- "L"

# GCL - left eye
OCT_matrices <- subset(OCT_matrices_longitud, Eye == "R")  # for R
#tmp <- subset(OCT_matrices, !is.na(OCT_matrices$GCL))
flme <- lmer(GCL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID), tmp)
slopeTable <- data.frame(GCL=coef(flme)$EPIC_ID)
slopeTable$EPICID <- row.names(slopeTable); rm(tmp)
tmp_R <- slopeTable
tmp_R$Eye <- "R"

# merge both eyes
tmp_LR <- rbind(tmp_L, tmp_R)
tmp_LR <- tmp_LR[order(tmp_LR$EPICID),]

###
ggplot(tmp_LR, aes(x=GCL.VisitSeq, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (with ON patients)")


###
# MAC - left eye
OCT_matrices <- subset(OCT_matrices_longitud, Eye == "L")  # for L
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$MAC))                           # for L
flme <- lmer(MAC ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID), tmp)
slopeTable <- data.frame(MAC=coef(flme)$EPIC_ID)
slopeTable$EPICID <- row.names(slopeTable); rm(tmp)
tmp_L <- slopeTable
tmp_L$Eye <- "L"

# MAC - right eye
OCT_matrices <- subset(OCT_matrices_longitud, Eye == "R")  # for R
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$MAC) & OCT_matrices$MAC < 400)  # for R
flme <- lmer(MAC ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID), tmp)
slopeTable <- data.frame(MAC=coef(flme)$EPIC_ID)
slopeTable$EPICID <- row.names(slopeTable); rm(tmp)
tmp_R <- slopeTable
tmp_R$Eye <- "R"

# merge both eyes
tmp_LR <- rbind(tmp_L, tmp_R)
tmp_LR <- tmp_LR[order(tmp_LR$EPICID),]

###
ggplot(tmp_LR, aes(x=MAC.VisitSeq, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (with ON patients)")


#####
ggplot(tmp_LR, aes(x=MAC.VisitSeq)) + 
  geom_histogram(alpha=0.5, position="identity", bins = 100) +
  ggtitle("adjusted slope with both eyes for OCT")



########
# left eye
OCT_matrices <- subset(OCT_matrices_longitud, Eye == "L" & ON == 0)  # for L
#OCT_matrices <- subset(OCT_matrices_longitud, Eye == "L")  # for L
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$pRNFL))                           # for L
flme <- lmer(pRNFL ~ VisitSeq + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID), tmp)
slopeTable <- data.frame(pRNFL=coef(flme)$EPIC_ID)
slopeTable$EPICID <- row.names(slopeTable); rm(tmp)
tmp_L <- slopeTable
tmp_L$Eye <- "L"

# right eye
OCT_matrices <- subset(OCT_matrices_longitud, Eye == "R" & ON == 0)  # for R
#OCT_matrices <- subset(OCT_matrices_longitud, Eye == "R")  # for R
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$pRNFL))  # for R
flme <- lmer(pRNFL ~ VisitSeq + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID), tmp)
slopeTable <- data.frame(pRNFL=coef(flme)$EPIC_ID)
slopeTable$EPICID <- row.names(slopeTable); rm(tmp)
tmp_R <- slopeTable
tmp_R$Eye <- "R"

# merge both eyes
tmp_LR <- rbind(tmp_L, tmp_R)
tmp_LR <- tmp_LR[order(tmp_LR$EPICID),]

##
ggplot(tmp_LR, aes(x=pRNFL.VisitSeq, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (without ON patients)")


#####
# MAC - both eyes
OCT_matrices <- OCT_matrices_longitud
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$MAC) & OCT_matrices$MAC < 400)
flme <- lmer(MAC ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq * Eye|EPIC_ID), tmp)
slopeTable <- data.frame(MAC=coef(flme)$EPIC_ID)
slopeTable$EPICID <- row.names(slopeTable)
head(slopeTable, 15)[1:8]
##




# GCL - both eyes
OCT_matrices <- OCT_matrices_longitud
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$GCL))
flme <- lmer(GCL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq + Eye|EPIC_ID), tmp)
slopeTable <- data.frame(GCL=coef(flme)$EPIC_ID); slopeTable$EPICID <- row.names(slopeTable)
head(slopeTable, 15)
##
ggplot(slopeTable, aes(x=GCL.VisitSeq)) + 
  geom_histogram(alpha=0.7, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (with ON patients)")


# GCL - both eyes
OCT_matrices <- OCT_matrices_longitud
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$GCL))
flme <- lmer(GCL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq * Eye|EPIC_ID), tmp)
slopeTable <- data.frame(GCL=coef(flme)$EPIC_ID); slopeTable$EPICID <- row.names(slopeTable)
head(slopeTable, 15)
##
ggplot(slopeTable, aes(x=GCL.VisitSeq)) + 
  geom_histogram(alpha=0.7, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (with ON patients)")




########################  NEW MODEL for slope calculation #####
# GCL - both eyes
OCT_matrices <- OCT_matrices_longitud
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$GCL))
flme <- lmer(GCL ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), tmp)
slopeTable <- data.frame(GCL=coef(flme)$EPIC_ID); slopeTable$EPICID <- row.names(slopeTable)
head(slopeTable, 12)
##
ggplot(slopeTable, aes(x=GCL.VisitSeq)) + 
  geom_histogram(alpha=0.7, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (with ON patients)")
###
df <- data.frame(matrix(unlist(strsplit(row.names(slopeTable), ":")), nrow=length(strsplit(row.names(slopeTable), ":")), byrow=T))
colnames(df) <- c("EPIC_ID","Eye")
slopeTable <- cbind(slopeTable, df)
ggplot(tmp_LR, aes(x=GCL.VisitSeq, color=X2)) + 
  geom_histogram(aes(fill=X2), alpha=0.3, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (with ON patients)")

# MAC - both eyes
OCT_matrices <- OCT_matrices_longitud
tmp <- subset(OCT_matrices, !is.na(OCT_matrices$MAC))
flme <- lmer(MAC ~ VisitSeq + ON + Gender + AgeAtExam + DiseaseDuration + (VisitSeq|EPIC_ID:Eye), tmp)
slopeTable <- data.frame(MAC=coef(flme)$EPIC_ID); slopeTable$EPICID <- row.names(slopeTable)
head(slopeTable, 12)
##
ggplot(slopeTable[slopeTable$MAC.VisitSeq < 50,], aes(x=MAC.VisitSeq)) + 
  geom_histogram(alpha=0.7, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (with ON patients)")
###
str <- strsplit(slopeTable$EPICID, ":")
df <- data.frame(matrix(unlist(str), nrow=length(str), byrow=T))
tmp_LR <- cbind(slopeTable, df)
ggplot(tmp_LR[tmp_LR$MAC.VisitSeq < 50,], aes(x=MAC.VisitSeq, color=Eye)) + 
  geom_histogram(aes(fill=Eye), alpha=0.3, position="identity", bins = 100) +
  ggtitle("adjusted slope for OCT (with ON patients)")



tmp <- merge(OCT_slope, slope_original, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"))
tmp <- merge(tmp_R, slope_original, by.x=c("EPICID","Eye"), by.y=c("EPICID","Eye"))
cor.test(tmp$GCL.VisitSeq, tmp$GCL_Slope)

