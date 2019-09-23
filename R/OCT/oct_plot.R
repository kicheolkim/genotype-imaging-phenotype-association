##### editing from code for MRI
library(SNPassoc)
library(ggplot2)

# MAC.VisitSeq, GCC.VisitSeq, GCL.VisitSeq, INL.VisitSeq, IPL.VisitSeq, pRNFL.VisitSeqp
EPIC_OCT_slope_snp

# rs1250568, chr11:118762073 (fdr < 0.05) - for MAC, Recessive
# rs35858108 (fdr < 0.1) - pRNFL, Recessive
# rs1399180 for MAC, addi, Codo1 // rs6435203 for GCL, Codo2 (p < 0.001)


##### plot for slope #####
snp <- "rs974666234"
geno <- dominant(snp(EPIC_OCT_slope_snp[,snp], sep=""))
geno <- EPIC_OCT_slope_snp[,snp]  # additive
geno <- recessive(snp(EPIC_OCT_slope_snp[,snp], sep=""))

tmp_plot <- data.frame(EPICID=EPIC_OCT_slope_snp$EPICID, Genotype=geno, 
                       Volume=EPIC_OCT_slope_snp$MAC.VisitSeq, Eye=EPIC_OCT_slope_snp$Eye)
tmp_plot <- tmp_plot[tmp_plot$Volume < 50,]
ggplot(tmp_plot, aes(x=Genotype, y=Volume, color=Genotype)) + geom_boxplot() + 
  geom_jitter(size=1.2, width=0.15, color="Blue", alpha=0.5) + 
  xlab("Genotype, recessive model") + ylab("Slope for MAC") + ggtitle(paste0("MAC: ",snp)) +
  theme_bw()   #+ theme(legend.position="none")
ggplot(tmp_plot, aes(x=Genotype, y=Volume, color=Eye, shape=Eye)) + geom_boxplot() + 
  geom_jitter(size=1.5, width=0.2, color="Blue", alpha=0.4) + 
  xlab("Genotype, recessive model") + ylab("Slope for MAC") + ggtitle(paste0("MAC: ",snp)) +
  theme_bw()   #+ theme(legend.position="none")
#

snp <- "rs35858108"
geno <- recessive(snp(EPIC_OCT_slope_snp[,snp], sep=""))

tmp_plot <- data.frame(EPICID=EPIC_OCT_slope_snp$EPICID, Genotype=geno, 
                       Volume=EPIC_OCT_slope_snp$pRNFL.VisitSeq, Eye=EPIC_OCT_slope_snp$Eye)
ggplot(tmp_plot, aes(x=Genotype, y=Volume, color=Genotype)) + geom_boxplot() + 
  geom_jitter(size=1.2, width=0.15, color="Blue", alpha=0.5) + 
  xlab("Genotype, recessive model") + ylab("Slope for pRNFL") + ggtitle(paste0("pRNFL: ",snp)) +
  theme_bw()   #+ theme(legend.position="none")
ggplot(tmp_plot, aes(x=Genotype, y=Volume, color=Eye, shape=Eye)) + geom_boxplot() + 
  geom_jitter(size=1.5, width=0.2, color="Blue", alpha=0.5) + 
  xlab("Genotype, recessive model") + ylab("Slope for pRNFL") + ggtitle(paste0("pRNFL: ",snp)) +
  theme_bw()   #+ theme(legend.position="none")
#

snp <- "rs6435203"
geno <- EPIC_OCT_slope_snp[,snp]  # additive
geno <- recessive(snp(EPIC_OCT_slope_snp[,snp], sep=""))

tmp_plot <- data.frame(EPICID=EPIC_OCT_slope_snp$EPICID, Genotype=geno, 
                       Volume=EPIC_OCT_slope_snp$GCL.VisitSeq, Eye=EPIC_OCT_slope_snp$Eye)
tmp_plot <- tmp_plot[tmp_plot$Volume > -0.2,]
ggplot(tmp_plot, aes(x=Genotype, y=Volume, color=Genotype)) + geom_boxplot() + 
  geom_jitter(size=1.2, width=0.15, color="Blue", alpha=0.5) + 
  xlab("Genotype, additive model") + ylab("Slope for GCL") + ggtitle(paste0("GCL: ",snp)) +
  theme_bw()   #+ theme(legend.position="none")
ggplot(tmp_plot, aes(x=Genotype, y=Volume, color=Eye, shape=Eye)) + geom_boxplot() + 
  geom_jitter(size=1.5, width=0.2, color="Blue", alpha=0.5) + 
  xlab("Genotype, additive model") + ylab("Slope for GCL") + ggtitle(paste0("GCL: ",snp)) +
  theme_bw()   #+ theme(legend.position="none")
#







##### plot for volume dataset #####
OCT_matrices <- OCT_matrices_longitud
OCT_matrices <- merge(OCT_matrices_longitud, EPIC_194snp, by.x="EPIC_ID", by.y="EPICID")
#

snp <- "chr11:118762073"
geno <- dominant(snp(OCT_matrices[,snp], sep=""))
geno <- OCT_matrices[,snp]  # additive
geno <- recessive(snp(OCT_matrices[,snp], sep=""))

tmp_plot <- data.frame(EPICID=OCT_matrices$EPIC_ID, Genotype=geno, Visit=factor(OCT_matrices$VisitSeq),
                       Volume=OCT_matrices$MAC, Eye=OCT_matrices$Eye)
tmp_plot <- tmp_plot[tmp_plot$Volume < 300,]

ggplot(tmp_plot, aes(x=Genotype, y=Volume, color=Visit)) + geom_boxplot() + 
  #geom_jitter(size=1.2, width=0.15, color="Blue", alpha=0.5) + 
  xlab("Genotype") + ylab("MAC thickness") + ggtitle(paste0("MAC: ",snp)) +
  theme_bw()
ggplot(tmp_plot, aes(x=Genotype, y=Volume, color=Eye, shape=Eye)) + geom_boxplot() + 
  geom_jitter(size=1.2, width=0.2, color="Blue", alpha=0.5) + 
  xlab("Genotype") + ylab("Slope for MAC") + ggtitle(paste0("MAC: ",snp)) +
  theme_bw()   #+ theme(legend.position="none")
#

