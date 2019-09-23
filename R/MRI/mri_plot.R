library(SNPassoc)
library(ggplot2)

# NBV_Slop, NGMV_Slop, NWMV_Slop, NT2Vlog_Slop, NVCSFlog_Slop, NCGM_Slop, PBVC_since_base_Slop
load(file="EPIC_MRI_slope+metadata+194genotypes.RData")
EPIC_MRI_lm <- EPIC_MRI
load(file="EPIC_MRI_slopeLME+metadata+194genotypes.RData")
EPIC_MRI_lme <- EPIC_MRI


#
library(SNPassoc)
snp <- "rs2364482"
geno <- dominant(snp(EPIC_MRI_snp[,snp], sep=""))
geno <- EPIC_MRI_snp[,snp]  # additive
geno <- recessive(snp(EPIC_MRI_snp[,snp], sep=""))

tmp_plot <- data.frame(Genotype=geno, Volume=EPIC_MRI_snp$NCGM.Time)
ggplot(tmp_plot, aes(Genotype, Volume)) + geom_boxplot() + geom_jitter(size=1.5, width=0.2, color="Blue", alpha=0.5) + 
  xlab("Genotype, recessive model") + ylab("Slope for CGM") + ggtitle(paste0("CGM: ",snp)) +
  theme_bw() + theme(legend.position="none")
#

snp <- "rs6435203"
geno <- dominant(snp(EPIC_MRI_snp[,snp], sep=""))
geno <- EPIC_MRI_snp[,snp]  # additive
geno <- recessive(snp(EPIC_MRI_snp[,snp], sep=""))

tmp_plot <- data.frame(Genotype=geno, Volume=EPIC_MRI_snp$NVCSFlog.Time)
ggplot(tmp_plot, aes(Genotype, Volume)) + geom_boxplot() + geom_jitter(size=1.5, width=0.2, color="Blue", alpha=0.5) + 
  xlab("Genotype, recessive model") + ylab("Log-transformed Slope for VCSF") + ggtitle(paste0("VCSF: ",snp)) +
  theme_bw() + theme(legend.position="none")
#





##### plot for volume dataset #####
MRIvol_meta_snp <- merge(EPIC_MRI_Vol_withMeta, EPIC_194snp, by.x="EPID", by.y="EPICID")

ggplot(MRIvol_meta_snp, aes(VisitType, NBV, color=rs6999228)) + geom_point() + geom_line(aes(group=EPICID)) +
  xlab("") + ylab("Slope") + ggtitle(paste0("slope using linear mixed effect model: ",snp)) +
  theme_bw()  ## + theme(legend.position="none")
#

geno <- dominant(snp(EPIC_MRI[,"rs2413436"], sep=""))
geno <- EPIC_MRI[,"rs2413436"]  # additive
geno <- recessive(snp(EPIC_MRI[,"rs2413436"], sep=""))
tmp_plot <- data.frame(Genotype=geno, Volume=EPIC_MRI$NCGM_Slop)
ggplot(tmp_plot, aes(Genotype, Volume, fill=Genotype)) + geom_boxplot() +   #geom_jitter(size=0.5, width=0.15) + 
  xlab("") + ylab("Slope") + ggtitle(paste0("")) +
  theme_bw() + theme(legend.position="none")
#

#
snp_list <- c("chr11:118762073","rs1051738","rs10877012","rs108990","rs10914539","rs10951042","rs11673987","rs11740512","rs12330493","rs13193887","rs137956","rs17030410","rs17051321","rs17066096","rs2104286","rs2413436","rs2445610","rs32658","rs35858108","rs4728142","rs4743150","rs61708525","rs6427540","rs6952809","rs706015","rs71624119","rs7550552","rs802725","rs879036","rs9282641","rs9913257")

for (i in 1:length(snp_list)){
  png(paste0("NBV_slope-",snp_list[i],".png"))
  tmp_plot <- data.frame(Genotype=EPIC_MRI[,snp_list[i]], Volume=EPIC_MRI$NBV_Slop)
  ggp <- ggplot(tmp_plot, aes(Genotype, Volume, fill=Genotype)) + geom_boxplot() +   #geom_jitter(size=0.5, width=0.15) + 
    xlab(snp_list[i]) + ylab("Slope of normalized brain volume (baseline ~ 2Y ~ 5Y") + ggtitle(paste0("NBV_Slop - ",snp_list[i])) +
    theme_bw() + theme(legend.position="none")
  print(ggp)
  dev.off()
  rm(tmp_plot)
}



#####  plot for baseline data  #####
library(ggplot2)
lmData <- EPIC_MRI_base_meta_snp
snp <- "rs35703946"   # rs9846396 (additive WM), rs735120 (recessive VCSF), rs35703946 (recessive VCSF), 

geno <- dominant(snp(lmData[,snp], sep=""))   # dominant
geno <- lmData[,snp]                          # additive
geno <- recessive(snp(lmData[,snp], sep=""))  # recessive

snp <- "rs35703946"
geno <- recessive(snp(lmData[,snp], sep=""))  # recessive
tmp_plot <- data.frame(Genotype=geno, Volume=lmData$NBV)
ggplot(tmp_plot, aes(Genotype, Volume)) + geom_boxplot() + geom_jitter(size=1.2, width=0.2, color="Blue", alpha=0.5) + 
  xlab("rs35703946 (recessive model)") + ylab("log VCSF") + ggtitle(paste0(snp," (FDR=0.097)")) +
  theme_bw() + theme(legend.position="none")

snp <- "rs735120"
geno <- recessive(snp(lmData[,snp], sep=""))  # recessive
tmp_plot <- data.frame(Genotype=geno, Volume=lmData$NVCSFlog)
ggplot(tmp_plot, aes(Genotype, Volume)) + geom_boxplot() + geom_jitter(size=1.2, width=0.2, color="Blue", alpha=0.5) + 
  xlab("rs735120 ") + ylab("log VCSF") + ggtitle(paste0(snp,"")) +
  theme_bw() + theme(legend.position="none")

snp <- "rs9846396"
geno <- lmData[,snp]
tmp_plot <- data.frame(Genotype=geno, Volume=lmData$NWMV)
ggplot(tmp_plot, aes(Genotype, Volume)) + geom_boxplot() + geom_jitter(size=1.2, width=0.2, color="Blue", alpha=0.5) + 
  xlab("rs9846396 (additive model)") + ylab("WM volume") + ggtitle(paste0(snp," (FDR = 0.046)")) +
  theme_bw() + theme(legend.position="none")


