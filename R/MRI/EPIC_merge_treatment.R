epic1 = read.csv("EPIC_MRI_Normalized_withMetadata_removed-NA&baseline.csv",head=T,as.is=T)
epic2 = read.csv("EPIC_MRI_treatments.csv",head=T,as.is=T)

# epic1 samples that 
if(!require("lubridate")){
  install.packages("lubridate")
}

naive = epic2[epic2$Treatment %in% "Treatment Naive",]
for(i in 1:nrow(epic1)){
  #match by EPICID
  matchsamples = epic2[epic2$EPICID %in% epic1$EPICID[i],]
  #match by date interval
  matchsamples = matchsamples[ymd(epic1$EPIC_scan_date[i]) %within% interval(ymd(matchsamples$Started), ymd(matchsamples$Ended)),]
  treatments = matchsamples$Treatment
  trtStart = matchsamples$Started
  trtEnd = matchsamples$Ended
  treatments = ifelse(length(treatments) > 1, paste(treatments, collapse= ","), treatments)
  trtStart = ifelse(length(trtStart) > 1, paste(trtStart, collapse= ","), trtStart)
  trtEnd = ifelse(length(trtEnd) > 1, paste(trtEnd, collapse= ","), trtEnd)
  epic1$Treatment[i] = treatments
  epic1$trtStarted[i] = trtStart
  epic1$trtEnd[i] = trtEnd
  if(epic1$EPICID[i] %in% naive$EPICID){
    epic1$Treatment[i] = "Treatment Naive"
  }
}
write.csv(epic1, "Merged_EPIC.csv",row.names=F)







