library(dplyr)

#############################          panIVT ENRICHED FILES GENERATION                    ############################################################
if(T){
#SPREADSHEET WITH THE KD DATA FOR PUS7 AND TRUB1 enzymes
SH <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PILEUP/init_gene_pileup/finalized_pileup/Merged_with_P_vals.csv", header = T)
SH_org <- SH

#panIVT dataset
pan_TRUB1KD<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/panIVT/TRUB1kd-panIVT.Merged_with_P_vals-COUNT.csv", header = T)
pan_TRUB1KD$ACP<-paste0(pan_TRUB1KD$Annotation,pan_TRUB1KD$chr,pan_TRUB1KD$position)
pan_PUS7KD<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/panIVT/PUS7kd-panIVT.Merged_with_P_vals-COUNT.csv", header = T)
pan_PUS7KD$ACP<-paste0(pan_PUS7KD$Annotation,pan_PUS7KD$chr,pan_PUS7KD$position)
panIVT_SCR<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/panIVT/SCRAMBLE-panIVT.Merged_with_P_vals-COUNT.csv", header = T)
panIVT_SCR$ACP<-paste0(panIVT_SCR$Annotation,panIVT_SCR$chr,panIVT_SCR$position)

# Merging dataframes and retaining all non-missing values
panIVT_TRUB1KD<- bind_rows(pan_TRUB1KD,panIVT_SCR) 
panIVT_TRUB1KD<- data.frame("ACP"=panIVT_TRUB1KD$ACP,"Annotation"=panIVT_TRUB1KD$Annotation,"chr"=panIVT_TRUB1KD$chr,"position"=panIVT_TRUB1KD$position,"N_reads_IVT"=panIVT_TRUB1KD$N_reads_IVT,"mm.IVT"=panIVT_TRUB1KD$mm.IVT,"T_IVT"=panIVT_TRUB1KD$T_IVT,"C_IVT"=panIVT_TRUB1KD$C_IVT)
panIVT_TRUB1KD<- distinct(panIVT_TRUB1KD)
panIVT_TRUB1KD<- panIVT_TRUB1KD[which(panIVT_TRUB1KD$N_reads_IVT>10),]
panIVT_TRUB1KD<- panIVT_TRUB1KD[which(panIVT_TRUB1KD$mm.IVT <=10),]

panIVT_PUS7KD<- bind_rows(pan_PUS7KD,panIVT_SCR)
panIVT_PUS7KD<- data.frame("ACP"=panIVT_PUS7KD$ACP,"Annotation"=panIVT_PUS7KD$Annotation,"chr"=panIVT_PUS7KD$chr,"position"=panIVT_PUS7KD$position,"N_reads_IVT"=panIVT_PUS7KD$N_reads_IVT,"mm.IVT"=panIVT_PUS7KD$mm.IVT,"T_IVT"=panIVT_PUS7KD$T_IVT,"C_IVT"=panIVT_PUS7KD$C_IVT)
panIVT_PUS7KD<- distinct(panIVT_PUS7KD)
panIVT_PUS7KD<- panIVT_PUS7KD[which(panIVT_PUS7KD$N_reads_IVT>10),]
panIVT_PUS7KD<- panIVT_PUS7KD[which(panIVT_PUS7KD$mm.IVT <=10),]

#data preparation, pan IVT enrichement and filtration (PSI-p-values, total reads IVT, total reads DRS)
if (T) {
################# For TRUB1_KD, SCR data
if (T) {
    ### Filter based on the p-value
    TRUB1_filter <- SH[which(SH$p.value.SCR < 0.001),]
    
    
    ### Add total_TRUB1_KD & total_SCR
    TRUB1_filter$total_TRUB1_KD <- 0
    TRUB1_filter$total_SCR <- 0
    TRUB1_filter$total_IVT <- 0
    for (i in c(1:nrow(TRUB1_filter))) {
      TRUB1_filter$total_TRUB1_KD[i] <- TRUB1_filter$T_TRUB1_KD[i] + TRUB1_filter$C_TRUB1_KD[i]
      TRUB1_filter$total_SCR[i] <- TRUB1_filter$T_SCR[i] + TRUB1_filter$C_SCR[i]
      TRUB1_filter$total_IVT[i] <- TRUB1_filter$T_IVT[i] + TRUB1_filter$C_IVT[i]
      
    }
    
    
    ### Add the mismatch difference between TRUB1_KD and IVT and SCR and IVT
    TRUB1_filter$mm.TRUB1_KDminusIVT <- 0
    TRUB1_filter$mm.SCRminusIVT <- 0
    for (i in c(1:nrow(TRUB1_filter))) {
      TRUB1_filter$mm.TRUB1_KDminusIVT[i] <- TRUB1_filter$mm.TRUB1_KD[i]-TRUB1_filter$mm.IVT[i]
      TRUB1_filter$mm.SCRminusIVT[i] <- TRUB1_filter$mm.SCR[i]-TRUB1_filter$mm.IVT[i]
    }
    
    ### Add the mismatch difference between TRUB1_KD and SCR
    TRUB1_filter$mm.TRUB1_KDminusSCR <- 0
    for (i in c(1:nrow(TRUB1_filter))) {
      TRUB1_filter$mm.TRUB1_KDminusSCR[i] <- TRUB1_filter$mm.TRUB1_KD[i]-TRUB1_filter$mm.SCR[i]
    }
    
    ### Filter to total_number_reads_TRUB1_KD & total_number_reads_SCR < 10
    TRUB1_filter <- TRUB1_filter[which(TRUB1_filter$total_TRUB1_KD > 10 & TRUB1_filter$total_SCR > 10),]
    
    
    ### Substitute mmIVT==0 with panIVT
    TRUB1_filter$ACP<-paste0(TRUB1_filter$Annotation,TRUB1_filter$chr,TRUB1_filter$position)
    ACP<-TRUB1_filter$ACP[which(TRUB1_filter$total_IVT < 10)]
    ACP_pan<-panIVT_TRUB1KD$ACP[which(panIVT_TRUB1KD$ACP %in% ACP)]
    ACP_pan_mm<-panIVT_TRUB1KD$mm.IVT[which(panIVT_TRUB1KD$ACP %in% ACP)]
    pan<-panIVT_TRUB1KD[which(panIVT_TRUB1KD$ACP %in% ACP),]
    for (i in c(1:length(ACP_pan))) {
      if (ACP_pan[i] %in% TRUB1_filter$ACP){
        idx<-which(TRUB1_filter$ACP==ACP_pan[i])
        TRUB1_filter$total_IVT[idx]<-panIVT_TRUB1KD$N_reads_IVT[which(panIVT_TRUB1KD$ACP %in% ACP_pan[i])]
        TRUB1_filter$T_IVT[idx]<-panIVT_TRUB1KD$T_IVT[which(panIVT_TRUB1KD$ACP %in% ACP_pan[i])]
        TRUB1_filter$C_IVT[idx]<-panIVT_TRUB1KD$C_IVT[which(panIVT_TRUB1KD$ACP %in% ACP_pan[i])]
        TRUB1_filter$mm.IVT[idx]<-panIVT_TRUB1KD$mm.IVT[which(panIVT_TRUB1KD$ACP %in% ACP_pan[i])]
      }
    }

    TRUB1_filter <- TRUB1_filter[which(TRUB1_filter$total_IVT > 10),]
    TRUB1_filter<-TRUB1_filter[which((TRUB1_filter$T_IVT+TRUB1_filter$C_IVT)!=0),]
   
    
    ### Add the color column based on the kmer
    TRUB1_filter$color <- "black"
    for (i in c(1:nrow(TRUB1_filter))) {
      if (TRUB1_filter$kmer[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
        TRUB1_filter$color[i] <- "blue"
      }
      if (TRUB1_filter$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) { #GUUCN
        TRUB1_filter$color[i] <- "red"
      }
      else {}
    }

  }
write.csv(TRUB1_filter,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TRUB1_filter_pan.csv",row.names = F)

################# For PUS7_KD, SCR data
if (T) {

  PUS7_filter <- SH[which(SH$p.value.SCR < 0.001),]
  # PUS7_filter<-SH[which(SH$mm.IVT<=10),]
  
  ### Add total_PUS7_KD & total_PUS7_KD
  PUS7_filter$total_PUS7_KD <- 0
  PUS7_filter$total_SCR <- 0
  for (i in c(1:nrow(PUS7_filter))) {
    PUS7_filter$total_PUS7_KD[i] <- PUS7_filter$T_PUS7_KD[i] + PUS7_filter$C_PUS7_KD[i]
    PUS7_filter$total_SCR[i] <- PUS7_filter$T_SCR[i] + PUS7_filter$C_SCR[i]
    PUS7_filter$total_IVT[i] <- PUS7_filter$T_IVT[i] + PUS7_filter$C_IVT[i]
  }
  
  
  ### Add the mismatch difference between TRUB1_KD and IVT and PUS7_KD and IVT
  PUS7_filter$mm.PUS7_KDminusIVT <- 0
  PUS7_filter$mm.SCRminusIVT <- 0
  for (i in c(1:nrow(PUS7_filter))) {
    PUS7_filter$mm.PUS7_KDminusIVT[i] <- PUS7_filter$mm.PUS7_KD[i]-PUS7_filter$mm.IVT[i]
    PUS7_filter$mm.SCRminusIVT[i] <- PUS7_filter$mm.SCR[i]-PUS7_filter$mm.IVT[i]
  }
  
  
  ### Add the mismatch difference between PUS7_KD and SCR
  PUS7_filter$mm.PUS7_KDminusSCR <- 0
  for (i in c(1:nrow(PUS7_filter))) {
    PUS7_filter$mm.PUS7_KDminusSCR[i] <- PUS7_filter$mm.PUS7_KD[i]-PUS7_filter$mm.SCR[i]
    
  }
  
  ### Filter to total_PUS7_KD & total_SCR < 10
  PUS7_filter <- PUS7_filter[which(PUS7_filter$total_PUS7_KD > 10 & PUS7_filter$total_SCR > 10),]
  
  ### Substitute mmIVT==0 with panIVT
  PUS7_filter$ACP<-paste0( PUS7_filter$Annotation, PUS7_filter$chr, PUS7_filter$position)
  ACP<- PUS7_filter$ACP[which( PUS7_filter$total_IVT < 10)]
  ACP_pan<-panIVT_PUS7KD$ACP[which(panIVT_PUS7KD$ACP %in% ACP)]
  ACP_pan_mm<-panIVT_PUS7KD$mm.IVT[which(panIVT_PUS7KD$ACP %in% ACP)]
  
  for (i in c(1:length(ACP_pan))) {
    if (ACP_pan[i] %in% PUS7_filter$ACP){
      idx<-which(PUS7_filter$ACP==ACP_pan[i])
      PUS7_filter$total_IVT[idx]<-panIVT_PUS7KD$N_reads_IVT[which(panIVT_PUS7KD$ACP %in% ACP_pan[i])]
      PUS7_filter$T_IVT[idx]<-panIVT_PUS7KD$T_IVT[which(panIVT_PUS7KD$ACP %in% ACP_pan[i])]
      PUS7_filter$C_IVT[idx]<-panIVT_PUS7KD$C_IVT[which(panIVT_PUS7KD$ACP %in% ACP_pan[i])]
      PUS7_filter$mm.IVT[idx]<-panIVT_PUS7KD$mm.IVT[which(panIVT_PUS7KD$ACP %in% ACP_pan[i])]
    }
  }
  

  # PUS7_filter$total_IVT[which( PUS7_filter$ACP %in% ACP_pan)]<-panIVT_PUS7KD$N_reads_IVT[which(panIVT_PUS7KD$ACP %in% ACP_pan)]
  # PUS7_filter$mm.IVT[which( PUS7_filter$ACP %in% ACP_pan)]<-panIVT_PUS7KD$mm.IVT[which(panIVT_PUS7KD$ACP %in% ACP_pan)]
  ## Remove the ones with less than 10 reads in IVT (in panIVT as we just substituted the pairedIVT with the panIVT)
  PUS7_filter <-  PUS7_filter[which( PUS7_filter$total_IVT > 10),]
  PUS7_filter<-PUS7_filter[which((PUS7_filter$T_IVT+PUS7_filter$C_IVT)!=0),]
  
  ### Add the color column based on the kmer
  PUS7_filter$color <- "black"
  for (i in c(1:nrow(PUS7_filter))) {
    if (PUS7_filter$kmer[i]  %in%  c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) { #UNUAR
      PUS7_filter$color[i] <- "blue"
    }
    if (PUS7_filter$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
      PUS7_filter$color[i] <- "red"
    }
    else {}
  }
  
  
}
write.csv(PUS7_filter,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PUS7_filter_pan.csv",row.names = F)

}

}


#############################    START HERE AFTER THE panIVT ENRICHED FILES HAVE BEEN GENERATED  ######################################################
#######################################################################################################################################################
library(tidyverse)
library(Cairo)
TRUB1_filter_orig<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TRUB1_filter_pan.csv", header = T) #with pan IVT enrichment
PUS7_filter_orig<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PUS7_filter_pan.csv", header = T) #with pan IVT enrichment

#IVT filtration to exclude sites with %IVT>10%. In this way we keep just sites with low IVT errors.
mmIVT=10
TRUB1_filter<-TRUB1_filter_orig[which(TRUB1_filter_orig$mm.IVT<=mmIVT),]
PUS7_filter<-PUS7_filter_orig[which(PUS7_filter_orig$mm.IVT<=mmIVT),]
#total number of C and T filtration
names_filtered<-paste0(TRUB1_filter$Annotation,TRUB1_filter$position)

#Supplementary1
if(F){
TRUB1_S1<-data.frame("Annotation"=TRUB1_filter$Annotation,"chr"=TRUB1_filter$chr,"position"=TRUB1_filter$position,"strand"=TRUB1_filter$strand,"kmer"=TRUB1_filter$kmer,
                     "T_TRUB1_KD"=TRUB1_filter$T_TRUB1_KD,"C_TRUB1_KD"=TRUB1_filter$C_TRUB1_KD,
                     "T_SCR"=TRUB1_filter$T_SCR,"C_SCR"=TRUB1_filter$C_SCR,
                     "T_IVT"=TRUB1_filter$T_IVT,"C_IVT"=TRUB1_filter$C_IVT,
                     "pvalue.SCR"=TRUB1_filter$p.value.SCR,
                     "mm.TRUB1KD"=TRUB1_filter$mm.TRUB1_KD,
                     "mm.SCR"=TRUB1_filter$mm.SCR,
                     "mm.IVT"=TRUB1_filter$mm.IVT)
write.csv(TRUB1_S1,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/TRUB1KD_s1.csv",row.names = F)
PUS7_S1<-data.frame("Annotation"=PUS7_filter$Annotation,"chr"=PUS7_filter$chr,"position"=PUS7_filter$position,"strand"=PUS7_filter$strand,"kmer"=PUS7_filter$kmer,
                    "T_PUS7_KD"=PUS7_filter$T_PUS7_KD,"C_PUS7_KD"=PUS7_filter$C_PUS7_KD,
                    "T_SCR"=PUS7_filter$T_SCR,"C_SCR"=PUS7_filter$C_SCR,
                    "T_IVT"=PUS7_filter$T_IVT,"C_IVT"=PUS7_filter$C_IVT,
                    "pvalue.SCR"=PUS7_filter$p.value.SCR,
                    "mm.PUS7KD"=PUS7_filter$mm.PUS7_KD,
                    "mm.SCR"=PUS7_filter$mm.SCR,
                    "mm.IVT"=PUS7_filter$mm.IVT)
write.csv(PUS7_S1,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/PUS7KD_s1.csv",row.names = F)
}

#TPM
TPM_Diff1 = read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/Diff1.tsv", sep = '\t', header = TRUE)
TPM_Diff2 = read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/Diff2.tsv", sep = '\t', header = TRUE)
TPM_Diff3 = read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/Diff3.tsv", sep = '\t', header = TRUE)
TPM_Undiff1 = read.table(file ="/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/Undiff1.tsv", sep = '\t', header = TRUE)
TPM_Undiff2 = read.table(file ="/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/Undiff2.tsv", sep = '\t', header = TRUE)
TPM_Undiff3 = read.table(file ="/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/Undiff3.tsv", sep = '\t', header = TRUE)
TPM_Undiff = read.table(file ="/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/Undiff_merged.tsv", sep = '\t', header = TRUE)
TPM_PB1 = read.table(file ="/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/PB1.tsv", sep = '\t', header = TRUE)
TPM_PB2 = read.table(file ="/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/PB2.tsv", sep = '\t', header = TRUE)
TPM_PB = read.table(file ="/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/TPM/PB_merged.tsv", sep = '\t', header = TRUE)

# table generation for SNR code. These tables will be the input for the SNR_marginal_likelihood.py script. 
# TRUB1_filter_Amr<-data.frame(TRUB1_filter$Annotation,TRUB1_filter$chr,TRUB1_filter$position,TRUB1_filter$T_TRUB1_KD,TRUB1_filter$C_TRUB1_KD,TRUB1_filter$T_IVT,TRUB1_filter$C_IVT)
# colnames(TRUB1_filter_Amr)<-c("Annotation","chr","position","DRS_U","DRS_C","IVT_U","IVT_C")
# write.csv(TRUB1_filter_Amr,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/TRUB1_KD_Stu_panIVT.csv",row.names = F)
# TRUB1_filter_Amr<-data.frame(TRUB1_filter$Annotation,TRUB1_filter$chr,TRUB1_filter$position,TRUB1_filter$T_SCR,TRUB1_filter$C_SCR,TRUB1_filter$T_IVT,TRUB1_filter$C_IVT)
# colnames(TRUB1_filter_Amr)<-c("Annotation","chr","position","DRS_U","DRS_C","IVT_U","IVT_C")
# write.csv(TRUB1_filter_Amr,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/TRUB1_SCR_Stu_panIVT.csv",row.names = F)
# 
# PUS7_filter_Amr<-data.frame(PUS7_filter$Annotation,PUS7_filter$chr,PUS7_filter$position,PUS7_filter$T_PUS7_KD,PUS7_filter$C_PUS7_KD,PUS7_filter$T_IVT,PUS7_filter$C_IVT)
# colnames(PUS7_filter_Amr)<-c("Annotation","chr","position","DRS_U","DRS_C","IVT_U","IVT_C")
# write.csv(PUS7_filter_Amr,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/PUS7_KD_Stu_panIVT.csv",row.names = F)
# PUS7_filter_Amr<-data.frame(PUS7_filter$Annotation,PUS7_filter$chr,PUS7_filter$position,PUS7_filter$T_SCR,PUS7_filter$C_SCR,PUS7_filter$T_IVT,PUS7_filter$C_IVT)
# colnames(PUS7_filter_Amr)<-c("Annotation","chr","position","DRS_U","DRS_C","IVT_U","IVT_C")
# write.csv(PUS7_filter_Amr,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/PUS7_SCR_Stu_panIVT.csv",row.names = F)

######################################################### Plots SCRminusKnockdowns vs N (Meni version)
sig=1
n_reads_high=30
n_reads_low=10
difference=0.2 #20%



######################################################### plot SCRminusTRUB1KD vs N and sigma>2
if (T){

  #Add snr information and filter SCR + TRUB1KD validated sites
  if(T){

  TRUB1_filter$mmSCRminusTRUB1KD = -TRUB1_filter$mm.TRUB1_KDminusSCR
  #SNR file load
  snrKD<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/snrTRUB1_KD_Stu_panIVT.csv", header=T)
  names(snrKD)[names(snrKD) == "SNR"] <- "SNR_TRUB1KD"
  snrSCR<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/snrTRUB1_SCR_Stu_panIVT.csv", header=T)
  names(snrSCR)[names(snrSCR) == "SNR"] <- "SNR_SCR"
  #SNR append
  selected_columns<-c("Annotation", "position","SNR_TRUB1KD")
  TRUB1_filter <- merge(TRUB1_filter, snrKD[selected_columns], by = c("Annotation", "position"), all.x = TRUE)
  selected_columns<-c("Annotation", "position","SNR_SCR")
  TRUB1_filter <- merge(TRUB1_filter, snrSCR[selected_columns], by = c("Annotation", "position"), all.x = TRUE)
  #N_value calculation
  TRUB1_filter$N_value<-pmax(TRUB1_filter$SNR_TRUB1KD,TRUB1_filter$SNR_SCR,na.rm = TRUE)
  
  #Supplementary Table2
  if(F){
    TRUB1_S1<-data.frame("Annotation"=TRUB1_filter$Annotation,"chr"=TRUB1_filter$chr,"position"=TRUB1_filter$position,"strand"=TRUB1_filter$strand,"kmer"=TRUB1_filter$kmer,
                         "T_TRUB1_KD"=TRUB1_filter$T_TRUB1_KD,"C_TRUB1_KD"=TRUB1_filter$C_TRUB1_KD,
                         "T_SCR"=TRUB1_filter$T_SCR,"C_SCR"=TRUB1_filter$C_SCR,
                         "T_IVT"=TRUB1_filter$T_IVT,"C_IVT"=TRUB1_filter$C_IVT,
                         "pvalue.SCR"=TRUB1_filter$p.value.SCR,
                         "mm.TRUB1KD"=TRUB1_filter$mm.TRUB1_KD,
                         "mm.SCR"=TRUB1_filter$mm.SCR,
                         "mm.IVT"=TRUB1_filter$mm.IVT,
                         "mm.SCRminusTRUB1KD"=-TRUB1_filter$mm.TRUB1_KDminusSCR,  
                         "SNR"=TRUB1_filter$N_value)
    write.csv(TRUB1_S1,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/TRUB1KD_s2.csv",row.names = F)

  }
  #Supplementary Table8 
  if(F){
    TRUB1_S1<-data.frame("<0"=length(which(TRUB1_filter$N_value<0)),
                         "0-0.5"=length(which(TRUB1_filter$N_value>=0 & TRUB1_filter$N_value<0.5)),
                         "0.5-1"=length(which(TRUB1_filter$N_value>=0.5 & TRUB1_filter$N_value<1)),
                         "1-2"=length(which(TRUB1_filter$N_value>=1 & TRUB1_filter$N_value<2)),
                         ">2"=length(which(TRUB1_filter$N_value>=2)))
    write.csv(TRUB1_S1,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/TRUB1KD_8",row.names = F)
    
  }
  
  #applying SNR threshold
  TRUB1_filter_sigma2<-TRUB1_filter[which(TRUB1_filter$N_value>sig),] #SNR>1
  
  length(which(TRUB1_filter_sigma2$N_value>sig))
  length(which(TRUB1_filter_sigma2$color=="red"))
  
  #filtering for number of reads in Scrambled
  TRUB1_filter_sigma2<-TRUB1_filter_sigma2[which( (TRUB1_filter_sigma2$total_TRUB1_KD>n_reads_high & TRUB1_filter_sigma2$total_SCR>n_reads_high)
                                                  | ( (TRUB1_filter_sigma2$total_TRUB1_KD>n_reads_low & TRUB1_filter_sigma2$total_SCR>n_reads_low) & (abs(TRUB1_filter_sigma2$mm.TRUB1_KDminusSCR) >= difference*TRUB1_filter_sigma2$mm.SCR)) ),]
  

  #add TPM
  TRUB1KD_TPM<-read.table(file="/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TPM/TRUB1_kd.tsv",sep='\t',header=T)
  
  TRUB1_filter_sigma2$TPM_KD <- 0
  for (i in c(1:nrow(TRUB1_filter_sigma2))) {
    row.on.Direct.TPM <- which(TRUB1KD_TPM$Gene.Name == TRUB1_filter_sigma2$Annotation[i])
    if (length(row.on.Direct.TPM == 1)) {
      TRUB1_filter_sigma2$TPM_KD[i] <- TRUB1KD_TPM$TPM[row.on.Direct.TPM][[1]]
    }
    else{}
  }
  TRUB1_filter_sigma2$TRUB1KD_TPM_log <- log10(TRUB1_filter_sigma2$TPM_KD)
  
  cutoff=15
  #sites in which we believe
  TRUB1KD_sites_all<-TRUB1_filter_sigma2[which(TRUB1_filter_sigma2$mmSCRminusTRUB1KD>=30 | 
                                              (TRUB1_filter_sigma2$mmSCRminusTRUB1KD>=cutoff & TRUB1_filter_sigma2$total_TRUB1_KD>=(30.5-TRUB1_filter_sigma2$mmSCRminusTRUB1KD)/0.5)),] #n_reads_high
  TRUB1KD_sites_grey<-TRUB1KD_sites_all[which(TRUB1KD_sites_all$color=="black"),]
  TRUB1KD_sites<-TRUB1KD_sites_all[which(TRUB1KD_sites_all$color=="red"),]
  TRUB1KD_sites$CP<-paste0(TRUB1KD_sites$chr,TRUB1KD_sites$position)
  for (i in c(1:nrow(TRUB1KD_sites))) {
    r <- which(TRUB1KD_sites$CP==TRUB1KD_sites$CP[i])
    if (length(r)>1) {
      
      TRUB1KD_sites$Annotation[r[1]]<-paste(TRUB1KD_sites$Annotation[r],collapse='-')
      row<-TRUB1KD_sites[r[1],]
      TRUB1KD_sites <-TRUB1KD_sites[-r,]
      TRUB1KD_sites[nrow(TRUB1KD_sites) + 1,] <- row
      
    }
    
  }
  
  
  #Supplementary Table3
  if(F){
    TRUB1_S1<-data.frame("Annotation"=TRUB1KD_sites_all$Annotation,"chr"=TRUB1KD_sites_all$chr,"position"=TRUB1KD_sites_all$position,"strand"=TRUB1KD_sites_all$strand,"kmer"=TRUB1KD_sites_all$kmer,
                         "T_TRUB1_KD"=TRUB1KD_sites_all$T_TRUB1_KD,"C_TRUB1_KD"=TRUB1KD_sites_all$C_TRUB1_KD,
                         "T_SCR"=TRUB1KD_sites_all$T_SCR,"C_SCR"=TRUB1KD_sites_all$C_SCR,
                         "T_IVT"=TRUB1KD_sites_all$T_IVT,"C_IVT"=TRUB1KD_sites_all$C_IVT,
                         "pvalue.SCR"=TRUB1KD_sites_all$p.value.SCR,
                         "mm.TRUB1KD"=TRUB1KD_sites_all$mm.TRUB1_KD,
                         "mm.SCR"=TRUB1KD_sites_all$mm.SCR,
                         "mm.IVT"=TRUB1KD_sites_all$mm.IVT,
                         "mm.SCRminusTRUB1KD"=TRUB1KD_sites_all$mmSCRminusTRUB1KD,
                         "SNR"=TRUB1KD_sites_all$N_value)
    write.csv(TRUB1_S1,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/TRUB1KD_s3.csv",row.names = F)

  }
  
  
  }
 
  plot_figs=T
  #plot KD minus SCR 
  if(plot_figs==T){
  min.x = -80
  max.x = 80
  min.y = -4
  max.y = 3
  line.thickness = 0.2
  label.size = 1
  text.size = 0.8
  point.size = 0.1
  tck.length = 0.01
  tick.thickness = 1
  transparency = 1
 
  file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Trub1_little_1_text_panIVT_sig",sig,"_stu_.pdf")
  
  Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
  
  plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = '',xlab='',ann = F,axes=T,type="l", 
       xaxt='n', yaxt='n', 
       bty='n', bg='transparent')
  
  
  points(x=-TRUB1_filter$mm.TRUB1_KDminusSCR[which(TRUB1_filter$color == "black")], y=log10(TRUB1_filter$N_value[which(TRUB1_filter$color=="black")]), pch=22, col="grey", bg=adjustcolor("grey",alpha=0.05), cex = log2(TRUB1_filter$total_TRUB1_KD[which(TRUB1_filter$color == "black")])*point.size)
  points(x=-TRUB1_filter$mm.TRUB1_KDminusSCR[which(TRUB1_filter$color == "red")], y=log10(TRUB1_filter$N_value[which(TRUB1_filter$color=="red")]), pch=22, col="orange", bg=adjustcolor("orange",alpha=0.6), cex = log2(TRUB1_filter$total_TRUB1_KD[which(TRUB1_filter$color == "red")])*point.size)
  points(x=-TRUB1_filter$mm.TRUB1_KDminusSCR[which(TRUB1_filter$color == "blue")], y=log10(TRUB1_filter$N_value[which(TRUB1_filter$color=="blue")]), pch=22, col=adjustcolor("blue",alpha=0.8), bg=adjustcolor("blue",alpha=0.2), cex = log2(TRUB1_filter$total_TRUB1_KD[which(TRUB1_filter$color == "blue")])*point.size)
  
  segments(x0 = 0, y0 = min.y, x1 = 0, y1 = max.y, lty=2, col = "grey3")
  segments(x0 = min.x, y0 = log10(sig), x1 = max.x, y1 =log10(sig), lty=2, col = "grey3")
  
  TRUB1KD_sites<-TRUB1KD_sites[order(TRUB1KD_sites$N_value),]
  
  rows <- nrow(TRUB1KD_sites) 
  # extracting odd rows  
  odd_rows <- seq_len(rows) %% 2 
  # getting data from odd data frame 
  TRUB1KD_sites_odd <- TRUB1KD_sites[odd_rows == 1, ] 

  rows <- nrow(TRUB1KD_sites) 
  # extracting odd rows  
  even_rows <- seq_len(rows) %% 2 
  # getting data from odd data frame 
  TRUB1KD_sites_even <- TRUB1KD_sites[even_rows == 0, ] 
  
  text(x=-TRUB1KD_sites_odd$mm.TRUB1_KDminusSCR+2, y=log10(TRUB1KD_sites_odd$N_value)+0.06,TRUB1KD_sites_odd$Annotation, col='grey3', cex=0.2)
  text(x=-TRUB1KD_sites_even$mm.TRUB1_KDminusSCR-2, y=log10(TRUB1KD_sites_even$N_value),TRUB1KD_sites_even$Annotation, col='grey3', cex=0.2)
  #grey
  text(x=-TRUB1KD_sites_grey$mm.TRUB1_KDminusSCR[which(TRUB1KD_sites_grey$mmSCRminusTRUB1KD>30 & TRUB1KD_sites_grey$N_value>10)], y=log10(TRUB1KD_sites_grey$N_value[which(TRUB1KD_sites_grey$mmSCRminusTRUB1KD>30 & TRUB1KD_sites_grey$N_value>10)]),TRUB1KD_sites_grey$Annotation[which(TRUB1KD_sites_grey$mmSCRminusTRUB1KD>30 & TRUB1KD_sites_grey$N_value>10)], col='dimgrey', cex=0.2)
  
   
  #for plot KD-SCR vs SCR
   axis(side = 2, at = seq(min.y, max.y, by = 1), labels = seq(min.y, max.y, by = 1),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  axis(side = 1, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  legend("topleft", c("TRUB1 motif","PUS7 motif", "Other"),     
         cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
         pch=c(16,16,16),lty=c(0,0,0), 
         lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
  mtext(text = expression(paste(log[10], SNR)),side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
  # mtext(text = "SNR",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
  mtext(text = "Scramble-TRUB1KD siRNA (U-C MM%)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
  ##ggsave(filename =file.name,width = 6,height = 6.6)
  alpha=0
  dev.off()
  }
    
  #KD minus SCR vs log10 for sigma>6  
  if(plot_figs==T){
      
      min.x = 1
      max.x = round(log10(max(TRUB1_filter_sigma2$total_TRUB1_KD)))
      min.y = 0
      max.y = 100
      line.thickness = 0.2
      label.size = 1
      text.size = 0.8
      point.size = 0.15
      tck.length = 0.01
      tick.thickness = 1
      transparency = 1
      
      file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Trub1_sigma2_reads_filterSCR_notext_sig",sig,"_panIVT_stu.pdf")
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
      
      plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = '',xlab='',ann = F,axes=T,type="l", 
           xaxt='n', yaxt='n', 
           bty='n', bg='transparent')
      
      
      points(x=log10(TRUB1_filter_sigma2$total_TRUB1_KD[which(TRUB1_filter_sigma2$color == "black")]), y=TRUB1_filter_sigma2$mmSCRminusTRUB1KD[which(TRUB1_filter_sigma2$color=="black")], pch=22, col="grey", bg=adjustcolor("grey",alpha=0.05), cex = log2(TRUB1_filter_sigma2$N_value[which(TRUB1_filter_sigma2$color == "black")])*point.size)
      points(x=log10(TRUB1_filter_sigma2$total_TRUB1_KD[which(TRUB1_filter_sigma2$color == "blue")]), y=TRUB1_filter_sigma2$mmSCRminusTRUB1KD[which(TRUB1_filter_sigma2$color=="blue")], pch=22, col="blue", bg=adjustcolor("blue",alpha=0.3), cex = log2(TRUB1_filter_sigma2$N_value[which(TRUB1_filter_sigma2$color == "blue")])*point.size)
      points(x=log10(TRUB1_filter_sigma2$total_TRUB1_KD[which(TRUB1_filter_sigma2$color == "red")]), y=TRUB1_filter_sigma2$mmSCRminusTRUB1KD[which(TRUB1_filter_sigma2$color=="red")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.6), cex = log2(TRUB1_filter_sigma2$N_value[which(TRUB1_filter_sigma2$color == "red")])*point.size)

      
      TRUB1KD_sites<-TRUB1KD_sites[order(TRUB1KD_sites$total_SCR),]
      
      rows <- nrow(TRUB1KD_sites) 
      # extracting odd rows  
      odd_rows <- seq_len(rows) %% 2 
      # getting data from odd data frame 
      TRUB1KD_sites_odd <- TRUB1KD_sites[odd_rows == 1, ] 
      
      rows <- nrow(TRUB1KD_sites) 
      # extracting odd rows  
      even_rows <- seq_len(rows) %% 2 
      # getting data from odd data frame 
      TRUB1KD_sites_even <- TRUB1KD_sites[even_rows == 0, ] 
      
      text(x=log10(TRUB1KD_sites_odd$total_TRUB1_KD)+0.03, y=TRUB1KD_sites_odd$mmSCRminusTRUB1KD,TRUB1KD_sites_odd$Annotation, col='grey3', cex=0.3)
      text(x=log10(TRUB1KD_sites_even$total_TRUB1_KD)+0.04, y=TRUB1KD_sites_even$mmSCRminusTRUB1KD+0.6,TRUB1KD_sites_even$Annotation, col='grey3', cex=0.3)
      #grey
      text(x=log10(TRUB1KD_sites_grey$total_TRUB1_KD[which(TRUB1KD_sites_grey$mmSCRminusTRUB1KD>30 & TRUB1KD_sites_grey$N_value>10)]), y=TRUB1KD_sites_grey$mmSCRminusTRUB1KD[which(TRUB1KD_sites_grey$mmSCRminusTRUB1KD>30 & TRUB1KD_sites_grey$N_value>10)],TRUB1KD_sites_grey$Annotation[which(TRUB1KD_sites_grey$mmSCRminusTRUB1KD>30 & TRUB1KD_sites_grey$N_value>10)], col='dimgrey', cex=0.2)
      
      segments(x0=1,y0=30,x1=log10(30),y1=15,lty=2,col="black")
      segments(x0=log10(30),y0=15,x1=max.x,y1=15,lty=2,col="black")
      
      #for plot KD-SCR vs SCR
      axis(side = 2, at = seq(min.y, max.y, by = 10), labels = seq(min.y, max.y, by = 10),
           lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
           las = 1, lwd = line.thickness, line = 0)
      axis(side = 1, at = seq(min.x,max.x, by = 0.5), labels = seq(min.x, max.x, by = 0.5),
           lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
           las = 1, lwd = line.thickness, line = 0)
      legend("topright", c("TRUB1 motif","PUS7 motif", "Other"),     
             cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
             pch=c(16,16,16),lty=c(0,0,0), 
             lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
      mtext(text = "Scramble - TRUB1 KD (U-to-C MM%) ",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
      mtext(text = "log10(Reads)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
      dev.off()
    }  

  ############################### seqLogos on the grey and orange motifs
  if (plot_figs==T){

   
    # kmers<-invisible(sapply(TRUB1KD_sites$kmer,function(x) str_replace_all(x,"T","U")))
    # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TRUB1logo_orange.pdf"

 
    # kmers<-sapply( TRUB1KD_sites_grey$kmer,function(x) str_replace_all(x,"T","U"))
    # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TRUB1logo_grey.pdf"

    
    kmers<-sapply( TRUB1KD_sites_all$kmer,function(x) str_replace_all(x,"T","U"))
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TRUB1logo_merged.pdf"

    
    my_list<-list("new motifs" = kmers)
    min.x = 0
    max.x = 200
    min.y = 0
    max.y = 20
    line.thickness = 1
    label.size = 1
    text.size = 0.8
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)

    ggseqlogo::ggseqlogo(my_list,ncol=1,method="prob")
    #same as
    #ggplot() + geom_logo(seqs_dna) + theme_logo() + facet_wrap(~seq_group, ncol=4, scales='free_x') 
    
    dev.off()
    
  }
  
  
}

######################################################### plot SCRminusPUS7_KD vs N and SNR>2
if(T){
  
  #Add snr info and filter SCR + PUS7KD validated sites
  if (T){
    
    PUS7_filter$mmSCRminusPUS7KD = -PUS7_filter$mm.PUS7_KDminusSCR

    snrKD<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/snrPUS7_KD_Stu_panIVT.csv", header=T)
    names(snrKD)[names(snrKD) == "SNR"] <- "SNR_PUS7KD"
    snrSCR<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/snrPUS7_SCR_Stu_panIVT.csv", header=T)
    names(snrSCR)[names(snrSCR) == "SNR"] <- "SNR_SCR"
    
    selected_columns<-c("Annotation", "position","SNR_PUS7KD")
    PUS7_filter <- merge(PUS7_filter, snrKD[selected_columns], by = c("Annotation", "position"), all.x = TRUE)
    selected_columns<-c("Annotation", "position","SNR_SCR")
    PUS7_filter <- merge(PUS7_filter, snrSCR[selected_columns], by = c("Annotation", "position"), all.x = TRUE)
    
    PUS7_filter$N_value<-pmax(PUS7_filter$SNR_PUS7KD,PUS7_filter$SNR_SCR,na.rm = TRUE)

    
    #Supplementary Table2
    if(F){
      PUS7_S1<-data.frame("Annotation"=PUS7_filter$Annotation,"chr"=PUS7_filter$chr,"position"=PUS7_filter$position,"strand"=PUS7_filter$strand,"kmer"=PUS7_filter$kmer,
                          "T_PUS7_KD"=PUS7_filter$T_PUS7_KD,"C_PUS7_KD"=PUS7_filter$C_PUS7_KD,
                          "T_SCR"=PUS7_filter$T_SCR,"C_SCR"=PUS7_filter$C_SCR,
                          "T_IVT"=PUS7_filter$T_IVT,"C_IVT"=PUS7_filter$C_IVT,
                          "pvalue.SCR"=PUS7_filter$p.value.SCR,
                          "mm.PUS7KD"=PUS7_filter$mm.PUS7_KD,
                          "mm.SCR"=PUS7_filter$mm.SCR,
                          "mm.IVT"=PUS7_filter$mm.IVT,
                          "SNR"=PUS7_filter$N_value)
      write.csv(PUS7_S1,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/PUS7KD_s2",row.names = F)
    }
    #Supplementary Table8 
    if(F){
      PUS7_s8<-data.frame("<0"=length(which(PUS7_filter$N_value<0)),
                           "0-0.5"=length(which(PUS7_filter$N_value>=0 & PUS7_filter$N_value<0.5)),
                           "0.5-1"=length(which(PUS7_filter$N_value>=0.5 & PUS7_filter$N_value<1)),
                           "1-2"=length(which(PUS7_filter$N_value>=1 & PUS7_filter$N_value<2)),
                           ">2"=length(which(PUS7_filter$N_value>=2)))
      write.csv(PUS7_s8,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/PUS7KD_8",row.names = F)
      
    }

    PUS7_filter_sigma2<-PUS7_filter[which(PUS7_filter$N_value>sig),]     #SNR>1
    
    length(which(PUS7_filter_sigma2$N_value>sig))
    length(which(PUS7_filter_sigma2$color=="blue"))
    
    #filtering for number of reads in PUS7KD with same parameters from TRUB1KD
    PUS7_filter_sigma2<-PUS7_filter_sigma2[which( (PUS7_filter_sigma2$total_PUS7_KD>n_reads_high & PUS7_filter_sigma2$total_SCR>n_reads_high)
                                      | ( (PUS7_filter_sigma2$total_PUS7_KD>n_reads_low & PUS7_filter_sigma2$total_SCR>n_reads_low) & (abs(PUS7_filter_sigma2$mm.PUS7_KDminusSCR) >= difference*PUS7_filter_sigma2$mm.SCR)) ),]

    #add TPM
    PUS7KD_TPM<-read.table(file="/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TPM/PUS7_kd.tsv",sep='\t',header=T)
    PUS7_filter_sigma2$TPM_KD <- 0
    for (i in c(1:nrow(PUS7_filter_sigma2))) {
      row.on.Direct.TPM <- which(PUS7KD_TPM$Gene.Name == PUS7_filter_sigma2$Annotation[i])
      if (length(row.on.Direct.TPM == 1)) {
        PUS7_filter_sigma2$TPM_KD[i] <- PUS7KD_TPM$TPM[row.on.Direct.TPM][[1]]
      }
      else{}
    }
    PUS7_filter_sigma2$PUS7KD_TPM_log <- log10(PUS7_filter_sigma2$TPM_KD)
    
    cutoff=15 #15
    #sites in which we believe
    #the equation of the line going from 30 to 15 delta is y=-0.5x+30.5 where y=mmSCRminusPUS7KD and x=reads
    PUS7KD_sites_all<-PUS7_filter_sigma2[which(PUS7_filter_sigma2$mmSCRminusPUS7KD>=30 | (PUS7_filter_sigma2$mmSCRminusPUS7KD>=cutoff & PUS7_filter_sigma2$total_PUS7_KD>=(30.5-PUS7_filter_sigma2$mmSCRminusPUS7KD)/0.5)),] #n_reads_high
    PUS7KD_sites_grey<-PUS7KD_sites_all[which(PUS7KD_sites_all$color=="black"),]
    PUS7KD_sites<-PUS7KD_sites_all[which(PUS7KD_sites_all$color=="blue"),]
    PUS7KD_sites$CP<-paste0(PUS7KD_sites$chr,PUS7KD_sites$position)
    for (i in c(1:nrow(PUS7KD_sites))) {
      r <- which(PUS7KD_sites$CP==PUS7KD_sites$CP[i])
      if (length(r)>1) {
        
        PUS7KD_sites$Annotation[r[1]]<-paste(PUS7KD_sites$Annotation[r],collapse='-')
        row<-PUS7KD_sites[r[1],]
        PUS7KD_sites <-PUS7KD_sites[-r,]
        PUS7KD_sites[nrow(PUS7KD_sites) + 1,] <- row
        
      }
      
    }
    
    #Supplementary Table3
    if(F){
      PUS7_S1<-data.frame("Annotation"=PUS7KD_sites_all$Annotation,"chr"=PUS7KD_sites_all$chr,"position"=PUS7KD_sites_all$position,"strand"=PUS7KD_sites_all$strand,"kmer"=PUS7KD_sites_all$kmer,
                          "T_PUS7_KD"=PUS7KD_sites_all$T_PUS7_KD,"C_PUS7_KD"=PUS7KD_sites_all$C_PUS7_KD,
                          "T_SCR"=PUS7KD_sites_all$T_SCR,"C_SCR"=PUS7KD_sites_all$C_SCR,
                          "T_IVT"=PUS7KD_sites_all$T_IVT,"C_IVT"=PUS7KD_sites_all$C_IVT,
                          "pvalue.SCR"=PUS7KD_sites_all$p.value.SCR,
                          "mm.PUS7KD"=PUS7KD_sites_all$mm.PUS7_KD,
                          "mm.SCR"=PUS7KD_sites_all$mm.SCR,
                          "mm.IVT"=PUS7KD_sites_all$mm.IVT,
                          "SNR"=PUS7KD_sites_all$N_value)
      write.csv(PUS7_S1,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/PUS7KD_s3.csv",row.names = F)
    }
    
    
    #creating even and odd labels
    PUS7KD_sites<-PUS7KD_sites[order(PUS7KD_sites$N_value),]
    
    rows <- nrow(PUS7KD_sites) 
    # extracting odd rows  
    odd_rows <- seq_len(rows) %% 2 
    # getting data from odd data frame 
    PUS7KD_sites_odd <- PUS7KD_sites[odd_rows == 1, ] 
    
    rows <- nrow(PUS7KD_sites) 
    # extracting odd rows  
    even_rows <- seq_len(rows) %% 2 
    # getting data from odd data frame 
    PUS7KD_sites_even <- PUS7KD_sites[even_rows == 0, ] 
    
    }
  
  plot_figs=T
  #plot KDminusSCR vs N_value 
  if(plot_figs==T){
  min.x = -80
  max.x = 80
  min.y = -4
  max.y = 3
  line.thickness = 1
  label.size = 1
  text.size = 0.8
  point.size = 0.1
  tck.length = 0.01
  tick.thickness = 1
  transparency = 1
 
  file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PUS7_little2_text_sigma",sig,".pdf")
  Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
  
  plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = '',xlab='',ann = F,axes=T,type="l", 
       xaxt='n', yaxt='n', 
       bty='n', bg='transparent')
  
  
  points(x=-PUS7_filter$mm.PUS7_KDminusSCR[which(PUS7_filter$color == "black")], y=log10(PUS7_filter$N_value[which(PUS7_filter$color=="black")]), pch=22, col="grey", bg=adjustcolor("grey",alpha=0.05), cex = log(PUS7_filter$total_PUS7_KD[which(PUS7_filter$color == "black")])*point.size)
  points(x=-PUS7_filter$mm.PUS7_KDminusSCR[which(PUS7_filter$color == "red")], y=log10(PUS7_filter$N_value[which(PUS7_filter$color=="red")]), pch=22, col="orange", bg=adjustcolor("orange",alpha=0.3), cex = log(PUS7_filter$total_PUS7_KD[which(PUS7_filter$color == "yellow")])*point.size)
  points(x=-PUS7_filter$mm.PUS7_KDminusSCR[which(PUS7_filter$color == "blue")], y=log10(PUS7_filter$N_value[which(PUS7_filter$color=="blue")]), pch=22, col="blue", bg=adjustcolor("blue",alpha=0.3), cex = log(PUS7_filter$total_PUS7_KD[which(PUS7_filter$color == "blue")])*point.size)
 
  segments(x0 = 0, y0 = min.y, x1 = 0, y1 = max.y, lty=2, col = "grey3")
  segments(x0 = -60, y0 = log10(sig), x1 = 60, y1 = log10(sig), lty=2, col = "grey3")
  

  text(x=-PUS7KD_sites_odd$mm.PUS7_KDminusSCR-2, y=log10(PUS7KD_sites_odd$N_value)+0.1,PUS7KD_sites_odd$Annotation, col='grey3', cex=0.4)
  text(x=-PUS7KD_sites_even$mm.PUS7_KDminusSCR+2, y=log10(PUS7KD_sites_even$N_value)-0.1,PUS7KD_sites_even$Annotation, col='grey3', cex=0.4)
  text(x=-PUS7KD_sites_grey$mm.PUS7_KDminusSCR[PUS7KD_sites_grey$mmSCRminusPUS7KD>30 & PUS7KD_sites_grey$N_value>10], y=log10(PUS7KD_sites_grey$N_value[PUS7KD_sites_grey$mmSCRminusPUS7KD>30 & PUS7KD_sites_grey$N_value>10]),PUS7KD_sites_grey$Annotation[PUS7KD_sites_grey$mmSCRminusPUS7KD>30 & PUS7KD_sites_grey$N_value>10], col='dimgrey', cex=0.4)
  
  axis(side = 2, at = seq(min.y, max.y, by = 1), labels = seq(min.y, max.y, by = 1),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  axis(side = 1, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  legend("topleft", c("PUS7 motif","TRUB1 motif", "Other"),     
         cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
         pch=c(16,16,16),lty=c(0,0,0), 
         lwd=c(2.5,2.5,2.5),col=c("blue", "orange", "grey"), bty = "n")
  mtext(text = expression(paste(log[10], SNR)),side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
  mtext(text = "Scramble-PUS7KD siRNA (U-C MM%)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
  ##ggsave(filename =file.name,width = 6,height = 6.6, device="eps")
  dev.off() 
  }

  #KD minus SCR vs log10 for sigma>2  
  if(plot_figs==T){
    
    min.x = 1
    max.x = round(log10(max(PUS7_filter_sigma2$total_PUS7_KD)))
    min.y = 0
    max.y = 100
    line.thickness = 1
    label.size = 1
    text.size = 0.8
    point.size = 0.2
    tck.length = 0.01
    tick.thickness = 1
    transparency = 1
    file.name = paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PUS7_mmIVT",mmIVT,"_sigma",sig,"_reads_filterSCR_",n_reads_high,"_text.pdf")
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    
    plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = '',xlab='',ann = F,axes=T,type="l", 
         xaxt='n', yaxt='n', 
         bty='n', bg='transparent')
    
    
    points(x=log10(PUS7_filter_sigma2$total_PUS7_KD[which(PUS7_filter_sigma2$color == "black")]), y=PUS7_filter_sigma2$mmSCRminusPUS7KD[which(PUS7_filter_sigma2$color=="black")], pch=22, col="grey", bg=adjustcolor("grey",alpha=0.05), cex = log2(PUS7_filter_sigma2$N_value[which(PUS7_filter_sigma2$color == "black")])*point.size)
    points(x=log10(PUS7_filter_sigma2$total_PUS7_KD[which(PUS7_filter_sigma2$color == "blue")]), y=PUS7_filter_sigma2$mmSCRminusPUS7KD[which(PUS7_filter_sigma2$color=="blue")], pch=22, col="blue", bg=adjustcolor("blue",alpha=0.3), cex = log2(PUS7_filter_sigma2$N_value[which(PUS7_filter_sigma2$color == "blue")])*point.size)
    points(x=log10(PUS7_filter_sigma2$total_PUS7_KD[which(PUS7_filter_sigma2$color == "red")]), y=PUS7_filter_sigma2$mmSCRminusPUS7KD[which(PUS7_filter_sigma2$color=="red")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.6), cex = log2(PUS7_filter_sigma2$N_value[which(PUS7_filter_sigma2$color == "red")])*point.size)
  
    text(x=log10(PUS7KD_sites_odd$total_PUS7_KD)+0.03, y=PUS7KD_sites_odd$mmSCRminusPUS7KD+0.5,PUS7KD_sites_odd$Annotation, col='grey3', cex=0.3)
    text(x=log10(PUS7KD_sites_even$total_PUS7_KD)+0.03, y=PUS7KD_sites_even$mmSCRminusPUS7KD,PUS7KD_sites_even$Annotation, col='grey3', cex=0.3)
    #grey
    text(x=log10(PUS7KD_sites_grey$total_PUS7_KD[PUS7KD_sites_grey$mmSCRminusPUS7KD>30 & PUS7KD_sites_grey$N_value>10]), y=PUS7KD_sites_grey$mmSCRminusPUS7KD[PUS7KD_sites_grey$mmSCRminusPUS7KD>30 & PUS7KD_sites_grey$N_value>10],PUS7KD_sites_grey$Annotation[PUS7KD_sites_grey$mmSCRminusPUS7KD>30 & PUS7KD_sites_grey$N_value>10], col='dimgrey', cex=0.4)

    segments(x0=log10(10),y0=30,x1=log10(30),y1=15,lty=2,col="black")
    segments(x0=log10(30),y0=15,x1=log10(3200),y1=15,lty=2,col="black")
    
    #for plot KD-SCR vs SCR
    axis(side = 2, at = seq(0, max.y, by = 10), labels = seq(0, max.y, by = 10),
         lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
         las = 1, lwd = line.thickness, line = 0)
    axis(side = 1, at = seq(min.x,max.x, by = 0.5), labels = seq(min.x, max.x, by = 0.5),
         lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
         las = 1, lwd = line.thickness, line = 0)
    legend("topright", c("TRUB1 motif","PUS7 motif", "Other"),     
           cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
           pch=c(16,16,16),lty=c(0,0,0), 
           lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
    mtext(text = "Scramble - PUS7 KD (U-to-C MM%) ",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    mtext(text = "log10(Reads)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
    

    dev.off()
  }  
  
  ############################### seqLogos on the grey and blue motifs
  if (plot_figs==T){
    
   
    # PUS7KD_sites_blue<-PUS7KD_sites$kmer[which(PUS7KD_sites$color=="blue")]
    # kmers<-sapply(PUS7KD_sites_blue,function(x) str_replace_all(x,"T","U"))
    # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PUS7logo_blue.pdf"

    # kmers<-sapply(PUS7KD_sites_grey$kmer,function(x) str_replace_all(x,"T","U"))
    # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PUS7logo_grey.pdf"

    kmers<-sapply(PUS7KD_sites_all$kmer,function(x) str_replace_all(x,"T","U"))
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PUS7logo_merged.pdf"

      
    my_list<-list("new motifs" = kmers)

    min.x = 0
    max.x = 200
    min.y = 0
    max.y = 20
    line.thickness = 1
    label.size = 1
    text.size = 0.8
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
   
    ggseqlogo::ggseqlogo(my_list,ncol=1,method="prob")
    #same as
    #ggplot() + geom_logo(seqs_dna) + theme_logo() + facet_wrap(~seq_group, ncol=4, scales='free_x') 
    
    dev.off()
    
  }

  }  

################################# Model for TRUB1KD and PUS7KD correlation
if(T){
  #real data
  ######################################################### plots for TRUB1KDminusSCR vs PUS7KDminusSCR
  if(T){
    #find TRUB1KD targets in the PUS7KD dataset, sites of interest SOI
    SOI<-c(TRUB1KD_sites_all$ACP,PUS7KD_sites_all$ACP)
    
    #all points with sigma>sig
    sig_sites<-c(TRUB1_filter_sigma2$ACP,PUS7_filter_sigma2$ACP)
    PUS7sig<-PUS7_filter[which(PUS7_filter$ACP %in% sig_sites),]
    TRUB1sig<-TRUB1_filter[which(TRUB1_filter$ACP %in% sig_sites),]
    commonsig<-intersect(PUS7sig$ACP,TRUB1sig$ACP)
    PUS7sigma<-PUS7sig[which(PUS7sig$ACP %in% commonsig),]
    TRUB1sigma<-TRUB1sig[which(TRUB1sig$ACP %in% commonsig),]

    ### dotplot TRUB1KD vs PUS7KD
    if(T){
      min.x = -100
      max.x = 100
      min.y = -100
      max.y = 100
      line.thickness = 1
      label.size = 1
      text.size = 0.8
      point.size = 0.1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 1

      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TRUB1KDvsPUS7KD_2.pdf"
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
  
      plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = '',xlab='',ann = F,axes=T,type="l",
           xaxt='n', yaxt='n',
           bty='n', bg='transparent')
  
      points(x=PUS7sigma$mm.PUS7_KDminusSCR[which(PUS7sigma$color=="black")], y=TRUB1sigma$mm.TRUB1_KDminusSCR[which(TRUB1sigma$color=="black")], pch=16,cex=0.5, col="grey")
      points(x=PUS7sigma$mm.PUS7_KDminusSCR[which(PUS7sigma$color=="red")], y=TRUB1sigma$mm.TRUB1_KDminusSCR[which(TRUB1sigma$color=="red")], pch=16,cex=0.5, col="orange")
      points(x=PUS7sigma$mm.PUS7_KDminusSCR[which(PUS7sigma$color=="blue")], y=TRUB1sigma$mm.TRUB1_KDminusSCR[which(TRUB1sigma$color=="blue")], pch=16,cex=0.5, col="blue")
  
  
      abline(h=0, col="black")
      abline(v=0,col= "black")
  
      
      # text(x=PUS7KDsoi$mm.PUS7_KDminusSCR[which(PUS7KDsoi$color=="red")], y=TRUB1KDsoi$mm.TRUB1_KDminusSCR[which(PUS7KDsoi$color=="red")]+2,TRUB1KDsoi$Annotation[which(PUS7KDsoi$color=="red")], col='grey3', cex=0.3)
      # text(x=PUS7KDsoi$mm.PUS7_KDminusSCR[which(PUS7KDsoi$color=="blue")], y=TRUB1KDsoi$mm.TRUB1_KDminusSCR[which(PUS7KDsoi$color=="blue")]+2,TRUB1KDsoi$Annotation[which(PUS7KDsoi$color=="blue")], col='grey3', cex=0.3)
      
  
      #for plot KD-SCR vs SCR
      axis(side = 2, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
           lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
           las = 1, lwd = line.thickness, line = 0)
      axis(side = 1, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
           lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
           las = 1, lwd = line.thickness, line = 0)
      legend("topright", c("TRUB1 motif","PUS7 motif", "Other"),
             cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
             pch=c(16,16,16),lty=c(0,0,0),
             lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
      mtext(text = "TRUB1 KD - Scramble siRNA (U-C MM%)",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
      mtext(text = "PUS7 KD - Scramble siRNA (U-C MM%)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
      ##ggsave(filename =file.name,width = 6,height = 6.6)
      alpha=0
      dev.off()
    }
  
  
  
  }
  
  #random
  ################################# Model for TRUB1KD and PUS7KD correlation
  if(T){
    SOI<-c(TRUB1KD_sites_all$ACP,PUS7KD_sites_all$ACP)
    #data creation, random model generation merged, null distribution, with mean and std of general kd libraries
    if(T){
      #working on all spots with sigma>2, this will variate according to mean and std of the general KD libraries
      SOI<-c(TRUB1_filter_sigma2$ACP,PUS7_filter_sigma2$ACP)
      #finding the TRUB1 sites in the unfiltered PUS7 dataset, pvalue=anything mmIVT<10
      PUS7soi_all<- PUS7_filter[which(PUS7_filter$ACP %in% SOI),]
      TRUB1soi_all <- TRUB1_filter[which(TRUB1_filter$ACP %in% SOI),]
      common<-intersect(PUS7soi_all$ACP,TRUB1soi_all$ACP)
      PUS7soi<-PUS7soi_all[which(PUS7soi_all$ACP %in% common),]
      TRUB1soi<-TRUB1soi_all[which(TRUB1soi_all$ACP %in% common),]
      
      #working on validated kd sites, this will variate between ranges
      SOI<-c(TRUB1KD_sites_all$ACP,PUS7KD_sites_all$ACP)
      PUS7KDsoi_all<- PUS7_filter[which(PUS7_filter$ACP %in% SOI),]
      TRUB1KDsoi_all <- TRUB1_filter[which(TRUB1_filter$ACP %in% SOI),]
      common<-intersect(PUS7KDsoi_all$ACP,TRUB1KDsoi_all$ACP)
      PUS7KDsoi<-PUS7KDsoi_all[which(PUS7KDsoi_all$ACP %in% common),]
      TRUB1KDsoi<-TRUB1KDsoi_all[which(TRUB1KDsoi_all$ACP %in% common),]
      
      #mean and std deviations of the trub1kd and pus7kd libraries
      mTRUB1=mean(TRUB1_filter$mm.TRUB1_KDminusSCR)
      sdTRUB1=sd(TRUB1_filter$mm.TRUB1_KDminusSCR)
      mPUS7=mean(PUS7_filter$mm.PUS7_KDminusSCR)
      sdPUS7=sd(PUS7_filter$mm.PUS7_KDminusSCR)
  
  
      #generate the min-max values for kd
      min_deltavalue_kd_trub1<-min(TRUB1KDsoi$mm.TRUB1_KDminusSCR)
      max_deltavalue_kd_trub1<-max(TRUB1KDsoi$mm.TRUB1_KDminusSCR)
      min_deltavalue_kd_pus7<-min(PUS7KDsoi$mm.PUS7_KDminusSCR)
      max_deltavalue_kd_pus7<-max(PUS7KDsoi$mm.PUS7_KDminusSCR)
  
      delta_onlyPUS7KD<-c()
      delta_onlyTRUB1KD<-c()
  
      loop=1
      set.seed(12)
      #no mean and std limits
      random_numbers_trub1kd<- runif(loop*(nrow(TRUB1KDsoi)),min=min_deltavalue_kd_trub1,max=max_deltavalue_kd_trub1)
  
      #with mean and std limits
      random_numbers_trub1 <- rnorm(loop*(nrow(TRUB1soi)), mean=mTRUB1,sd=sdTRUB1)
      #random_numbers_trub1 <- pmax(pmin(random_numbers_trub1,max_deltavalue_kd_trub1 ), min_deltavalue_kd_trub1) 
      
      set.seed(123)
      #no mean and std limits
      random_numbers_pus7kd <- runif(loop*(nrow(PUS7KDsoi)), min=min_deltavalue_kd_pus7, max=max_deltavalue_kd_pus7)
  
      #with mean and std limits
      random_numbers_pus7<- rnorm(loop*(nrow(PUS7soi)), mean = mPUS7, sd = sdPUS7)
      #random_numbers_pus7 <- pmin(pmax(random_numbers_pus7, min_deltavalue_kd_pus7), max_deltavalue_kd_pus7)
  
      delta_onlyTRUB1KD<-c(random_numbers_trub1,random_numbers_trub1kd)
      delta_onlyPUS7KD<-c(random_numbers_pus7,random_numbers_pus7kd)
      model<-data.frame("mm.TRUB1minusSCR"=delta_onlyTRUB1KD, "mm.PUS7minusSCR"=delta_onlyPUS7KD)
      # model<-model[-which(model$mm.TRUB1minusSCR>0 & model$mm.PUS7minusSCR>0),]
      # model<-model[-which(model$mm.TRUB1minusSCR>-10 & model$mm.PUS7minusSCR>-10),]
      write.csv(model,paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/model2sigma2_merged_st_mean_ranges_x",loop,".csv"),row.names = F)
  
  
    }
  }

  #comparison    
  if(T){
    model_merged_sigma2_ranges_mean_sd<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/model2sigma2_merged_st_mean_ranges_x1.csv",header=T)
    model<-model_merged_sigma2_ranges_mean_sd
    
    #plot
    if(T){
  
       #plot models separately
       if(T){
        min.x = -100
        max.x = 100
        min.y = -100
        max.y = 100
        line.thickness = 1
        label.size = 1
        text.size = 0.8
        point.size = 0.1
        tck.length = 0.01
        tick.thickness = 1
        transparency = 1

        file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/TRUB1KDvsPUS7KDmodel_merged_sd_mean_sigma2_x1.pdf"
        Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
        CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
        par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
  
        plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
             ylab = '',xlab='',ann = F,axes=T,type="l",
             xaxt='n', yaxt='n',
             bty='n', bg='transparent')
  
        points(x=model$mm.PUS7minusSCR, y=model$mm.TRUB1minusSCR, pch=16,cex=0.5, col="grey")
  
        abline(h=0, col="black")
        abline(v=0,col= "black")
  
        #for plot KD-SCR vs SCR
        axis(side = 2, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
             lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
             las = 1, lwd = line.thickness, line = 0)
        axis(side = 1, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
             lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
             las = 1, lwd = line.thickness, line = 0)
        mtext(text = "TRUB1 KD - Scramble siRNA (U-C MM%)",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
        mtext(text = "PUS7 KD - Scramble siRNA (U-C MM%)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
        ##ggsave(filename =file.name,width = 6,height = 6.6)
        alpha=0
        dev.off()
      }
    }
    
    #tests
    if(T){
      #wilkoxon test
      if(T){
         #wilcox to check if my data are statistically different from the random data
         v=-15
         trub1_bottomleft<-TRUB1sigma$mm.TRUB1_KDminusSCR[which(TRUB1sigma$mm.TRUB1_KDminusSCR<v & (TRUB1sigma$mm.PUS7_KD-TRUB1sigma$mm.SCR)<v)]
         trub1_model_bottomleft<-model$mm.TRUB1minusSCR[which(model$mm.TRUB1minusSCR<v & model$mm.PUS7minusSCR<v)]
         wTRUB1<-wilcox.test(trub1_bottomleft,trub1_model_bottomleft)
         pTRUB1<-wTRUB1$p.value
         pus7_bottomleft<-PUS7sigma$mm.PUS7_KDminusSCR[which(PUS7sigma$mm.PUS7_KDminusSCR<v & (PUS7sigma$mm.TRUB1_KD-PUS7sigma$mm.SCR)<v)]
         pus7_model_bottomleft<-model$mm.PUS7minusSCR[which(model$mm.PUS7minusSCR<v & model$mm.TRUB1minusSCR<v)]
         wPUS<-wilcox.test(pus7_bottomleft,pus7_model_bottomleft)
         pPUS<-wPUS$p.value
         print(paste("W pValue trub1:",pTRUB1))
         print(paste("W pValue pus7:",pPUS))
       }

      #density plots
      if(T){ 
       v=0
       #real TRUB1 
       den<-density(TRUB1sigma$mm.TRUB1_KDminusSCR[which(TRUB1sigma$mm.TRUB1_KDminusSCR<v & (TRUB1sigma$mm.PUS7_KD-TRUB1sigma$mm.SCR)<v)])
       plot(den, lwd=2,ylim=c(0,0.25),col="dodgerblue2",main="KD data",xlab = "TRUB1 KD - Scrambled (control) (% U-to_C error)")
       #random TRUB1
       den2<-density(model$mm.TRUB1minusSCR[which(model$mm.TRUB1minusSCR<v & model$mm.PUS7minusSCR<v)])
       plot(den2,lwd=2, ylim=c(0,0.25),col="grey",main="Random model",xlab= "TRUB1 KD - Scrambled (control) (% U-to_C error)")
       #real PUS7
       den3<-density(PUS7sigma$mm.PUS7_KDminusSCR[which(PUS7sigma$mm.PUS7_KDminusSCR<v & (PUS7sigma$mm.TRUB1_KD-PUS7sigma$mm.SCR)<v)])
       plot(den3, lwd=2,xlim=c(-40,3),ylim=c(0,0.25),col="dodgerblue4",main="KD data",xlab= "PUS7 KD - Scrambled (control) (% U-to_C error)")
       #random PUS7
       den4<-density(model$mm.PUS7minusSCR[which(model$mm.PUS7minusSCR<v & model$mm.TRUB1minusSCR<v)])
       plot(den4, lwd=2,xlim=c(-40,3), ylim=c(0,0.25),col="grey",main="Random model",xlab= "PUS7 KD - Scrambled (control) (% U-to_C error)")
      }
       
      #tailed test
      if(T){
         data <-TRUB1KDsoi$mm.TRUB1_KDminusSCR[which(TRUB1KDsoi$color=="red")]
  
         # Step 2: Perform the Wilcoxon signed-rank test
         # Null hypothesis: The data comes from a distribution with median equal to 0 (standard Gaussian distribution)
         # Alternative hypothesis: The data does not come from a distribution with median equal to 0
  
         # Perform the Wilcoxon signed-rank test
         wT=wilcox.test(data, mu = 0, alternative = "less")
  
         # Step 3: Analyze the direction of the data
         # The test statistic from the Wilcoxon signed-rank test gives information about the direction of the data.
         # If the test statistic is positive, it indicates that the data is shifted to the right (higher values).
         # If the test statistic is negative, it indicates that the data is shifted to the left (lower values).
  
         data <-PUS7KDsoi$mm.PUS7_KDminusSCR[which(PUS7KDsoi$color=="blue")]
         wP=wilcox.test(data, mu = 0, alternative = "less")
         
         print(paste("pValue trub1 tailed, real model is less than a gaussian distribution with mean=0:",wT$p.value))
         print(paste("pValue pus7 tailed, real model is less than a gaussian distribution with mean=0",wP$p.value))
         
       }
    }
    
  }
}


###################### Deseq2 for TRUB1_KD and SCR, gene expression differential analysis
if(T) {
  ###### Deseq2 
  ### Load Libraries
  library("DESeq2")
  library(ggplot2)
  library(Rsamtools)

  ### Read Inputs
  TRUB1_KD_nanocount_Rep1 <- read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/nanocount/random_TRUB1_KD_SH.Rep0.tsv", sep = '\t', header = TRUE)
  TRUB1_KD_nanocount_Rep2 <- read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/nanocount/random_TRUB1_KD_SH.Rep1.tsv", sep = '\t', header = TRUE)
  SCR_nanocount_Rep1 <- read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/nanocount/random_SCRAMBLED_SH_guppy.Rep0.tsv", sep = '\t', header = TRUE)
  SCR_nanocount_Rep2 <- read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/nanocount/random_SCRAMBLED_SH_guppy.Rep1.tsv", sep = '\t', header = TRUE)

  ### Merge all counts
  all.transcripts = unique(c(TRUB1_KD_nanocount_Rep1$transcript_name,
                             TRUB1_KD_nanocount_Rep2$transcript_name,
                             SCR_nanocount_Rep1$transcript_name,
                             SCR_nanocount_Rep2$transcript_name))
  
  countData.mine = data.frame("ensgene" = all.transcripts,
                              "TRUB1_KD_Rep1"=0,
                              "TRUB1_KD_Rep2"=0,
                              "SCR_Rep1"=0,
                              "SCR_Rep2"=0
  )
  countData.mine = countData.mine[order(countData.mine$ensgene),]
  
  SCR_nanocount_Rep1 = SCR_nanocount_Rep1[order(SCR_nanocount_Rep1$transcript_name),]
  SCR_nanocount_Rep2 = SCR_nanocount_Rep2[order(SCR_nanocount_Rep2$transcript_name),]
  TRUB1_KD_nanocount_Rep1 = TRUB1_KD_nanocount_Rep1[order(TRUB1_KD_nanocount_Rep1$transcript_name),]
  TRUB1_KD_nanocount_Rep2 = TRUB1_KD_nanocount_Rep2[order(TRUB1_KD_nanocount_Rep2$transcript_name),]
  
  
  
  # assign the values for each replicate to the right column
  countData.mine$SCR_Rep1[which(countData.mine$ensgene%in%SCR_nanocount_Rep1$transcript_name)] = round(SCR_nanocount_Rep1$est_count)
  countData.mine$SCR_Rep2[which(countData.mine$ensgene%in%SCR_nanocount_Rep2$transcript_name)] = round(SCR_nanocount_Rep2$est_count)
  countData.mine$TRUB1_KD_Rep1[which(countData.mine$ensgene%in%TRUB1_KD_nanocount_Rep1$transcript_name)] = round(TRUB1_KD_nanocount_Rep1$est_count)
  countData.mine$TRUB1_KD_Rep2[which(countData.mine$ensgene%in%TRUB1_KD_nanocount_Rep2$transcript_name)] = round(TRUB1_KD_nanocount_Rep2$est_count)
  
  
  #Filter to the number of reads
  countData.mine = countData.mine[which( countData.mine$TRUB1_KD_Rep1 > 5000 & countData.mine$TRUB1_KD_Rep2 > 5000 &
                                         countData.mine$SCR_Rep1 > 5000 & countData.mine$SCR_Rep2 > 5000),]
  
  # Create metadata
  metaData.mine = data.frame("id"=c("SCR_Rep1","SCR_Rep2","TRUB1_KD_Rep1","TRUB1_KD_Rep2"),
                             "dex"=c("control","control","treated","treated"),   #these are the levels
                             "celltype"=c("SH_SY5Y","SH_SY5Y","SH_SY5Y","SH_SY5Y"),
                             "geo_id"=c("GSM1275862","GSM1275862","GSM1275862","GSM1275862"))
  
#levels order is important  A positive gene fold change means that the gene is upregulated in the treated condition relatively to the control condition.
#so the level you see first is the reference and the level you see second is the condition you are comapring to the reference
#We are performing a differential expression analysis between control(1st) and treated(2nd).
#The first control condition is used as reference level.
#A positive log2 fold change for a gene would mean that gene is more abundant in the 2nd condition (treated) than in the 1st condition (control)  
#What is the biological meaning of a log2  fold change equal to 1 for gene X?
#What is the biological meaning of a log2  fold change equal to -1?
#A log2  equal to 1 means that gene X has a higher expression (x2, two-fold) in the treated condition compared to the control condition.
#A log2  equal to -1 means that gene X has a smaller expression (1/2) in the treated condition.  
    
      
  # create dds
  dds_TRUB1 <- DESeqDataSetFromMatrix(countData=countData.mine, 
                                colData=metaData.mine, 
                                design=~dex, tidy = TRUE)
  dds_TRUB1 <- DESeq(dds_TRUB1)
  res_TRUB1 <- results(dds_TRUB1)
  # Filter to significant ones
  row.on.data <- which(res_TRUB1$pvalue < 0.001)
  
  plotMA(res_TRUB1, ylim=c(-10,10))
  
}
###################### PUS7_KD and SCR
if(T) {
  ###### Deseq2 
  ### Load Libraries
  library("DESeq2")
  library(ggplot2)
  library(Rsamtools)
  
  countData.mine<-0
  
  
  ### Read Inputs
  PUS7_KD_nanocount_Rep1 <- read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/nanocount/random_PUS7_KD_SH_guppy.Rep0.tsv", sep = '\t', header = TRUE)
  PUS7_KD_nanocount_Rep2 <- read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/nanocount/random_PUS7_KD_SH_guppy.Rep1.tsv", sep = '\t', header = TRUE)
  SCR_nanocount_Rep1 <- read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/nanocount/random_SCRAMBLED_SH_guppy.Rep0.tsv", sep = '\t', header = TRUE)
  SCR_nanocount_Rep2 <- read.table(file = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/nanocount/random_SCRAMBLED_SH_guppy.Rep1.tsv", sep = '\t', header = TRUE)
  
  ### Merge all counts
  all.transcripts = unique(c(PUS7_KD_nanocount_Rep1$transcript_name,
                             PUS7_KD_nanocount_Rep2$transcript_name,
                             SCR_nanocount_Rep1$transcript_name,
                             SCR_nanocount_Rep2$transcript_name))
  
  countData.mine = data.frame("ensgene" = all.transcripts,
                              "PUS7_KD_Rep1"=0,
                              "PUS7_KD_Rep2"=0,
                              "SCR_Rep1"=0,
                              "SCR_Rep2"=0
  )
  countData.mine = countData.mine[order(countData.mine$ensgene),]
  
  SCR_nanocount_Rep1 = SCR_nanocount_Rep1[order(SCR_nanocount_Rep1$transcript_name),]
  SCR_nanocount_Rep2 = SCR_nanocount_Rep2[order(SCR_nanocount_Rep2$transcript_name),]
  PUS7_KD_nanocount_Rep1 = PUS7_KD_nanocount_Rep1[order(PUS7_KD_nanocount_Rep1$transcript_name),]
  PUS7_KD_nanocount_Rep2 = PUS7_KD_nanocount_Rep2[order(PUS7_KD_nanocount_Rep2$transcript_name),]
  
  
  
  # assign the values for each replicate to the right column
  countData.mine$SCR_Rep1[which(countData.mine$ensgene%in%SCR_nanocount_Rep1$transcript_name)] = round(SCR_nanocount_Rep1$est_count)
  countData.mine$SCR_Rep2[which(countData.mine$ensgene%in%SCR_nanocount_Rep2$transcript_name)] = round(SCR_nanocount_Rep2$est_count)
  countData.mine$PUS7_KD_Rep1[which(countData.mine$ensgene%in%PUS7_KD_nanocount_Rep1$transcript_name)] = round(PUS7_KD_nanocount_Rep1$est_count)
  countData.mine$PUS7_KD_Rep2[which(countData.mine$ensgene%in%PUS7_KD_nanocount_Rep2$transcript_name)] = round(PUS7_KD_nanocount_Rep2$est_count)
  
  
  #Filter to the number of reads
  # countData.mine = countData.mine[which( countData.mine$PUS7_KD_Rep1 > 150 & countData.mine$PUS7_KD_Rep2 > 150 &
  #                                        countData.mine$SCR_Rep1 > 150 & countData.mine$SCR_Rep2 > 150),]
  
  # Create metadata
  metaData.mine = data.frame("id"=c("SCR_Rep1","SCR_Rep2","PUS7_KD_Rep1","PUS7_KD_Rep2"),
                             "dex"=c("control","control","treated","treated"),
                             "celltype"=c("SH_SY5Y","SH_SY5Y","SH_SY5Y","SH_SY5Y"),
                             "geo_id"=c("GSM1275862","GSM1275862","GSM1275862","GSM1275862"))
  
  
  # create dds
  dds_PUS7 <- DESeqDataSetFromMatrix(countData=countData.mine, 
                                colData=metaData.mine, 
                                design=~dex, tidy = TRUE)
  dds_PUS7 <- DESeq(dds_PUS7)
  res_PUS7 <- results(dds_PUS7)
  # Filter to significant ones
  row.on.data <- which(res_PUS7$pvalue < 0.01)
  
  plotMA(res_PUS7, ylim=c(-10,10))
  
}

###################### Pheatmap from youtube: https://www.youtube.com/watch?v=S2_FTg9kaZU&ab_channel=Sanbomics
if(T) {
  library("pheatmap")
  
  PB_Undiff_deseq2<-res_TRUB1 #change to PUS7 of TRUB1
  dds<-dds_TRUB1
  
  res_filter1 <- PB_Undiff_deseq2[which(PB_Undiff_deseq2$pvalue < 0.001),]
  sigs.df <- as.data.frame(res_filter1)
  mat <- counts(dds, normalized = T)[rownames(sigs.df),]
  ## Z score
  mat.z <- t(apply(mat,1, scale))  #THIS IS THE TABLE WE WANT TO PLOT
  colnames(mat.z) <- rownames(colData)
  res_filter1$transcript_name <- row.names(res_filter1)
  res_filter1$gene <- ""
  res_filter1$gene_ENST <- ""
  for (i in c(1:nrow(res_filter1))) {
    res_filter1$gene[i] <- strsplit(res_filter1$transcript_name[i], "\\|") [[1]][6]
    res_filter1$ENST[i] <- strsplit(res_filter1$transcript_name[i], "\\|") [[1]][1]
    res_filter1$gene_ENST[i] <- paste0(res_filter1$gene[i],"|",res_filter1$ENST[i])
    
  }
  
  rownames(mat.z) <- res_filter1$gene_ENST
  
  pheatmap(mat.z, cluster_rows = T , cluster_columns = T, fontsize = 3,
           column_labels = colnames(mat.z),name = "Z-score",
           row_labels = sigs.df[rownames(mat.z),]$gene)  
  
}


################################## Plot: Barplot of the PUS enzymes expressions for KDs 
if (T) {

  if(T){    
    library(stringr)
    
    #### List of the PUS enzyme genes:
    pus_enzymes <- c("PUS1", "PUSL1", "PUS3", "TRUB1", "TRUB2",
                     "DKC1", "PUS7", "PUS7L", "RPUSD1", "RPUSD2",
                     "RPUSD3", "RPUSD4", "PUS10","GAPDH")
    ###  Housekeeping genes 
    #pus_enzymes <- c("GAPDH", "HPRT", "ACTB","UBC")
    #### Make the data frame of TPM of the pus_enzymes 
    pus_df <- data.frame("gene" = pus_enzymes)
    
    files <- list.files(path="/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/stringtie_10reps", pattern="*.tsv", full.names=TRUE, recursive=FALSE)
    
    for (filename in files){

      if ( grepl("TRUB1", filename )) {
        rep=str_split_i(filename, "/",-1) #take last field
        rep=str_split_i(rep,"[.]",-2) #take the repN
        nam<-paste0("TRUB1",rep)
        f<-read.table(file =filename, sep = '\t', header = TRUE)
        assign(nam,f)
        nam_tpm<-paste(nam,"TPM",sep="_")
        pus_df[,nam_tpm]<-0
        #For the PUS enzymes, add TPMs for each replicate
        for (i in c(1:nrow(pus_df))) {
          pus_df[,nam_tpm][i]<-f$TPM[which(f$Gene.Name == pus_df$gene[i])][1]
         
        }
      }
      
      if ( grepl("PUS7", filename )) {
        rep=str_split_i(filename, "/",-1) #take last field
        rep=str_split_i(rep,"[.]",-2) #take the repN
        nam<-paste0("PUS7",rep)
        f<-read.table(file =filename, sep = '\t', header = TRUE)
        assign(nam,f)
        nam_tpm<-paste(nam,"TPM",sep="_")
        pus_df[,nam_tpm]<-0
        for (i in c(1:nrow(pus_df))) {
          pus_df[,nam_tpm][i]<-f$TPM[which(f$Gene.Name == pus_df$gene[i])][1]
        }
        
      }
      
      if ( grepl("SCR", filename )) {
        rep=str_split_i(filename, "/",-1) #take last field
        rep=str_split_i(rep,"[.]",-2) #take the repN
        nam<-paste0("SCR",rep)
        f<-read.table(file =filename, sep = '\t', header = TRUE)
        assign(nam,f)
        nam_tpm<-paste(nam,"TPM",sep="_")
        pus_df[,nam_tpm]<-0
        for (i in c(1:nrow(pus_df))) {
          pus_df[,nam_tpm][i]<-f$TPM[which(f$Gene.Name == pus_df$gene[i])][1]
        }
        
      }
      
    }
    
    
    #### Add the standard deviation to the pus_df
    pus_df$TRUB1_SD <- 0
    pus_df$PUS7_SD <- 0
    pus_df$Scram_SD <- 0
    for (i in c(1:nrow(pus_df))) {
      pus_df$TRUB1_SD[i] <- sd(c(pus_df$TRUB1Rep0_TPM[i],pus_df$TRUB1Rep1_TPM[i],pus_df$TRUB1Rep2_TPM[i],pus_df$TRUB1Rep3_TPM[i],pus_df$TRUB1Rep4_TPM[i],pus_df$TRUB1Rep5_TPM[i],pus_df$TRUB1Rep6_TPM[i],pus_df$TRUB1Rep7_TPM[i],pus_df$TRUB1Rep8_TPM[i], pus_df$TRUB1Rep9_TPM[i]))
      pus_df$PUS7_SD[i] <- sd(c(pus_df$PUS7Rep0_TPM[i],pus_df$PUS7Rep1_TPM[i],pus_df$PUS7Rep2_TPM[i],pus_df$PUS7Rep3_TPM[i],pus_df$PUS7Rep4_TPM[i],pus_df$PUS7Rep5_TPM[i],pus_df$PUS7Rep6_TPM[i],pus_df$PUS7Rep7_TPM[i],pus_df$PUS7Rep8_TPM[i], pus_df$PUS7Rep9_TPM[i]))
      pus_df$Scram_SD[i] <- sd(c(pus_df$SCRRep0_TPM[i],pus_df$SCRRep1_TPM[i],pus_df$SCRRep2_TPM[i],pus_df$SCRRep3_TPM[i],pus_df$SCRRep4_TPM[i],pus_df$SCRRep5_TPM[i],pus_df$SCRRep6_TPM[i],pus_df$SCRRep7_TPM[i],pus_df$SCRRep8_TPM[i]))                       
    }
    
    #### Add the mean to the pus_df 
    pus_df$TRUB1_mean <- 0
    pus_df$PUS7_mean <- 0
    pus_df$Scram_mean <- 0
    for (i in c(1:nrow(pus_df))) {
      pus_df$TRUB1_mean[i] <- mean(c(pus_df$TRUB1Rep0_TPM[i],pus_df$TRUB1Rep1_TPM[i],pus_df$TRUB1Rep2_TPM[i],pus_df$TRUB1Rep3_TPM[i],pus_df$TRUB1Rep4_TPM[i],pus_df$TRUB1Rep5_TPM[i],pus_df$TRUB1Rep6_TPM[i],pus_df$TRUB1Rep7_TPM[i],pus_df$TRUB1Rep8_TPM[i], pus_df$TRUB1Rep9_TPM[i]))
      pus_df$PUS7_mean[i] <- mean(c(pus_df$PUS7Rep0_TPM[i],pus_df$PUS7Rep1_TPM[i],pus_df$PUS7Rep2_TPM[i],pus_df$PUS7Rep3_TPM[i],pus_df$PUS7Rep4_TPM[i],pus_df$PUS7Rep5_TPM[i],pus_df$PUS7Rep6_TPM[i],pus_df$PUS7Rep7_TPM[i],pus_df$PUS7Rep8_TPM[i], pus_df$PUS7Rep9_TPM[i]))
      pus_df$Scram_mean[i] <- mean(c(pus_df$SCRRep0_TPM[i],pus_df$SCRRep1_TPM[i],pus_df$SCRRep2_TPM[i],pus_df$SCRRep3_TPM[i],pus_df$SCRRep4_TPM[i],pus_df$SCRRep5_TPM[i],pus_df$SCRRep6_TPM[i],pus_df$SCRRep7_TPM[i],pus_df$SCRRep8_TPM[i]))
    }

    
    #ttest for each gene between TRUB1KD and SCR and between PUS7KD and SCR
    for (i in c(1:nrow(pus_df))) {
     
      tT<-t.test(c(pus_df$TRUB1Rep0_TPM[i],pus_df$TRUB1Rep1_TPM[i],pus_df$TRUB1Rep2_TPM[i],pus_df$TRUB1Rep3_TPM[i],pus_df$TRUB1Rep4_TPM[i],pus_df$TRUB1Rep5_TPM[i],pus_df$TRUB1Rep6_TPM[i],pus_df$TRUB1Rep7_TPM[i],pus_df$TRUB1Rep8_TPM[i], pus_df$TRUB1Rep9_TPM[i]),c(pus_df$SCRRep0_TPM[i],pus_df$SCRRep1_TPM[i],pus_df$SCRRep2_TPM[i],pus_df$SCRRep3_TPM[i],pus_df$SCRRep4_TPM[i],pus_df$SCRRep5_TPM[i],pus_df$SCRRep6_TPM[i],pus_df$SCRRep7_TPM[i],pus_df$SCRRep8_TPM[i]))
      pus_df$TRUB1_SCR_ttest[i] <- tT$p.value
     
      tP<- t.test(c(pus_df$PUS7Rep0_TPM[i],pus_df$PUS7Rep1_TPM[i],pus_df$PUS7Rep2_TPM[i],pus_df$PUS7Rep3_TPM[i],pus_df$PUS7Rep4_TPM[i],pus_df$PUS7Rep5_TPM[i],pus_df$PUS7Rep6_TPM[i],pus_df$PUS7Rep7_TPM[i],pus_df$PUS7Rep8_TPM[i], pus_df$PUS7Rep9_TPM[i]),c(pus_df$SCRRep0_TPM[i],pus_df$SCRRep1_TPM[i],pus_df$SCRRep2_TPM[i],pus_df$SCRRep3_TPM[i],pus_df$SCRRep4_TPM[i],pus_df$SCRRep5_TPM[i],pus_df$SCRRep6_TPM[i],pus_df$SCRRep7_TPM[i],pus_df$SCRRep8_TPM[i]))
      pus_df$PUS7_SCR_ttest[i] <- tP$p.value
      
    
    }
    
   
    #fold increase/decrease
    #reference=siRNA 
    #sample to assess= PUS7 or TRUB1KD
    if(T){
    nreplicates=8
    for (i in c(0:nreplicates)){
      namTRUB1<-paste0("TRUB1Rep",i,"_foldchange")
      TRUB1_col_name<-paste0("TRUB1Rep",i,"_TPM")
      namPUS7<-paste0("PUS7Rep",i,"_foldchange")
      PUS7_col_name<-paste0("PUS7Rep",i,"_TPM")
      SCR_col_name<-paste0("SCRRep",i,"_TPM")
      pus_df[,namTRUB1]<-pus_df[,TRUB1_col_name]/pus_df[,SCR_col_name]
      pus_df[,namPUS7]<-pus_df[,PUS7_col_name]/pus_df[,SCR_col_name]
    }
    
    #Add mean and STD
    pus_df$TRUB1_fold_mean<-0
    pus_df$PUS7_fold_mean<-0
    pus_df$TRUB1_fold_std<-0
    pus_df$PUS7_fold_std<-0
    
    for (i in c(1:nrow(pus_df))) {
      # v<-c(pus_df$TRUB1Rep0_foldchange[i],pus_df$TRUB1Rep1_foldchange[i],pus_df$TRUB1Rep2_foldchange[i],pus_df$TRUB1Rep3_foldchange[i],pus_df$TRUB1Rep4_foldchange[i],pus_df$TRUB1Rep5_foldchange[i],pus_df$TRUB1Rep6_foldchange[i],pus_df$TRUB1Rep7_foldchange[i],pus_df$TRUB1Rep8_foldchange[i])
      # v<-v[is.finite(v)]
      # pus_df$TRUB1_fold_mean[i] <- mean(v, na.rm = TRUE)
      # v<-c(pus_df$PUS7Rep0_foldchange[i],pus_df$PUS7Rep1_foldchange[i],pus_df$PUS7Rep2_foldchange[i],pus_df$PUS7Rep3_foldchange[i],pus_df$PUS7Rep4_foldchange[i],pus_df$PUS7Rep5_foldchange[i],pus_df$PUS7Rep6_foldchange[i],pus_df$PUS7Rep7_foldchange[i],pus_df$PUS7Rep8_foldchange[i])
      # v<-v[is.finite(v)]
      # pus_df$PUS7_fold_mean[i] <- mean(v, na.rm=TRUE)
      v<-c(pus_df$TRUB1Rep0_foldchange[i],pus_df$TRUB1Rep1_foldchange[i],pus_df$TRUB1Rep2_foldchange[i],pus_df$TRUB1Rep3_foldchange[i],pus_df$TRUB1Rep4_foldchange[i],pus_df$TRUB1Rep5_foldchange[i],pus_df$TRUB1Rep6_foldchange[i],pus_df$TRUB1Rep7_foldchange[i],pus_df$TRUB1Rep8_foldchange[i])
      v<-v[is.finite(v)]
      pus_df$TRUB1_fold_std[i] <- sd(v,na.rm = TRUE)
      v<-c(pus_df$PUS7Rep0_foldchange[i],pus_df$PUS7Rep1_foldchange[i],pus_df$PUS7Rep2_foldchange[i],pus_df$PUS7Rep3_foldchange[i],pus_df$PUS7Rep4_foldchange[i],pus_df$PUS7Rep5_foldchange[i],pus_df$PUS7Rep6_foldchange[i],pus_df$PUS7Rep7_foldchange[i],pus_df$PUS7Rep8_foldchange[i])
      v<-v[is.finite(v)]
      pus_df$PUS7_fold_std[i] <- sd(v,na.rm = TRUE)
      
      pus_df$TRUB1_fold_mean[i]<-pus_df$TRUB1_mean[i]/pus_df$Scram_mean[i]
      pus_df$PUS7_fold_mean[i]<-pus_df$PUS7_mean[i]/pus_df$Scram_mean[i]
    }
    
    pus_df$TRUB1_fold_mean[which(pus_df$TRUB1_fold_mean<1)]=-1/(pus_df$TRUB1_fold_mean[which(pus_df$TRUB1_fold_mean<1)])
    pus_df$PUS7_fold_mean[which(pus_df$PUS7_fold_mean<1)]=-1/(pus_df$PUS7_fold_mean[which(pus_df$PUS7_fold_mean<1)])
    }
}
    
  #### barplot PUS enzymes Plot 
  if(T){ 
      min.x = 0
      max.x = 14
      min.y = 0
      max.y = 110
      line.thickness = 0.5
      label.size = 1
      text.size = 1
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 0.2

      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/BarplotKD_Bootsrap_stringtie10_points.pdf"
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
      plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = NA,xlab=NA,ann = F,axes=F)
      #box(lwd = line.thickness)
      #### Add the lines
      
      #for TRUB1
      for (i in c(1:nrow(pus_df))) {
        lines(x=c(i-.2,i-.2),y=c(0,pus_df$TRUB1_mean[i]),
              lwd=5 ,lend=1, lty=1, col = "dodgerblue3")
        pts=c(pus_df$TRUB1Rep0_TPM[i],pus_df$TRUB1Rep1_TPM[i],pus_df$TRUB1Rep2_TPM[i],pus_df$TRUB1Rep3_TPM[i],pus_df$TRUB1Rep4_TPM[i],pus_df$TRUB1Rep5_TPM[i],pus_df$TRUB1Rep6_TPM[i],pus_df$TRUB1Rep7_TPM[i],pus_df$TRUB1Rep8_TPM[i], pus_df$TRUB1Rep9_TPM[i])
        points(x=rep(i-.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
      }
      #for PUS7
      for (i in c(1:nrow(pus_df))) {
        lines(x=c(i,i),y=c(0,pus_df$PUS7_mean[i]),
              lwd=5 ,lend=1, lty=1, col = "dodgerblue4")
        pts=c(pus_df$PUS7Rep0_TPM[i],pus_df$PUS7Rep1_TPM[i],pus_df$PUS7Rep2_TPM[i],pus_df$PUS7Rep3_TPM[i],pus_df$PUS7Rep4_TPM[i],pus_df$PUS7Rep5_TPM[i],pus_df$PUS7Rep6_TPM[i],pus_df$PUS7Rep7_TPM[i],pus_df$PUS7Rep8_TPM[i], pus_df$PUS7Rep9_TPM[i])
        points(x=rep(i,length(pts)), y=pts, pch=20,col="black",cex=0.1)
      }
      #for Scramble
      for (i in c(1:nrow(pus_df))) {
        lines(x=c(i+.2,i+.2),y=c(0,pus_df$Scram_mean[i]),
              lwd=5 ,lend=1, lty=1, col = "grey")
        pts=c(pus_df$SCRRep0_TPM[i],pus_df$SCRRep1_TPM[i],pus_df$SCRRep2_TPM[i],pus_df$SCRRep3_TPM[i],pus_df$SCRRep4_TPM[i],pus_df$SCRRep5_TPM[i],pus_df$SCRRep6_TPM[i],pus_df$SCRRep7_TPM[i],pus_df$SCRRep8_TPM[i])
        points(x=rep(i+.2,length(pts)), y=pts, pch=20,col="black",cex=0.1)
      }
      
      
      #### Add the SD line
      for (i in c(1:nrow(pus_df))) {
        lines(x=c(i-.2,i-0.2),
              y=c(pus_df$TRUB1_mean[i]-pus_df$TRUB1_SD[i],pus_df$TRUB1_mean[i]+pus_df$TRUB1_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }
      for (i in c(1:nrow(pus_df))) {
        lines(x=c(i,i),
              y=c(pus_df$PUS7_mean[i]-pus_df$PUS7_SD[i],pus_df$PUS7_mean[i]+pus_df$PUS7_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }
      
      for (i in c(1:nrow(pus_df))) {
        lines(x=c(i+.2,i+0.2),
              y=c(pus_df$Scram_mean[i]-pus_df$Scram_SD[i],pus_df$Scram_mean[i]+pus_df$Scram_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }
      
      
      #### Add Asterisks for p-values
      # *   < 0.05
      # **  < 0.01
      # *** < 0.001
      
      for (i in c(1:nrow(pus_df))){
        

      
        if (pus_df$TRUB1_SCR_ttest[i]<0.001) {
          text(i-0.1,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+10, labels = "***", col = "dodgerblue3", cex = 0.6)
          lines(x=c(i-0.2,i+0.2),y=c(pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+8,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
        } else if (pus_df$TRUB1_SCR_ttest[i]<0.01) {
          text(i-0.1,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+10, labels = "**", col = "dodgerblue3", cex = 0.6)
          lines(x=c(i-0.2,i+0.2),y=c(pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+8,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
        } else if (pus_df$TRUB1_SCR_ttest[i]<0.05) {
          text(i-0.1,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+10, labels = "*", col = "dodgerblue3", cex = 0.6)
          lines(x=c(i-0.2,i+0.2),y=c(pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+8,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
        } 
        
        
        
        if (pus_df$PUS7_SCR_ttest[i]<0.001) {
          text(i+0.1,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+15, labels = "***", col = "dodgerblue4", cex = 0.6)
          lines(x=c(i,i+0.2),y=c(pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+13,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
        } else if (pus_df$PUS7_SCR_ttest[i]<0.01) {
          text(i+0.1,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+15, labels = "**", col = "dodgerblue4", cex = 0.6)
          lines(x=c(i,i+0.2),y=c(pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+13,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
        } else if (pus_df$PUS7_SCR_ttest[i]<0.05) {
          text(i+0.1,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+15, labels = "*", col = "dodgerblue4", cex = 0.6)
          lines(x=c(i,i+0.2),y=c(pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+13,pus_df$Scram_mean[i]+pus_df$Scram_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
        } 
        
          
      }
        
        
      
      
      
      axis(side = 2,at = seq(0, max.y, by = 10) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=0)
      axis(side = 2,at = seq(0, max.y, by = 10),labels = seq(0, max.y, by = 10),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
      axis(side = 1,at = c(0:max.x) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
      axis(side = 1,at = c(1:length(pus_df$gene)),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
      mtext(text = "TPM",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
      legend("topright", c("TRUB1 KD","PUS7 KD","SiRNA control"),
             cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
             pch=c(15,15,15),lty=c(0,0,0),
             lwd=c(2.5,2.5,2.5),col=c("dodgerblue3", "dodgerblue4","grey"), bty = "n")
      
      dev.off() 
    }
    
    
  }

########################## Kinetic model
#([mRNAu][mRNApus])a -> ([mRNAu][mRNApus])b
#a=KD
#b=SCR
if(T){
  #trub1 and pus7 sites validated with KD libraries
  confirmed_sites_names<-c(TRUB1KD_sites$ACP,PUS7KD_sites$ACP)
  #plot
  if(T){
    totPUS7KD=4486086
    totTRUB1KD=9924610
    totSCR=8093571
    scale=1000000
    
    #run the pus_df from the Barplot of the PUS enzymes that is above line 1531
    
    #decide if you want to have trub1 or pus7
    trub1=T
    namefile="/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Kinetic_modelTRUB1_text.pdf"
    # trub1=F
    # namefile="/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Kinetic_modelPUS7_text.pdf"
    # 
    #TRUB1_KD
    if(trub1==T){
      kinetic_sites<-TRUB1_filter[which(TRUB1_filter$ACP %in% confirmed_sites_names),]
      
      M=totTRUB1KD #normalization
      kinetic_sites$ConditionA[which(kinetic_sites$color=="red")] <- (kinetic_sites$total_TRUB1_KD[which(kinetic_sites$color=="red")]/M)*scale*pus_df$TRUB1_mean[pus_df$gene=="TRUB1"]
      #kinetic_sites$ConditionA[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_TRUB1_KD[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$TRUB1_mean[pus_df$gene=="PUS7"]
      kinetic_sites$ConditionA[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_TRUB1_KD[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$TRUB1_mean[pus_df$gene=="TRUB1"]
      kinetic_sites$yConditionA <- kinetic_sites$C_TRUB1_KD/M*scale 
      M=totSCR #normalization
      kinetic_sites$ConditionB[which(kinetic_sites$color=="red")] <- (kinetic_sites$total_SCR[which(kinetic_sites$color=="red")]/M)*scale*pus_df$Scram_mean[pus_df$gene=="TRUB1"]
      #kinetic_sites$ConditionB[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_SCR[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$Scram_mean[pus_df$gene=="PUS7"] 
      kinetic_sites$ConditionB[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_SCR[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$Scram_mean[pus_df$gene=="TRUB1"] #
      kinetic_sites$yConditionB <- kinetic_sites$C_SCR/M*scale 
      
      kinetic_sites$BminusA  <- kinetic_sites$ConditionB - kinetic_sites$ConditionA
      kinetic_sites$yBminusA <- kinetic_sites$yConditionB - kinetic_sites$yConditionA
    }
    
    #PUS7_KD
    if(trub1==F){
      kinetic_sites<-PUS7_filter[which( PUS7_filter$ACP %in% confirmed_sites_names ),]
      
      M=totPUS7KD #normalization
      kinetic_sites$ConditionA<-0
      kinetic_sites$ConditionA[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_PUS7_KD[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$PUS7_mean[pus_df$gene=="PUS7"]
      #kinetic_sites$ConditionA[which(kinetic_sites$color=="red")] <- (kinetic_sites$total_TRUB1_KD[which(kinetic_sites$color=="red")]/M)*scale*pus_df$TRUB1_mean[pus_df$gene=="TRUB1"]
      kinetic_sites$ConditionA[which(kinetic_sites$color=="red")] <- (kinetic_sites$total_PUS7_KD[which(kinetic_sites$color=="red")]/M)*scale*pus_df$PUS7_mean[pus_df$gene=="PUS7"]
      kinetic_sites$yConditionA <- kinetic_sites$C_PUS7_KD/M*scale 
      M=totSCR #normalization
      kinetic_sites$ConditionB<-0
      kinetic_sites$ConditionB[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_SCR[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$Scram_mean[pus_df$gene=="PUS7"]
      #kinetic_sites$ConditionB[which(kinetic_sites$color=="red")] <- (kinetic_sites$total_SCR[which(kinetic_sites$color=="red")]/M)*scale*pus_df$Scram_mean[pus_df$gene=="TRUB1"] 
      kinetic_sites$ConditionB[which(kinetic_sites$color=="red")] <- (kinetic_sites$total_SCR[which(kinetic_sites$color=="red")]/M)*scale*pus_df$Scram_mean[pus_df$gene=="PUS7"] #
      kinetic_sites$yConditionB <- kinetic_sites$C_SCR/M*scale 
      
      kinetic_sites$BminusA  <- kinetic_sites$ConditionB - kinetic_sites$ConditionA
      kinetic_sites$yBminusA <- kinetic_sites$yConditionB - kinetic_sites$yConditionA
    }
    
    
    #plot conditionB-conditionA
    if(T){
      min.x = -100 #round(min(kinetic_sites$BminusA))-100 
      max.x = 2000 #round(max(kinetic_sites$BminusA))+500 
      min.y = -10
      max.y = 65
      line.thickness = 1
      label.size = 1
      text.size = 0.8
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 1

      file.name = namefile
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
      
      plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = '',xlab='',ann = F,axes=T,type="l", 
           xaxt='n', yaxt='n', 
           bty='n', bg='transparent')
      
      
      points(x=kinetic_sites$BminusA[which(kinetic_sites$color == "red")], y=kinetic_sites$yBminusA[which(kinetic_sites$color=="red")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.3), cex = point.size)
      
      points(x=kinetic_sites$BminusA[which(kinetic_sites$color == "blue")], y=kinetic_sites$yBminusA[which(kinetic_sites$color=="blue")], pch=22, col="blue", bg=adjustcolor("blue",alpha=0.6), cex=point.size)
      
      # Correlation TRUB1 sites only 
      x <- kinetic_sites$BminusA[which(kinetic_sites$color == "red")]
      y <- kinetic_sites$yBminusA[which(kinetic_sites$color=="red")]
      
      # Calculate correlation coefficient
      correlation_coefficient <- cor(x, y)
      # Calculate R-squared value
      r_squared <- correlation_coefficient^2
      
      # Print the correlation coefficient and R-squared value
      print(paste("Correlation Coefficient trub1:", correlation_coefficient))
      print(paste("R-squared Value:", r_squared))
      
      # Add correlation line
      abline(lm(y ~ x), col = "orange")
      
      
      x<-kinetic_sites$BminusA[which(kinetic_sites$color == "blue")]
      y<-kinetic_sites$yBminusA[which(kinetic_sites$color=="blue")]
      
      # Calculate correlation coefficient
      correlation_coefficient <- cor(x, y)
      # Calculate R-squared value
      r_squared <- correlation_coefficient^2
      
      # Add correlation line
      abline(lm(y ~ x), col = "blue")
      
      # Print the correlation coefficient and R-squared value
      print(paste("Correlation Coefficient pus7:", correlation_coefficient))
      print(paste("R-squared Value:", r_squared))
      
      abline(v=0, col = "black")
      
      
      text(x=kinetic_sites$BminusA, y=kinetic_sites$yBminusA,kinetic_sites$Annotation, col='grey3', cex=0.4)
      
      
      #for plot KD-SCR vs SCR
      axis(side = 2, at = seq(min.y, max.y, by = 10), labels = seq(min.y, max.y, by = 10),
           lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
           las = 1, lwd = line.thickness, line = 0)
      axis(side = 1, at = seq(min.x,max.x, by = 100), labels = seq(min.x, max.x, by = 100),
           lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
           las = 1, lwd = line.thickness, line = 0)
      legend("topright", c("TRUB1 validated sites","PUS7 validated sites", "Other"),     
             cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
             pch=c(16,16,16),lty=c(0,0,0), 
             lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
      mtext(text = expression(" ψScr - ψKD per million "),side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
      mtext(text = "([mRNA U][mRNA PUS])Scr - ([mRNA U][mRNA PUS])KD",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
      ##ggsave(filename =file.name,width = 6,height = 6.6)
      alpha=0
      dev.off()
      
      
    }
    
  }
  
  #CDS UTR kinetic model
  if(T){
  #read once per session
    library(rtracklayer)
    g = readGFF("/home/sasha/Rouhanifard lab Dropbox/Sasha/CellLinesPVal/gencode/gencode.v27.annotation.gff3")
    gff3=as.data.frame(g)
    gff3 <- gff3[which(gff3$gene_type == "protein_coding" & ((gff3$type == "CDS") | (gff3$type == "three_prime_UTR") | (gff3$type == "five_prime_UTR") | (gff3$type == "start_codon") | (gff3$type == "stop_codon") | (gff3$type == "stop_codon_redefined_as_selenocysteine"))), c("type", "seqid","start", "end", "gene_name")]
    
    #Confirmed sites TRUB1
    if(T){
      library(rtracklayer)
      TRUB1KD_sites$CDS = F
      TRUB1KD_sites$Three_prime_UTR = F
      TRUB1KD_sites$Five_Prime_UTR = F
      TRUB1KD_sites$start_codon = F
      TRUB1KD_sites$stop_codon = F
      TRUB1KD_sites$stop_codon_redefined_as_selenocysteine = F
      
      for (row in 1:nrow(TRUB1KD_sites)){
    
        type = gff3$type[which(gff3$start <= TRUB1KD_sites$position[row] & 
                                 gff3$end >= TRUB1KD_sites$position[row] &
                                 gff3$seqid == TRUB1KD_sites$chr[row])]
        
        if (length(type)>0){
          if ("CDS" %in% type) {
            TRUB1KD_sites$CDS[row] = T
          }
          if ("three_prime_UTR" %in% type) {
            TRUB1KD_sites$Three_prime_UTR[row] = T
          }
          if ("five_prime_UTR" %in% type) {
            TRUB1KD_sites$Five_Prime_UTR[row] = T
          }
          if ("start_codon" %in% type) {
            TRUB1KD_sites$start_codon[row] = T
          }
          if ("stop_codon" %in% type) {
            TRUB1KD_sites$stop_codon[row] = T
          }
          if ("stop_codon_redefined_as_selenocysteine" %in% type) {
            TRUB1KD_sites$stop_codon_redefined_as_selenocysteine[row] = T
          }
        }
      }
      
      
      start.codon <- TRUB1KD_sites[which(TRUB1KD_sites$start_codon == 1),]
      
      all.false <- TRUB1KD_sites[which(TRUB1KD_sites$CDS == 0 &
                                 TRUB1KD_sites$Three_prime_UTR == 0 &
                                 TRUB1KD_sites$Five_Prime_UTR == 0 &
                                 TRUB1KD_sites$start_codon == 0 &
                                 TRUB1KD_sites$stop_codon == 0),]
    }
    #Confirmed sites PUS7
    if(T){
      library(rtracklayer)
      PUS7KD_sites$CDS = F
      PUS7KD_sites$Three_prime_UTR = F
      PUS7KD_sites$Five_Prime_UTR = F
      PUS7KD_sites$start_codon = F
      PUS7KD_sites$stop_codon = F
      PUS7KD_sites$stop_codon_redefined_as_selenocysteine = F
      
      for (row in 1:nrow(PUS7KD_sites)){
        
        type = gff3$type[which(gff3$start <= PUS7KD_sites$position[row] & 
                                 gff3$end >= PUS7KD_sites$position[row] &
                                 gff3$seqid == PUS7KD_sites$chr[row])]
        
        if (length(type)>0){
          if ("CDS" %in% type) {
            PUS7KD_sites$CDS[row] = T
          }
          if ("three_prime_UTR" %in% type) {
            PUS7KD_sites$Three_prime_UTR[row] = T
          }
          if ("five_prime_UTR" %in% type) {
            PUS7KD_sites$Five_Prime_UTR[row] = T
          }
          if ("start_codon" %in% type) {
            PUS7KD_sites$start_codon[row] = T
          }
          if ("stop_codon" %in% type) {
            PUS7KD_sites$stop_codon[row] = T
          }
          if ("stop_codon_redefined_as_selenocysteine" %in% type) {
            PUS7KD_sites$stop_codon_redefined_as_selenocysteine[row] = T
          }
        }
      }
      
      
      start.codon <- PUS7KD_sites[which(PUS7KD_sites$start_codon == 1),]
      
      all.false <- PUS7KD_sites[which(PUS7KD_sites$CDS == 0 &
                                 PUS7KD_sites$Three_prime_UTR == 0 &
                                 PUS7KD_sites$Five_Prime_UTR == 0 &
                                 PUS7KD_sites$start_codon == 0 &
                                 PUS7KD_sites$stop_codon == 0),]
    }
    
    #plot
    if(T){
      min.x = 0
      max.x = 4
      min.y = 0
      max.y = 50
      label.size = 1
      text.size = 0.8
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 0.2
      file.name ="/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/UTR-CDS.pdf"
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
      plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = NA,xlab=NA,ann = F,axes=F)
      #box(lwd = line.thickness)
      #### Add the lines
      
      #normalizing for number of reads
      max_reads= 1#max(c(sum(TRUB1KD_sites$total_TRUB1_KD),sum(PUS7KD_sites$total_PUS7_KD)))
      factor_norm=1#000
      spacer=0.1
      
      lines(x=c(3-spacer,3-spacer), y=c(0,length(which(PUS7KD_sites$Three_prime_UTR==T))/max_reads*factor_norm),lwd=10 ,lend=1, lty=1, col = "blue")
      lines(x=c(3+spacer,3+spacer), y=c(0,length(which(TRUB1KD_sites$Three_prime_UTR==T))/max_reads*factor_norm),lwd=10 ,lend=1, lty=1, col = "orange")
  
      lines(x=c(2-spacer,2-spacer), y=c(0,length(which(PUS7KD_sites$CDS==T))/max_reads*factor_norm),lwd=10 ,lend=1, lty=1, col ="blue")
      lines(x=c(2+spacer,2+spacer), y=c(0,length(which(TRUB1KD_sites$CDS==T))/max_reads*factor_norm),lwd=10 ,lend=1, lty=1, col = "orange")
  
      lines(x=c(1-spacer,1-spacer), y=c(0,length(which(PUS7KD_sites$Five_Prime_UTR==T)*factor_norm)/max_reads),lwd=10 ,lend=1, lty=1, col ="blue")
      lines(x=c(1+spacer,1+spacer), y=c(0,length(which(TRUB1KD_sites$Five_Prime_UTR==T))/max_reads*factor_norm),lwd=10 ,lend=1, lty=1, col = "orange")
  
      
      axis(side = 2,at = seq(min.y, max.y, by = 1),labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=0)
      axis(side = 2,at = seq(min.y, max.y, by = 1),seq(min.y, max.y, by = 1),lwd.ticks = 0 ,cex.axis=label.size-0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-0.5 )
      axis(side = 1,at = c(min.x:max.x) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
           tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
      axis(side = 1,at = c(1:3),labels = c("5'UTR","CDS", "3'UTR"),lwd.ticks = 0 ,cex.axis=text.size ,tck= 0 ,las=2 ,lwd = 0, line=-1)
      legend("topright", c("TRUB1","PUS7"),
             cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
             pch=c(16,16),lty=c(0,0),
             lwd=c(2.5,2.5),col=c("orange","blue"), bty = "n")
      mtext(text = "Number of validated sites",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
      dev.off() 
    }
  }
}


###########################################################################################################################################################################
###########################################################################################################################################################################
####################################################### Figure 3 ##########################################################################################################
#panIVT enrichment data creation, run once
if(T){
panIVT_PB<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/panIVT/PB-panIVT.Merged_with_P_vals-COUNT.csv", header = T)
panIVT_PB$ACP<-paste0(panIVT_PB$Annotation,panIVT_PB$chr,panIVT_PB$position)
panIVT_Diff<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/panIVT/DIFF-panIVT.Merged_with_P_vals-COUNT.csv", header = T)
panIVT_Diff$ACP<-paste0(panIVT_Diff$Annotation,panIVT_Diff$chr,panIVT_Diff$position)
panIVT_Undiff<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/panIVT/UNDIFF-panIVT.Merged_with_P_vals-COUNT.csv", header = T)
panIVT_Undiff$ACP<-paste0(panIVT_Undiff$Annotation,panIVT_Undiff$chr,panIVT_Undiff$position)

#mergind Diff with Undiff 
panIVT_SH_Diff<- bind_rows(panIVT_Diff,panIVT_Undiff) 
panIVT_SH_Diff<- data.frame("ACP"=panIVT_SH_Diff$ACP,"Annotation"=panIVT_SH_Diff$Annotation,"chr"=panIVT_SH_Diff$chr,"position"=panIVT_SH_Diff$position,"N_reads_IVT"=panIVT_SH_Diff$N_reads_IVT,"mm.IVT"=panIVT_SH_Diff$mm.IVT,"T_IVT"=panIVT_SH_Diff$T_IVT,"C_IVT"==panIVT_SH_Diff$C_IVT)
panIVT_SH_Diff<- distinct(panIVT_SH_Diff)
panIVT_SH_Diff<- panIVT_SH_Diff[which(panIVT_SH_Diff$N_reads_IVT>10),]
panIVT_SH_Diff<- panIVT_SH_Diff[which(panIVT_SH_Diff$mm.IVT <=10),]
#merging PB and Undiff
panIVT_SH_PB<- bind_rows(panIVT_PB,panIVT_Undiff) 
panIVT_SH_PB<- data.frame("ACP"=panIVT_SH_PB$ACP,"Annotation"=panIVT_SH_PB$Annotation,"chr"=panIVT_SH_PB$chr,"position"=panIVT_SH_PB$position,"N_reads_IVT"=panIVT_SH_PB$N_reads_IVT,"mm.IVT"=panIVT_SH_PB$mm.IVT,"T_IVT"=panIVT_SH_PB$T_IVT,"C_IVT"==panIVT_SH_PB$C_IVT)
panIVT_SH_PB<- distinct(panIVT_SH_PB)
panIVT_SH_PB<- panIVT_SH_PB[which(panIVT_SH_PB$N_reads_IVT>10),]
panIVT_SH_PB<- panIVT_SH_PB[which(panIVT_SH_PB$mm.IVT <=10),]

#Diff vs Undiff
if(T){
  ### Diff and Undiff data
  SH <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/SH.csv", header = T)
  
  ### Filter based on the p-value < 0.001
  SH_filter <- SH[which( ((SH$p.value.Diff1 < 0.001 & SH$p.value.Diff2 < 0.001) | 
                            (SH$p.value.Diff1 < 0.001 & SH$p.value.Diff3 < 0.001) |
                            (SH$p.value.Diff2 < 0.001 & SH$p.value.Diff3 < 0.001)) |
                           ((SH$p.value.Undiff1 < 0.001 & SH$p.value.Undiff2 < 0.001) | 
                              (SH$p.value.Undiff1 < 0.001 & SH$p.value.Undiff3 < 0.001) |
                              (SH$p.value.Undiff2 < 0.001 & SH$p.value.Undiff3 < 0.001)) ),] 
  
  ### Add the total for C_Diff, T_Diff , C_Undiff, T_Undiff
  SH_filter$T_Diff <- 0
  SH_filter$C_Diff <- 0
  SH_filter$T_Undiff <- 0
  SH_filter$C_Undiff <- 0
  for (i in c(1:nrow(SH_filter))) {
    SH_filter$T_Diff[i] <- SH_filter$T_Diff1[i] + SH_filter$T_Diff2[i] + SH_filter$T_Diff3[i]
    SH_filter$C_Diff[i] <- SH_filter$C_Diff1[i] + SH_filter$C_Diff2[i] + SH_filter$C_Diff3[i]
    SH_filter$T_Undiff[i] <- SH_filter$T_Undiff1[i] + SH_filter$T_Undiff2[i] + SH_filter$T_Undiff3[i]
    SH_filter$C_Undiff[i] <- SH_filter$C_Undiff1[i] + SH_filter$C_Undiff2[i] + SH_filter$C_Undiff3[i]
  }
  
  ### Add total_Diff & total_Undiff
  SH_filter$total_Diff <- 0
  SH_filter$total_Undiff <- 0
  SH_filter$total_IVT <- 0
  for (i in c(1:nrow(SH_filter))) {
    SH_filter$total_Diff[i] <- SH_filter$T_Diff[i] + SH_filter$C_Diff[i]
    SH_filter$total_Undiff[i] <- SH_filter$T_Undiff[i] + SH_filter$C_Undiff[i]
    SH_filter$total_IVT[i] <- SH_filter$T_IVT[i] + SH_filter$C_IVT[i]
    
  }
  
  ### Add totalUndiff1 2 3
  SH_filter$total_Undiff1 <- 0
  SH_filter$total_Undiff2 <- 0
  SH_filter$total_Undiff3 <- 0
  for (i in c(1:nrow(SH_filter))) {
    SH_filter$total_Undiff1[i] <-  SH_filter$T_Undiff1[i] + SH_filter$C_Undiff1[i]
    SH_filter$total_Undiff2[i] <-  SH_filter$T_Undiff2[i] + SH_filter$C_Undiff2[i]
    SH_filter$total_Undiff3[i] <-  SH_filter$T_Undiff3[i] + SH_filter$C_Undiff3[i]
  }
  
  ### Add totalPB1 2 3
  SH_filter$total_Diff1 <- 0
  SH_filter$total_Diff2 <- 0
  SH_filter$total_Diff3 <- 0
  for (i in c(1:nrow(SH_filter))) {
    SH_filter$total_Diff1[i] <- SH_filter$T_Diff1[i] + SH_filter$C_Diff1[i]
    SH_filter$total_Diff2[i] <- SH_filter$T_Diff2[i] + SH_filter$C_Diff2[i]
    SH_filter$total_Diff3[i] <- SH_filter$T_Diff3[i] + SH_filter$C_Diff3[i] 
  }
  
  ### Add the mismatches mm.Diff & mm.Undiff
  SH_filter$mm.Diff <- 0
  SH_filter$mm.Undiff <- 0
  for (i in c(1:nrow(SH_filter))) {
    SH_filter$mm.Diff[i] <- (SH_filter$C_Diff[i]/ (SH_filter$C_Diff[i] + SH_filter$T_Diff[i]))*100
    SH_filter$mm.Undiff[i] <- (SH_filter$C_Undiff[i] / (SH_filter$C_Undiff[i] + SH_filter$T_Undiff[i]))*100
    
  }
  
  ### Add the mismatch differences between Diff and IVT and Undiff and IVT
  SH_filter$mm.DiffMINUSmm.IVT <- 0
  SH_filter$mm.UndiffMINUSmm.IVT <- 0
  for (i in c(1:nrow(SH_filter))) {
    SH_filter$mm.DiffMINUSmm.IVT[i] <- SH_filter$mm.Diff[i]-SH_filter$mm.IVT[i]
    SH_filter$mm.UndiffMINUSmm.IVT[i] <- SH_filter$mm.Undiff[i]-SH_filter$mm.IVT[i]
  }
  
  ### Add the mismatch differences between Diff and undiff
  SH_filter$mm.DiffMINUSmm.Undiff <- 0
  for (i in c(1:nrow(SH_filter))) {
    SH_filter$mm.DiffMINUSmm.Undiff[i] <- SH_filter$mm.Diff[i]-SH_filter$mm.Undiff[i]
  }
  
  ### Filter number of reads to total_diff & total_Undiff > 10
  SH_filter <- SH_filter[which(SH_filter$total_Diff > 10 & SH_filter$total_Undiff > 10),]
  
  ### Substitute mmIVT==0 with panIVT
  SH_filter$ACP<-paste0(SH_filter$Annotation,SH_filter$chr,SH_filter$position)
  ACP<-SH_filter$ACP[which(SH_filter$total_IVT < 10)]
  ACP_pan<-panIVT_SH_Diff$ACP[which(panIVT_SH_Diff$ACP %in% ACP)]
  ACP_pan_mm<-panIVT_SH_Diff$mm.IVT[which(panIVT_SH_Diff$ACP %in% ACP)]
  pan<-panIVT_SH_Diff[which(panIVT_SH_Diff$ACP %in% ACP),]
  for (i in c(1:length(ACP_pan))) {
    if (ACP_pan[i] %in% SH_filter$ACP){
      idx<-which(SH_filter$ACP==ACP_pan[i])
      SH_filter$total_IVT[idx]<-panIVT_SH_Diff$N_reads_IVT[which(panIVT_SH_Diff$ACP %in% ACP_pan[i])]
      SH_filter$T_IVT[idx]<-panIVT_SH_Diff$T_IVT[which(panIVT_SH_Diff$ACP %in% ACP_pan[i])]
      SH_filter$C_IVT[idx]<-panIVT_SH_Diff$C_IVT[which(panIVT_SH_Diff$ACP %in% ACP_pan[i])]
      SH_filter$mm.IVT[idx]<-panIVT_SH_Diff$mm.IVT[which(panIVT_SH_Diff$ACP %in% ACP_pan[i])]
    }
  }
  
  SH_filter <- SH_filter[which(SH_filter$total_IVT > 10),]
  SH_filter <- SH_filter[which((SH_filter$T_IVT+SH_filter$C_IVT)!=0),]
  ##### Filter to the ones that has mmIVT < 10 
  SH_filter <- SH_filter[which(SH_filter$mm.IVT <10),]
  
  ### Add the color column based on the kmer
  SH_filter$color <- "black"
  for (i in c(1:nrow(SH_filter))) {
    if (SH_filter$kmer[i] %in% c("TGTAG", "TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) {
      SH_filter$color[i] <- "blue"
    }
    if (SH_filter$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
      SH_filter$color[i] <- "red"
    }
    else {}
  }
  
}
write.csv(SH_filter,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Diff_filter_pan.csv",row.names = F)

#PB vs Undiff
if(T){
  
  PB <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Data/PB.csv",header = T)
  
 
  PB_filter <- PB[which(PB$mm.IVT <=10),]
  PB_filter <- PB[which(((PB$p.value.PB1 < 0.001 & PB$p.value.PB2 < 0.001)) |
                          ((PB$p.value.Undiff1 < 0.001 & PB$p.value.Undiff2 < 0.001) | 
                             (PB$p.value.Undiff1 < 0.001 & PB$p.value.Undiff3 < 0.001) |
                             (PB$p.value.Undiff2 < 0.001 & PB$p.value.Undiff3 < 0.001)) ),] 
  
  
  ### Add the total for C_Diff, T_Diff , C_Undiff, T_Undiff
  PB_filter$T_PB <- 0
  PB_filter$C_PB <- 0
  PB_filter$T_Undiff <- 0
  PB_filter$C_Undiff <- 0
  for (i in c(1:nrow(PB_filter))) {
    PB_filter$T_PB[i] <- PB_filter$T_PB1[i] + PB_filter$T_PB2[i] 
    PB_filter$C_PB[i] <- PB_filter$C_PB1[i] + PB_filter$C_PB2[i] 
    PB_filter$T_Undiff[i] <- PB_filter$T_Undiff1[i] + PB_filter$T_Undiff2[i] + PB_filter$T_Undiff3[i]
    PB_filter$C_Undiff[i] <- PB_filter$C_Undiff1[i] + PB_filter$C_Undiff2[i] + PB_filter$C_Undiff3[i]
  }
  ### Add total_PB & total_Undiff
  PB_filter$total_PB <- 0
  PB_filter$total_Undiff <- 0
  PB_filter$total_IVT[i]
  for (i in c(1:nrow(PB_filter))) {
    PB_filter$total_PB[i] <- PB_filter$T_PB[i] + PB_filter$C_PB[i]
    PB_filter$total_Undiff[i] <- PB_filter$T_Undiff[i] + PB_filter$C_Undiff[i]
    PB_filter$total_IVT[i] <- PB_filter$T_IVT[i] + PB_filter$C_IVT[i]
  }
  
  ### Add totalUndiff1 2 3
  PB_filter$total_Undiff1 <- 0
  PB_filter$total_Undiff2 <- 0
  PB_filter$total_Undiff3 <- 0
  for (i in c(1:nrow(PB_filter))) {
    PB_filter$total_Undiff1[i] <-  PB_filter$T_Undiff1[i] + PB_filter$C_Undiff1[i]
    PB_filter$total_Undiff2[i] <-  PB_filter$T_Undiff2[i] + PB_filter$C_Undiff2[i]
    PB_filter$total_Undiff3[i] <-  PB_filter$T_Undiff3[i] + PB_filter$C_Undiff3[i]
  }
  
  ### Add totalPB1 2 3
  PB_filter$total_PB1 <- 0
  PB_filter$total_PB2 <- 0
  for (i in c(1:nrow(PB_filter))) {
    PB_filter$total_PB1[i] <- PB_filter$T_PB1[i] + PB_filter$C_PB1[i]
    PB_filter$total_PB2[i] <- PB_filter$T_PB2[i] + PB_filter$C_PB2[i]
  }
  
  
  ### Add the mismatches mm.PB & mm.Undiff
  PB_filter$mm.PB <- 0
  PB_filter$mm.Undiff <- 0
  for (i in c(1:nrow(PB_filter))) {
    PB_filter$mm.PB[i] <- (PB_filter$C_PB[i]/ (PB_filter$C_PB[i] + PB_filter$T_PB[i]))*100
    PB_filter$mm.Undiff[i] <- (PB_filter$C_Undiff[i] / (PB_filter$C_Undiff[i] + PB_filter$T_Undiff[i]))*100
    
  }
  
  ### Add the mismatch difference between Diff and IVT and Undiff and IVT
  PB_filter$mm.PBMINUSmm.IVT <- 0
  PB_filter$mm.UndiffMINUSmm.IVT <- 0
  for (i in c(1:nrow(PB_filter))) {
    PB_filter$mm.PBMINUSmm.IVT[i] <- PB_filter$mm.PB[i]-PB_filter$mm.IVT[i]
    PB_filter$mm.UndiffMINUSmm.IVT[i] <- PB_filter$mm.Undiff[i]-PB_filter$mm.IVT[i]
  }
  
  
  ### Add the mismatch difference between PB and undiff
  PB_filter$mm.PBMINUSmm.Undiff <- 0
  for (i in c(1:nrow(PB_filter))) {
    PB_filter$mm.PBMINUSmm.Undiff[i] <- PB_filter$mm.PB[i]-PB_filter$mm.Undiff[i]
    
  }
  ### Filter to total_PB & total_Undiff > 10
  PB_filter <- PB_filter[which(PB_filter$total_PB > 10 & PB_filter$total_Undiff > 10),]
  
  ### Substitute mmIVT==0 with panIVT
  PB_filter$ACP<-paste0(PB_filter$Annotation,PB_filter$chr,PB_filter$position)
  ACP<-PB_filter$ACP[which(PB_filter$total_IVT < 10)]
  ACP_pan<-panIVT_SH_PB$ACP[which(panIVT_SH_PB$ACP %in% ACP)]
  ACP_pan_mm<-panIVT_SH_PB$mm.IVT[which(panIVT_SH_PB$ACP %in% ACP)]
  pan<-panIVT_SH_PB[which(panIVT_SH_PB$ACP %in% ACP),]
  for (i in c(1:length(ACP_pan))) {
    if (ACP_pan[i] %in% PB_filter$ACP){
      idx<-which(PB_filter$ACP==ACP_pan[i])
      PB_filter$total_IVT[idx]<-panIVT_SH_PB$N_reads_IVT[which(panIVT_SH_PB$ACP %in% ACP_pan[i])]
      PB_filter$T_IVT[idx]<-panIVT_SH_PB$T_IVT[which(panIVT_SH_PB$ACP %in% ACP_pan[i])]
      PB_filter$C_IVT[idx]<-panIVT_SH_PB$C_IVT[which(panIVT_SH_PB$ACP %in% ACP_pan[i])]
      PB_filter$mm.IVT[idx]<-panIVT_SH_PB$mm.IVT[which(panIVT_SH_PB$ACP %in% ACP_pan[i])]
    }
  }
  
  PB_filter <- PB_filter[which(PB_filter$total_IVT > 10),]
  PB_filter <- PB_filter[which((PB_filter$T_IVT+PB_filter$C_IVT)!=0),]
  PB_filter <- PB_filter[which(PB_filter$mm.IVT<10),]
  
  
  ### Add the color column based on the kmer
  PB_filter$color <- "black"
  for (i in c(1:nrow(PB_filter))) {
    if (PB_filter$kmer[i] %in% c("TGTAG", "TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")) {
      PB_filter$color[i] <- "blue"
    }
    if (PB_filter$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
      PB_filter$color[i] <- "red"
    }
    else {}
  }
}
write.csv(PB_filter,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PB_filter_pan.csv",row.names = F)

#Undiff
if(T){
  SHSY5Y<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PILEUP/init_gene_pileup/finalized_pileup/SH.csv", header = T)
  SHSY5Y_filter <- SHSY5Y[which( ((SHSY5Y$p.value.Undiff1 < 0.001 & SHSY5Y$p.value.Undiff2 < 0.001) | 
                                    (SHSY5Y$p.value.Undiff1 < 0.001 & SHSY5Y$p.value.Undiff3 < 0.001) |
                                    (SHSY5Y$p.value.Undiff2 < 0.001 & SHSY5Y$p.value.Undiff3 < 0.001)) ),] 
  
  ### Add the total for C_Diff, T_Diff , C_Undiff, T_Undiff
  SHSY5Y_filter$T_Undiff <- 0
  SHSY5Y_filter$C_Undiff <- 0
  for (i in c(1:nrow(SHSY5Y_filter))) {
    SHSY5Y_filter$T_Undiff[i] <- SHSY5Y_filter$T_Undiff1[i] + SHSY5Y_filter$T_Undiff2[i] + SHSY5Y_filter$T_Undiff3[i]
    SHSY5Y_filter$C_Undiff[i] <- SHSY5Y_filter$C_Undiff1[i] + SHSY5Y_filter$C_Undiff2[i] + SHSY5Y_filter$C_Undiff3[i]
  }
  
  ### Add total_Diff & total_Undiff
  SHSY5Y_filter$total_Undiff <- 0
  SHSY5Y_filter$total_IVT <- 0
  for (i in c(1:nrow(SHSY5Y_filter))) {
    SHSY5Y_filter$total_Undiff[i] <- SHSY5Y_filter$T_Undiff[i] + SHSY5Y_filter$C_Undiff[i]
    SHSY5Y_filter$total_IVT[i] <- SHSY5Y_filter$T_IVT[i] + SHSY5Y_filter$C_IVT[i]
    
  }
  
  ### Add the mismatches mm.Diff & mm.Undiff
  SHSY5Y_filter$mm.Undiff <- 0
  for (i in c(1:nrow(SHSY5Y_filter))) {
    SHSY5Y_filter$mm.Undiff[i] <- (SHSY5Y_filter$C_Undiff[i] / (SHSY5Y_filter$C_Undiff[i] + SHSY5Y_filter$T_Undiff[i]))*100
    
  }
  
  ### Add the mismatch difference between Diff and IVT and Undiff and IVT
  SHSY5Y_filter$mm.UndiffMINUSmm.IVT <- 0
  for (i in c(1:nrow(SHSY5Y_filter))) {
    SHSY5Y_filter$mm.UndiffMINUSmm.IVT[i] <- SHSY5Y_filter$mm.Undiff[i]-SHSY5Y_filter$mm.IVT[i]
  }
  
  
  ### Filter to total_Undiff < 10
  SHSY5Y_filter <- SHSY5Y_filter[which(SHSY5Y_filter$total_Undiff > 10),]
  
  
  ### Substitute mmIVT==0 with panIVT
  SHSY5Y_filter$ACP<-paste0(SHSY5Y_filter$Annotation,SHSY5Y_filter$chr,SHSY5Y_filter$position)
  ACP<-SHSY5Y_filter$ACP[which(SHSY5Y_filter$total_IVT < 10)]
  ACP_pan<-panIVT_Undiff$ACP[which(panIVT_Undiff$ACP %in% ACP)]
  ACP_pan_mm<-panIVT_Undiff$mm.IVT[which(panIVT_Undiff$ACP %in% ACP)]
  pan<-panIVT_Undiff[which(panIVT_Undiff$ACP %in% ACP),]
  for (i in c(1:length(ACP_pan))) {
    if (ACP_pan[i] %in% SHSY5Y_filter$ACP){
      idx<-which(SHSY5Y_filter$ACP==ACP_pan[i])
      SHSY5Y_filter$total_IVT[idx]<-panIVT_Undiff$N_reads_IVT[which(panIVT_Undiff$ACP %in% ACP_pan[i])]
      SHSY5Y_filter$T_IVT[idx]<-panIVT_Undiff$T_IVT[which(panIVT_Undiff$ACP %in% ACP_pan[i])]
      SHSY5Y_filter$C_IVT[idx]<-panIVT_Undiff$C_IVT[which(panIVT_Undiff$ACP %in% ACP_pan[i])]
      SHSY5Y_filter$mm.IVT[idx]<-panIVT_Undiff$mm.IVT[which(panIVT_Undiff$ACP %in% ACP_pan[i])]
    }
  }
  
  SHSY5Y_filter <- SHSY5Y_filter[which(SHSY5Y_filter$total_IVT > 10),]
  SHSY5Y_filter <- SHSY5Y_filter[which((SHSY5Y_filter$T_IVT+SHSY5Y_filter$C_IVT)!=0),]
  SHSY5Y_filter <- SHSY5Y_filter[which(SHSY5Y_filter$mm.IVT<10),]
  
  
  ### Remove the ones with no reads in IVT
  SHSY5Y_filter <- SHSY5Y_filter[which(SHSY5Y_filter$total_IVT != 0),]
  
  ### Add the color column based on the kmer
  SHSY5Y_filter$color <- "black"
  for (i in c(1:nrow(SHSY5Y_filter))) {
    if (SHSY5Y_filter$kmer[i] == "TGTAG") {
      SHSY5Y_filter$color[i] <- "blue"
    }
    if (SHSY5Y_filter$kmer[i] %in%  c("GTTCA", "GTTCT", "GTTCC", "GTTCG")) {
      SHSY5Y_filter$color[i] <- "red"
    }
    else {}
  }
  
  ### Add the TPM of undiff
  SHSY5Y_filter$TPM_Undiff <- 0
  for (i in c(1:nrow(SHSY5Y_filter))) {
    row.on.Undiff.TPM <- which(TPM_Undiff$Gene.Name == SHSY5Y_filter$Annotation[i])
    if (length(row.on.Undiff.TPM == 1)) {
      SHSY5Y_filter$TPM_Undiff[i] <- TPM_Undiff$TPM[row.on.Undiff.TPM][[1]]
    }
    else{}
  }
  SHSY5Y_filter$TPM_Undiff_log <- log10(SHSY5Y_filter$TPM_Undiff)

}
write.csv(SHSY5Y_filter,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/SHSY5_Undiff_panIVTfiltered.csv",row.names=F)


}

################################################ start here after panIVT enrichment is done ##############################################################################
library(Cairo)
sig=1
Stu=T
totDiff=sum(c(2497734,2421424,1396031))
totUndiff=sum(c(3035434,3035434,873025))
totPB=sum(c(616787,1022734))

SH <- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Diff_filter_pan.csv", header = T) #diff
SH_filter<-SH
PB<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PB_filter_pan.csv",header=T) #pb
PB_filter<-PB
Undiff<- read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/SHSY5_Undiff_panIVTfiltered.csv", header=T) #undiff

#Think about this, is this correct?
#number of psi sites found in Undiff library normalized for the total number of reads in Undiff library
totPsiUndiff=(nrow(Undiff)/totUndiff)*100
#number of psi sites found in Pb library normalized for the total number of reads in Pb library
totPsiPB=(nrow(PB)/totPB)*100
#number of psi sites found in Diff library normalized for the total number of reads in Diff library
totPsiDiff=(nrow(SH)/totDiff)*100

################ Filtration: Diff, Undiff data
if(T){
  
  ### Calculate N using bootstrap SD information (YES)
  if (T) {
    ##### Read the snr deviation data
    SH_filter$mm.DiffminusUndiff = SH_filter$mm.Diff-SH_filter$mm.Undiff
    snrDiff<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/snrD_Diff.csv", header=T)
    names(snrDiff)[names(snrDiff) == "SNR"] <- "SNR_Diff"
    snrUndiff<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/snrD_Undiff.csv", header=T)
    names(snrUndiff)[names(snrUndiff) == "SNR"] <- "SNR_Undiff"
    
    selected_columns<-c("Annotation", "position","SNR_Diff")
    SH_filter <- merge(SH_filter, snrDiff[selected_columns], by = c("Annotation", "position"), all.x = TRUE)
    selected_columns<-c("Annotation", "position","SNR_Undiff")
    SH_filter <- merge(SH_filter, snrUndiff[selected_columns], by = c("Annotation", "position"), all.x = TRUE)
    
    SH_filter$N_value<-pmax(SH_filter$SNR_Undiff,SH_filter$SNR_Diff,na.rm = TRUE)

    
    #supplementary table 8
    if(F){
      PUS7_s8<-data.frame("<0"=length(which(SH_filter$N_value<0)),
                          "0-0.5"=length(which(SH_filter$N_value>=0 & SH_filter$N_value<0.5)),
                          "0.5-1"=length(which(SH_filter$N_value>=0.5 & SH_filter$N_value<1)),
                          "1-2"=length(which(SH_filter$N_value>=1 & SH_filter$N_value<2)),
                          ">2"=length(which(SH_filter$N_value>=2)))
      write.csv(PUS7_s8,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/Diff_8",row.names = F)
      
    }

  }
  
  #sites in which we believe from PUS7KD and TRUB1_KD
  TRUB1sites<-TRUB1KD_sites
  TRUB1sites$ACP<-paste0(TRUB1sites$Annotation,TRUB1sites$CP)
  SH_filter$ACP<-paste0(SH_filter$Annotation,SH_filter$chr,SH_filter$position)
  common<-intersect(SH_filter$ACP,TRUB1sites$ACP)
  TRUB1sites<-SH_filter[which(SH_filter$ACP %in% common),]
  SH_filter$color[which(SH_filter$ACP %in% common)]<-"orange"
  
  PUS7sites<-PUS7KD_sites
  PUS7sites$ACP<-paste0(PUS7sites$Annotation,PUS7sites$CP)
  common<-intersect(SH_filter$ACP, PUS7sites$ACP)
  PUS7sites<-SH_filter[which(SH_filter$ACP %in% common),]
  SH_filter$color[which(SH_filter$ACP %in% common)]<-"blue2"
  
  #calculating sigma2 after the new color assignment 
  SH_filter_sigma2<-SH_filter[which(SH_filter$N_value>sig),] #sigma >6
  TRUB1sites_sigma2<-TRUB1sites[which(TRUB1sites$N_value>sig),]
  PUS7sites_sigma2<-PUS7sites[which(PUS7sites$N_value>sig),]
  

  
  #Supplemtary Table 4
  if(F){
    
    SH_tab<-data.frame("Annotation"=SH_filter_sigma2$Annotation,"chr"=SH_filter_sigma2$chr,"position"=SH_filter_sigma2$position,
                       "strand"=SH_filter_sigma2$strand,"kmer"=SH_filter_sigma2$kmer,
                       "T_Differentiated"=SH_filter_sigma2$T_Diff,"C_Differentiated"=SH_filter_sigma2$C_Diff,
                       "T_Untreated"=SH_filter_sigma2$T_Undiff,"C_Untreated"=SH_filter_sigma2$C_Undiff,
                       "T_IVT"=SH_filter_sigma2$T_IVT,"C_IVT"=SH_filter_sigma2$C_IVT,
                       "mm.IVT"=SH_filter_sigma2$mm.IVT,
                       "mm.Differentiated_minus_Untreated"=SH_filter_sigma2$mm.DiffminusUndiff,
                       "SNR"=SH_filter_sigma2$N_value)
    write.csv(SH_tab,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/SupplementaryTable4.csv",row.names = F)
    
  }
  
  #Supplementary Table 5
  if(F){
    SH_filter_validated<-SH_filter_sigma2[which(SH_filter_sigma2$color=="orange" | SH_filter_sigma2$color=="blue2" ),]
    SH_tab<-data.frame("Annotation"=SH_filter_validated$Annotation,"chr"=SH_filter_validated$chr,"position"=SH_filter_validated$position,
                       "strand"=SH_filter_validated$strand,"kmer"=SH_filter_validated$kmer,
                       "T_Differentiated"=SH_filter_validated$T_Diff,"C_Differentiated"=SH_filter_validated$C_Diff,
                       "T_Untreated"=SH_filter_validated$T_Undiff,"C_Untreated"=SH_filter_validated$C_Undiff,
                       "T_IVT"=SH_filter_validated$T_IVT,"C_IVT"=SH_filter_validated$C_IVT,
                       "mm.IVT"=SH_filter_validated$mm.IVT,
                       "mm.Differentiated_minus_Untreated"=SH_filter_validated$mm.DiffminusUndiff,
                       "SNR"=SH_filter_validated$N_value)
    write.csv(SH_tab,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/SupplementaryTable5.csv",row.names = F)
    
  }
}

#SH Diff plots
if(T){
################################## Plot: DotPlot U-to-C mismatch Diff-Undiff VS Undiff
if(T){  
  min.x = -60
  max.x = 60
  min.y = -4
  max.y = 3
  line.thickness = 1
  label.size = 1
  text.size = 0.8
  point.size = 0.1
  tck.length = 0.01
  tick.thickness = 1
  transparency = 1
  file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/DiffvsUndiff_little_bootstrap.pdf"
  Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
  
  plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = '',xlab='',ann = F,axes=T,type="l", 
       xaxt='n', yaxt='n', 
       bty='n', bg='transparent')
  
  
  points(x=SH_filter$mm.DiffminusUndiff[which(SH_filter$color == "black")], y=log10(SH_filter$N_value[which(SH_filter$color=="black")]), pch=22, col="grey", bg=adjustcolor("grey",alpha=0.1), cex = log2(SH_filter$total_Diff[which(SH_filter$color == "black")])*point.size)
  points(x=SH_filter$mm.DiffminusUndiff[which(SH_filter$color == "blue")], y=log10(SH_filter$N_value[which(SH_filter$color=="blue")]), pch=22, col="grey", bg=adjustcolor("grey",alpha=0.3), cex = log2(SH_filter$total_Diff[which(SH_filter$color == "blue")])*point.size)
  points(x=SH_filter$mm.DiffminusUndiff[which(SH_filter$color == "red")], y=log10(SH_filter$N_value[which(SH_filter$color=="red")]), pch=22, col="grey", bg=adjustcolor("grey",alpha=0.5), cex = log2(SH_filter$total_Diff[which(SH_filter$color == "red")])*point.size)
  points(x=SH_filter$mm.DiffminusUndiff[which(SH_filter$color == "orange")], y=log10(SH_filter$N_value[which(SH_filter$color=="orange")]), pch=22, col="orange", bg=adjustcolor("orange",alpha=0.8), cex = log2(SH_filter$total_Diff[which(SH_filter$color == "orange")])*point.size)
  points(x=SH_filter$mm.DiffminusUndiff[which(SH_filter$color == "blue2")], y=log10(SH_filter$N_value[which(SH_filter$color=="blue2")]), pch=22, col="blue2", bg=adjustcolor("blue2",alpha=0.3), cex = log2(SH_filter$total_Diff[which(SH_filter$color == "blue2")])*point.size)
  
  segments(x0 = 0, y0 = -16, x1 = 0, y1 = 10, lty=2, col = "grey3")
  segments(x0 = -60, y0 = log10(sig), x1 = 60, y1 = log10(sig), lty=2, col = "grey3")
  
  
  SH_sites_odd <- TRUB1sites
  SH_sites_even <- PUS7sites
  
  text(x=SH_sites_odd$mm.DiffminusUndiff, y=log10(SH_sites_odd$N_value),SH_sites_odd$Annotation, col='grey3', cex=0.3)
  text(x=SH_sites_even$mm.DiffminusUndiff-2, y=log10(SH_sites_even$N_value)-0.1,SH_sites_even$Annotation, col='grey3', cex=0.3)
  
  
  SH_filter$N_value[which(SH_filter$N_value==0)]<-NA
  
  axis(side = 2,at = seq(min.y, max.y, by=1) ,labels = seq(min.y,max.y , by=1) ,lwd.ticks = tick.thickness ,cex.axis= label.size,
       tck=-tck.length,las=1,lwd=line.thickness, line=0)
  #axis(side = 2,at = seq(0,max(log2(SH_filter_sigma2$SH.standard.dev.....)), by=1),labels = seq(0, 10, by = 1),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
  axis(side = 1,at = seq(-60,60,by=20) ,labels = seq(-60,60,by=20)  ,lwd.ticks = tick.thickness ,cex.axis= label.size,
       tck=-tck.length,las=1,lwd=line.thickness, line=0)
  #axis(side = 1,at = c(-6:6)*10,labels = c(-6:6)*10,lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
  legend("topright", c("TRUB1 sites","PUS7 sites", "Other"),     
         cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
         pch=c(16,16,16),lty=c(0,0,0), 
         lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
  mtext(text = expression(paste(log[10], SNR)),side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
  mtext(text = "Differentiated - Untreated (U-C MM%)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
  alpha=0
  dev.off() 
}
################################## Plot: sigma2
if(T){
  min.x = 1
  max.x = round(log10(max(SH_filter_sigma2$total_Diff)))+0.5
  min.y = -60
  max.y = 60
  line.thickness = 1
  label.size = 1
  text.size = 0.8
  point.size = 0.2
  tck.length = 0.01
  tick.thickness = 1
  transparency = 1
  file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/DiffUndiff_log2_bootstrap_text.pdf"
  Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
  
  plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = '',xlab='',ann = F,axes=T,type="l", 
       xaxt='n', yaxt='n', 
       bty='n', bg='transparent')
  
  
  points(x=log10(SH_filter_sigma2$total_Diff[which(SH_filter_sigma2$color == "black")]), y=-SH_filter_sigma2$mm.DiffminusUndiff[which(SH_filter_sigma2$color=="black")], pch=22, col="grey", bg=adjustcolor("grey",alpha=0.05), cex = log2(SH_filter_sigma2$N_value[which(SH_filter_sigma2$color == "black")])*point.size)
  points(x=log10(SH_filter_sigma2$total_Diff[which(SH_filter_sigma2$color == "blue")]), y=-SH_filter_sigma2$mm.DiffminusUndiff[which(SH_filter_sigma2$color=="blue")], pch=22, col="grey", bg=adjustcolor("grey",alpha=0.3), cex = log2(SH_filter_sigma2$N_value[which(SH_filter_sigma2$color == "blue")])*point.size)
  points(x=log10(SH_filter_sigma2$total_Diff[which(SH_filter_sigma2$color == "red")]), y=-SH_filter_sigma2$mm.DiffminusUndiff[which(SH_filter_sigma2$color=="red")], pch=22, col="grey", bg=adjustcolor("grey",alpha=0.5), cex = log2(SH_filter_sigma2$N_value[which(SH_filter_sigma2$color == "red")])*point.size)
  points(x=log10(SH_filter_sigma2$total_Diff[which(SH_filter_sigma2$color == "orange")]), y=-SH_filter_sigma2$mm.DiffminusUndiff[which(SH_filter_sigma2$color=="orange")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.8), cex = log2(SH_filter_sigma2$N_value[which(SH_filter_sigma2$color == "orange")])*point.size)
  points(x=log10(SH_filter_sigma2$total_Diff[which(SH_filter_sigma2$color == "blue2")]), y=-SH_filter_sigma2$mm.DiffminusUndiff[which(SH_filter_sigma2$color=="blue2")], pch=22, col="blue2", bg=adjustcolor("blue2",alpha=0.3), cex = log2(SH_filter_sigma2$N_value[which(SH_filter_sigma2$color == "blue2")])*point.size)
  
  SH_sites_odd_sigma2<-TRUB1sites[which(TRUB1sites$N_value>sig),]
  SH_sites_even_sigma2<-PUS7sites[which(PUS7sites$N_value>sig),]
  
  text(x=log10(SH_sites_odd_sigma2$total_Diff), y=-SH_sites_odd_sigma2$mm.DiffminusUndiff+0.5,SH_sites_odd_sigma2$Annotation, col='grey3', cex=0.3)
  text(x=log10(SH_sites_even_sigma2$total_Diff), y=-SH_sites_even_sigma2$mm.DiffminusUndiff,SH_sites_even_sigma2$Annotation, col='grey3', cex=0.3)

  
  #for plot KD-SCR vs SCR
  axis(side = 2, at = seq(min.y, max.y, by = 10), labels = seq(min.y, max.y, by = 10),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  axis(side = 1, at = seq(1,max.x, by = 0.5), labels = seq(1, max.x, by = 0.5),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
       las = 1, lwd = line.thickness, line = 0)
  legend("topright", c("TRUB1 sites","PUS7 sites", "Other"),     
         cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
         pch=c(16,16,16),lty=c(0,0,0), 
         lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
  mtext(text = "Untreated - Differetiated (U-to-C MM%) ",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
  mtext(text = "log10(Reads)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
  
  
  dev.off()
}  

################################## Barplot for 9 sites
if (T) {
  
  library(stringr)
  #### List of the TRUB1 sites with sigma>2 and with a delta>5 or <-5
  TRUB1sites5<-TRUB1sites_sigma2[which(abs(TRUB1sites_sigma2$mm.DiffminusUndiff)>5),]
  PUS7sites5<-PUS7sites_sigma2[which(abs(PUS7sites_sigma2$mm.DiffminusUndiff)>5),]
  pus_enzymes <-c(TRUB1sites5$Annotation,PUS7sites5$Annotation)

  pus_df <- data.frame("gene" = pus_enzymes, "Diff1_TPM" = 0,"Diff2_TPM" = 0,"Diff3_TPM" = 0,
                       "Undiff1_TPM" = 0,"Undiff2_TPM" = 0, "Undiff3_TPM" = 0)
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff1_TPM[i] <- TPM_Diff1$TPM[which(TPM_Diff1$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff2_TPM[i] <- TPM_Diff2$TPM[which(TPM_Diff2$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff3_TPM[i] <- TPM_Diff3$TPM[which(TPM_Diff3$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff1_TPM[i] <- TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff2_TPM[i] <- TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff3_TPM[i] <- TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == pus_df$gene[i])][1]
  }
  
  #### Add the standard deviation to the pus_df
  pus_df$Diff_SD <- 0
  pus_df$Undiff_SD <- 0
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff_SD[i] <- sd(c(pus_df$Diff1_TPM[i], pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
    pus_df$Undiff_SD[i] <- sd(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
  }
  #### Add the mean to the pus_df 
  pus_df$Diff_mean <- 0
  pus_df$Undiff_mean <- 0
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff_mean[i] <- mean(c(pus_df$Diff1_TPM[i], pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
    pus_df$Undiff_mean[i] <- mean(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
  }

  #pvalues, is Diff significantly different compared to undiff?
  for (i in c(1:nrow(pus_df))) {
    tDiff<-t.test(c(pus_df$Diff1_TPM[i],pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]),c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
    pus_df$pvalue[i] <- tDiff$p.value
  }

  #### Plot with segment on y axis because the rpll22 is way higher than other genes
  break_plot=F
  if(break_plot==T){ 
    library(plotrix)
    
    min.x = 0
    max.x = 12
    min.y = -1
    realmax.y=round(max(pus_df[,2:length(pus_df)]))+20
    max.y = realmax.y-100
    line.thickness = 1
    label.size = 1
    text.size = 1.4
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/DiffvsUndiff_barplot_pvals.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
    #box(lwd = line.thickness)
    #### Add the lines
  
    off=100
    
    #for Undiff
    for (i in c(1:nrow(pus_df))) {
      if (pus_df$gene[i]=="RPL22"){
        lines(x=c(i-.1,i-.1),y=c(0,pus_df$Undiff_mean[i]-off),
              lwd=5 ,lend=1, lty=1, col = "dodgerblue4")
      }else{
        lines(x=c(i-.1,i-.1),y=c(0,pus_df$Undiff_mean[i]),
              lwd=5 ,lend=1, lty=1, col = "dodgerblue4")
      }
    }
    
    
    #for Diff
    for (i in c(1:nrow(pus_df))) {
      
      if (pus_df$gene[i]=="RPL22"){
        lines(x=c(i+.1,i+.1),y=c(0,pus_df$Diff_mean[i]-off),
              lwd=5 ,lend=1, lty=1, col = "dodgerblue3")
      }else{
        lines(x=c(i+.1,i+.1),y=c(0,pus_df$Diff_mean[i]),
              lwd=5 ,lend=1, lty=1, col = "dodgerblue3")
      }
    }
    
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      if (pus_df$gene[i]=="RPL22"){
        lines(x=c(i-.1,i-.1),
              y=c(pus_df$Undiff_mean[i]-off-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]-off+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }else{
        lines(x=c(i-.1,i-.1),
              y=c(pus_df$Undiff_mean[i]-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }
    }
    
    
    for (i in c(1:nrow(pus_df))) {
      if (pus_df$gene[i]=="RPL22"){
        lines(x=c(i+.1,i+.1),
              y=c(pus_df$Diff_mean[i]-off-pus_df$Diff_SD[i],pus_df$Diff_mean[i]-off+pus_df$Diff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }else{
        lines(x=c(i+.1,i+.1),
              y=c(pus_df$Diff_mean[i]-pus_df$Diff_SD[i],pus_df$Diff_mean[i]+pus_df$Diff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }
    }

    
    
    #plot pvalue asterisk
    #### Add Asterisks for p-values
    # *   < 0.05
    # **  < 0.01
    # *** < 0.001
    
    for (i in c(1:nrow(pus_df))){
     if (pus_df$pvalue[i]<0.001) {
        text(i,max.y-5, labels = "***", col = "dodgerblue3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(realmax.y-8,realmax.y-8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalue[i]<0.01) {
        text(i,max.y-5, labels = "**", col = "dodgerblue3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(realmax.y-8,realmax.y-8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalue[i]<0.05) {
        text(i,max.y-5, labels = "*", col = "dodgerblue3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(realmax.y-8,realmax.y-8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } 
    }
    
    

    ngenes=nrow(pus_df)
    
    axis(side = 2,at = c(-1,50,100,max.y) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    axis(side = 2,at =  c(0,50,100,max.y),labels = c(0,50,100,realmax.y),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
    axis.break(2, 100, style = "gap") 
    axis(side = 1,at = c(0:max.x) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
    axis(side = 1,at = c(1:ngenes),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
    mtext(text = "TPM",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    legend("topleft", c("Untreated","Differentiated"),
           cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
           pch=c(15,15,15),lty=c(0,0,0),
           lwd=c(2.5,2.5,2.5),col=c("dodgerblue4", "dodgerblue3"), bty = "n")
    
    
    dev.off() 
 
  }
  
  #### Plot log(TPM)
  if(T){ 
    library(plotrix)
    
    min.x = 0
    max.x = nrow(pus_df)+1
    min.y = 0
    max.y = round(max(log10(pus_df[,2:length(pus_df)])))+2
    line.thickness = 1
    label.size = 1
    text.size = 1.4
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/DiffvsUndiff_barplot_pvals_log.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 8, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
    #box(lwd = line.thickness)
    #### Add the lines
    
    
    #for Undiff
    for (i in c(1:nrow(pus_df))) {
        lines(x=c(i-.1,i-.1),y=c(0.1,log10(pus_df$Undiff_mean[i])),lwd=3 ,lend=1, lty=1, col = "dodgerblue4")
        pts=log10(c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
        points(x=rep(i-.1,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    
    #for Diff
    for (i in c(1:nrow(pus_df))) {
        lines(x=c(i+.1,i+.1),y=c(0.1,log10(pus_df$Diff_mean[i])),lwd=3 ,lend=1, lty=1, col = "dodgerblue3")
        pts=log10(c(pus_df$Diff1_TPM[i],pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
        points(x=rep(i+.1,length(pts)), y=pts, pch=20,col="black",cex=0.01)

    }
    
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
        lines(x=c(i-.1,i-.1),
              y=c(log10(pus_df$Undiff_mean[i])-log10(pus_df$Undiff_SD[i]),log10(pus_df$Undiff_mean[i])+log10(pus_df$Undiff_SD[i])),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    
    for (i in c(1:nrow(pus_df))) {

        lines(x=c(i+.1,i+.1),
              y=c(log10(pus_df$Diff_mean[i])-log10(pus_df$Diff_SD[i]),log10(pus_df$Diff_mean[i])+log10(pus_df$Diff_SD[i])),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    
    
    #plot pvalue asterisk
    #### Add Asterisks for p-values
    # *   < 0.05
    # **  < 0.01
    # *** < 0.001
    for (i in c(1:nrow(pus_df))){
      if (pus_df$pvalue[i]<0.001) {
        text(i,max.y-1.8, labels = "***", col = "dodgerblue3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-2.,max.y-2),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalue[i]<0.01) {
        text(i,max.y-0.5, labels = "**", col = "dodgerblue3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-1,max.y-1),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalue[i]<0.05) {
        text(i,max.y-0.5, labels = "*", col = "dodgerblue3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-1,max.y-1),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } 
    }
    
    
    
    ngenes=nrow(pus_df)
    
    axis(side = 2,at = c(0:max.y) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    axis(side = 2,at =  c(0:max.y),labels = c(0:max.y),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
    axis(side = 1,at = c(0:max.x) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
    axis(side = 1,at = c(1:ngenes),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
    mtext(text = "log10TPM",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    legend("topleft", c("Untreated","Differentiated"),
           cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
           pch=c(15,15,15),lty=c(0,0,0),
           lwd=c(2.5,2.5,2.5),col=c("dodgerblue4", "dodgerblue3"), bty = "n")
    
    
    dev.off() 
    
  }
}

#####write table fig2 for Diff Undiff
if(T){
sites_validated<-SH_filter_sigma2[which(SH_filter_sigma2$color=="orange" | SH_filter_sigma2$color=="blue2"),]
sites_validated<-sites_validated[which(abs(sites_validated$mm.DiffMINUSmm.Undiff)>5),]  
SitesTable<-data.frame(Annotation=sites_validated$Annotation, chr=sites_validated$chr,
                         position=sites_validated$position, kmer=sites_validated$kmer,
                         mm.UntreatedminusIVT=sites_validated$mm.UndiffMINUSmm.IVT,
                         mm.DiffminusIVT=sites_validated$mm.DiffMINUSmm.IVT)
  
  
  
write.csv(SitesTable,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Tablefig2.csv",row.names = F)
}

############################### seqLogos on the 5 different intervals
if(T){
####int1: >20
####int2: >=5 & <=20
####int3: >-5 & <5
####int4: <=-5 & >=-20
####int5: <-20
SH_filter_sigma2$mm.UndiffminusDiff=-SH_filter_sigma2$mm.DiffminusUndiff
int1<-SH_filter_sigma2[which(SH_filter_sigma2$mm.UndiffminusDiff >   20 ),]
int2<-SH_filter_sigma2[which(SH_filter_sigma2$mm.UndiffminusDiff >=  5 & SH_filter_sigma2$mm.UndiffminusDiff <=  20 ),]
int3<-SH_filter_sigma2[which(SH_filter_sigma2$mm.UndiffminusDiff >  -5 & SH_filter_sigma2$mm.UndiffminusDiff <   5  ),]
int4<-SH_filter_sigma2[which(SH_filter_sigma2$mm.UndiffminusDiff <= -5 & SH_filter_sigma2$mm.UndiffminusDiff >= -20 ),]
int5<-SH_filter_sigma2[which(SH_filter_sigma2$mm.UndiffminusDiff <  -20),]

if (T){
  
  # kmers<-sapply(int1$kmer,function(x) str_replace_all(x,"T","U"))
  # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Logo>20.pdf"
  
  # kmers<-sapply(int2$kmer,function(x) str_replace_all(x,"T","U"))
  # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Logo20to5.pdf"
   
  # kmers<-sapply(int3$kmer,function(x) str_replace_all(x,"T","U"))
  # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Logo5to-5.pdf"

  # kmers<-sapply(int4$kmer,function(x) str_replace_all(x,"T","U"))
  # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Logo-5to-20.pdf"
  
  # kmers<-sapply(int5$kmer,function(x) str_replace_all(x,"T","U"))
  # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Logo<-20.pdf"
  
  my_list<-list("new motifs" = kmers)

  min.x = 0
  max.x = 200
  min.y = 0
  max.y = 20
  line.thickness = 1
  label.size = 1
  text.size = 0.8
  point.size = 1
  tck.length = 0.01
  tick.thickness = 1
  transparency = 0.2
  Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
  plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = NA,xlab=NA,ann = F,axes=F)
  ggseqlogo::ggseqlogo(my_list,ncol=1,method="prob")
  #same as
  #ggplot() + geom_logo(seqs_dna) + theme_logo() + facet_wrap(~seq_group, ncol=4, scales='free_x') 
  
  dev.off()
  
}
}
  
########################## Kinetic model Diff Undiff
if(T){
#([mRNAu][mRNApus])a -> ([mRNAu][mRNApus])b
#a=Diff
#b=Undiff
#plot BminusA

#trub1 and pus7 sites validated with KD libraries
confirmed_sites_names<-c(TRUB1sites_sigma2$ACP,PUS7sites_sigma2$ACP)

#plot kinetic


  scale=1000000
  
  diff=T
  #run the pus_df from the Barplot of the PUS enzymes that is above 
  
  #Diff
  if(diff==T){
    kinetic_sites<-SH_filter[which(SH_filter$ACP %in% confirmed_sites_names),]
    
    Diff_TPM_TRUB1 <- mean(c(TPM_Diff1$TPM[which(TPM_Diff1$Gene.Name == "TRUB1")],TPM_Diff2$TPM[which(TPM_Diff2$Gene.Name == "TRUB1")],TPM_Diff3$TPM[which(TPM_Diff3$Gene.Name == "TRUB1")]))
    Undiff_TPM_TRUB1 <-mean(c(TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == "TRUB1")],TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == "TRUB1")],TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == "TRUB1")]))
    Diff_TPM_PUS7 <- mean(c(TPM_Diff1$TPM[which(TPM_Diff1$Gene.Name == "PUS7")],TPM_Diff2$TPM[which(TPM_Diff2$Gene.Name == "PUS7")],TPM_Diff3$TPM[which(TPM_Diff3$Gene.Name == "PUS7")]))
    Undiff_TPM_PUS7 <-mean(c(TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == "PUS7")],TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == "PUS7")],TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == "PUS7")]))
    
    M=totDiff #normalization
    kinetic_sites$ConditionA<-0
    kinetic_sites$ConditionA[which(kinetic_sites$color=="orange")] <- (kinetic_sites$total_Diff[which(kinetic_sites$color=="orange")]/M)*scale*Diff_TPM_TRUB1
    #kinetic_sites$ConditionA[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_TRUB1_KD[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$TRUB1_mean[pus_df$gene=="PUS7"]
    kinetic_sites$ConditionA[which(kinetic_sites$color=="blue2")] <- (kinetic_sites$total_Diff[which(kinetic_sites$color=="blue2")]/M)*scale*Diff_TPM_PUS7
    kinetic_sites$yConditionA <- kinetic_sites$C_Diff/M*scale 
    M=totUndiff #normalization
    kinetic_sites$ConditionB<-0
    kinetic_sites$ConditionB[which(kinetic_sites$color=="orange")] <- (kinetic_sites$total_Undiff[which(kinetic_sites$color=="orange")]/M)*scale*Undiff_TPM_TRUB1
    #kinetic_sites$ConditionB[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_SCR[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$Scram_mean[pus_df$gene=="PUS7"] 
    kinetic_sites$ConditionB[which(kinetic_sites$color=="blue2")] <- (kinetic_sites$total_Undiff[which(kinetic_sites$color=="blue2")]/M)*scale*Undiff_TPM_PUS7 #
    kinetic_sites$yConditionB <- kinetic_sites$C_Undiff/M*scale 
    
    kinetic_sites$BminusA  <- kinetic_sites$ConditionB - kinetic_sites$ConditionA
    kinetic_sites$yBminusA <- kinetic_sites$yConditionB - kinetic_sites$yConditionA
  }
  
  #plot conditionUndiff-conditionDiff
  if(T){
    min.x = -300 #round(min(kinetic_sites$BminusA))-49
    max.x = 300#round(max(kinetic_sites$BminusA))+50
    min.y = -40
    max.y = 40
    line.thickness = 1
    label.size = 1
    text.size = 0.8
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 1
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Kinetic_model_Diff_annotated.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    
    plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = '',xlab='',ann = F,axes=T,type="l", 
         xaxt='n', yaxt='n', 
         bty='n', bg='transparent')
    
    
    points(x=kinetic_sites$BminusA[which(kinetic_sites$color == "orange")], y=kinetic_sites$yBminusA[which(kinetic_sites$color=="orange")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.3), cex = point.size)
    
    points(x=kinetic_sites$BminusA[which(kinetic_sites$color == "blue2")], y=kinetic_sites$yBminusA[which(kinetic_sites$color=="blue2")], pch=22, col="blue", bg=adjustcolor("blue",alpha=0.6), cex=point.size)
    
    # Correlation TRUB1 sites only 
    x <- kinetic_sites$BminusA[which(kinetic_sites$color == "orange")]
    y <- kinetic_sites$yBminusA[which(kinetic_sites$color=="orange")]
    
    # Calculate correlation coefficient
    correlation_coefficient <- cor(x, y)
    # Calculate R-squared value
    r_squared <- correlation_coefficient^2
    
    # Print the correlation coefficient and R-squared value
    print(paste("Correlation Coefficient trub1:", correlation_coefficient))
    print(paste("R-squared Value:", r_squared))
    
    # Add correlation line
    abline(lm(y ~ x), col = "orange")
    
    
    x<-kinetic_sites$BminusA[which(kinetic_sites$color == "blue2")]
    y<-kinetic_sites$yBminusA[which(kinetic_sites$color=="blue2")]
    
    # Calculate correlation coefficient
    correlation_coefficient <- cor(x, y)
    # Calculate R-squared value
    r_squared <- correlation_coefficient^2
    
    # Add correlation line
    abline(lm(y ~ x), col = "blue")
    
    # Print the correlation coefficient and R-squared value
    print(paste("Correlation Coefficient pus7:", correlation_coefficient))
    print(paste("R-squared Value:", r_squared))
    
    abline(v=0, col = "black")
    
    
    text(x=kinetic_sites$BminusA, y=kinetic_sites$yBminusA,kinetic_sites$Annotation, col='grey3', cex=0.4)
    
    
    #for plot KD-SCR vs SCR
    axis(side = 2, at = seq(min.y, max.y, by = 10), labels = seq(min.y, max.y, by = 10),
         lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
         las = 1, lwd = line.thickness, line = 0)
    axis(side = 1, at = seq(min.x,max.x, by = round(max.x/10)), labels = seq(min.x, max.x, by = round(max.x/10)),
         lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
         las = 1, lwd = line.thickness, line = 0)
    legend("topright", c("TRUB1 validated sites","PUS7 validated sites"),     
           cex= text.size-0.1, pt.cex=1, y.intersp=1,x.intersp=0.01,
           pch=c(16,16,16),lty=c(0,0,0), 
           lwd=c(2.5,2.5,2.5),col=c("orange","blue"), bty = "n")
    mtext(text = expression(" ψUndiff - ψDiff per million "),side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    mtext(text = "([mRNA U][mRNA PUS])Undiff - ([mRNA U][mRNA PUS])Diff",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
    ##ggsave(filename =file.name,width = 6,height = 6.6)
    alpha=0
    dev.off()
    
    
  }
  
  
}
  
########################## Barplot Diff Undiff PUS ezymes
if (T) {
  #### List of the PUS enzyme genes:
  pus_enzymes <- c("PUS1", "PUSL1", "PUS3", "TRUB1", "TRUB2",
                   "DKC1", "PUS7", "PUS7L", "RPUSD1", "RPUSD2",
                   "PUS10", "PRUSD3", "PRUSD4")
  #### Make the data frame of TPM of the pus_enzymes 
  pus_df <- data.frame("gene" = pus_enzymes, "Diff1_TPM" = 0,"Diff2_TPM" = 0,"Diff3_TPM" = 0,
                       "Undiff1_TPM" = 0,"Undiff2_TPM" = 0, "Undiff3_TPM" = 0)
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff1_TPM[i] <- TPM_Diff1$TPM[which(TPM_Diff1$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff2_TPM[i] <- TPM_Diff2$TPM[which(TPM_Diff2$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff3_TPM[i] <- TPM_Diff3$TPM[which(TPM_Diff3$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff1_TPM[i] <- TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff2_TPM[i] <- TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff3_TPM[i] <- TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == pus_df$gene[i])][1]
  }
  
  #### Add the standard deviation to the pus_df
  pus_df$Diff_SD <- 0
  pus_df$Undiff_SD <- 0
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff_SD[i] <- sd(c(pus_df$Diff1_TPM[i], pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
    pus_df$Undiff_SD[i] <- sd(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
  }
  #### Add the mean to the pus_df 
  pus_df$Diff_mean <- 0
  pus_df$Undiff_mean <- 0
  for (i in c(1:nrow(pus_df))) {
    pus_df$Diff_mean[i] <- mean(c(pus_df$Diff1_TPM[i], pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
    pus_df$Undiff_mean[i] <- mean(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
  }
  
  #### Remove the unwanted rows
  pus_df = pus_df[-c(12,13),]
  
  #pvalues, is Diff significantly different compared to undiff?
  for (i in c(1:nrow(pus_df))) {
    tDiff<-t.test(c(pus_df$Diff1_TPM[i],pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]),c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
    pus_df$pvalueDiff[i] <- tDiff$p.value
  }
  
  

  #### Plot 
  if(T){ 
    min.x = 0
    max.x = 15
    min.y = 0
    max.y = 50
    line.thickness = 1
    label.size = 1
    text.size = 1.4
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/DiffvsUndiff_barplot_PUSenz_pvals_log.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
    #box(lwd = line.thickness)
    #### Add the lines
    
    #for Undiff
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.1,i-.1),y=c(0,pus_df$Undiff_mean[i]),lwd=5 ,lend=1, lty=1, col = "dodgerblue4")
      pts=c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i])
      points(x=rep(i-.1,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    #for Diff
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i+.1,i+.1),y=c(0,pus_df$Diff_mean[i]),lwd=5 ,lend=1, lty=1, col = "dodgerblue3")
      pts=c(pus_df$Diff1_TPM[i],pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i])
      points(x=rep(i+.1,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i+.1,i+.1),
            y=c(pus_df$Diff_mean[i]-pus_df$Diff_SD[i],pus_df$Diff_mean[i]+pus_df$Diff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.1,i-.1),
            y=c(pus_df$Undiff_mean[i]-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }

    for (i in c(1:nrow(pus_df))) {
      if (pus_df$pvalueDiff[i]<0.001) {
        text(i+.1,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "***", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalueDiff[i]<0.01) {
        text(i+.1,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "**", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalueDiff[i]<0.05) {
        text(i+.1,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "*", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } 
      
    }  
    
    axis(side = 2,at = c(0:5)*10 ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    axis(side = 2,at = c(0:5)*10,labels = c(0:5)*10,lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
    axis(side = 1,at = c(0:12) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
    axis(side = 1,at = c(1:11),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
    legend("topright", c("Undifferentiated","Differentiated"),
         cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
         pch=c(16,16),lty=c(0,0),
         lwd=c(2.5,2.5,2.5),col=c("dodgerblue4", "dodgerblue3"), bty = "n")
    mtext(text = "TPM",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    dev.off() 
  }
  
  
  }
  
}



### PB and Undiff data
### PB vs Undiff
################# Filtration: PB, Undiff data
if (T) {
  ### Calculate N using bootstrap SD information (NO)
  if (T) {
    ##### Read the snr deviation data
    PB_filter$mm.PBminusUndiff = PB_filter$mm.PB-PB_filter$mm.Undiff

   
      snrPB<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/snrP_Pb.csv", header=T)
      names(snrPB)[names(snrPB) == "SNR"] <- "SNR_PB"
      snrUndiff<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Stu/snrP_Undiff.csv", header=T)
      names(snrUndiff)[names(snrUndiff) == "SNR"] <- "SNR_Undiff"
      
      selected_columns<-c("Annotation", "position","SNR_PB")
      PB_filter <- merge(PB_filter, snrPB[selected_columns], by = c("Annotation", "position"), all.x = TRUE)
      selected_columns<-c("Annotation", "position","SNR_Undiff")
      PB_filter <- merge(PB_filter, snrUndiff[selected_columns], by = c("Annotation", "position"), all.x = TRUE)
      
      PB_filter$N_value<-pmax(PB_filter$SNR_PB,PB_filter$SNR_Undiff,na.rm = TRUE)

    #Supplementary Table8
    if(F){
      PUS7_s8<-data.frame("<0"=length(which(PB_filter$N_value<0)),
                          "0-0.5"=length(which(PB_filter$N_value>=0 & PB_filter$N_value<0.5)),
                          "0.5-1"=length(which(PB_filter$N_value>=0.5 & PB_filter$N_value<1)),
                          "1-2"=length(which(PB_filter$N_value>=1 & PB_filter$N_value<2)),
                          ">2"=length(which(PB_filter$N_value>=2)))
      write.csv(PUS7_s8,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/Pb_8",row.names = F)
      
    }

  }

  #sites in which we believe from PUS7KD and TRUB1_KD
  TRUB1sites<-TRUB1KD_sites
  TRUB1sites$ACP<-paste0(TRUB1sites$Annotation,TRUB1sites$CP)
  PB_filter$ACP<-paste0(PB_filter$Annotation,PB_filter$chr,PB_filter$position)
  common<-intersect(PB_filter$ACP,TRUB1sites$ACP)
  TRUB1sites<-PB_filter[which(PB_filter$ACP %in% common),]
  PB_filter$color[which(PB_filter$ACP %in% common)]<-"orange"
  
  PUS7sites<-PUS7KD_sites
  PUS7sites$ACP<-paste0(PUS7sites$Annotation,PUS7sites$CP)
  common<-intersect(PB_filter$ACP, PUS7sites$ACP)
  PUS7sites<-PB_filter[which(PB_filter$ACP %in% common),]
  PB_filter$color[which(PB_filter$ACP %in% common)]<-"blue2"

  PB_filter_sigma2<-PB_filter[which(PB_filter$N_value>sig),] #sigma >6
  TRUB1sites_sigma2<-TRUB1sites[which(TRUB1sites$N_value>sig),]
  PUS7sites_sigma2<-PUS7sites[which(PUS7sites$N_value>sig),]
  
  

  
  
  #Supplemtary Table 6
  if(F){
    
    PB_tab<-data.frame("Annotation"=PB_filter_sigma2$Annotation,"chr"=PB_filter_sigma2$chr,"position"=PB_filter_sigma2$position,
                       "strand"=PB_filter_sigma2$strand,"kmer"=PB_filter_sigma2$kmer,
                       "T_PB"=PB_filter_sigma2$T_PB,"C_PB"=PB_filter_sigma2$T_PB,
                       "T_Untreated"=PB_filter_sigma2$T_Undiff,"C_Untreated"=PB_filter_sigma2$C_Undiff,
                       "T_IVT"=PB_filter_sigma2$T_IVT,"C_IVT"=PB_filter_sigma2$C_IVT,
                       "mm.IVT"=PB_filter_sigma2$mm.IVT,
                       "mm.Pb_minus_Untreated"=PB_filter_sigma2$mm.PBminusUndiff,
                       "SNR"=PB_filter_sigma2$N_value)
    write.csv(PB_tab,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/SupplementaryTable6.csv",row.names = F)
    
  }
  #Supplemtary Table 7
  if(F){
    PB_filter_sigma2_validated<-PB_filter_sigma2[which(PB_filter_sigma2$color=="orange" | PB_filter_sigma2$color=="blue2"),]
    PB_tab<-data.frame("Annotation"=PB_filter_sigma2_validated$Annotation,"chr"=PB_filter_sigma2_validated$chr,"position"=PB_filter_sigma2_validated$position,
                       "strand"=PB_filter_sigma2_validated$strand,"kmer"=PB_filter_sigma2_validated$kmer,
                       "T_PB"=PB_filter_sigma2_validated$T_PB,"C_PB"=PB_filter_sigma2_validated$T_PB,
                       "T_Untreated"=PB_filter_sigma2_validated$T_Undiff,"C_Untreated"=PB_filter_sigma2_validated$C_Undiff,
                       "T_IVT"=PB_filter_sigma2_validated$T_IVT,"C_IVT"=PB_filter_sigma2_validated$C_IVT,
                       "mm.IVT"=PB_filter_sigma2_validated$mm.IVT,
                       "mm.Pb_minus_Untreated"=PB_filter_sigma2_validated$mm.PBminusUndiff,
                       "SNR"=PB_filter_sigma2_validated$N_value)
    write.csv(PB_tab,"/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Figures/SupplementaryTables/SupplementaryTable7.csv",row.names = F)
    
  }
  
  
}

### Pb plots
if(T){
################################## Plot: DotPlot U-to-C mismatch Pb-Undiff VS Undiff
if(T){  

    min.x = -60
    max.x = 60
    min.y = -3
    max.y =  3
    line.thickness = 1
    label.size = 1
    text.size = 0.8
    point.size = 0.2
    tck.length = 0.01
    tick.thickness = 1
    transparency = 1
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBvsUndiff_little_bootstrap_text.pdf"
    # Set up fonts
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")

    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    
    plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = '',xlab='',ann = F,axes=T,type="l", 
         xaxt='n', yaxt='n', 
         bty='n', bg='transparent')
    
    
    points(x=PB_filter$mm.PBminusUndiff[which(PB_filter$color == "black")], y=log10(PB_filter$N_value[which(PB_filter$color=="black")]), pch=22, col="grey", bg=adjustcolor("grey",alpha=0.05), cex = log2(PB_filter$total_PB[which(PB_filter$color == "black")])*point.size)
    points(x=PB_filter$mm.PBminusUndiff[which(PB_filter$color == "blue")], y=log10(PB_filter$N_value[which(PB_filter$color=="blue")]), pch=22, col="grey", bg=adjustcolor("grey",alpha=0.3), cex = log2(PB_filter$total_PB[which(PB_filter$color == "blue")])*point.size)
    points(x=PB_filter$mm.PBminusUndiff[which(PB_filter$color == "red")], y=log10(PB_filter$N_value[which(PB_filter$color=="red")]), pch=22, col="grey", bg=adjustcolor("grey",alpha=0.5), cex = log2(PB_filter$total_PB[which(PB_filter$color == "red")])*point.size)
    points(x=PB_filter$mm.PBminusUndiff[which(PB_filter$color == "orange")], y=log10(PB_filter$N_value[which(PB_filter$color=="orange")]), pch=22, col="orange", bg=adjustcolor("orange",alpha=0.8), cex = log2(PB_filter$total_PB[which(PB_filter$color == "orange")])*point.size)
    points(x=PB_filter$mm.PBminusUndiff[which(PB_filter$color == "blue2")], y=log10(PB_filter$N_value[which(PB_filter$color=="blue2")]), pch=22, col="blue2", bg=adjustcolor("blue2",alpha=0.3), cex = log2(PB_filter$total_PB[which(PB_filter$color == "blue2")])*point.size)
    
    segments(x0 = 0, y0 = -12, x1 = 0, y1 = 10, lty=2, col = "grey3")
    segments(x0 = -60, y0 = log10(sig), x1 = 60, y1 = log10(sig), lty=2, col = "grey3") 
    

    PB_sites_odd <- TRUB1sites
    text(x=PB_sites_odd$mm.PBminusUndiff, y=log10(PB_sites_odd$N_value),PB_sites_odd$Annotation, col='grey3', cex=0.3)
    
    PB_sites_even <- PUS7sites
    text(x=PB_sites_even$mm.PBminusUndiff-2, y=log10(PB_sites_even$N_value)-0.1,PB_sites_even$Annotation, col='grey3', cex=0.3)

    
    axis(side = 2,at = seq(min.y,max.y, by=1) ,labels = seq(min.y,max.y, by=1),lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    #axis(side = 2,at = seq(0,max(log2(PB_filter_sigma2$PB.standard.dev.....)), by=1),labels = seq(0, 10, by = 1),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
    axis(side = 1,at = seq(min.x,max.x,by=20) ,labels = seq(min.x,max.x,by=20)  ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    #axis(side = 1,at = c(-6:6)*10,labels = c(-6:6)*10,lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line= -0.5 )
    legend("topright", c("TRUB1 sites","PUS7 sites", "Other"),     
           cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
           pch=c(16,16,16),lty=c(0,0,0), 
           lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
    mtext(text = "log2(SNR)",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    mtext(text = "Pb treated - Untreated (U-C MM%)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
    alpha=0
    dev.off() 
  }
################################## Plot: sigma2 
if(T){
  min.x = 1
  max.x = round(log10(max(PB_filter_sigma2$total_PB)))
  min.y = -50
  max.y = 50
  line.thickness = 1
  label.size = 1
  text.size = 0.8
  point.size = 0.2
  tck.length = 0.01
  tick.thickness = 1
  transparency = 1

  file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBUndiff_log2_bootstrap_text.pdf"
  Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  
  plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = '',xlab='',ann = F,axes=T,type="l", 
       xaxt='n', yaxt='n', 
       bty='n', bg='transparent')
  
  
  points(x=log10(PB_filter_sigma2$total_PB[which(PB_filter_sigma2$color == "black")]), y=-PB_filter_sigma2$mm.PBminusUndiff[which(PB_filter_sigma2$color=="black")], pch=22, col="grey", bg=adjustcolor("grey",alpha=0.05), cex = log2(PB_filter_sigma2$mm.Undiff[which(PB_filter_sigma2$color == "black")])*point.size)
  points(x=log10(PB_filter_sigma2$total_PB[which(PB_filter_sigma2$color == "blue")]), y=-PB_filter_sigma2$mm.PBminusUndiff[which(PB_filter_sigma2$color=="blue")], pch=22, col="grey", bg=adjustcolor("grey",alpha=0.3), cex = log2(PB_filter_sigma2$mm.Undiff[which(PB_filter_sigma2$color == "blue")])*point.size)
  points(x=log10(PB_filter_sigma2$total_PB[which(PB_filter_sigma2$color == "red")]), y=-PB_filter_sigma2$mm.PBminusUndiff[which(PB_filter_sigma2$color=="red")], pch=22, col="grey", bg=adjustcolor("grey",alpha=0.5), cex = log2(PB_filter_sigma2$mm.Undiff[which(PB_filter_sigma2$color == "red")])*point.size)
  points(x=log10(PB_filter_sigma2$total_PB[which(PB_filter_sigma2$color == "orange")]), y=-PB_filter_sigma2$mm.PBminusUndiff[which(PB_filter_sigma2$color=="orange")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.8), cex = log2(PB_filter_sigma2$mm.Undiff[which(PB_filter_sigma2$color == "orange")])*point.size)
  points(x=log10(PB_filter_sigma2$total_PB[which(PB_filter_sigma2$color == "blue2")]), y=-PB_filter_sigma2$mm.PBminusUndiff[which(PB_filter_sigma2$color=="blue2")], pch=22, col="blue2", bg=adjustcolor("blue2",alpha=0.3), cex = log2(PB_filter_sigma2$mm.Undiff[which(PB_filter_sigma2$color == "blue2")])*point.size)
  
  PB_sites_odd_sigma2<-TRUB1sites[which(TRUB1sites$N_value>sig),]
  text(x=log10(PB_sites_odd_sigma2$total_PB), y=-PB_sites_odd_sigma2$mm.PBminusUndiff,PB_sites_odd_sigma2$Annotation, col='grey3', cex=0.5)
  
  PB_sites_even_sigma2<-PUS7sites[which(PUS7sites$N_value>sig),]
  text(x=log10(PB_sites_even_sigma2$total_PB), y=-PB_sites_even_sigma2$mm.PBminusUndiff+1,PB_sites_even_sigma2$Annotation, col='grey3', cex=0.5)
  
  
  #text(x=PBKD_sites$total_PB_KD, y=PBKD_sites$mmSCRminusPBKD,PBKD_sites$Annotation, col='grey3', cex=0.4)
  
  
  #for plot KD-SCR vs SCR
  axis(side = 2, at = seq(min.y, max.y, by = 10), labels = seq(min.y, max.y, by = 10),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  axis(side = 1, at = seq(min.x,max.x, by = 0.5), labels = seq(min.x,max.x, by = 0.5),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
       las = 1, lwd = line.thickness, line = 0)
  legend("topright", c("TRUB1 sites","PUS7 sites", "Other"),     
         cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
         pch=c(16,16,16),lty=c(0,0,0), 
         lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
  mtext(text = "Untreated - Pb treated (U-to-C MM%) ",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
  mtext(text = "log10(Reads)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
  
  
  dev.off()
} 
  
################################## Barplot interesting sites  
if(T){
  library(stringr)
  #### List of the TRUB1 sites with sigma>2 and with a delta>5 or <-5
  if(T){
    TRUB1sites5<-TRUB1sites_sigma2[which(abs(TRUB1sites_sigma2$mm.PBminusUndiff)>5),]
    PUS7sites5<-PUS7sites_sigma2[which(abs(PUS7sites_sigma2$mm.PBminusUndiff)>5),]
    pus_enzymes <-c(TRUB1sites5$Annotation,PUS7sites5$Annotation)
    
    pus_df <- data.frame("gene" = pus_enzymes, "PB1_TPM" = 0,"PB2_TPM" = 0,
                         "Undiff1_TPM" = 0,"Undiff2_TPM" = 0, "Undiff3_TPM" = 0)
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB1_TPM[i] <- TPM_PB1$TPM[which(TPM_PB1$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB2_TPM[i] <- TPM_PB2$TPM[which(TPM_PB2$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff1_TPM[i] <- TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff2_TPM[i] <- TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff3_TPM[i] <- TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == pus_df$gene[i])][1]
    }
    
    #### Add the standard deviation to the pus_df
    pus_df$PD_SD <- 0
    pus_df$Undiff_SD <- 0
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB_SD[i] <- sd(c(pus_df$PB1_TPM[i], pus_df$PB2_TPM[i]))
      pus_df$Undiff_SD[i] <- sd(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
    }
    #### Add the mean to the pus_df 
    pus_df$PB_mean <- 0
    pus_df$Undiff_mean <- 0
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB_mean[i] <- mean(c(pus_df$PB1_TPM[i], pus_df$PB2_TPM[i]))
      pus_df$Undiff_mean[i] <- mean(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
    }
    
    #pvalues, is Diff significantly different compared to undiff?
    for (i in c(1:nrow(pus_df))) {
      tPB<-t.test(c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i]),c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
      pus_df$pvalue[i] <- tPB$p.value
    }
  }
  
  #### Plot with segment on y axis because the rpll22 is way higher than other genes
  break_plot=T
  if(break_plot==T){ 
    library(plotrix)
    off=160
    min.x = 0
    max.x = nrow(pus_df)+2
    min.y = -1
    realmax.y=round(max(pus_df[,2:length(pus_df)]))+15
    max.y = realmax.y-off
    line.thickness = 1
    label.size = 1
    text.size = 1.4
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBvsUndiff_barplot_pvals.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 7, height = 6, family = "Helvetica")
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
    #box(lwd = line.thickness)
    #### Add the lines
   
    
    #for Undiff
    for (i in c(1:nrow(pus_df))) {
      if (pus_df$gene[i]=="RPL22"){
        lines(x=c(i-.2,i-.2),y=c(0,pus_df$Undiff_mean[i]-off), lwd=5 ,lend=1, lty=1, col = adjustcolor("black", alpha=0.7))
        pts=c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i])-off      
      }else if (pus_df$gene[i]=="ERH"){
        lines(x=c(i-.2,i-.2), y=c(0, pus_df$Undiff_mean[i]-(off/2)),lwd=5 ,lend=1, lty=1, col = adjustcolor("black", alpha=0.7))
        pts=c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i])-(off/2)     
      }else{
        lines(x=c(i-.2,i-.2),y=c(0,pus_df$Undiff_mean[i]),lwd=5 ,lend=1, lty=1, col = adjustcolor("black", alpha=0.7))
        pts=c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i])
      }
      
     
      points(x=rep(i-.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    
    #for PB
    for (i in c(1:nrow(pus_df))) {
      
      if (pus_df$gene[i]=="RPL22"){
        lines(x=c(i+.2,i+.2),y=c(0,pus_df$PB_mean[i]-off), lwd=5 ,lend=1, lty=1, col = "red")
        pts=c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i])-off
      }else if (pus_df$gene[i]=="ERH"){
        lines(x=c(i+.2,i+.2),y=c(0,pus_df$PB_mean[i]-(off/2)), lwd=5 ,lend=1, lty=1, col = "red")
        pts=c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i])-(off/2)
      }else{
        lines(x=c(i+.2,i+.2),y=c(0,pus_df$PB_mean[i]),lwd=5 ,lend=1, lty=1, col = "red")
        pts=c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i])
      }
      
      points(x=rep(i+.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
      
    }
    
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      if (pus_df$gene[i]=="RPL22"){
        lines(x=c(i-.2,i-.2), y=c(pus_df$Undiff_mean[i]-off-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]-off+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }else if (pus_df$gene[i]=="ERH"){
        lines(x=c(i-.2,i-.2),y=c(pus_df$Undiff_mean[i]-(off/2)-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]-(off/2)+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
        }else{
        lines(x=c(i-.2,i-.2),
              y=c(pus_df$Undiff_mean[i]-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }
    }
    
    
    for (i in c(1:nrow(pus_df))) {
      if (pus_df$gene[i]=="RPL22"){
        lines(x=c(i+.2,i+.2),y=c(pus_df$PB_mean[i]-off-pus_df$PB_SD[i],pus_df$PB_mean[i]-off+pus_df$PB_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }else if (pus_df$gene[i]=="ERH"){
        lines(x=c(i+.2,i+.2), y=c(pus_df$PB_mean[i]-(off/2)-pus_df$PB_SD[i],pus_df$PB_mean[i]-(off/2)+pus_df$PB_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }else{
        lines(x=c(i+.2,i+.2),y=c(pus_df$PB_mean[i]-pus_df$PB_SD[i],pus_df$PB_mean[i]+pus_df$PB_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }
    }
    
    
    
    #plot pvalue asterisk
    #### Add Asterisks for p-values
    # *   < 0.05
    # **  < 0.01
    # *** < 0.001
    
    for (i in c(1:nrow(pus_df))){
      if (pus_df$pvalue[i]<0.001) {
        text(i,max.y-30, labels = "***", col = "black", cex = 0.6)
        
      } else if (pus_df$pvalue[i]<0.01) {
        text(i,max.y-30, labels = "**", col = "black", cex = 0.6)
       
      } else if (pus_df$pvalue[i]<0.05) {
        text(i,max.y-30, labels = "*", col = "black", cex = 0.6)
      
      } 
    }
    
    
    
    ngenes=nrow(pus_df)
    
    axis(side = 2,at = c(-1,50,100,max.y) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    axis(side = 2,at =  c(0,50,100,110,max.y),labels = c(0,50,100,max.y,realmax.y),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
    axis.break(2, 100, style = "gap") 
    axis(side = 1,at = c(0:max.x) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
    axis(side = 1,at = c(1:ngenes),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
    mtext(text = "TPM",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    legend("topleft", c("Untreated","Pb"),
           cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
           pch=c(15,15,15),lty=c(0,0,0),
           lwd=c(2.5,2.5,2.5),col=c("black", "red"), bty = "n")
    
    
    dev.off() 
    
  }  
  
  #### Plot log(TPM)
  if(F){ 
    library(plotrix)
    min.x = 0
    max.x = nrow(pus_df)+1
    min.y = -1
    max.y = round(log10(max(pus_df[,2:length(pus_df)])))+2
    line.thickness = 1
    label.size = 1
    text.size = 1.4
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
  
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBvsUndiff_barplot_pvals_log.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
    #box(lwd = line.thickness)
    #### Add the lines
    
    
    #for Undiff
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.2,i-.2),y=c(0,log10(pus_df$Undiff_mean[i])),lwd=5 ,lend=1, lty=1, col = adjustcolor("black",alpha=0.7))
      pts=log10(c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
      points(x=rep(i-.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    
    #for PB
    for (i in c(1:nrow(pus_df))) {
      
      lines(x=c(i+.2,i+.2),y=c(0,log10(pus_df$PB_mean[i])),lwd=5 ,lend=1, lty=1, col = "red3")
      pts=log10(c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i]))
      points(x=rep(i+.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.2,i-.2),
            y=c(log10(pus_df$Undiff_mean[i])-log10(pus_df$Undiff_SD[i]),log10(pus_df$Undiff_mean[i])+log10(pus_df$Undiff_SD[i])),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    
    for (i in c(1:nrow(pus_df))) {
      
      lines(x=c(i+.2,i+.2),
            y=c(log10(pus_df$PB_mean[i])-log10(pus_df$PB_SD[i]),log10(pus_df$PB_mean[i])+log10(pus_df$PB_SD[i])),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    abline(h=0)
    
    #plot pvalue asterisk
    #### Add Asterisks for p-values
    # *   < 0.05
    # **  < 0.01
    # *** < 0.001
    
    for (i in c(1:nrow(pus_df))){
      if (pus_df$pvalue[i]<0.001) {
        text(i,max.y-0.5, labels = "***", col = "red3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-1.,max.y-1),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalue[i]<0.01) {
        text(i,max.y-0.5, labels = "**", col = "red3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-1,max.y-1),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalue[i]<0.05) {
        text(i,max.y-0.5, labels = "*", col = "red3", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-1,max.y-1),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } 
    }
    
    
    
    ngenes=nrow(pus_df)
    
    axis(side = 2,at = c(min.y:max.y) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    axis(side = 2,at =  c(min.y:max.y),labels = c(-1:max.y),lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
    axis(side = 1,at = c(0:max.x) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
    axis(side = 1,at = c(1:ngenes),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
    mtext(text = "log10(TPM)",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    legend("topleft", c("Untreated","PB treated"),
           cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
           pch=c(15,15,15),lty=c(0,0,0),
           lwd=c(2.5,2.5,2.5),col=c("black", "red3"), bty = "n")
    
    
    dev.off() 
    
  }
  
}

#####write table fig3 for Pb Undiff
if(T){
  # sites_minusplus20<-0
  # sites_minusplus20<-PB_filter_sigma2[which(abs(PB_filter_sigma2$mm.PBminusUndiff)>20),]
  # SitesTable<-0
  # SitesTable<-data.frame(Annotation=c(TRUB1sites5$Annotation,PUS7sites5$Annotation,sites_minusplus20$Annotation), chr=c(TRUB1sites5$chr,PUS7sites5$chr,sites_minusplus20$chr),
  #                        position=c(TRUB1sites5$position,PUS7sites5$position,sites_minusplus20$position), kmer=c(TRUB1sites5$kmer,PUS7sites5$kmer,sites_minusplus20$kmer),
  #                        mm.UntreatedminusIVT=c(TRUB1sites5$mm.UndiffMINUSmm.IVT,PUS7sites5$mm.UndiffMINUSmm.IVT,sites_minusplus20$mm.UndiffMINUSmm.IVT),
  #                        mm.PBminusIVT=c(TRUB1sites5$mm.PBMINUSmm.IVT,PUS7sites5$mm.PBMINUSmm.IVT,sites_minusplus20$mm.PBMINUSmm.IVT))
  
  sites_validated<-PB_filter_sigma2[which(PB_filter_sigma2$color=="orange" | PB_filter_sigma2$color=="blue2"),]
  sites_validated<-sites_validated[which(abs(sites_validated$mm.PBMINUSmm.Undiff)>5),]  
  SitesTable<-data.frame(Annotation=sites_validated$Annotation, chr=sites_validated$chr,
                         position=sites_validated$position, kmer=sites_validated$kmer,
                         mm.UntreatedminusIVT=sites_validated$mm.UndiffMINUSmm.IVT,
                         mm.PBminusIVT=sites_validated$mm.PBMINUSmm.IVT)
  
  write.csv(SitesTable,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Tablefig3_PB.csv",row.names = F)
  
}    
  
############################### seqLogos on the 5 different intervals
if(T){
    ####int1: >20
    ####int2: >=5 & <=20
    ####int3: >-5 & <5
    ####int4: <=-5 & >=-20
    ####int5: <-20
    PB_filter_sigma2$mm.UndiffminusPB=-PB_filter_sigma2$mm.PBminusUndiff
    int1<-PB_filter_sigma2[which(PB_filter_sigma2$mm.UndiffminusPB >   20 ),]
    int2<-PB_filter_sigma2[which(PB_filter_sigma2$mm.UndiffminusPB >=  5 & PB_filter_sigma2$mm.UndiffminusPB <=  20 ),]
    int3<-PB_filter_sigma2[which(PB_filter_sigma2$mm.UndiffminusPB >  -5 & PB_filter_sigma2$mm.UndiffminusPB <   5  ),]
    int4<-PB_filter_sigma2[which(PB_filter_sigma2$mm.UndiffminusPB <= -5 & PB_filter_sigma2$mm.UndiffminusPB >= -20 ),]
    int5<-PB_filter_sigma2[which(PB_filter_sigma2$mm.UndiffminusPB <  -20),]
    
    if (T){
      
      # kmers<-sapply(int1$kmer,function(x) str_replace_all(x,"T","U"))
      # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBLogo>20.pdf"
      # 
      # kmers<-sapply(int2$kmer,function(x) str_replace_all(x,"T","U"))
      # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBLogo20to5.pdf"
      # 
      # kmers<-sapply(int3$kmer,function(x) str_replace_all(x,"T","U"))
      # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBLogo5to-5.pdf"
      # 
      # kmers<-sapply(int4$kmer,function(x) str_replace_all(x,"T","U"))
      # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBLogo-5to-20.pdf"
      # 
      # kmers<-sapply(int5$kmer,function(x) str_replace_all(x,"T","U"))
      # file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBLogo<-20.pdf"

      my_list<-list("new motifs" = kmers)
      min.x = 0
      max.x = 200
      min.y = 0
      max.y = 20
      line.thickness = 1
      label.size = 1
      text.size = 0.8
      point.size = 1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 0.2
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 

      ggseqlogo::ggseqlogo(my_list,ncol=1,method="prob")
      #same as
      #ggplot() + geom_logo(seqs_dna) + theme_logo() + facet_wrap(~seq_group, ncol=4, scales='free_x') 
      
      dev.off()
      
    }
  }

  
########################## Kinetic model Diff Undiff
if(T){
#([mRNAu][mRNApus])a -> ([mRNAu][mRNApus])b
#a=Diff
#b=Undiff
#plot BminusA

#trub1 and pus7 sites validated with KD libraries
confirmed_sites_names<-c(TRUB1sites_sigma2$ACP,PUS7sites_sigma2$ACP)

#plot
if(T){
  
  scale=1000000
  
  pb=T
  #run the pus_df from the Barplot of the PUS enzymes that is above 
  
  #PB
  if(pb==T){
 
    kinetic_sites<-PB_filter[which(PB_filter$ACP %in% confirmed_sites_names),]
    
    PB_TPM_TRUB1 <- mean(c(TPM_PB1$TPM[which(TPM_PB1$Gene.Name == "TRUB1")],TPM_PB2$TPM[which(TPM_PB2$Gene.Name == "TRUB1")]))
    Undiff_TPM_TRUB1 <-mean(c(TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == "TRUB1")],TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == "TRUB1")],TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == "TRUB1")]))
    
    PB_TPM_PUS7 <- mean(c(TPM_PB1$TPM[which(TPM_PB1$Gene.Name == "PUS7")],TPM_PB2$TPM[which(TPM_PB2$Gene.Name == "PUS7")]))
    Undiff_TPM_PUS7 <-mean(c(TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == "PUS7")],TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == "PUS7")],TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == "PUS7")]))
    
    M=totPB #normalization
    kinetic_sites$ConditionA<-0
    kinetic_sites$ConditionA[which(kinetic_sites$color=="orange")] <- (kinetic_sites$total_PB[which(kinetic_sites$color=="orange")]/M)*scale*PB_TPM_TRUB1
    #kinetic_sites$ConditionA[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_TRUB1_KD[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$TRUB1_mean[pus_df$gene=="PUS7"]
    kinetic_sites$ConditionA[which(kinetic_sites$color=="blue2")] <- (kinetic_sites$total_PB[which(kinetic_sites$color=="blue2")]/M)*scale*PB_TPM_PUS7
    kinetic_sites$yConditionA <- kinetic_sites$C_PB/M*scale 
    M=totUndiff #normalization
    kinetic_sites$ConditionB<-0
    kinetic_sites$ConditionB[which(kinetic_sites$color=="orange")] <- (kinetic_sites$total_Undiff[which(kinetic_sites$color=="orange")]/M)*scale*Undiff_TPM_TRUB1
    #kinetic_sites$ConditionB[which(kinetic_sites$color=="blue")] <- (kinetic_sites$total_SCR[which(kinetic_sites$color=="blue")]/M)*scale*pus_df$Scram_mean[pus_df$gene=="PUS7"] 
    kinetic_sites$ConditionB[which(kinetic_sites$color=="blue2")] <- (kinetic_sites$total_Undiff[which(kinetic_sites$color=="blue2")]/M)*scale*Undiff_TPM_PUS7 
    kinetic_sites$yConditionB <- kinetic_sites$C_Undiff/M*scale 
    
    kinetic_sites$BminusA  <- kinetic_sites$ConditionB - kinetic_sites$ConditionA
    kinetic_sites$yBminusA <- kinetic_sites$yConditionB - kinetic_sites$yConditionA
  }
  
  #plot conditionUndiff-conditionPB
  if(T){
    min.x = -200 #round(min(kinetic_sites$BminusA))-49
    max.x = 200#round(max(kinetic_sites$BminusA))+50
    min.y = -20
    max.y = 20
    line.thickness = 1
    label.size = 1
    text.size = 0.8
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 1
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/Kinetic_model_PB_text.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    
    plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = '',xlab='',ann = F,axes=T,type="l", 
         xaxt='n', yaxt='n', 
         bty='n', bg='transparent')
    
    
    points(x=kinetic_sites$BminusA[which(kinetic_sites$color == "orange")], y=kinetic_sites$yBminusA[which(kinetic_sites$color=="orange")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.3), cex = point.size)
    points(x=kinetic_sites$BminusA[which(kinetic_sites$color == "blue2")], y=kinetic_sites$yBminusA[which(kinetic_sites$color=="blue2")], pch=22, col="blue", bg=adjustcolor("blue",alpha=0.6), cex=point.size)
    
    # Correlation TRUB1 sites only 
    x <- kinetic_sites$BminusA[which(kinetic_sites$color == "orange")]
    y <- kinetic_sites$yBminusA[which(kinetic_sites$color=="orange")]
    
    # Calculate correlation coefficient
    correlation_coefficient <- cor(x, y)
    # Calculate R-squared value
    r_squared <- correlation_coefficient^2
    
    # Print the correlation coefficient and R-squared value
    print(paste("Correlation Coefficient trub1:", correlation_coefficient))
    print(paste("R-squared Value:", r_squared))
    
    # Add correlation line
    abline(lm(y ~ x), col = "orange")
    
    
    x<-kinetic_sites$BminusA[which(kinetic_sites$color == "blue2")]
    y<-kinetic_sites$yBminusA[which(kinetic_sites$color=="blue2")]
    
    # Calculate correlation coefficient
    correlation_coefficient <- cor(x, y)
    # Calculate R-squared value
    r_squared <- correlation_coefficient^2
    
    # Add correlation line
    abline(lm(y ~ x), col = "blue")
    
    # Print the correlation coefficient and R-squared value
    print(paste("Correlation Coefficient pus7:", correlation_coefficient))
    print(paste("R-squared Value:", r_squared))
    
    abline(v=0, col = "black")
    
    
    text(x=kinetic_sites$BminusA, y=kinetic_sites$yBminusA,kinetic_sites$Annotation, col='grey3', cex=0.4)
    
    
    #for plot KD-SCR vs SCR
    axis(side = 2, at = seq(min.y, max.y, by = 10), labels = seq(min.y, max.y, by = 10),
         lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
         las = 1, lwd = line.thickness, line = 0)
    axis(side = 1, at = seq(min.x,max.x, by = 50), labels = seq(min.x, max.x, by = 50),
         lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
         las = 1, lwd = line.thickness, line = 0)
    legend("topright", c("TRUB1 validated sites","PUS7 validated sites"),     
           cex= text.size-0.1, pt.cex=1, y.intersp=1,x.intersp=0.01,
           pch=c(16,16,16),lty=c(0,0,0), 
           lwd=c(2.5,2.5,2.5),col=c("orange","blue"), bty = "n")
    mtext(text = expression(" ψUndiff - ψPb per million "),side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    mtext(text = "([mRNA U][mRNA PUS])Undiff - ([mRNA U][mRNA PUS])Pb",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
    ##ggsave(filename =file.name,width = 6,height = 6.6)
    alpha=0
    dev.off()
    
    
  }
  
  
}    

}


########################## Barplot Pb Undiff PUS ezymes
if (T) {
  #### List of the PUS enzyme genes:
  pus_enzymes <- c("PUS1", "PUSL1", "PUS3", "TRUB1", "TRUB2",
                   "DKC1", "PUS7", "PUS7L", "RPUSD1", "RPUSD2",
                   "PUS10", "PRUSD3", "PRUSD4")
  #### Make the data frame of TPM of the pus_enzymes 
  pus_df <- data.frame("gene" = pus_enzymes, "Diff1_TPM" = 0,"Diff2_TPM" = 0,"Diff3_TPM" = 0,
                       "Undiff1_TPM" = 0,"Undiff2_TPM" = 0, "Undiff3_TPM" = 0)
  for (i in c(1:nrow(pus_df))) {
    pus_df$PB1_TPM[i] <- TPM_PB1$TPM[which(TPM_PB1$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$PB2_TPM[i] <- TPM_PB2$TPM[which(TPM_PB2$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff1_TPM[i] <- TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff2_TPM[i] <- TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == pus_df$gene[i])][1]
  }
  for (i in c(1:nrow(pus_df))) {
    pus_df$Undiff3_TPM[i] <- TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == pus_df$gene[i])][1]
  }
  
  #### Add the standard deviation to the pus_df
  pus_df$PB_SD <- 0
  pus_df$Undiff_SD <- 0
  for (i in c(1:nrow(pus_df))) {
    pus_df$PB_SD[i] <- sd(c(pus_df$PB1_TPM[i], pus_df$PB2_TPM[i]))
    pus_df$Undiff_SD[i] <- sd(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
  }
  #### Add the mean to the pus_df 
  pus_df$PB_mean <- 0
  pus_df$Undiff_mean <- 0
  for (i in c(1:nrow(pus_df))) {
    pus_df$PB_mean[i] <- mean(c(pus_df$PB1_TPM[i], pus_df$PB2_TPM[i]))
    pus_df$Undiff_mean[i] <- mean(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
  }
  
  #pvalues, is Diff significantly different compared to undiff?
  for (i in c(1:nrow(pus_df))) {
    tPB<-t.test(c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i]),c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
    pus_df$pvalue[i] <- tPB$p.value
  }
  
  #### Remove the unwanted rows
  pus_df = pus_df[-c(12,13),]
  

  
  #### Plot 
  if(T){ 
    min.x = 0
    max.x = 15
    min.y = 0
    max.y = 45
    line.thickness = 1
    label.size = 1
    text.size = 1.4
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PbvsUndiff_barplot_PUSenz_pvals_log.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
    #box(lwd = line.thickness)
    #### Add the lines
    
    #for Undiff
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.2,i-.2),y=c(0,pus_df$Undiff_mean[i]),lwd=5 ,lend=1, lty=1, col = adjustcolor("black", alpha=0.7))
      pts=c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i])
      points(x=rep(i-.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    #for PB
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i+.2,i+.2),y=c(0,pus_df$PB_mean[i]),lwd=5 ,lend=1, lty=1, col = "red")
      pts=c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i])
      points(x=rep(i+.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i+.2,i+.2),
            y=c(pus_df$PB_mean[i]-pus_df$PB_SD[i],pus_df$PB_mean[i]+pus_df$PB_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.2,i-.2),
            y=c(pus_df$Undiff_mean[i]-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    for (i in c(1:nrow(pus_df))){
      if (pus_df$pvalue[i]<0.001) {
        text(i,max.y-0.5, labels = "***", col = "red", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-1.,max.y-1),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalue[i]<0.01) {
        text(i,max.y-0.5, labels = "**", col = "red", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-1,max.y-1),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalue[i]<0.05) {
        text(i,max.y-0.5, labels = "*", col = "red", cex = 0.6)
        lines(x=c(i-0.2,i+0.2),y=c(max.y-1,max.y-1),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } 
    }
    
    
    axis(side = 2,at = c(0:5)*10 ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    axis(side = 2,at = c(0:5)*10,labels = c(0:5)*10,lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
    axis(side = 1,at = c(0:12) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
    axis(side = 1,at = c(1:11),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
    legend("topright", c("Undifferentiated","Lead treated"),
           cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
           pch=c(16,16),lty=c(0,0),
           lwd=c(2.5,2.5,2.5),col=c("black", "red"), bty = "n")
    mtext(text = "TPM",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    dev.off() 
  }
  
  
}  
  
########################## Gene Onthology
#Load deseq2 data
PB_Undiff_deseq2<-read.csv(file="/home/sasha/Rouhanifard lab Dropbox/Papers/Pseudouridine neuronal differentiation/Supplementary_tables/PBUndiff_Deseq2.csv",header=T,row.names = F)

### for gene experession differences in PB and Undiff
if (T) {
  ## Filter to significantly different positions
  res_filter1_Up <- PB_Undiff_deseq2[which(PB_Undiff_deseq2$pvalue < 0.01 & PB_Undiff_deseq2$log2FoldChange >1),]
  res_filter1_Down <- PB_Undiff_deseq2[which(PB_Undiff_deseq2$pvalue < 0.01 & PB_Undiff_deseq2$log2FoldChange < -1),]
  sigs.df_Up <- as.data.frame(res_filter1_Up)
  sigs.df_Down <- as.data.frame(res_filter1_Down)
  ## Add the separate gene name
  res_filter1_Up$transcript_name <- res_filter1_Up$X
  res_filter1_Up$gene <- ""
  res_filter1_Up$gene_ENST <- ""
  for (i in c(1:nrow(res_filter1_Up))) {
    res_filter1_Up$gene[i] <- strsplit(res_filter1_Up$transcript_name[i], "\\|") [[1]][6]
    res_filter1_Up$ENST[i] <- strsplit(res_filter1_Up$transcript_name[i], "\\|") [[1]][1]
    res_filter1_Up$gene_ENST[i] <- paste0(res_filter1_Up$gene[i],"|",res_filter1_Up$ENST[i])
  }
  res_filter1_Down$transcript_name <- res_filter1_Down$X
  res_filter1_Down$gene <- ""
  res_filter1_Down$gene_ENST <- ""
  for (i in c(1:nrow(res_filter1_Down))) {
    res_filter1_Down$gene[i] <- strsplit(res_filter1_Down$transcript_name[i], "\\|") [[1]][6]
    res_filter1_Down$ENST[i] <- strsplit(res_filter1_Down$transcript_name[i], "\\|") [[1]][1]
    res_filter1_Down$gene_ENST[i] <- paste0(res_filter1_Down$gene[i],"|",res_filter1_Down$ENST[i])
  }
  ## Extract the gene names
  write.csv(res_filter1_Up$gene, "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PbUndiff_Up.csv",row.names = F)
  write.csv(res_filter1_Down$gene, "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PbUndiff_Down.csv",row.names = F)
  ## Input the GO
  
  
  
  
}

}
 






##### Figure 5 #################################################################

############################# merged PB and SH
if(T){
  PB_filter$ACP<-paste(PB_filter$Annotation,PB_filter$chr,PB_filter$position,sep='')
  SH_filter$ACP<-paste(SH_filter$Annotation,SH_filter$chr,SH_filter$position,sep='')
  common <- intersect(PB_filter$ACP, SH_filter$ACP)  
  ALL_PB<-PB_filter[which(PB_filter$ACP %in% common),]
  ALL_SH<-SH_filter[which(SH_filter$ACP %in% common),]
  ALL_PB_SH<-cbind(ALL_PB, ALL_SH$total_Diff,ALL_SH$mm.Diff)
  ALL_PB_SH$color[which(ALL_PB_SH$ACP %in% TRUB1KD_sites$ACP )]<-"orange"
  ALL_PB_SH$color[which(ALL_PB_SH$ACP %in% PUS7KD_sites$ACP)]<-"blue2"
  file_to_save<-data.frame("Annotation"=ALL_PB_SH$Annotation,"chr"=ALL_PB_SH$chr,"pos"=ALL_PB_SH$position,"kmer"=ALL_PB_SH$kmer,"color"=ALL_PB_SH$color,"mm.IVT"=ALL_PB_SH$mm.IVT,"mm.Undiff"=ALL_PB_SH$mm.Undiff, "mm.Diff"=ALL_PB_SH$`ALL_SH$mm.Diff`,"mm.PB"=ALL_PB_SH$mm.PB,"N_value"=ALL_PB_SH$N_value)
  file_to_save_less10ivt<-file_to_save[which(file_to_save$mm.IVT<=10),]
  write.csv(file_to_save_less10ivt,"/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PB_SH_merged.less10ivt.csv",row.names = F)
}

#plot Diff-Lead vs Undiff
PB_SH<-read.csv(file="/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PB_SH_merged.less10ivt.csv", header=T)
PB_SH$mm.UndiffMINUSmm.IVT<-PB_SH$mm.Undiff-PB_SH$mm.IVT
PB_SH$mm.DiffMINUSmm.IVT<-PB_SH$mm.Diff-PB_SH$mm.IVT
PB_SH$mm.PBMINUSmm.IVT<-PB_SH$mm.PB-PB_SH$mm.IVT

#for the heatmap
PB_SH2 <- PB_SH[which( ((PB_SH$mm.UndiffMINUSmm.IVT>40) | 
                          (PB_SH$mm.DiffMINUSmm.IVT >40) |
                          (PB_SH$mm.PBMINUSmm.IVT >40))),]
### Add SD to data_matrix
PB_SH2$SD <- 0
for (i in c(1:nrow(PB_SH2))) {
  PB_SH2$SD[i] <- sd(c(PB_SH2$mm.UndiffMINUSmm.IVT[i], PB_SH2$mm.DiffMINUSmm.IVT[i], PB_SH2$mmPBMINUSmm.IVT[i]))
}
PB_SH2 = PB_SH2[order(PB_SH2$SD, decreasing = F),]

#Pb-Diff vs Undiff
if(T){
  min.x = 0 #round(min(kinetic_sites$BminusA))-49
  max.x = 100#round(max(kinetic_sites$BminusA))+50
  min.y = -50
  max.y = 50
  line.thickness = 1
  label.size = 1
  text.size = 0.8
  point.size = 1
  tck.length = 0.01
  tick.thickness = 1
  transparency = 1
  
  file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/fig5_DiffLeadvsUntreated.pdf"
  Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
  
  plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = '',xlab='',ann = F,axes=T,type="l", 
       xaxt='n', yaxt='n', 
       bty='n', bg='transparent')
  
  points(y=PB_SH$mm.Diff-PB_SH$mm.PB, x=PB_SH$mm.Undiff, pch=22, col="grey", bg=adjustcolor("grey",alpha=0.3), cex = point.size)
  points(y=PB_SH$mm.Diff[which(PB_SH$color == "orange")]-PB_SH$mm.PB[which(PB_SH$color == "orange")], x=PB_SH$mm.Undiff[which(PB_SH$color == "orange")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.3), cex = point.size)
  points(y=PB_SH$mm.Diff[which(PB_SH$color == "blue2")]-PB_SH$mm.PB[which(PB_SH$color == "blue2")], x=PB_SH$mm.Undiff[which(PB_SH$color == "blue2")], pch=22, col="blue2", bg=adjustcolor("blue2",alpha=0.3), cex = point.size)
  
  abline(h=0,col="black",lty=2)
  
  text(x=PB_SH$mm.Undiff[which(PB_SH$color=="orange")], y=PB_SH$mm.Diff[which(PB_SH$color == "orange")]-PB_SH$mm.PB[which(PB_SH$color == "orange")]+2,PB_SH$Annotation[which(PB_SH$color == "orange")], col='grey3', cex=0.4)
  text(x=PB_SH$mm.Undiff[which(PB_SH$color=="blue2")], y=PB_SH$mm.Diff[which(PB_SH$color=="blue2")]-PB_SH$mm.PB[which( PB_SH$color=="blue2")],PB_SH$Annotation[which( PB_SH$color=="blue2")], col='grey3', cex=0.4)
  
  #for plot KD-SCR vs SCR
  axis(side = 2, at = seq(min.y, max.y, by = 10), labels = seq(min.y, max.y, by = 10),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  axis(side = 1, at = seq(min.x,max.x, by = round(max.x/10)), labels = seq(min.x, max.x, by = round(max.x/10)),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  legend("topright", c("TRUB1 validated sites","PUS7 validated sites","Other"),     
         cex= text.size-0.1, pt.cex=1, y.intersp=1,x.intersp=0.01,
         pch=c(16,16,16),lty=c(0,0,0), 
         lwd=c(2.5,2.5,2.5),col=c("orange","blue","grey"), bty = "n")
  mtext(text = " Differentiated - Lead (U-C %mm)  ",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
  mtext(text = " Undifferentiated U-C %mm ",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
  ##ggsave(filename =file.name,width = 6,height = 6.6)
  alpha=0
  dev.off()
  
  
  
    
}

#Pb-Diff vs Undiff
if(T){
  min.x = -50 #round(min(kinetic_sites$BminusA))-49
  max.x = 50#round(max(kinetic_sites$BminusA))+50
  min.y = -50
  max.y = 50
  line.thickness = 1
  label.size = 1
  text.size = 0.8
  point.size = 1
  tck.length = 0.01
  tick.thickness = 1
  transparency = 1
  
  file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/fig5_DiffUntr_vs_LeadUntr.pdf"
  Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
  
  plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = '',xlab='',ann = F,axes=T,type="l", 
       xaxt='n', yaxt='n', 
       bty='n', bg='transparent')
  
  points(y=PB_SH$mm.PB-PB_SH$mm.Undiff, x=PB_SH$mm.Diff-PB_SH$mm.Undiff, pch=22, col="grey", bg=adjustcolor("grey",alpha=0.3), cex = point.size)
  points(y=PB_SH$mm.PB[which(PB_SH$color == "orange")]-PB_SH$mm.Undiff[which(PB_SH$color == "orange")], x=PB_SH$mm.Diff[which(PB_SH$color == "orange")]-PB_SH$mm.Undiff[which(PB_SH$color == "orange")], pch=22, col="orange", bg=adjustcolor("orange",alpha=0.3), cex = point.size)
  points(y=PB_SH$mm.PB[which(PB_SH$color == "blue2")]-PB_SH$mm.Undiff[which(PB_SH$color == "blue2")], x=PB_SH$mm.Diff[which(PB_SH$color == "blue2")]-PB_SH$mm.Undiff[which(PB_SH$color == "blue2")], pch=22, col="blue2", bg=adjustcolor("blue2",alpha=0.3), cex = point.size)
  
  # abline(h=0,col="black",lty=2)
  
  text(x=PB_SH$mm.Diff[which(PB_SH$color == "orange")]-PB_SH$mm.Undiff[which(PB_SH$color=="orange")], y=PB_SH$mm.PB[which(PB_SH$color == "orange")]-PB_SH$mm.Undiff[which(PB_SH$color=="orange")],PB_SH$Annotation[which(PB_SH$color == "orange")], col='grey3', cex=0.3)
  text(x=PB_SH$mm.Diff[which(PB_SH$color=="blue2")]-PB_SH$mm.Undiff[which(PB_SH$color=="blue2")], y=PB_SH$mm.PB[which(PB_SH$color=="blue2")]-PB_SH$mm.Undiff[which(PB_SH$color=="blue2")],PB_SH$Annotation[which( PB_SH$color=="blue2")], col='grey3', cex=0.3)
  
  #for plot KD-SCR vs SCR
  axis(side = 2, at = seq(min.y, max.y, by = 10), labels = seq(min.y, max.y, by = 10),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  axis(side = 1, at = seq(min.x,max.x, by = round(max.x/10)), labels = seq(min.x, max.x, by = round(max.x/10)),
       lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length, 
       las = 1, lwd = line.thickness, line = 0)
  legend("topright", c("TRUB1 validated sites","PUS7 validated sites","Other"),     
         cex= text.size-0.1, pt.cex=1, y.intersp=1,x.intersp=0.01,
         pch=c(16,16,16),lty=c(0,0,0), 
         lwd=c(2.5,2.5,2.5),col=c("orange","blue","grey"), bty = "n")
  mtext(text = " Lead - Untreated (U-C %mm)  ",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
  mtext(text = " Differentiated - Untreated (U-C %mm) ",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
  ##ggsave(filename =file.name,width = 6,height = 6.6)
  alpha=0
  dev.off()

  }

### Heatmpap 3 columns
if (T) {
  ################################# Heatmapt Jurkat and Naive to see different level of mismatch for the common genes
  library(reshape2)
  library(ggplot2)
  library(gplots)
  library(dplyr)
  library(ggseqlogo)
  

  
  data_matrix <- PB_SH2[,c(11,12,13,14)]
  colnames(data_matrix) <- c("Untreated", "Diff", "Lead", "SD")
  AC<-paste(PB_SH2$Annotation,PB_SH2$chr,sep = " - ")
  ACP<-paste(AC,PB_SH2$pos,sep = " : ")
  ACPS<-paste(sprintf("%.2f",PB_SH2$SD),ACP,sep = "  ")
  rownames(data_matrix) <- ACPS
  
  row_name_colors <- PB_SH2$color
  row_name_colors[which(row_name_colors=='black')]<-"white"
  row_name_colors[which(row_name_colors=="blue")]<-"lightblue"

  colors <- c("#047FE1","#047FE1", "#1184A7","#15A2A2","#E1AAFF", "#FFAFF0", "#B4418E", "#D94A8C", "#EA515F",
              "#FE7434", "#FEE009")
  breaks <- c(0,9.99, 19.99, 29.99, 39.99, 49.99, 59.99, 69.99, 79.99, 89.99, 99.99, 100)
  
  labels <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "100")
  
  #heatmap plot
  if(T){
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/fig5_heatmap.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 5, height = 6, family = "Helvetica")
    min.x = 0
    max.x = 100
    min.y = 0
    max.y = 800
    line.thickness = 1
    label.size = 1
    text.size = 0.8
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    #par(mfrow=c(1,1),mar = c(1,1,6,1), lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
   
    gplots::heatmap.2(as.matrix(data_matrix[,c(1,2,3)]),
                      Rowv = FALSE,
                      Colv = FALSE,
                      col = colorRampPalette(colors)(length(breaks)-1),
                      #col = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(length(breaks) - 1),
                      breaks = breaks,
                      dendrogram = "none",
                      key = FALSE,
                      key.title = "Color Key",
                      density.info = "none",
                      trace ="none",
                      margins=c(9,15), #0 starts from the bottom right, the first value controls if it's up/down. The second value left/right
                      cexRow=0.4,
                      cexCol = 0.8,
                      lhei = c(1,8),
                      lwid = c(0.5,4),
                      colsep=1:nrow(data_matrix), # Add vertical grid lines
                      # rowsep=1:nrow(data_matrix), # Add horizontal grid lines
                      sepcolor = "white", #modify color of separating grid lines
                      RowSideColors = row_name_colors)
    # Add the color key manually at the bottom center
    legend("top", fill = colors, xpd = TRUE,
           legend = labels,
           bty = "y", horiz = TRUE, inset = c(-0.5, -0.15), x.intersp = 0.5,cex=0.5)
   
    dev.off()
  }
}

#table creation with pnas dataframe
if(T){ 
  
  pnas_filter=PB_SH2
  pnas_filter$color[which(pnas_filter$color=="blue")]<-"dodgerblue2"
  pnas_filter$kmer<-sapply(pnas_filter$kmer,function(x) str_replace_all(x,"T","U")) #tidyverse necessary
  
  ### Criteria
  min.x = -50
  max.x = 70
  min.y = -10
  max.y = nrow(pnas_filter)+ 5
  line.thickness = 0.8
  txt.size = 0.2
  point.size = 0.3
  label.size = 0.5
  anot.pos = -20
  file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/fig5_heatmaplist.pdf"
pdf(file=file.name)
    # Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
  # CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
  
  par(mfrow=c(1,1),lend=1)
  plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
       ylab = NA,xlab=NA,ann = F,axes=F)
  ### Add the gene texts on Y axis and dot for each methods
  for (i in c(1:nrow(pnas_filter))) {
    text(x = anot.pos ,y = nrow(pnas_filter)-i, labels = pnas_filter$Annotation[i],
         srt = 0,cex = 1.2*txt.size, pos = 4, col = pnas_filter$color[i], family = "Helvetica", font = 3)
    ### Add the chr
    text(x = anot.pos+10 ,y = nrow(pnas_filter)-i,
         labels = pnas_filter$chr[i] ,cex = 1.2*txt.size, pos = 4, srt = 0, col = pnas_filter$color[i]) 
    ### Add the pos 
    text(x = anot.pos+30 ,y = nrow(pnas_filter)-i,
         labels = pnas_filter$pos[i] ,cex = 1.2*txt.size, pos = 2, srt = 0) 
    ### Add the kmer 
    text(x=anot.pos+40 ,y = nrow(pnas_filter)-i,
         labels = pnas_filter$kmer[i] ,cex = 1.2*txt.size, pos = 2, srt = 0)
    ### Add similarity score
    text(x = anot.pos+50 ,y = nrow(pnas_filter)-i,
         labels = round(pnas_filter$SD[i],digits = 1) ,cex = 1.2*txt.size, pos = 2, srt = 0)
    
  }
  ### Add the axis
  text(x = c(anot.pos,anot.pos+10,anot.pos+20,anot.pos+30,anot.pos+40), y = (nrow(pnas_filter)+2)*c(1,1,1,1,1), labels = c("Gene","chr","position","kmer","Similarity score"),
       srt = 0,cex = 1.2*txt.size, pos = 4, col = pnas_filter$color[i], family = "Helvetica", font = 3)
  
  dev.off()
}

#barplot genes with a similarity score >10
if(T){

  if(T){  
    similarity=10
    
    pus_enzymes <-PB_SH2$Annotation[which(PB_SH2$SD>10)]
    
    pus_df <- data.frame("gene" = pus_enzymes, "PB1_TPM" = 0,"PB2_TPM" = 0,
                         "Diff1_TPM" =0, "Diff2_TPM" = 0,"Diff3_TPM" = 0,
                         "Undiff1_TPM" = 0,"Undiff2_TPM" = 0, "Undiff3_TPM" = 0)
    
    
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB1_TPM[i] <- TPM_PB1$TPM[which(TPM_PB1$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB2_TPM[i] <- TPM_PB2$TPM[which(TPM_PB2$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Diff1_TPM[i] <- TPM_Diff1$TPM[which(TPM_Diff1$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Diff2_TPM[i] <- TPM_Diff2$TPM[which(TPM_Diff2$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Diff3_TPM[i] <- TPM_Diff3$TPM[which(TPM_Diff3$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff1_TPM[i] <- TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff2_TPM[i] <- TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff3_TPM[i] <- TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == pus_df$gene[i])][1]
    }
    
    #### Add the standard deviation to the pus_df
    pus_df$Diff_SD <- 0
    pus_df$Undiff_SD <- 0
    pus_df$PB_SD <- 0
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB_SD[i] <- sd(c(pus_df$PB1_TPM[i], pus_df$PB2_TPM[i]))
      pus_df$Diff_SD[i] <- sd(c(pus_df$Diff1_TPM[i], pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
      pus_df$Undiff_SD[i] <- sd(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
    }
    
    #### Add the mean to the pus_df 
    pus_df$PB_mean <- 0
    pus_df$Diff_mean <- 0
    pus_df$Undiff_mean <- 0
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB_mean[i] <- mean(c(pus_df$PB1_TPM[i], pus_df$PB2_TPM[i]))
      pus_df$Diff_mean[i] <- mean(c(pus_df$Diff1_TPM[i], pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
      pus_df$Undiff_mean[i] <- mean(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
    }
    
    #### Remove the unwanted rows
    pus_df = pus_df[-c(12,13),]
    
    #pvalues, is Diff significantly different compared to undiff?
    for (i in c(1:nrow(pus_df))) {
      tPB<-t.test(c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i]),c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
      tDiff<-t.test(c(pus_df$Diff1_TPM[i],pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]),c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
      pus_df$pvaluePB[i] <- tPB$p.value
      pus_df$pvalueDiff[i] <- tDiff$p.value
    }
  }
  
  #### Plot 
  if(T){ 
    min.x = 0
    max.x = 15
    min.y = 0
    max.y = 80
    line.thickness = 1
    label.size = 1
    text.size = 0.8
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBvsDiffvsUndiff_barplot_fromHeatmap.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 7, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
    #box(lwd = line.thickness)
    #### Add the lines
    
    #for Undiff
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.2,i-.2),y=c(0,pus_df$Undiff_mean[i]),lwd=4 ,lend=1, lty=1, col = adjustcolor("black", alpha=0.7))
      pts=c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i])
      points(x=rep(i-.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    #for Diff
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i,i),y=c(0,pus_df$Diff_mean[i]),lwd=4 ,lend=1, lty=1, col = "dodgerblue3")
      pts=c(pus_df$Diff1_TPM[i],pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i])
      points(x=rep(i,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    #for PB
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i+.2,i+.2),y=c(0,pus_df$PB_mean[i]),lwd=4 ,lend=1, lty=1, col = "red")
      pts=c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i])
      points(x=rep(i+.2,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i,i),
            y=c(pus_df$Diff_mean[i]-pus_df$Diff_SD[i],pus_df$Diff_mean[i]+pus_df$Diff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.2,i-.2),
            y=c(pus_df$Undiff_mean[i]-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i+.2,i+.2),
            y=c(pus_df$PB_mean[i]-pus_df$PB_SD[i],pus_df$PB_mean[i]+pus_df$PB_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    
    #### Add Asterisks for p-values
    # *   < 0.05
    # **  < 0.01
    # *** < 0.001
    
    for (i in c(1:nrow(pus_df))){
      
      if (pus_df$pvaluePB[i]<0.001) {
        text(i+.2,pus_df$PB_mean[i]+pus_df$PB_SD[i]+10, labels = "***", col = "red3", cex = 0.6)
        # lines(x=c(i-0.2,i+0.2),y=c(pus_df$PB_mean[i]+pus_df$PB_SD[i]+8,pus_df$PB_mean[i]+pus_df$PB_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvaluePB[i]<0.01) {
        text(i+.2,pus_df$PB_mean[i]+pus_df$PB_SD[i]+10, labels = "**", col = "red3", cex = 0.6)
        # lines(x=c(i-0.2,i+0.2),y=c(pus_df$PB_mean[i]+pus_df$PB_SD[i]+8,pus_df$PB_mean[i]+pus_df$PB_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvaluePB[i]<0.05) {
        text(i+.2,pus_df$PB_mean[i]+pus_df$PB_SD[i]+10, labels = "*", col = "red3", cex = 0.6)
        # lines(x=c(i-0.2,i+0.2),y=c(pus_df$PB_mean[i]+pus_df$PB_SD[i]+8,pus_df$PB_mean[i]+pus_df$PB_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } 
      
      if (pus_df$pvalueDiff[i]<0.001) {
        text(i-.2,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "***", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalueDiff[i]<0.01) {
        text(i-.2,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "**", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalueDiff[i]<0.05) {
        text(i-.2,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "*", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } 
      
      
    }
    
    
    
    
    
    
    axis(side = 2,at = c(0:8)*10 ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    axis(side = 2,at = c(0:8)*10,labels = c(0:8)*10,lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
    axis(side = 1,at = c(-1:nrow(pus_df)+1) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
    axis(side = 1,at = c(1:nrow(pus_df)),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
    legend("topright", c("Undifferentiated","Differentiated","Lead treated"),
           cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
           pch=c(16,16),lty=c(0,0),
           lwd=c(2.5,2.5,2.5),col=c("black", "dodgerblue3","red"), bty = "n")
    mtext(text = "TPM",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    dev.off() 
  }
  
}

#barplot PUS enzymes
if (T) {
  
  if(T){
    #### List of the PUS enzyme genes:
    pus_enzymes <- c("PUS1", "PUSL1", "PUS3", "TRUB1", "TRUB2",
                     "DKC1", "PUS7", "PUS7L", "RPUSD1", "RPUSD2",
                     "PUS10", "PRUSD3", "PRUSD4")
    #### Make the data frame of TPM of the pus_enzymes 
    pus_df <- data.frame("gene" = pus_enzymes,"PB1_TPM" = 0,"PB2_TPM" = 0,
                         "Diff1_TPM" = 0,"Diff2_TPM" = 0,"Diff3_TPM" = 0,
                         "Undiff1_TPM" = 0,"Undiff2_TPM" = 0, "Undiff3_TPM" = 0)
    
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB1_TPM[i] <- TPM_PB1$TPM[which(TPM_PB1$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB2_TPM[i] <- TPM_PB2$TPM[which(TPM_PB2$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Diff1_TPM[i] <- TPM_Diff1$TPM[which(TPM_Diff1$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Diff2_TPM[i] <- TPM_Diff2$TPM[which(TPM_Diff2$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Diff3_TPM[i] <- TPM_Diff3$TPM[which(TPM_Diff3$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff1_TPM[i] <- TPM_Undiff1$TPM[which(TPM_Undiff1$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff2_TPM[i] <- TPM_Undiff2$TPM[which(TPM_Undiff2$Gene.Name == pus_df$gene[i])][1]
    }
    for (i in c(1:nrow(pus_df))) {
      pus_df$Undiff3_TPM[i] <- TPM_Undiff3$TPM[which(TPM_Undiff3$Gene.Name == pus_df$gene[i])][1]
    }
    
    #### Add the standard deviation to the pus_df
    pus_df$Diff_SD <- 0
    pus_df$Undiff_SD <- 0
    pus_df$PB_SD <- 0
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB_SD[i] <- sd(c(pus_df$PB1_TPM[i], pus_df$PB2_TPM[i]))
      pus_df$Diff_SD[i] <- sd(c(pus_df$Diff1_TPM[i], pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
      pus_df$Undiff_SD[i] <- sd(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
    }
    
    #### Add the mean to the pus_df 
    pus_df$PB_mean <- 0
    pus_df$Diff_mean <- 0
    pus_df$Undiff_mean <- 0
    for (i in c(1:nrow(pus_df))) {
      pus_df$PB_mean[i] <- mean(c(pus_df$PB1_TPM[i], pus_df$PB2_TPM[i]))
      pus_df$Diff_mean[i] <- mean(c(pus_df$Diff1_TPM[i], pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]))
      pus_df$Undiff_mean[i] <- mean(c(pus_df$Undiff1_TPM[i], pus_df$Undiff2_TPM[i], pus_df$Undiff3_TPM[i]))
    }
    
    #### Remove the unwanted rows
    pus_df = pus_df[-c(12,13),]
    
    #pvalues, is Diff significantly different compared to undiff?
    for (i in c(1:nrow(pus_df))) {
      tPB<-t.test(c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i]),c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
      tDiff<-t.test(c(pus_df$Diff1_TPM[i],pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i]),c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i]))
      pus_df$pvaluePB[i] <- tPB$p.value
      pus_df$pvalueDiff[i] <- tDiff$p.value
    }
  
  }  
  
  
  
  #### Plot 
  if(T){ 
    min.x = 0
    max.x = 15
    min.y = 0
    max.y = 50
    line.thickness = 1
    label.size = 1
    text.size = 1.4
    point.size = 1
    tck.length = 0.01
    tick.thickness = 1
    transparency = 0.2
    file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PBvsDiffvsUndiff_barplot_PUSenz_pvals_log.pdf"
    Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
    CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
    par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2) 
    plot(-100,-100,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
         ylab = NA,xlab=NA,ann = F,axes=F)
    #box(lwd = line.thickness)
    #### Add the lines
    
    #for Undiff
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.3,i-.3),y=c(0,pus_df$Undiff_mean[i]),lwd=4 ,lend=1, lty=1, col = adjustcolor("black", alpha=0.7))
      pts=c(pus_df$Undiff1_TPM[i],pus_df$Undiff2_TPM[i],pus_df$Undiff3_TPM[i])
      points(x=rep(i-.3,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    #for Diff
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i,i),y=c(0,pus_df$Diff_mean[i]),lwd=4 ,lend=1, lty=1, col = "dodgerblue3")
      pts=c(pus_df$Diff1_TPM[i],pus_df$Diff2_TPM[i],pus_df$Diff3_TPM[i])
      points(x=rep(i,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    #for PB
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i+.3,i+.3),y=c(0,pus_df$PB_mean[i]),lwd=4 ,lend=1, lty=1, col = "red")
      pts=c(pus_df$PB1_TPM[i],pus_df$PB2_TPM[i])
      points(x=rep(i+.3,length(pts)), y=pts, pch=20,col="black",cex=0.01)
    }
    
    

    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i,i),
            y=c(pus_df$Diff_mean[i]-pus_df$Diff_SD[i],pus_df$Diff_mean[i]+pus_df$Diff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i-.3,i-.3),
            y=c(pus_df$Undiff_mean[i]-pus_df$Undiff_SD[i],pus_df$Undiff_mean[i]+pus_df$Undiff_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    #### Add the SD line
    for (i in c(1:nrow(pus_df))) {
      lines(x=c(i+.3,i+.3),
            y=c(pus_df$PB_mean[i]-pus_df$PB_SD[i],pus_df$PB_mean[i]+pus_df$PB_SD[i]),lwd=line.thickness ,lend=1, lty=1, col = "black")
    }
    
    
    #### Add Asterisks for p-values
    # *   < 0.05
    # **  < 0.01
    # *** < 0.001
    
    for (i in c(1:nrow(pus_df))){
      
      if (pus_df$pvaluePB[i]<0.001) {
        text(i+.3,pus_df$PB_mean[i]+pus_df$PB_SD[i]+10, labels = "***", col = "red3", cex = 0.6)
        # lines(x=c(i-0.2,i+0.2),y=c(pus_df$PB_mean[i]+pus_df$PB_SD[i]+8,pus_df$PB_mean[i]+pus_df$PB_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvaluePB[i]<0.01) {
        text(i+.3,pus_df$PB_mean[i]+pus_df$PB_SD[i]+10, labels = "**", col = "red3", cex = 0.6)
        # lines(x=c(i-0.2,i+0.2),y=c(pus_df$PB_mean[i]+pus_df$PB_SD[i]+8,pus_df$PB_mean[i]+pus_df$PB_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvaluePB[i]<0.05) {
        text(i+.3,pus_df$PB_mean[i]+pus_df$PB_SD[i]+10, labels = "*", col = "red3", cex = 0.6)
        # lines(x=c(i-0.2,i+0.2),y=c(pus_df$PB_mean[i]+pus_df$PB_SD[i]+8,pus_df$PB_mean[i]+pus_df$PB_SD[i]+8),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }

      if (pus_df$pvalueDiff[i]<0.001) {
        text(i-.3,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "***", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalueDiff[i]<0.01) {
        text(i-.3,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "**", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      } else if (pus_df$pvalueDiff[i]<0.05) {
        text(i-.3,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+15, labels = "*", col = "dodgerblue3", cex = 0.6)
        # lines(x=c(i,i+0.2),y=c(pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13,pus_df$Diff_mean[i]+pus_df$Diff_SD[i]+13),lwd=line.thickness ,lend=1, lty=1, col = "black")
      }

      
    }
    
    
    
    
    
    
    axis(side = 2,at = c(0:5)*10 ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=0)
    axis(side = 2,at = c(0:5)*10,labels = c(0:5)*10,lwd.ticks = 0 ,cex.axis=label.size ,tck= 0 ,las=1 ,lwd = 0, line=-0.5 )
    axis(side = 1,at = c(-1:nrow(pus_df)+1) ,labels = NA ,lwd.ticks = tick.thickness ,cex.axis= label.size,
         tck=-tck.length,las=1,lwd=line.thickness, line=-0.8)
    axis(side = 1,at = c(1:nrow(pus_df)),labels = pus_df$gene,lwd.ticks = 0 ,cex.axis=0.5 ,tck= 0 ,las=2 ,lwd = 0, line=-1)
    legend("topright", c("Undifferentiated","Differentiated","Lead treated"),
           cex= text.size, pt.cex=point.size, y.intersp=1,x.intersp=0.01,
           pch=c(16,16),lty=c(0,0),
           lwd=c(2.5,2.5,2.5),col=c("black", "dodgerblue3","red"), bty = "n")
    mtext(text = "TPM",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
    dev.off() 
  }
  
  
}

#random model
if(T){
  
  PB_SH_sig<-PB_SH[which(PB_SH$N_value>sig),]
  PB_SH_sig$mm.PBminusUndiff <-   PB_SH_sig$mm.PB     - PB_SH_sig$mm.Undiff
  PB_SH_sig$mm.DiffminusUndiff <- PB_SH_sig$mm.Diff - PB_SH_sig$mm.Undiff
  
  #real data
  ######################################################### plots for TRUB1KDminusSCR vs PUS7KDminusSCR
  if(T){

    ### dotplot Diff vs PB
    if(T){
      min.x = -50
      max.x = 50
      min.y = -50
      max.y = 50
      line.thickness = 1
      label.size = 1
      text.size = 0.8
      point.size = 0.1
      tck.length = 0.01
      tick.thickness = 1
      transparency = 1
      
      file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/fig5LeadDiff_realmodel.pdf"
      Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
      CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
      par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
      
      plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           ylab = '',xlab='',ann = F,axes=T,type="l",
           xaxt='n', yaxt='n',
           bty='n', bg='transparent')
      
      points(y=PB_SH_sig$mm.PBminusUndiff[which(PB_SH_sig$color=="black")], x=PB_SH_sig$mm.DiffminusUndiff[which(PB_SH_sig$color=="black")], pch=16,cex=0.5, col="grey")
      points(y=PB_SH_sig$mm.PBminusUndiff[which(PB_SH_sig$color=="red" | PB_SH_sig$color=="orange") ], x=PB_SH_sig$mm.DiffminusUndiff[which(PB_SH_sig$color=="red" | PB_SH_sig$color=="orange")], pch=16,cex=0.5, col="orange")
      points(y=PB_SH_sig$mm.PBminusUndiff[which(PB_SH_sig$color=="blue" | PB_SH_sig$color=="blue2")], x=PB_SH_sig$mm.DiffminusUndiff[which(PB_SH_sig$color=="blue" | PB_SH_sig$color=="blue2")], pch=16,cex=0.5, col=adjustcolor("blue",alpha=0.5))
      
      
      abline(h=0, col="black")
      abline(v=0,col= "black")
      
      text(y=PB_SH_sig$mm.PBminusUndiff[which(PB_SH_sig$color=="orange")], x=PB_SH_sig$mm.DiffminusUndiff[which(PB_SH_sig$color=="orange")]+2,PB_SH_sig$Annotation[which(PB_SH_sig$color=="orange")],col='grey3', cex=0.3)
      text(y=PB_SH_sig$mm.PBminusUndiff[which(PB_SH_sig$color=="blue2")], x=PB_SH_sig$mm.DiffminusUndiff[which(PB_SH_sig$color=="blue2")]+2, PB_SH_sig$Annotation[which(PB_SH_sig$color=="blue2")],col='grey3', cex=0.3)
      
      #for plot KD-SCR vs SCR
      axis(side = 2, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
           lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
           las = 1, lwd = line.thickness, line = 0)
      axis(side = 1, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
           lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
           las = 1, lwd = line.thickness, line = 0)
      legend("topright", c("TRUB1 motif","PUS7 motif", "Other"),
             cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
             pch=c(16,16,16),lty=c(0,0,0),
             lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
      mtext(text = "Lead - Untreated (U-C MM%)",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
      mtext(text = "Differentiated - Untreated (U-C MM%)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
      ##ggsave(filename =file.name,width = 6,height = 6.6)
      alpha=0
      dev.off()
    }
    
    
    
  }
  
  #random
  ################################# Model for TRUB1KD and PUS7KD correlation
  if(T){

    #data creation, random model generation merged, null distribution, with mean and std of general kd libraries
    if(T){
      #working on all spots with sigma>2, this will variate according to mean and std of the general KD libraries
      SOI<-PB_SH_sig
      
      #mean and std deviations of the trub1kd and pus7kd libraries
      mPB=mean(PB_filter$mm.PBMINUSmm.Undiff)
      sdPB=sd(PB_filter$mm.PBMINUSmm.Undiff)
      mDiff=mean(SH_filter$mm.DiffminusUndiff)
      sdDiff=sd(SH_filter$mm.DiffminusUndiff)
      
      #generate the min-max values for kd
      min_deltavalue_PB<-min(SOI$mm.PBminusUndiff)
      max_deltavalue_PB<-max(SOI$mm.PBminusUndiff)
      min_deltavalue_Diff<-min(SOI$mm.DiffminusUndiff)
      max_deltavalue_Diff<-max(SOI$mm.DiffminusUndiff)
      
      delta_onlyPUS7KD<-c()
      delta_onlyTRUB1KD<-c()
      
      loop=100
      set.seed(12)
      #no mean and std limits
      random_numbers_PB_range<- runif(loop*round(nrow(SOI)/2),min=min_deltavalue_PB,max=max_deltavalue_PB)
      
      #with mean and std limits
      random_numbers_PB <- rnorm(loop*round(nrow(SOI)/2), mean=mPB,sd=sdPB)
     
      
      set.seed(123)
      #no mean and std limits
      random_numbers_Diff_range <- runif(loop*round(nrow(SOI)/2), min=min_deltavalue_Diff, max=max_deltavalue_Diff)
      
      #with mean and std limits
      random_numbers_Diff<- rnorm(loop*round(nrow(SOI)/2), mean = mDiff, sd = sdDiff)
     
      
      # random_numbers_PB<-0
      # random_numbers_Diff<-0
      delta_onlyPB<-c(random_numbers_PB,random_numbers_PB_range)
      delta_onlyDiff<-c(random_numbers_Diff,random_numbers_Diff_range)
      model<-data.frame("mm.PBminusUndiff"=delta_onlyPB,"mm.DiffminusUndiff"=delta_onlyDiff)
      # model<-model[-which(model$mm.TRUB1minusSCR>0 & model$mm.PUS7minusSCR>0),]
      # model<-model[-which(model$mm.TRUB1minusSCR>-10 & model$mm.PUS7minusSCR>-10),]
      write.csv(model,paste0("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PB_Diff_model_x",loop,".csv"),row.names = F)

    }
    
    #plot ttest 
    if(T){
      
      model_merged_sigma2_ranges_mean_sd<-read.csv("/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/PB_Diff_model_x1.csv",header=T)
      model<-model_merged_sigma2_ranges_mean_sd
      #plot models separately
      if(T){
        min.x = -50
        max.x = 50
        min.y = -50
        max.y = 50
        line.thickness = 1
        label.size = 1
        text.size = 0.8
        point.size = 0.1
        tck.length = 0.01
        tick.thickness = 1
        transparency = 1
        
        file.name = "/home/sasha/Rouhanifard lab Dropbox/Sasha/Figure1_Sepi/fig5PbvsDiffmodel_merged_sd_mean_sigma2_x1.pdf"
        Cairo::CairoFonts(regular = "Helvetica", bold = "Helvetica Bold", italic = "Helvetica Italic", bolditalic = "Helvetica Bold Italic")
        CairoPDF(file.name, width = 6, height = 6, family = "Helvetica")
        par(mfrow=c(1,1),lend=1)   #if you want for exampm 2X2 plot(4plots) in one mfrow-c(2,2)
        
        plot(0,0,xlim=c(min.x,max.x),ylim=c(min.y,max.y),
             ylab = '',xlab='',ann = F,axes=T,type="l",
             xaxt='n', yaxt='n',
             bty='n', bg='transparent')
        
        points(x=model$mm.DiffminusUndiff, y=model$mm.PBminusUndiff, pch=16,cex=0.5, col="grey")
        
        
        
        abline(h=0, col="black")
        abline(v=0,col= "black")
        
        
        #for plot KD-SCR vs SCR
        axis(side = 2, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
             lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
             las = 1, lwd = line.thickness, line = 0)
        axis(side = 1, at = seq(min.x, max.x, by = 20), labels = seq(min.x, max.x, by = 20),
             lwd.ticks = tick.thickness, cex.axis = label.size, tck = -tck.length,
             las = 1, lwd = line.thickness, line = 0)
        # legend("topright", c("TRUB1 motif","PUS7 motif", "Other"),
        #        cex= text.size, pt.cex=1, y.intersp=1,x.intersp=0.01,
        #        pch=c(16,16,16),lty=c(0,0,0),
        #        lwd=c(2.5,2.5,2.5),col=c("orange","blue", "grey"), bty = "n")
        mtext(text = "Lead - Untreated (U-C MM%)",side = 2,cex = text.size ,line = 2.3, family = "Helvetica")
        mtext(text = "Differentiated - Untreated (U-C MM%)",side = 1,cex = text.size ,line = 2.3, family = "Helvetica")
        ##ggsave(filename =file.name,width = 6,height = 6.6)
        alpha=0
        dev.off()
      }
    }
    
    #test
    if(T){
      #wilkoxon test
      if(T){
        #wilcox to check if my data are statistically different from the random data (ALL quadrants)
        wPB<-wilcox.test(PB_SH_sig$mm.PBminusUndiff,model$mm.PBminusUndiff)
        pPB<-wPB$p.value
        wDiff<-wilcox.test(PB_SH_sig$mm.DiffminusUndiff,model$mm.DiffminusUndiff)
        pDiff<-wDiff$p.value
        print(paste("whole dataset pValue PB:",pPB))
        print(paste("whole dataset pValue Diff:",pDiff))
        
        wPB<-wilcox.test(PB_SH_sig$mm.PBminusUndiff[which(PB_SH_sig$mm.PBminusUndiff<0 & PB_SH_sig$mm.DiffminusUndiff<0)],model$mm.PBminusUndiff[which(model$mm.PBminusUndiff<0 & model$mm.DiffminusUndiff<0)])
        pPB<-wPB$p.value
        wDiff<-wilcox.test(PB_SH_sig$mm.DiffminusUndiff[which(PB_SH_sig$mm.PBminusUndiff<0 & PB_SH_sig$mm.DiffminusUndiff<0)],model$mm.DiffminusUndiff[which(model$mm.PBminusUndiff<0 & model$mm.DiffminusUndiff<0)])
        pDiff<-wDiff$p.value
        print(paste("bottom-left pValue PB:",pPB))
        print(paste("bottom-left pValue Diff:",pDiff))
        
      }
      
      #density plots
      #real TRUB1 
      den<-density(PB_SH_sig$mm.TRUB1_KDminusSCR[which(PB_SH_sig$mm.TRUB1_KDminusSCR<0 & (PB_SH_sig$mm.PUS7_KD-PB_SH_sig$mm.SCR)<0)])
      plot(den, lwd=2,ylim=c(0,0.12),col="dodgerblue2",main="KD data",xlab = "TRUB1 KD - Scrambled (control) (% U-to_C error)")
      #random TRUB1
      den2<-density(model$mm.TRUB1minusSCR[which(model$mm.TRUB1minusSCR<0 & model$mm.PUS7minusSCR<0)])
      plot(den2,lwd=2, ylim=c(0,0.12),col="grey",main="Random model",xlab= "TRUB1 KD - Scrambled (control) (% U-to_C error)")
      #real PUS7
      den3<-density(PB_SH_sig$mm.PUS7_KDminusSCR[which(PB_SH_sig$mm.PUS7_KDminusSCR<0 & (PB_SH_sig$mm.TRUB1_KD-PB_SH_sig$mm.SCR)<0)])
      plot(den3, lwd=2,xlim=c(-40,3),ylim=c(0,0.12),col="dodgerblue4",main="KD data",xlab= "PUS7 KD - Scrambled (control) (% U-to_C error)")
      #random PUS7
      den4<-density(model$mm.PUS7minusSCR[which(model$mm.PUS7minusSCR<0 & model$mm.TRUB1minusSCR<0)])
      plot(den4, lwd=2,xlim=c(-40,3), col="grey",main="Random model",xlab= "PUS7 KD - Scrambled (control) (% U-to_C error)")
      
      
      # Example data (replace with your actual data)
      real_errors <- c(0.1, 0.2, 0.3, 0.4, 0.5)
      random_errors <- c(0.05, 0.15, 0.25, 0.35, 0.45)
      
      # Plot density of errors in TRUB1 KD for real dataset
      dens_real <- density(real_errors)
      plot(dens_real, xlim = c(0, max(max(real_errors), max(random_errors))), main = "Density Plot of TRUB1 KD Errors", xlab = "Error", ylab = "Density", col = "blue", lwd = 2)
      
      # Add density of errors in TRUB1 KD for random dataset
      dens_random <- density(random_errors)
      lines(dens_random, col = "red", lwd = 2)
      
      # Shade the area under the curve for the real dataset
      polygon(c(dens_real$x, max(dens_real$x)), c(dens_real$y, 0), col = "lightblue", border = NA)
      
      # Shade the area under the curve for the random dataset
      polygon(c(dens_random$x, max(dens_random$x)), c(dens_random$y, 0), col = "pink", border = NA)
      
      # Add a legend
      legend("topright", legend = c("Real Data", "Random Data"), col = c("blue", "red"), lty = 1, lwd = 2)
      
      
      
      #tailed test
      if(T){
        data <-TRUB1KDsoi$mm.TRUB1_KDminusSCR[which(TRUB1KDsoi$color=="red")]
        
        # Step 2: Perform the Wilcoxon signed-rank test
        # Null hypothesis: The data comes from a distribution with median equal to 0 (standard Gaussian distribution)
        # Alternative hypothesis: The data does not come from a distribution with median equal to 0
        
        # Perform the Wilcoxon signed-rank test
        wT=wilcox.test(data, mu = 0, alternative = "less")
        
        # Step 3: Analyze the direction of the data
        # The test statistic from the Wilcoxon signed-rank test gives information about the direction of the data.
        # If the test statistic is positive, it indicates that the data is shifted to the right (higher values).
        # If the test statistic is negative, it indicates that the data is shifted to the left (lower values).
        
        data <-PUS7KDsoi$mm.PUS7_KDminusSCR[which(PUS7KDsoi$color=="blue")]
        wP=wilcox.test(data, mu = 0, alternative = "less")
        
        print(paste("pValue trub1 tailed, real model is less than a gaussian distribution with mean=0:",wT$p.value))
        print(paste("pValue pus7 tailed, real model is less than a gaussian distribution with mean=0",wP$p.value))
        
      }
    }
  }
  
}








