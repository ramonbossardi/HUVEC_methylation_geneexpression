#Experiment 96 hrs recover

library("ChAMP")

#Local Mac
setwd("~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ")

myLoad <- champ.load(directory = "~/Analises R/r_array_huvec_120121/all_comparison_72hrs_96hrs_champ",arraytype = "EPIC", filterSNPs= FALSE, filterXY= FALSE)


CpG.GUI( arraytype="EPIC")

#Save specific data
#save(Anno, EPIC.manifest.hg19, multi.hit, myLoad, probe.features, bloodCtl, MatchGeneName, myCNA,
#     myDMP, myNorm, PathwayList, probeInfoALL.lv, file = "champ_huvec_120121_72hrs.RData")

#Save all 

save.image(file = "champ_huvec_022322_72hr_96wash_allsamples.RData")
load("champ_huvec_120121_72hrs_recover_96hrs_wash.RData")

#QC.GUI(beta=myLoad$beta)
#NOrmalization method of BMIQ
myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC")

head(myNorm)


#QC.GUI for BMIQ
QC.GUI(beta=myNorm,arraytype = "EPIC")


champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)
champ.SVD(beta=myNorm, resultsDir="./CHAMP_SVDimages_withoutnormalization/")

############ Batch Effect Correction
#This function implements the ComBat normalization method for BMIQ
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))
#Slide is pd file is numeric, which is not correct, so we firstly manually change it into charater
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
#Transform array as a factor
myLoad$pd$Array <- as.factor(myLoad$pd$Array)
champ.SVD(beta=myCombat, resultsDir=paste(getwd(), "resultsChamp", sep = "/"))
#principal components (PC-1) correlated to the covariates information provided
QC.GUI(beta=myCombat,arraytype = "EPIC")

write.csv(myCombat,file="./raw_result_72hrs_il6_96hrs_wash.csv")
#,quote=F,row.names = F


#DMP BMIQ
# To adjust the p value to identify = , adjPVal = 0.1 #compare.group=c("pPBS_72hrs", "pIL6_72hrs")
myDMP <- champ.DMP(beta=myCombat,adjPVal = 0.2,  arraytype = "EPIC", pheno=myLoad$pd$Sample_Group)
DMP.GUI()

DMP.GUI(DMP=myDMP[[1]],beta=myCombat,pheno=myLoad$pd$Sample_Group)
## myDMP is a list now, each data frame is stored as myDMP[[1]], myDMP[[2]], myDMP[[3]]
write_xls(myDMP[[1]], "DMP_IL6_72hrs.xlsx_withoutSNP.xlsx")
write.csv(myDMP[[1]],file="./DMP_IL6_72hrs.xlsx_withoutSNP.csv",quote=F,row.names = F)
write.csv(myDMP[[6]],file="./DMP_IL6_72hrs_vs_pbs_72hrs.csv")

head(myDMP[[6]])

#DMP for the PBC normalization results
myDMP_2 <- champ.DMP(beta = myCombat2,
                     pheno = myLoad$pd$Sample_Group,
                     compare.group = NULL,
                     adjPVal = 0.05,
                     adjust.method = "BH",
                     arraytype = "EPIC")
myDMP_2 <- champ.DMP(beta=myCombat2, arraytype = "EPIC")

head(myDMP_2[[1]])
write.csv2(myDMP[[1]],"/Users/pbanerjee/Documents/Payal_Scripts/R/Methylation_Champ/Methylation_idat_files/Group1vsGroup2/Results/DMP_850_grp1vsgrp2.csv" )
DMP.GUI(DMP=myDMP_2[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)

# DMP for SWAN normalization results
myDMP_3 <- champ.DMP(beta=myCombat3, arraytype = "EPIC")
DMP.GUI(DMP=myDMP_3[[1]],beta=myCombat3,pheno=myLoad$pd$Sample_Group)


#DMR
myDMR <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="Bumphunter", arraytype = "EPIC")
DMR.GUI(DMR=myDMR,arraytype="EPIC")
#,compare.group=c("PrEC_cells","LNCaP_cells"))

#Need to fix
myDMR_2 <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="DMRcate",arraytype="EPIC")

#<- champ.DMR(arraytype = "EPIC",method="DMRcate",cores=1)
#DMR.GUI(DMR=myDMR2,arraytype="EPIC")


#GSEA
myGSEA <- champ.GSEA(beta=myNorm,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
# myDMP and myDMR could (not must) be used directly.

head(myGSEA$DMP)
head(myGSEA$DMR)

##Block 
myBlock <- champ.Block(arraytype = "EPIC")
Block.GUI(arraytype="EPIC")

##myEpiMod <- champ.EpiMod(beta=myNorm,arraytype = "EPIC", resultsDir="./CHAMP_EpiMod/", PDFplot=TRUE)

myEpiMod <- champ.EpiMod(arraytype="EPIC")