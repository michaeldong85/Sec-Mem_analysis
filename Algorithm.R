####This R script is designed to screen for Sec-Mem interactions between
####one subgroup from Group 1 and another subgroup from Group 2
####In this script, 6 subgroups from Epithelial group and 4 subgroups from 
####CAF group are used as examples
####Sec-Mem matrix between secreted proteins from subgroups of epithelial cells
####and membrane proteins from subgroups of CAF cells is generated and stored 
####in RESULT and output as output.csv

getwd()

####Data input (Epithelial/CAF Expression Matrix and Subgroup information)
EXP_Epi <- read.csv("Epithelial_expression_matrix.csv", header = T) ####Epithelial expression matrix
rownames(EXP_Epi) = EXP_Epi[,1]
EXP_Epi <- EXP_Epi[,2:ncol(EXP_Epi)]
colnames(EXP_Epi) <- gsub("\\.","-",colnames(EXP_Epi))
INFO_Epi <- read.csv("Epithelial_subgroup_info.csv", header = T) ####Epithelial subgroup information
rownames(INFO_Epi) = INFO_Epi[,1]

EXP_CAF <- read.csv("CAF_expression_matrix.csv", header = T) ####CAF expression matrix
rownames(EXP_CAF) = EXP_CAF[,1]
EXP_CAF <- EXP_CAF[,2:ncol(EXP_CAF)]
colnames(EXP_CAF) <- gsub("\\.","-",colnames(EXP_CAF))
INFO_CAF <- read.csv("CAF_subgroup_info.csv", header = T) ####CAF subgroup information
rownames(INFO_CAF) = INFO_CAF[,1]

####Sec-Mem information input
Mem <- read.csv("Membrane_gene_list.csv", header = T) ####Proteins that locates on plasma membranes
Mem <- Mem$x
Sec <- read.csv("Secreted_gene_list.csv", header = T) ####Proteins that are secreted into blood
Sec <- Sec$x
PPI <- read.csv("Secreted_Membrane_PPI_list.csv", header = T) ####PPIs in the order of "Secreted, Membrane"

####Cutoff input
Nexp_cutoff = 0.6;               ###Expression cutoff     (detailed in methods)
pos_cutoff = 0.2;         ###positive ratio cutoff (detailed in methods)


####Step 1: normalize both expression matrix
NEXP_Epi <- scale(EXP_Epi)
NEXP_CAF <- scale(EXP_CAF)

####Step 2: generate "Sec-Mem" Matrix between Subgroups of Epithelial cells and Subgroups of CAF cells
RESULT <- data.frame()   ####Epithelial as secreted
                         ####CAF as membrane receptors    
                         ####Results output file
for (i in 1:nrow(PPI)){
  sec_gene = PPI$Secreted[i]; mem_gene = PPI$Membrane[i];
  ####Epithelail subgroups are "Pro"   "Neck"  "Pari"  "Pit"   "Endo"  "Chief"
  Sub1 = c("Pro","Neck","Pari","Pit","Endo","Chief")
  ####CAF subgroups are "iCAF" "mCAF" "SMC"  "cCAF"
  Sub2 = c("iCAF","mCAF","SMC","cCAF")
  
  Value <- list();
  for (num1 in 1:length(Sub1)){
    Value[num1] = mean(NEXP_Epi[sec_gene,rownames(INFO_Epi[which(INFO_Epi$Subgroup == Sub1[num1]),])])
  }
  Value <- unlist(Value)
  
  for (num1 in 1:length(Sub1)){
    ####Secreted gene normalized expression value in all Subgroups of Epithelial cells
    Value1 = mean(NEXP_Epi[sec_gene,rownames(INFO_Epi[which(INFO_Epi$Subgroup == Sub1[num1]),])])
    ####Ratio of cells expressing sec_gene in Subgroup1 of Epithelial cells
    Ratio1 = sum(EXP_Epi[sec_gene,rownames(INFO_Epi[which(INFO_Epi$Subgroup == Sub1[num1]),])] > 0)/length(rownames(INFO_Epi[which(INFO_Epi$Subgroup == Sub1[num1]),]))
    
    for (num2 in 1:length(Sub2)){
      ####Membrane gene normalized expression value in Subgroup2 of CAF cells
      Value2 = mean(NEXP_CAF[mem_gene,rownames(INFO_CAF[which(INFO_CAF$Subgroup == Sub2[num2]),])])
      ####Ratio of cells expressing mem_gene in Subgroup2 of CAF cells
      Ratio2 = sum(NEXP_CAF[mem_gene,rownames(INFO_CAF[which(INFO_CAF$Subgroup == Sub2[num2]),])] > 0)/length(rownames(INFO_CAF[which(INFO_CAF$Subgroup == Sub2[num2]),]))
      
      if ((Value1 > Nexp_cutoff)  ####Expression value in Subgroup1 of Epithelial cells over 0.6
          & (Value1 == max(Value)) #####Expression value in Subgroup1 is the highest among all subgroups of Epithelial cells
          & (Value2 > Nexp_cutoff) 
          & (Ratio1 > pos_cutoff) 
          & (Ratio2 > pos_cutoff)){
        RESULT <- rbind(RESULT, c(sec_gene, mem_gene,Sub1[num1],"Epithelial",Sub2[num2],"CAF"))
      }
      
    }
  }
}

colnames(RESULT) = c("Secreted","Membrane","Subgroup1","Epithelial","Subgroup2","CAF")

write.csv(RESULT, file = "output.csv", quote = F)












