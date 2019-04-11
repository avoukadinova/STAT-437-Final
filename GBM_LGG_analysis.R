# STAT 437
# April 9th, 2019
# https://www.abta.org/tumor_types/glioblastoma-gbm/

require(reshape2)
library(gdata)
library(varhandle)

raw_data1 <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/data_subset.xlsx")
data1 <- raw_data1[with(raw_data1, order(project_id)),]
data2 = as.data.frame.matrix(data1) 

GBM <- subset(data2, subset=(project_id =='TCGA-GBM'))
LGG <- subset(data2, subset=(project_id =='TCGA-LGG'))

TP53 <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "TP53")
TP53_GBM = subset(TP53, subset=(project_id =='TCGA-GBM'))
TP53_LGG = subset(TP53, subset=(project_id =='TCGA-LGG'))

TTN <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "TTN")
TTN_GBM = subset(TTN, subset=(project_id =='TCGA-GBM'))
TTN_LGG = subset(TTN, subset=(project_id =='TCGA-LGG'))

PTEN <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "PTEN")
PTEN_GBM = subset(PTEN, subset=(project_id =='TCGA-GBM'))
PTEN_LGG = subset(PTEN, subset=(project_id =='TCGA-LGG'))

EGFR <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "EGFR")
EGFR_GBM = subset(EGFR, subset=(project_id =='TCGA-GBM'))
EGFR_LGG = subset(EGFR, subset=(project_id =='TCGA-LGG'))

IDH1 <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "IDH1")
IDH1_GBM = subset(IDH1, subset=(project_id =='TCGA-GBM'))
IDH1_LGG = subset(IDH1, subset=(project_id =='TCGA-LGG'))

FLG <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "FLG")
FLG_GBM = subset(FLG, subset=(project_id =='TCGA-GBM'))
FLG_LGG = subset(FLG, subset=(project_id =='TCGA-LGG'))

ATRX <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "ATRX")
ATRX_GBM = subset(ATRX, subset=(project_id =='TCGA-GBM'))
ATRX_LGG = subset(ATRX, subset=(project_id =='TCGA-LGG'))

MUC16 <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "MUC16")
MUC16_GBM = subset(MUC16, subset=(project_id =='TCGA-GBM'))
MUC16_LGG = subset(MUC16, subset=(project_id =='TCGA-LGG'))

CIC <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "CIC")
CIC_GBM = subset(CIC, subset=(project_id =='TCGA-GBM'))
CIC_LGG = subset(CIC, subset=(project_id =='TCGA-LGG'))

genes_GBM = rbind(TP53_GBM, TTN_GBM, PTEN_GBM, EGFR_GBM, IDH1_GBM, FLG_GBM, ATRX_GBM, MUC16_GBM, CIC_GBM)
genes_LGG = rbind(TP53_LGG, TTN_LGG, PTEN_LGG, EGFR_LGG, IDH1_LGG, FLG_LGG, ATRX_LGG, MUC16_LGG, CIC_LGG)

index_GBM = c(rep(1,length(TP53_GBM$days_to_death)), rep(2,length(TTN_GBM$days_to_death)), rep(3,length(PTEN_GBM$days_to_death)), rep(4,length(EGFR_GBM$days_to_death)), rep(5,length(IDH1_GBM$days_to_death)), rep(6,length(FLG_GBM$days_to_death)), rep(7,length(ATRX_GBM$days_to_death)), rep(8,length(MUC16_GBM$days_to_death)), rep(9,length(CIC_GBM$days_to_death)))
index_LGG = c(rep(1,length(TP53_LGG$days_to_death)), rep(2,length(TTN_LGG$days_to_death)), rep(3,length(PTEN_LGG$days_to_death)), rep(4,length(EGFR_LGG$days_to_death)), rep(5,length(IDH1_LGG$days_to_death)), rep(6,length(FLG_LGG$days_to_death)), rep(7,length(ATRX_LGG$days_to_death)), rep(8,length(MUC16_LGG$days_to_death)), rep(9,length(CIC_LGG$days_to_death)))

genes <- unique(genes_GBM$gene_name)

days_to_death_GBM <- genes_GBM$days_to_death
days_to_death_LGG <- genes_LGG$days_to_death

boxplot(days_to_death_GBM~index_GBM, names=genes, xlab = "Gene Name", ylab = "Patient Survival Time (days)", main = "GBM")
abline(h = 450, col = "red")

boxplot(days_to_death_LGG~index_LGG, names=genes, xlab = "Gene Name", ylab = "Patient Survival Time (days)", main = "LGG")
abline(h = 2700, col = "red")
