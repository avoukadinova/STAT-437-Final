# STAT 437
# April 9th, 2019
# https://www.abta.org/tumor_types/glioblastoma-gbm/

require(reshape2)
library(gdata)
library(varhandle)
library(ggbiplot)
library(survival)
library(survminer)

raw_data1 <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/data_subset.xlsx")
data1 <- raw_data1[with(raw_data1, order(project_id)),]
data2 = as.data.frame.matrix(data1) 

GBM <- subset(data2, subset=(project_id =='TCGA-GBM'))
LGG <- subset(data2, subset=(project_id =='TCGA-LGG'))

TP53 <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "TP53")
TP53_GBM = subset(TP53, subset=(project_id =='TCGA-GBM'))
TP53_LGG = subset(TP53, subset=(project_id =='TCGA-LGG'))
wilcox.test(TP53_GBM$days_to_death,mu=381)
wilcox.test(TP53_LGG$days_to_death,mu=902.5)

TTN <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "TTN")
TTN_GBM = subset(TTN, subset=(project_id =='TCGA-GBM'))
TTN_LGG = subset(TTN, subset=(project_id =='TCGA-LGG'))
wilcox.test(TTN_GBM$days_to_death,mu=381)
wilcox.test(TTN_LGG$days_to_death,mu=902.5)

PTEN <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "PTEN")
PTEN_GBM = subset(PTEN, subset=(project_id =='TCGA-GBM'))
PTEN_LGG = subset(PTEN, subset=(project_id =='TCGA-LGG'))
wilcox.test(PTEN_GBM$days_to_death,mu=381)
wilcox.test(PTEN_LGG$days_to_death,mu=902.5)

EGFR <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "EGFR")
EGFR_GBM = subset(EGFR, subset=(project_id =='TCGA-GBM'))
EGFR_LGG = subset(EGFR, subset=(project_id =='TCGA-LGG'))
wilcox.test(EGFR_GBM$days_to_death,mu=381)
wilcox.test(EGFR_LGG$days_to_death,mu=902.5)

IDH1 <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "IDH1")
IDH1_GBM = subset(IDH1, subset=(project_id =='TCGA-GBM'))
IDH1_LGG = subset(IDH1, subset=(project_id =='TCGA-LGG'))
wilcox.test(IDH1_GBM$days_to_death,mu=381)
wilcox.test(IDH1_LGG$days_to_death,mu=902.5)

FLG <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "FLG")
FLG_GBM = subset(FLG, subset=(project_id =='TCGA-GBM'))
FLG_LGG = subset(FLG, subset=(project_id =='TCGA-LGG'))
wilcox.test(FLG_GBM$days_to_death,mu=381)
wilcox.test(FLG_LGG$days_to_death,mu=902.5)

ATRX <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "ATRX")
ATRX_GBM = subset(ATRX, subset=(project_id =='TCGA-GBM'))
ATRX_LGG = subset(ATRX, subset=(project_id =='TCGA-LGG'))
wilcox.test(ATRX_GBM$days_to_death,mu=381)
wilcox.test(ATRX_LGG$days_to_death,mu=902.5)

MUC16 <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "MUC16")
MUC16_GBM = subset(MUC16, subset=(project_id =='TCGA-GBM'))
MUC16_LGG = subset(MUC16, subset=(project_id =='TCGA-LGG'))
wilcox.test(MUC16_GBM$days_to_death,mu=381)
wilcox.test(MUC16_LGG$days_to_death,mu=902.5)

CIC <- read.xls("/Users/adivoukadinova/Desktop/Adi's Folder/Loyola Semester 2/STAT 437/Final/all_data.xls", sheet = "CIC")
CIC_GBM = subset(CIC, subset=(project_id =='TCGA-GBM'))
CIC_LGG = subset(CIC, subset=(project_id =='TCGA-LGG'))
wilcox.test(CIC_GBM$days_to_death,mu=381)
wilcox.test(CIC_LGG$days_to_death,mu=902.5)

genes_GBM = rbind(TP53_GBM, TTN_GBM, PTEN_GBM, EGFR_GBM, IDH1_GBM, FLG_GBM, ATRX_GBM, MUC16_GBM, CIC_GBM)
genes_LGG = rbind(TP53_LGG, TTN_LGG, PTEN_LGG, EGFR_LGG, IDH1_LGG, FLG_LGG, ATRX_LGG, MUC16_LGG, CIC_LGG)

index_GBM = c(rep(1,length(TP53_GBM$days_to_death)), rep(2,length(TTN_GBM$days_to_death)), rep(3,length(PTEN_GBM$days_to_death)), rep(4,length(EGFR_GBM$days_to_death)), rep(5,length(IDH1_GBM$days_to_death)), rep(6,length(FLG_GBM$days_to_death)), rep(7,length(ATRX_GBM$days_to_death)), rep(8,length(MUC16_GBM$days_to_death)), rep(9,length(CIC_GBM$days_to_death)))
index_LGG = c(rep(1,length(TP53_LGG$days_to_death)), rep(2,length(TTN_LGG$days_to_death)), rep(3,length(PTEN_LGG$days_to_death)), rep(4,length(EGFR_LGG$days_to_death)), rep(5,length(IDH1_LGG$days_to_death)), rep(6,length(FLG_LGG$days_to_death)), rep(7,length(ATRX_LGG$days_to_death)), rep(8,length(MUC16_LGG$days_to_death)), rep(9,length(CIC_LGG$days_to_death)))

gene_names <- unique(genes_GBM$gene_name)

days_to_death_GBM <- genes_GBM$days_to_death
days_to_death_LGG <- genes_LGG$days_to_death

boxplot(days_to_death_GBM~index_GBM, names=gene_names, xlab = "Gene Name", ylab = "Patient Survival Time (days)", main = "GBM")
abline(h = 381, col = "red",lwd = 2)

boxplot(days_to_death_LGG~index_LGG, names=gene_names, xlab = "Gene Name", ylab = "Patient Survival Time (days)", main = "LGG")
abline(h = 902.5, col = "red",lwd = 2)

pca <- prcomp(GBM[,10:18], center = TRUE,scale. = TRUE)
summary(pca)

print(ggbiplot(pca, obs.scale = 1, var.scale = 1,ellipse = TRUE, circle = TRUE, choices = c(1,2)))
print(ggbiplot(pca, obs.scale = 1, var.scale = 1,ellipse = TRUE, circle = TRUE, choices = c(1,3)))
print(ggbiplot(pca, obs.scale = 1, var.scale = 1,ellipse = TRUE, circle = TRUE, choices = c(1,4)))
print(ggbiplot(pca, obs.scale = 1, var.scale = 1,ellipse = TRUE, circle = TRUE, choices = c(1,5)))
print(ggbiplot(pca, obs.scale = 1, var.scale = 1,ellipse = TRUE, circle = TRUE, choices = c(1,6)))
print(ggbiplot(pca, obs.scale = 1, var.scale = 1,ellipse = TRUE, circle = TRUE, choices = c(1,7)))

(avg_GBM_IDH1 = mean(IDH1_GBM$days_to_death))
(avg_LGG_IDH1 = mean(IDH1_LGG$days_to_death))
(sd_GBM_IDH1 = sd(IDH1_GBM$days_to_death))
(sd_GBM_IDH1 = sd(IDH1_LGG$days_to_death))

(avg_GBM_ATRX = mean(ATRX_GBM$days_to_death))
(avg_LGG_ATRX = mean(ATRX_LGG$days_to_death))
(sd_GBM_ATRX = sd(ATRX_GBM$days_to_death))
(sd_LGG_ATRX = sd(ATRX_LGG$days_to_death))

(avg_GBM_CIC = mean(CIC_GBM$days_to_death))
(avg_LGG_CIC =  mean(CIC_LGG$days_to_death))
(sd_GBM_CIC = sd(CIC_GBM$days_to_death))
(sd_LGG_CIC = sd(CIC_LGG$days_to_death))

(avg_GBM_EGFR = mean(EGFR_GBM$days_to_death))
(avg_LGG_EGFR = mean(EGFR_LGG$days_to_death))
(sd_GBM_EGFR = sd(EGFR_GBM$days_to_death))
(sd_LGG_EGFR = sd(EGFR_LGG$days_to_death))

(avg_GBM_FLG = mean(FLG_GBM$days_to_death))
(avg_LGG_FLG = mean(FLG_LGG$days_to_death))
(sd_GBM_FLG = sd(FLG_GBM$days_to_death))
(sd_LGG_FLG = sd(FLG_LGG$days_to_death))

(avg_GBM_MUC16 = mean(MUC16_GBM$days_to_death))
(avg_LGG_MUC16 = mean(MUC16_LGG$days_to_death))
(sd_GBM_MUC16 = sd(MUC16_GBM$days_to_death))
(sd_LGG_MUC16 = sd(MUC16_LGG$days_to_death))

(avg_GBM_PTEN = mean(PTEN_GBM$days_to_death))
(avg_LGG_PTEN = mean(PTEN_LGG$days_to_death))
(sd_GBM_PTEN = sd(PTEN_GBM$days_to_death))
(sd_LGG_PTEN = sd(PTEN_LGG$days_to_death))

(avg_GBM_TP53 = mean(TP53_GBM$days_to_death))
(avg_LGG_TP53 = mean(TP53_LGG$days_to_death))
(sd_GBM_TP53 = sd(TP53_GBM$days_to_death))
(sd_LGG_TP53 = sd(TP53_LGG$days_to_death))

(avg_GBM_TTN = mean(TTN_GBM$days_to_death))
(avg_LGG_TTN = mean(TTN_LGG$days_to_death))
(sd_GBM_TTN = sd(TTN_GBM$days_to_death))
(sd_LGG_TTN = sd(TTN_LGG$days_to_death))

fit <- coxph(Surv(days_to_death) ~ IDH1 + ATRX + CIC + EGFR + FLG + MUC16 + TP53 + PTEN + TTN , data = GBM)
summary(fit)

fit2 <- step(fit, data = GBM)
summary(fit2)

ggsurvplot(fit = survfit(fit2), data = GBM, palette = "#2E9FDF", 
           ggtheme = theme_minimal(), title = "GBM",
           font.x = 18, font.y = 18, font.ticks = 18)

fit2$formula
fit2$coefficients

fit3 <- coxph(Surv(days_to_death) ~ IDH1 + ATRX + CIC + EGFR + FLG + MUC16 + TP53 + PTEN + TTN , data = LGG)
summary(fit3)

fit4 <- step(fit3, data = LGG)
summary(fit4)

ggsurvplot(fit = survfit(fit4), data = LGG, palette = "#2E9FDF", 
           ggtheme = theme_minimal(), title = "LGG",
           font.x = 18, font.y = 18, font.ticks = 18)
fit4$formula
fit4$coefficients
