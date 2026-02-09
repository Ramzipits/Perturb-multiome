datalnc_BRAIN_RCC<-read.table(file='lnc_matrix_brain_RCC.txt',header=TRUE)
dp <- datalnc_BRAIN_RCC[complete.cases(datalnc_BRAIN_RCC),]
dp<- datalnc_BRAIN_RCC
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)
p <- ggplot(data=dp, aes(x=BRAIN_AGING_level), y=-log10(BRAIN_AGING_pvalue_RCC)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=2, col="red") + geom_hline(yintercept=-log10(0.05), col="red")
dp$correlation <- "NO"
dp$correlation[dp$BRAIN_AGING_level >2 & dp$BRAIN_AGING_pvalue_RCC<0.05] <- "Related"
dp$correlation[dp$BRAIN_AGING_level <2 | dp$BRAIN_AGING_pvalue_RCC>0.05] <- "Unrelated"

list_brain <- subset(dp,dp$BRAIN_AGING_level > 2 & dp$BRAIN_AGING_pvalue_RCC < 0.05)
write.table(list_brain,"brain_RCC_lnclist.txt")



dp$dplabel <- NA
dp$dplabel[dp$correlation != "Unrelated"] <- dp$BRAIN_AGING_genename[dp$correlation != "Unrelated"]
ggplot(data=dp, aes(x=BRAIN_AGING_level, y=-log10(BRAIN_AGING_pvalue_RCC), col=correlation, label=dplabel)) + geom_point() + theme_minimal() + geom_text()


brain_fdr = p.adjust(dp$BRAIN_AGING_pvalue_RCC, method ="BH")
dp =cbind(dp,brain_fdr)

ggplot(data=dp, aes(x=BRAIN_AGING_level, y=-log10(BRAIN_AGING_pvalue_RCC), col=correlation, label=dplabel2)) +
  scale_x_continuous(breaks =c(0,3,6,9,12,15))+
  geom_point(aes(shape=correlation, color=correlation, size=correlation))+
  scale_shape_manual(values=c(16, 16))+
  scale_color_manual(values=c('#999999','#E69F00'))+
  scale_size_manual(values=c(1.4,1))+
  scale_color_gradient('set1', low = "#cfd3ed", high = "#4954c1") + 
  theme_minimal() +
  geom_text_repel(segment.color='black',segment.size=0.3,
                  segment.alpha=0.4,max.overlaps =30,color='black',size=2.4) +
  scale_color_manual(values=c("orange", "gray")) + 
  geom_vline(xintercept=2, linetype=2,col="black") +
  geom_hline(yintercept=-log10(0.05),linetype=2,col="black") +
  labs(x="log2(fpkm+1)",
       y="-log10 (p-value of RCC)", title="lnc RNA's correlation with aging in BRAIN") + 
  theme(axis.title.x = element_text(size = 13,face='bold'),
        axis.title.y = element_text(size = 13,face='bold'),
        axis.text=element_text(size=11,face='bold'),
        legend.title = element_text(color = "black", size = 12,face='bold'),
        legend.text = element_text(color = "black", size = 10,face='bold'),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(panel.grid=element_blank())


###blood
datalnc_BLOOD_RCC<-read.table(file='lnc_matrix_blood_RCC.txt',header=TRUE)
dp <- datalnc_BLOOD_RCC[complete.cases(datalnc_BLOOD_RCC),]
dp<- datalnc_BLOOD_RCC
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)
p <- ggplot(data=dp, aes(x=BLOOD_AGING_level), y=-log10(BLOOD_AGING_pvalue_RCC)) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=2, col="red") + geom_hline(yintercept=-log10(0.05), col="red")
dp$correlation <- "NO"
dp$correlation[dp$BLOOD_AGING_level >2 & dp$BLOOD_AGING_pvalue_RCC<0.05] <- "Related"
dp$correlation[dp$BLOOD_AGING_level <2 | dp$BLOOD_AGING_pvalue_RCC>0.05] <- "Unrelated"

list_BLOOD <- subset(dp,dp$BLOOD_AGING_level > 2 & dp$BLOOD_AGING_pvalue_RCC < 0.05)

dp$dplabel <- NA
dp$dplabel[dp$correlation != "Unrelated"] <- dp$BLOOD_AGING_genename[dp$correlation != "Unrelated"]
ggplot(data=dp, aes(x=BLOOD_AGING_level, y=-log10(BLOOD_AGING_pvalue_RCC), col=correlation, label=dplabel)) + geom_point() + theme_minimal() + geom_text()


BLOOD_fdr = p.adjust(dp$BLOOD_AGING_pvalue_RCC, method ="BH")
dp =cbind(dp,BLOOD_fdr)

ggplot(data=dp, aes(x=BLOOD_AGING_level, y=-log10(BLOOD_AGING_pvalue_RCC), col=correlation, label=dplabel2)) +
  scale_x_continuous(breaks =c(0,3,6,9,12,15))+
  geom_point(aes(shape=correlation, color=correlation, size=correlation))+
  scale_shape_manual(values=c(16, 16))+
  scale_color_manual(values=c('#999999','#E69F00'))+
  scale_size_manual(values=c(1.4,1))+
  scale_color_gradient('set1', low = "#cfd3ed", high = "#4954c1") + 
  theme_minimal() +
  geom_text_repel(segment.color='black',segment.size=0.3,
                  segment.alpha=0.4,max.overlaps =30,color='black',size=2.4) +
  scale_color_manual(values=c("orange", "gray")) + 
  geom_vline(xintercept=2, linetype=2,col="black") +
  geom_hline(yintercept=-log10(0.05),linetype=2,col="black") +
  labs(x="log2(fpkm+1)",
       y="-log10 (p-value of RCC)", title="lnc RNA's correlation with aging in BLOOD") + 
  theme(axis.title.x = element_text(size = 13,face='bold'),
        axis.title.y = element_text(size = 13,face='bold'),
        axis.text=element_text(size=11,face='bold'),
        legend.title = element_text(color = "black", size = 12,face='bold'),
        legend.text = element_text(color = "black", size = 10,face='bold'),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(panel.grid=element_blank())


  