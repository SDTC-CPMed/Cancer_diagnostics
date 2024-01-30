library(ggplot2)
library(ggpubr)


######################## Boxplot ###################################

Boxplot_data = read.csv("boxplot_data_cancer_diagnosis.csv", sep=",", header=T, stringsAsFactors = F, row.names = 1) 

Boxplot_data_genomics = read.csv("graph_cancer_diagnosis_genomics.csv", sep=",", header=T, stringsAsFactors = F, row.names = 1)

Boxplot_data = rbind(Boxplot_data, Boxplot_data_genomics)

Boxplot_data = Boxplot_data[Boxplot_data$Model == 'prevalent',]
Boxplot_data = Boxplot_data[Boxplot_data$control_data == 'HC',]

Boxplot_data_CPTAC = read.csv("boxplot_data_cancer_diagnosis_CPTAC.csv", sep=",", header=T, stringsAsFactors = F, row.names = 1) 
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Primary Tumor"]<- 'Tumor'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Solid Tissue Normal"]<- 'Normal'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Breast PDC000120"]<- 'Breast'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Colon PDC000116"]<- 'Colon'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="HeadNeck PDC000221"]<- 'HeadNeck'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Kidney PDC000127"]<- 'Kidney'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Liver PDC000198"]<- 'Liver'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Lung PDC000153"]<- 'Lung'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Ovary PDC000110"]<- 'Ovarian'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="Pancreas PDC000270"]<- 'Pancreas'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="UCEC PDC000125"]<- 'Uterine'
Boxplot_data_CPTAC[Boxplot_data_CPTAC=="tissue"]<- 'CPTAC'

Boxplot_data_CPTAC$prob = 1 - Boxplot_data_CPTAC$prob
colnames(Boxplot_data_CPTAC)[2] = 'Legend'

#Boxplot_data = rbind(Boxplot_data, Boxplot_data_CPTAC)


Boxplot_data[Boxplot_data=="incident"]<- 'Incident'

Boxplot_data[Boxplot_data=="prevalent"]<- 'Prevalent'
Boxplot_data[Boxplot_data=="obesity"]<- 'Obesity'
Boxplot_data[Boxplot_data=="met"]<- 'Metabolomics'
Boxplot_data[Boxplot_data=="prot"]<- 'Proteomics'
Boxplot_data[Boxplot_data=="all"]<- 'Genomics'
Boxplot_data[Boxplot_data=="Malignant_melanoma"]<- 'Melanoma'
Boxplot_data[Boxplot_data=="AllCancers"]<- 'All cancers'
Boxplot_data[Boxplot_data=="Leukaemia"]<- 'Leukemia'

#Boxplot_data[Boxplot_data=="HC"]<- 'Non-Cancer'


colnames(Boxplot_data)[2] = 'Legend'

Boxplot_data = Boxplot_data[Boxplot_data$disease != c('All cancers'),]
Boxplot_data = Boxplot_data[Boxplot_data$data_type %in% 'Genomics',]


pd <- position_dodge(0.5)
TH2_COLOR <- c("#A9CD75","orange", '#005292', 'yellow')

my_comparisons <- list( c('HC', 'Prevalent') )



labels_AUC<-  paste(round(as.numeric(Boxplot_data$AUC),2))
Ploter_X <- ggplot(data=Boxplot_data, aes(x=Legend, y=prob, fill=Legend)) +
  #geom_boxplot(outlier.size = 0.2, width = 0.5) +
  geom_violin(width = 0.5, draw_quantiles = c(0.5)) +
  theme(text=element_text(family="Arial Narrow", size=12))+
  scale_fill_manual(values= TH2_COLOR)+
  labs( x = 'Disease', y="Probability to develop the disease")+
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.99, 0.99), method = "t.test", label = "p.signif")+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme_bw()+  
  ylim(0.0, 1.25) 
Ploter_X + facet_grid(data_type ~ disease ) +
  theme(strip.text = element_text(face = 'bold')) 




labels_AUC<-  paste(round(as.numeric(Boxplot_data_CPTAC$AUC),2))
Ploter_X <- ggplot(data=Boxplot_data_CPTAC, aes(x=Legend, y=prob, fill=Legend)) +
  geom_boxplot(outlier.size = 0.2, width = 0.5) +
  #geom_violin(width = 0.5, draw_quantiles = c(0.5)) +
  theme(text=element_text(family="Arial Narrow", size=12))+
  scale_fill_manual(values= TH2_COLOR)+
  labs( x = 'Disease', y="Probability to develop the disease")+
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.99, 0.99), method = "t.test", label = "p.signif")+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme_bw()+
  ylim(0.0, 1.25) 


Ploter_X + facet_grid(data_type ~ disease) +
  theme(strip.text = element_text(face = 'bold'))




######################## line plot #################################################

AUC_data = read.csv('ShinyApp/ShinyData_cancer_diagnosis.csv', row.names = 1)


AUC_data[AUC_data=="obesity"]<- 'Obesity'
#AUC_data[AUC_data=="genomics"]<- 'Genomics'
AUC_data[AUC_data=="metabolomics"]<- 'Metabolomics'
AUC_data[AUC_data=="proteomics"]<- 'Proteomics'
AUC_data[AUC_data=="met"]<- 'Metabolomics'
AUC_data[AUC_data=="prot"]<- 'Proteomics'
AUC_data[AUC_data=="Atherosclerosis"]<- 'ASVD'
AUC_data[AUC_data=="Malignant_melanoma"]<- 'Melanoma'
AUC_data[AUC_data=="AllCancers"]<- 'All cancers'


AUC_data = AUC_data[AUC_data$N < 16,]

AUC_data_HC = AUC_data[AUC_data['Model'] == 'prevalent' & AUC_data['control_data'] == 'HC',]
colnames(AUC_data_HC)[6] = 'AUC_train_HC'
colnames(AUC_data_HC)[7] = 'AUC_test_HC'

AUC_data_NonCancer = AUC_data[AUC_data['Model'] == 'prevalent' & AUC_data['control_data'] == 'NonCancer',]
colnames(AUC_data_NonCancer)[6] = 'AUC_train_NonCancer'
colnames(AUC_data_NonCancer)[7] = 'AUC_test_NonCancer'

AUC_data = cbind(AUC_data_HC, AUC_data_NonCancer)
AUC_data = AUC_data[,colnames(AUC_data)[c(1:11, 17,18)]]

AUC_data = AUC_data[!AUC_data$Disease %in% c('All cancers'),]


pd <- position_dodge(0.5)
Ploter_X <- ggplot(data=AUC_data) +
  geom_line(aes(x = N, y = AUC_test_HC, color = "HC"), linetype = 1)+
  geom_point(aes(x = N, y = AUC_test_HC, color = "HC"), size = 1)+
  geom_line(aes(x = N, y = AUC_test_NonCancer, color = "NonCancer"), linetype = 1)+
  geom_point(aes(x = N, y = AUC_test_NonCancer, color = "NonCancer"), size = 1)+
  

  theme(text=element_text(family="Arial Narrow", size=12))+
  scale_color_manual(name = "Model", values = c("HC" = "orange", 'NonCancer' = 'blue'))+
  labs(x="# Features", y = "AUC")+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  scale_x_continuous(breaks = c(3, 6, 9, 12, 15)) +
  
  theme_bw()


Model_Comparision = Ploter_X + facet_grid(Data ~ Disease) +
  theme(legend.position="top")+
  theme(strip.text = element_text(face = 'bold')) 
Model_Comparision
