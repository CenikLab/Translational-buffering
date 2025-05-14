library(data.table)
library(dplyr)
library(ggplot2)
library(MatchIt)
library(tidyr)
library(readxl)
library(ggrepel)
library(ggpubr)
library(seqinr)
library(cubar)
library(Biostrings)
library(cowplot)
library(ribor)
library(stringr)
library(limma)
library(edgeR)
library(treemap)
library(patchwork)
library(zCompositions)
library(foreach)
library(doParallel)


# here TB_score=1/ TB1 refers to "TB high"
#TB_score=2/ TB2 refers to "TB low"
#TB_core=3/TB3 refers to "Others"

# table 
# to find correlation convert 


RNA_clr=read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/RNA_clr.csv")
Ribo_clr=read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/Ribo_clr.csv")
TE_GBM_clr=read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/TE_GBM_clr.csv")
dim(RNA_clr)
dim(Ribo_clr)
dim(TE_GBM_clr)

RNA_clr[1:5,1:5]
RNA_clr=RNA_clr[,-1]
Ribo_clr=Ribo_clr[,-1]
TE_GBM_clr= TE_GBM_clr[,-1]

RNA_clr[1:5,1:5]
Ribo_clr[1:5,1:5]
TE_GBM_clr[1:5,1:5]
RNA_clr=as.data.table(RNA_clr)
Ribo_clr=as.data.table(Ribo_clr)
TE_GBM_clr=as.data.table(TE_GBM_clr)


#calculate the correlation


numeric_RNA_nond<-RNA_clr[, apply(.SD,1, as.numeric),.SDcols =-"transcript"]
numeric_RNA_nond[1:5,1:5]

numeric_TE_nond <-TE_GBM_clr[, apply(.SD,1, as.numeric),.SDcols =-"transcript"]
numeric_TE_nond[1:5,1:5]
# in the below code each column is corelated to each other column such that the corleation between each geneis stored as diagonal
TE_RNA_cor_nond= cor(numeric_TE_nond, numeric_RNA_nond, method = "spearman")
TE_RNA_cor_value_nond= diag(TE_RNA_cor_nond)
TE_RNA_sample_corl_nond= data.table(RNA_clr$transcript, TE_RNA_cor_value_nond)
colnames(TE_RNA_sample_corl_nond)= c("transcript","TE_RNA_cor_value_nond")
median_TE_RNA_sample_corl_nond = median(TE_RNA_sample_corl_nond$TE_RNA_cor_value_nond)


TE_RNA_nond=ggplot(TE_RNA_sample_corl_nond, aes(x = TE_RNA_cor_value_nond)) +
  geom_histogram(binwidth = 0.02, fill = "#B7E6A5",alpha= 0.5)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines")) +
  geom_vline(xintercept = c(median_TE_RNA_sample_corl_nond), linetype = "dashed", color = c("blue"), linewidth = 1)+
  scale_x_continuous(name = "Spearman Correlation  coeffcient", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(name = "Frequency", breaks = scales::pretty_breaks(n = 10),expand = c(0, 0))+coord_cartesian(ylim = c(0, 400), xlim=c(-0.9
                                                                                                                                           , 0.65))
TE_RNA_nond


numeric_Ribo_nond <-Ribo_clr[, apply(.SD,1, as.numeric),.SDcols =-"transcript"]
Ribo_RNA_cor_value_nond= cor(numeric_Ribo_nond, numeric_RNA_nond, method = "spearman")
Ribo_RNA_cor_value_nond= diag(Ribo_RNA_cor_value_nond)
length(Ribo_RNA_cor_value_nond)
Ribo_RNA_cor_sample_nond=as.data.table(Ribo_RNA_cor_value_nond)
Ribo_RNA_cor_sample_nond$transcript= RNA_clr$transcript
colnames(Ribo_RNA_cor_sample_nond)= c("Ribo_RNA_cor_value_nond","transcript")

median_Ribo_RNA_cor_sample_nond = median(Ribo_RNA_cor_sample_nond$Ribo_RNA_cor_value_nond)
Ribo_RNA_nond_plot=ggplot(Ribo_RNA_cor_sample_nond, aes(x = Ribo_RNA_cor_value_nond)) +
  geom_histogram(binwidth = 0.02, fill = "#B7E6A5",alpha= 0.5)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines")) +
  geom_vline(xintercept = c(median_Ribo_RNA_cor_sample_nond), linetype = "dashed", color = c("blue"), linewidth = 1)+
  scale_x_continuous(name = "Spearman Correlation  coeffcient", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(name = "Frequency", breaks = scales::pretty_breaks(n = 10),expand = c(0, 0)) +coord_cartesian(ylim = c(0, 450), xlim=c(-0.1
                                                                                                                                                              , 1))



Ribo_RNA_nond_plot


MAD_RNA1= apply(numeric_RNA_nond, 2, mad)
MAD_Ribo1= apply(numeric_Ribo_nond, 2, mad)
length(MAD_RNA1)
length(MAD_Ribo1)
MAD_RNA_Ribo1= cbind(MAD_RNA1, MAD_Ribo1, RNA_clr$transcript)
MAD_RNA_Ribo1=as.data.table(MAD_RNA_Ribo1)
colnames(MAD_RNA_Ribo1)= c("MAD_RNA","MAD_Ribo","transcript")
#convert to numeric
library(sp)
MAD_RNA_Ribo1=MAD_RNA_Ribo1%>%
  mutate_at(vars(-all_of(c("transcript"))), as.numeric)
dim(MAD_RNA_Ribo1)
tail(MAD_RNA_Ribo1,20)


MAD_RNA_Ribo_ratio= mutate(MAD_RNA_Ribo1, MAD_ratio= MAD_Ribo/MAD_RNA)
MAD_RNA_Ribo_ratio=MAD_RNA_Ribo_ratio[order(MAD_RNA_Ribo_ratio$MAD_ratio)]
head(MAD_RNA_Ribo_ratio[1200:1300,],100)
dim(MAD_RNA_Ribo_ratio)
#merge other parameters. 
Buffering_human= merge(TE_RNA_sample_corl_nond, MAD_RNA_Ribo_ratio, by.x= "transcript", by.y= "transcript")
Buffering_human



Buffering_human=merge(Buffering_human, Ribo_RNA_cor_sample_nond, by.x= "transcript", by.y= "transcript")
# 
# the idea is to rank each criteria
# Give weight to each criteria. Multiple the rank with the weight and sum it up . Rank according to score again (lowest to highest)
#remove dummy genes
dim(Buffering_human)

# rank according to the parameter $# Multiply each parameter by 2,1, 0. 5 and  0.5 then take the sum
Buffering_human_rank= Buffering_human %>% mutate(
  TE_RNA_rank = rank(TE_RNA_cor_value_nond, ties.method = "first"),
  MAD_rank = rank(MAD_ratio, ties.method = "first")) %>%
  mutate(
    TE_RNA_rank_adj = TE_RNA_rank * 1,
    MAD_rank_adj = MAD_rank * 1) %>%
  rowwise() %>%
  mutate(Row_Sum = sum(c(TE_RNA_rank_adj,  MAD_rank_adj))) %>%
  ungroup() %>%
  mutate(TB_rank=rank(Row_Sum, ties.method="first")) %>%
  arrange(TB_rank)

as.list(Buffering_human_rank[251:500,1])

# 
Top2_buffered_human= subset(Buffering_human_rank, Buffering_human_rank$TB_rank > 250 & Buffering_human_rank$TB_rank <= 500) %>% arrange(TB_rank)
Top_buffered_human= subset(Buffering_human_rank, Buffering_human_rank$TB_rank > 0 & Buffering_human_rank$TB_rank <= 250) %>% arrange(TB_rank)
dim(Top2_buffered_human)
# add TB score if present in Top_buffered genes
Buffering_human_rank <- Buffering_human_rank %>%
  mutate(Buffering = if_else(transcript %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(transcript %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))




head(Top_buffered_human,30)


# Distribution of RNA Ribo MAD Ratio (Human) . 

MAD_ratio_plot_human= ggplot(Buffering_human, aes(x= MAD_ratio))+geom_histogram(binwidth = 0.05, fill = "#B7E6A5",alpha= 0.5)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines")) +
  scale_x_continuous(name = "MAD Ribo/MAD RNA ratio", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(name = "Frequency", breaks = scales::pretty_breaks(n = 10),expand = c(0, 0))+ coord_cartesian(ylim = c(0, 700), xlim=c(0 , 3.5))
                                                                                                                                                            

MAD_ratio_plot_human



# Make the MAD plot with Buffering score


MAD_table_plot_human=ggplot(Buffering_human_rank, aes(x = MAD_RNA, y = MAD_Ribo, shape = Buffering, color = as.factor(Buffering))) +
  # Outer layer for border effect with a larger size
  geom_point(color = NA, stroke = 1.5) +
  # Inner layer for filled points
  geom_point(size = 1, aes(fill = Buffering, alpha = ifelse(Buffering== "TB_score=1", 0.4, ifelse(Buffering == "TB_score=2", 0.4, 0.1)))) +
  geom_abline(intercept = 0, slope = 1, color = "#003147", linetype= "dashed")+
  scale_shape_manual(values = c(24, 23, 21)) +  # Use shapes with borders
  scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5")) +  # Set fill colors
  scale_color_manual(values = c("#03334a", "#067179", "#a0de87")) +  # Set border colors
  labs(shape = "Condition", fill = "Fill Color", color = "Border Color")+guides(alpha = "none", color="none", fill="none")+theme(axis.text = element_text(size = 12),
                                                                                                                                 axis.title = element_text(size = 12),
                                                                                                                                 plot.title = element_text(hjust = 0.5,face="plain",size= 12),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
                                                                                                                                 panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "Median Absolute Deviation (mRNA abundance)") +scale_y_continuous(name = "Median Absolute Deviation (Ribosome Occupancy)")+theme(
    axis.ticks = element_line(color = "black"),  # Customize tick marks
    axis.ticks.length = unit(0.2, "cm"),  # Adjust length of tick marks
    axis.text.x = element_text(hjust = 0.5),  # Rotate x-axis text
    axis.text.y = element_text(hjust = 0.5)  # Customize y-axis text size
  )+labs(shape= "Buffering_score")+theme(legend.position ="none")+coord_cartesian(ylim = c(0, 3), xlim=c(0 , 3))
MAD_table_plot_human


#to make the function work
RNA_human_nond= RNA_clr
human_TE_nond= TE_GBM_clr
Ribo_human_nond= Ribo_clr





# Examples of buffering 


FUS_RNA_nond <- RNA_human_nond[transcript =="FUS"] 
FUS_RNA_nond[,1:5]
FUS_TE_nond<- human_TE_nond[transcript =="FUS"]
FUS_TE_nond[,1:5]
FUS_Ribo_nond = Ribo_human_nond[transcript=="FUS"]


# corelation of FUS RNA and TE across samples
FUS_cor_nond= cor.test(as.numeric(FUS_TE_nond[,2:ncol(FUS_TE_nond)]), as.numeric(FUS_RNA_nond[,2:ncol(FUS_TE_nond)]), method = "spearman")
FUS_cor_nond
FUS_RNA_TE_nond = data.table(as.numeric(FUS_RNA_nond[,2:ncol(FUS_TE_nond)]), as.numeric(FUS_TE_nond[,2:ncol(FUS_TE_nond)]))

colnames(FUS_RNA_TE_nond) = c("FUS_RNA_nond", "FUS_TE_nond")

cell_line_info= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/infor_filter.csv", header = TRUE)
dim(cell_line_info)
head(cell_line_info, 5)

# each row is sample, therefore add sample name. 
sample_name=colnames(RNA_human_nond)
sample_name= subset(sample_name, sample_name != "transcript") 
FUS_RNA_TE_nond$sample_name=sample_name

cell_line_info_inuse= subset(cell_line_info, cell_line_info$experiment_alias %in% sample_name)

merged_df <- merge(FUS_RNA_TE_nond,cell_line_info_inuse[, c("experiment_alias","cell_line")], by.x = "sample_name",by.y= "experiment_alias", all.x = TRUE)


merged_df=na.omit(merged_df)

merged_df$color_group= ifelse(merged_df$cell_line =="HEK293T", "HEK293T","Other")
merged_df$color_group= factor(merged_df$color_group,levels=c("HEK293T","Other"))
FUS_RNA_TE_nond_plot= ggplot(merged_df, aes(x=FUS_RNA_nond, y=FUS_TE_nond, color=color_group, fill=color_group))+geom_point(alpha=0.5, size= 1, shape=24)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Translation efficiency", breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values= c("HEK293T"= "#03334a", Other= "#03334a"))+labs(color= "Cell_line")+theme(legend.position =  "none")+
  scale_fill_manual(values= c("HEK293T"= "#045275", Other= "#045275"))

FUS_RNA_TE_nond_plot

#For RNA Ribo FUS 
FUS_cor_nond_Ribo= cor.test(as.numeric(FUS_Ribo_nond[,2:ncol(FUS_TE_nond)]), as.numeric(FUS_RNA_nond[,2:ncol(FUS_TE_nond)]), method = "spearman")
FUS_cor_nond_Ribo
FUS_RNA_Ribo_nond = data.table(as.numeric(FUS_Ribo_nond[,2:ncol(FUS_TE_nond)]), as.numeric(FUS_RNA_nond[,2:ncol(FUS_TE_nond)]))

colnames(FUS_RNA_Ribo_nond) = c("FUS_Ribo_nond", "FUS_RNA_nond")
FUS_RNA_Ribo_nond$sample_name=sample_name
merged_df <- merge(FUS_RNA_Ribo_nond,cell_line_info_inuse[, c("experiment_alias","cell_line")], by.x = "sample_name",by.y= "experiment_alias", all.x = TRUE)
merged_df=na.omit(merged_df)

merged_df$color_group= ifelse(merged_df$cell_line =="HEK293T", "HEK293T","Other")
merged_df$color_group= factor(merged_df$color_group,levels=c("HEK293T","Other"))
FUS_RNA_Ribo_nond_plot= ggplot(merged_df, aes(x=FUS_RNA_nond, y=FUS_Ribo_nond, color=color_group,fill=color_group))+geom_point(alpha=0.5, size= 1, shape=24)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Ribosome occupancy", breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values= c("HEK293T"= "#03334a", Other= "#03334a"))+labs(color= "Cell_line")+theme(legend.position =  "none")+
  scale_fill_manual(values= c("HEK293T"= "#045275", Other= "#045275"))
FUS_RNA_Ribo_nond_plot+geom_text(aes(x=3,y=6.5, label = "ρ= 0.27", color="black"))
FUS_RNA_Ribo_nond_plot



#G6PD
# Selected G6PD as this is  -.05 corleation between in RNA_TE corelation
G6PD_RNA_nond <- RNA_human_nond[transcript =="G6PD"] 
G6PD_RNA_nond[,1:5]
G6PD_TE_nond<- human_TE_nond[transcript =="G6PD"]
G6PD_TE_nond[,1:5]
G6PD_Ribo_nond = Ribo_human_nond[transcript=="G6PD"]


G6PD_cor_nond= cor.test(as.numeric(G6PD_TE_nond[,2:ncol(G6PD_TE_nond)]), as.numeric(G6PD_RNA_nond[,2:ncol(G6PD_TE_nond)]), method = "spearman")
G6PD_cor_nond
G6PD_RNA_TE_nond = data.table(as.numeric(G6PD_RNA_nond[,2:ncol(G6PD_TE_nond)]), as.numeric(G6PD_TE_nond[,2:ncol(G6PD_TE_nond)]))

colnames(G6PD_RNA_TE_nond) = c("G6PD_RNA_nond", "G6PD_TE_nond")
cell_line_info= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/infor_filter.csv", header = TRUE)
dim(cell_line_info)
head(cell_line_info, 5)


# each row is sample, therefore add sample name. 
sample_name=colnames(RNA_human_nond)
sample_name= subset(sample_name, sample_name != "transcript") 
G6PD_RNA_TE_nond$sample_name=sample_name

cell_line_info_inuse= subset(cell_line_info, cell_line_info$experiment_alias %in% sample_name)

merged_df <- merge(G6PD_RNA_TE_nond,cell_line_info_inuse[, c("experiment_alias","cell_line")], by.x = "sample_name",by.y= "experiment_alias", all.x = TRUE)


merged_df=na.omit(merged_df)

merged_df$color_group= ifelse(merged_df$cell_line =="HEK293T", "HEK293T","Other")
merged_df$color_group= factor(merged_df$color_group,levels=c("HEK293T","Other"))
G6PD_RNA_TE_nond_plot= ggplot(merged_df, aes(x=G6PD_RNA_nond, y=G6PD_TE_nond, color=color_group,fill=color_group))+geom_point(alpha=0.5,size= 1, shape=24)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Translation efficiency", breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values= c("HEK293T"= "#a0de87", Other= "#a0de87"))+labs(color= "Cell_line")+theme(legend.position =  "none")+
  scale_fill_manual(values= c("HEK293T"= "#70ce4b", Other= "#70ce4b"))
G6PD_RNA_TE_nond_plot



#G6PD_RNA Ribo
G6PD_cor_nond_Ribo= cor.test(as.numeric(G6PD_Ribo_nond[,2:ncol(G6PD_TE_nond)]), as.numeric(G6PD_RNA_nond[,2:ncol(G6PD_TE_nond)]), method = "spearman")
G6PD_cor_nond_Ribo
G6PD_RNA_Ribo_nond = data.table(as.numeric(G6PD_Ribo_nond[,2:ncol(G6PD_TE_nond)]), as.numeric(G6PD_RNA_nond[,2:ncol(G6PD_TE_nond)]))

colnames(G6PD_RNA_Ribo_nond) = c("G6PD_Ribo_nond", "G6PD_RNA_nond")
G6PD_RNA_Ribo_nond$sample_name=sample_name
merged_df <- merge(G6PD_RNA_Ribo_nond,cell_line_info_inuse[, c("experiment_alias","cell_line")], by.x = "sample_name",by.y= "experiment_alias", all.x = TRUE)
merged_df=na.omit(merged_df)

merged_df$color_group= ifelse(merged_df$cell_line =="HEK293T", "HEK293T","Other")
merged_df$color_group= factor(merged_df$color_group,levels=c("HEK293T","Other"))
G6PD_RNA_Ribo_nond_plot= ggplot(merged_df, aes(x=G6PD_RNA_nond, y=G6PD_Ribo_nond, color=color_group,fill=color_group))+geom_point(alpha=0.5,size= 1, shape=24)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Ribosome occupancy", breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values= c("HEK293T"= "#a0de87", Other= "#a0de87"))+labs(color= "Cell_line")+theme(legend.position = "none")+
  scale_fill_manual(values= c("HEK293T"= "#70ce4b", Other= "#70ce4b"))
G6PD_RNA_Ribo_nond_plot


# other examples for supplementary
##RPS11
# Selected G6PD as this is  0.1 corleation between in RNA_TE corelation
PPIA_RNA_nond <- RNA_human_nond[transcript =="PPIA"] 
PPIA_RNA_nond[,1:5]
PPIA_TE_nond<- human_TE_nond[transcript =="PPIA"]
PPIA_TE_nond[,1:5]
PPIA_Ribo_nond = Ribo_human_nond[transcript=="PPIA"]


PPIA_cor_nond= cor.test(as.numeric(PPIA_TE_nond[,2:ncol(PPIA_TE_nond)]), as.numeric(PPIA_RNA_nond[,2:ncol(PPIA_TE_nond)]), method = "spearman")
PPIA_cor_nond
PPIA_RNA_TE_nond = data.table(as.numeric(PPIA_RNA_nond[,2:ncol(PPIA_TE_nond)]), as.numeric(PPIA_TE_nond[,2:ncol(PPIA_TE_nond)]))

colnames(PPIA_RNA_TE_nond) = c("PPIA_RNA_nond", "PPIA_TE_nond")
# each row is sample, therefore add sample name. 
sample_name=colnames(RNA_human_nond)
sample_name= subset(sample_name, sample_name != "transcript") 
PPIA_RNA_TE_nond$sample_name=sample_name

cell_line_info_inuse= subset(cell_line_info, cell_line_info$experiment_alias %in% sample_name)

merged_df <- merge(PPIA_RNA_TE_nond,cell_line_info_inuse[, c("experiment_alias","cell_line")], by.x = "sample_name",by.y= "experiment_alias", all.x = TRUE)


merged_df=na.omit(merged_df)

merged_df$color_group= ifelse(merged_df$cell_line =="HEK293T", "HEK293T","Other")
merged_df$color_group= factor(merged_df$color_group,levels=c("HEK293T","Other"))
PPIA_RNA_TE_nond_plot= ggplot(merged_df, aes(x=PPIA_RNA_nond, y=PPIA_TE_nond, color=color_group,fill=color_group))+geom_point(alpha=0.5,size= 1, shape=24)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Translation efficiency", breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values= c("HEK293T"= "#03334a", Other= "#03334a"))+labs(color= "Cell_line")+theme(legend.position = "none")+
  scale_fill_manual(values= c("HEK293T"= "#045275", Other= "#045275"))
PPIA_RNA_TE_nond_plot



#PPIA_RNA Ribo
PPIA_cor_nond_Ribo= cor.test(as.numeric(PPIA_Ribo_nond[,2:ncol(PPIA_TE_nond)]), as.numeric(PPIA_RNA_nond[,2:ncol(PPIA_TE_nond)]), method = "spearman")
PPIA_cor_nond_Ribo
PPIA_RNA_Ribo_nond = data.table(as.numeric(PPIA_Ribo_nond[,2:ncol(PPIA_TE_nond)]), as.numeric(PPIA_RNA_nond[,2:ncol(PPIA_TE_nond)]))

colnames(PPIA_RNA_Ribo_nond) = c("PPIA_Ribo_nond", "PPIA_RNA_nond")
PPIA_RNA_Ribo_nond$sample_name=sample_name
merged_df <- merge(PPIA_RNA_Ribo_nond,cell_line_info_inuse[, c("experiment_alias","cell_line")], by.x = "sample_name",by.y= "experiment_alias", all.x = TRUE)
merged_df=na.omit(merged_df)

merged_df$color_group= ifelse(merged_df$cell_line =="HEK293T", "HEK293T","Other")
merged_df$color_group= factor(merged_df$color_group,levels=c("HEK293T","Other"))
PPIA_RNA_Ribo_nond_plot= ggplot(merged_df, aes(x=PPIA_RNA_nond, y=PPIA_Ribo_nond, color=color_group,fill=color_group))+geom_point(alpha=0.5,size= 1, shape=24)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Ribosome occupancy", breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values= c("HEK293T"= "#03334a", Other= "#03334a"))+labs(color= "Cell_line")+theme(legend.position = "none")+
  scale_fill_manual(values= c("HEK293T"= "#045275", Other= "#045275"))
PPIA_RNA_Ribo_nond_plot



## non buffered genes example 2
##FKBP11 no longer exist in the list . Add another example with RNA _Ribo is >0.8 and TE _RNA >0
head(subset(Buffering_human_rank, Buffering_human_rank$TE_RNA_cor_value_nond >0 & Buffering_human_rank$Ribo_RNA_cor_value_nond >0.8),30)

# 
JUNB_RNA_nond <- RNA_human_nond[transcript =="JUNB"] 
JUNB_RNA_nond[,1:5]
JUNB_TE_nond<- human_TE_nond[transcript =="JUNB"]
JUNB_TE_nond[,1:5]
JUNB_Ribo_nond = Ribo_human_nond[transcript=="JUNB"]


JUNB_cor_nond= cor.test(as.numeric(JUNB_TE_nond[,2:ncol(JUNB_TE_nond)]), as.numeric(JUNB_RNA_nond[,2:ncol(JUNB_TE_nond)]), method = "spearman")
JUNB_cor_nond
JUNB_RNA_TE_nond = data.table(as.numeric(JUNB_RNA_nond[,2:ncol(JUNB_TE_nond)]), as.numeric(JUNB_TE_nond[,2:ncol(JUNB_TE_nond)]))

colnames(JUNB_RNA_TE_nond) = c("JUNB_RNA_nond", "JUNB_TE_nond")
# each row is sample, therefore add sample name. 
sample_name=colnames(RNA_human_nond)
sample_name= subset(sample_name, sample_name != "transcript") 
JUNB_RNA_TE_nond$sample_name=sample_name

cell_line_info_inuse= subset(cell_line_info, cell_line_info$experiment_alias %in% sample_name)

merged_df <- merge(JUNB_RNA_TE_nond,cell_line_info_inuse[, c("experiment_alias","cell_line")], by.x = "sample_name",by.y= "experiment_alias", all.x = TRUE)


merged_df=na.omit(merged_df)

merged_df$color_group= ifelse(merged_df$cell_line =="HEK293T", "HEK293T","Other")
merged_df$color_group= factor(merged_df$color_group,levels=c("HEK293T","Other"))
JUNB_RNA_TE_nond_plot= ggplot(merged_df, aes(x=JUNB_RNA_nond, y=JUNB_TE_nond, color=color_group,fill=color_group))+geom_point(alpha=0.5,size= 1, shape=24)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Translation efficiency", breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values= c("HEK293T"= "#03334a", Other= "#03334a"))+labs(color= "Cell_line")+theme(legend.position = "none")+
  scale_fill_manual(values= c("HEK293T"= "#045275", Other= "#045275"))
JUNB_RNA_TE_nond_plot



#FKBP11_RNA Ribo
JUNB_cor_nond_Ribo= cor.test(as.numeric(JUNB_Ribo_nond[,2:ncol(JUNB_TE_nond)]), as.numeric(JUNB_RNA_nond[,2:ncol(JUNB_TE_nond)]), method = "spearman")
JUNB_cor_nond_Ribo
JUNB_RNA_Ribo_nond = data.table(as.numeric(JUNB_Ribo_nond[,2:ncol(JUNB_TE_nond)]), as.numeric(JUNB_RNA_nond[,2:ncol(JUNB_TE_nond)]))

colnames(JUNB_RNA_Ribo_nond) = c("JUNB_Ribo_nond", "JUNB_RNA_nond")
JUNB_RNA_Ribo_nond$sample_name=sample_name
merged_df <- merge(JUNB_RNA_Ribo_nond,cell_line_info_inuse[, c("experiment_alias","cell_line")], by.x = "sample_name",by.y= "experiment_alias", all.x = TRUE)
merged_df=na.omit(merged_df)

merged_df$color_group= ifelse(merged_df$cell_line =="HEK293T", "HEK293T","Other")
merged_df$color_group= factor(merged_df$color_group,levels=c("HEK293T","Other"))
JUNB_RNA_Ribo_nond_plot= ggplot(merged_df, aes(x=JUNB_RNA_nond, y=JUNB_Ribo_nond, color=color_group,fill=color_group))+geom_point(alpha=0.5,size= 1, shape=24)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Ribosome occupancy", breaks = scales::pretty_breaks(n = 5))+
  scale_color_manual(values= c("HEK293T"= "#03334a", Other= "#03334a"))+labs(color= "Cell_line")+theme(legend.position = "none")+
  scale_fill_manual(values= c("HEK293T"= "#045275", Other= "#045275"))
JUNB_RNA_Ribo_nond_plot





#2 C. human sequence
human_stat=read.table("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/human_stats.txt", header = TRUE)
dim(human_stat)
human_stat=as.data.table(human_stat)
# there are many duplicates so collapse to one.
human_stat <- unique(human_stat, by= "Gene")

rep="CDS_length"

human_stat <- human_stat %>%
  mutate(Buffering = if_else(Gene %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(Gene %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))

# for Matchit
human_stat[, treat := ifelse(Buffering == "TB_score=1", 1,
                             ifelse(Buffering == "TB_score=2", 1, 0))]
#remove NA
human_stat <- na.omit(human_stat)
sum(is.na(human_stat))

#set seed
set.seed(123)
m.out <- matchit(treat ~ UTR3_len+UTR5_len, human_stat, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched


# make a box plot of GC content of the data matched
library(ggpubr)
data_matched_long <- pivot_longer(data_matched, cols = c(CDS_len), names_to = "Variable", values_to = "Value")                   
data_matched_long  
data_matched_long_CDS_length= data_matched_long
CDS_length_plot=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot(width= 0.3)+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs( y = "CDS length (log10)")+theme(axis.text = element_text(size = 12),
                                        axis.title = element_text(size = 12),
                                        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.1, 0.1))+coord_cartesian(ylim = c(0,4.5))

CDS_length_plot


my_comparison = list( c("TB_score=1", "TB_score=2"), c("TB_score=1", "TB_score=3"), c("TB_score=2", "TB_score=3"))


significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )

significance_data


#5UTR length
rep="5UTR_length"

set.seed(123)
m.out <- matchit(treat ~ CDS_len, human_stat, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched


# make a box plot of GC content of the data matched
library(ggpubr)
data_matched_long <- pivot_longer(data_matched, cols = c(UTR5_len), names_to = "Variable", values_to = "Value")                   
data_matched_long
data_matched_long_5UTR_length= data_matched_long 
UTR5_length_plot=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs( y = "5UTR length")+theme(axis.text = element_text(size = 8),
                                 axis.title = element_text(size = 8),
                                 plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0,4.5))

UTR5_length_plot
data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )


# #3 UTR length
rep="3UTR_length"
set.seed(123)

m.out <- matchit(treat ~ CDS_len, human_stat, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched


# make a box plot of GC content of the data matched
library(ggpubr)
data_matched_long <- pivot_longer(data_matched, cols = c(UTR3_len), names_to = "Variable", values_to = "Value")                   
data_matched_long 
data_matched_long_3UTR_length= data_matched_long
UTR3_length_plot=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs( y = "3UTR length")+theme(axis.text = element_text(size = 12),
                                 axis.title = element_text(size = 12),
                                 plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0,4.5))

UTR3_length_plot
data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )


# to make a single plot of all feature lengths. 
# to combine all variable, value, buffering for each 
length_new= rbind(data_matched_long_CDS_length[,c(1, 7,12,13)],data_matched_long_5UTR_length[,c(1, 7,12,13)], data_matched_long_3UTR_length[,c(1, 7,12,13)])
length_new$Variable <- factor(length_new$Variable, levels = c("CDS_len", "UTR5_len", "UTR3_len"))
Combined_length= ggplot(length_new, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs( y = "length (log10)")+theme(axis.text = element_text(size = 12),
                                    axis.title = element_text(size = 12),
                                    plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0.5,4.3))




Combined_length


# For supplementary


# # 2A human CDS GC
# 2C human UTR5 GC
# 2E human UTR 5 GC 

m.out <- matchit(treat ~ UTR3_GC+UTR5_GC, human_stat, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched


# make a box plot of GC content of the data matched

data_matched_long <- pivot_longer(data_matched, cols = c(CDS_GC), names_to = "Variable", values_to = "Value")                   
data_matched_long 
data_matched_long_CDS_GC= data_matched_long
CDS_GC_plot=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs(y = "CDS GC content")+theme(axis.text = element_text(size = 12),
                                   axis.title = element_text(size = 12),
                                   plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0,0.9))

CDS_GC_plot


my_comparison = list( c("TB_score=1", "TB_score=2"), c("TB_score=1", "TB_score=3"), c("TB_score=2", "TB_score=3"))


significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )

significance_data


#5UTR GC
rep="5UTR_GC"

set.seed(123)
m.out <- matchit(treat ~ CDS_GC, human_stat, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched


# make a box plot of GC content of the data matched
library(ggpubr)
data_matched_long <- pivot_longer(data_matched, cols = c(UTR5_GC), names_to = "Variable", values_to = "Value")                   
data_matched_long  
data_matched_long_5UTR_GC= data_matched_long
UTR5_GC_plot=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs( y = "5 UTR GC content")+theme(axis.text = element_text(size = 12),
                                      axis.title = element_text(size = 12),
                                      plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0,0.9))

UTR5_GC_plot
significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )

significance_data



# #3 UTR GC
rep="3UTR_GC"
set.seed(123)

m.out <- matchit(treat ~ CDS_GC, human_stat, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched


# make a box plot of GC content of the data matched
library(ggpubr)
data_matched_long <- pivot_longer(data_matched, cols = c(UTR3_GC), names_to = "Variable", values_to = "Value")                   
data_matched_long 
data_matched_long_3UTR_GC= data_matched_long
UTR3_GC_plot=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs( y = "3 UTR GC content ")+theme(axis.text = element_text(size = 12),
                                       axis.title = element_text(size = 12),
                                       plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0,0.9))

UTR3_GC_plot

significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data_matched_long = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )

significance_data


length_new= rbind(data_matched_long_CDS_GC[,c(1, 7,12,13)],data_matched_long_5UTR_GC[,c(1, 7,12,13)], data_matched_long_3UTR_GC[,c(1, 7,12,13)])
length_new$Variable <- factor(length_new$Variable, levels = c("CDS_GC", "UTR5_GC", "UTR3_GC"))
Combined_length_GC= ggplot(length_new, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs( y = "GC content")+theme(axis.text = element_text(size = 12),
                                    axis.title = element_text(size = 12),
                                    plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0.2,0.9))




Combined_length_GC



## # Plot the median mRNA expression levels  of all genes
numeric_RNA=RNA_human_nond%>%
  mutate_at(vars(-all_of(c("transcript"))), as.numeric)
# Add the medians as a new column to the data table
medians <- numeric_RNA[, apply(.SD,1, median),.SDcols =-"transcript"]
# 
RNA_median<- numeric_RNA %>%
  mutate(median = medians)
# also add if buffering or non buffering
RNA_median <-RNA_median %>%
  mutate(Buffering = if_else(transcript %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(transcript %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))


# Effect of expression levels on buffering.
#with all genes in x axis and y xis mRNA median levels
RNA_median= cbind(RNA_median$transcript, RNA_median$median, RNA_median$Buffering)
colnames(RNA_median)= c("transcript","median", "Buffering")
RNA_median= as.data.table(RNA_median)
RNA_median=RNA_median%>%
  mutate_at(vars(-all_of(c("transcript", "Buffering"))), as.numeric)
#plot the median 
Median_RNA_expression_human=ggplot(RNA_median, aes(x= Buffering, y= median, fill= Buffering))+geom_boxplot()+theme(axis.text = element_text(size = 12),
                                                                                                                   axis.title = element_text(size = 12),
                                                                                                                   plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                   panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_fill_manual(values = c( "#045275","#089099", "#B7E6A5"))+labs(
                                                                                                                     x = "Buffering",y = "Median RNA expression")+ylim(-1,5)


Median_RNA_expression_human=ggplot(RNA_median, aes(x= median, color= Buffering))+geom_density(size = 1.0)+theme(axis.text = element_text(size = 12),
                                                                                                                axis.title = element_text(size = 12),
                                                                                                                plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_color_manual(values = c( "#045275","#089099", "#B7E6A5"))+labs(
                                                                                                                  x = "Median RNA expression",y = "Frequency")+ylim(0,0.5)
Median_RNA_expression_human


#calculate the significance


significance_data <- RNA_median %>%
  summarise(
    TB1_TB2 = wilcox.test(median ~ Buffering, RNA_median = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(median ~ Buffering, RNA_median = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(median ~ Buffering, RNA_median = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )


# Figure 3

#Are buffered genes more variable or less compared to non buffered set?                                                                                                                                                                                                                 
# plot MAD values as function of buffering
head(Buffering_human)
Buffering_human <-Buffering_human %>%
  mutate(Buffering = if_else(transcript %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(transcript %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))

#match by similar RNA expression by joing the median to the table

RNA_MAD_human= ggplot(Buffering_human, aes(x= Buffering, y= MAD_RNA, fill= Buffering))+geom_violin()+stat_summary(fun = median, geom = "crossbar", 
                                                                                                                   width = 0.2, fill = "black", color = "black")+theme(axis.text = element_text(size = 12),
                                                                                                           axis.title = element_text(size = 12),
                                                                                                           plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                           panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs(
                                                                                                             x = "Buffering",y = "MAD RNA expression")+coord_cartesian(ylim=c(0,2.6))
RNA_MAD_human








#codon bias

# clean the FASTA sequnce
fasta_sequences <- readDNAStringSet("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/non deduplicated data/final_TE_8_2_23/CDS_human.fa")
fasta_sequences
#calculate teh freequncy of each codon in each gene
cf_all <- count_codons(fasta_sequences)
cf_all

cf_buffering= subset(cf_all, row.names(cf_all) %in% Top_buffered_human$transcript)
dim(cf_buffering)
# 322 genes present . 
score1= Buffering_human_rank[Buffering_human_rank$Buffering =="TB_score=1", "transcript"]
score1=as.list(score1)
cf_buffering2= subset(cf_all, row.names(cf_all) %in% Top2_buffered_human$transcript)
dim(cf_buffering2)
# for a random set
#first remove those that are buffered 
cf_random= subset(cf_all, ! row.names(cf_all) %in% row.names(cf_buffering))
sampled_indices <- sample(nrow(cf_random), 250, replace=FALSE)
cf_random1=cf_random[sampled_indices, ]
dim(cf_random1)
#est_rscu returns the RSCU value of codons=he observed frequency of a particular codon to the expected frequency under equal usage of synonymous codons.
#wCAI is a measure that evaluates the adaptation of a gene's codon usage to the tRNA pool in an organism.
rscu_buffering=
  est_rscu(cf_buffering, weight = 1, pseudo_cnt = 1, codon_table = get_codon_table())
rscu_random=est_rscu(cf_random1, weight = 1, pseudo_cnt = 1, codon_table = get_codon_table())

#how does it look if it were to be compared to non buffered gene set 
Non_buffered_transcript= Buffering_human_rank[Buffering_human_rank$Buffering =="TB_score=3", "transcript"]
Non_buffered_transcript=as.list(Non_buffered_transcript)
cf_non_buffered= subset(cf_all, row.names(cf_all) %in%  Non_buffered_transcript$transcript)
dim(cf_non_buffered)


# can we calculate the RSCu of highly expressed gnes

quantiles <- quantile(RNA_median$median, probs = c(0.95))
highly_expressed_genes= RNA_median[median > quantiles]
#remove dummy gene
dim( highly_expressed_genes)
highly_expressed_genes= highly_expressed_genes[-nrow(highly_expressed_genes), ]
#calculate Cf
cf_HEG= subset(cf_all, row.names(cf_all) %in% highly_expressed_genes$transcript)
dim(cf_HEG)
rscu_HEG=
  est_rscu(cf_HEG, weight = 1, pseudo_cnt = 1, codon_table = get_codon_table())




#how does the cai work. Codon adaptation index
library(data.table)
cai_random= get_cai(cf_random1, rscu_HEG)
cai_random= as.data.table(cai_random)
cai_random$transcript= row.names(cf_random1)
dim(cai_random)
cai_random=melt(cai_random)

cai_buffered= get_cai(cf_buffering, rscu_HEG)
cai_buffered= as.data.table(cai_buffered)
cai_buffered$transcript= row.names(cf_buffering)
cai_buffered=melt(cai_buffered)

cai_non_buffered= get_cai(cf_non_buffered, rscu_HEG)
cai_non_buffered= as.data.table(cai_non_buffered)
cai_non_buffered$transcript= row.names(cf_non_buffered)
cai_non_buffered=melt(cai_non_buffered)

cai_bufferd2= get_cai(cf_buffering2, rscu_HEG)
cai_bufferd2= as.data.table(cai_bufferd2)
cai_bufferd2$transcript= row.names(cf_buffering2)
cai_bufferd2=melt(cai_bufferd2)


#Combine both together
cai_buf_ran= rbind(cai_buffered ,cai_bufferd2,cai_non_buffered)

cai_buf_ran_long_plot= ggplot(cai_buf_ran, aes(x = variable, y = value, fill= variable)) +
  geom_boxplot()+
  labs(x = "Condition", y = "Codon adaptation index ")+theme(axis.text = element_text(size = 12),
                                                             axis.title = element_text(size = 12),
                                                             plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                             panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))

cai_buf_ran_long_plot  
# make a table with genes included.


# # Find if there is any association with CSC and buffering.  
#make a plot of the codon stabilization score (CSC) *** the Pearson correlation coefficient between mRNA stability and codon occurrence.
CSC_table= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/non deduplicated data/final_TE_8_2_23/elife-45396-fig1-data2-v2.csv")

head(CSC_table)
# we need to plot between the CSC and the number of each codon. Since each gene will be of differnt length
# so convert each codon into percentage of total codons in each each gene. We need a table of freequency of codons for buffered, non buffered set
# count the number of each codon percentage for all genes and plot against CSC 
# we need to find the freequency of each codon in a gene set.
# so total number of each codon(colsums) and divide by the total number of codons in the gene set . See if there is any  correlation between codon frequency and CSC bteween the three gene sets  
calculate_colsum_ratio <- function(matrix_input) {
  # Calculate column sums
  col_sums <- colSums(matrix_input)
  
  # Calculate total sum of column sums
  total_sum <- sum(col_sums)
  
  # Calculate ratio of each column sum to total sum
  ratios <- col_sums / total_sum
  
  # Create a data frame with column names and ratios
  result <- data.table(codon = names(col_sums), Ratios = ratios)
  
  # Set the name of the result
  
  return(result)
}
cf_buffering_fre= calculate_colsum_ratio(cf_buffering)
cf_buffering_fre=merge(cf_buffering_fre, CSC_table, by="codon")
CSC_buffering_cor= cor(cf_buffering_fre$X293T_ORFome, cf_buffering_fre$Ratios, method= "pearson")
CSC_buffering_cor_endo= cor(cf_buffering_fre$X293T_endo, cf_buffering_fre$Ratios, method= "pearson")
#buffered 2
cf_buffering2_fre= calculate_colsum_ratio(cf_buffering2)
cf_buffering2_fre=merge(cf_buffering2_fre, CSC_table, by="codon")
CSC_buffered_TB_Score2_cor= cor(cf_buffering2_fre$X293T_ORFome, cf_buffering_fre$Ratios, method= "pearson")
CSC_buffered_TB_Score2_cor_endo= cor(cf_buffering2_fre$X293T_endo, cf_buffering_fre$Ratios, method= "pearson")
# non _buffered
cf_nonbuffered_fre= calculate_colsum_ratio(cf_non_buffered)
cf_nonbuffered_fre=merge(cf_nonbuffered_fre, CSC_table, by="codon")
CSC_nonbuffered_cor=cor(cf_nonbuffered_fre$Ratios, cf_nonbuffered_fre$X293T_ORFome, method= "pearson")
CSC_nonbuffered_cor_endo=cor(cf_nonbuffered_fre$Ratios, cf_nonbuffered_fre$X293T_endo, method= "pearson")

cf_random_fre= calculate_colsum_ratio(cf_random1)
cf_random_fre=merge(cf_random_fre, CSC_table, by="codon")
CSC_random_cor=cor(cf_random_fre$Ratios, cf_random_fre$X293T_ORFome, method= "pearson")
CSC_random_cor_endo=cor(cf_random_fre$Ratios, cf_random_fre$X293T_endo, method= "pearson")
CSC_cor_all= cbind(CSC_buffering_cor, CSC_buffered_TB_Score2_cor, CSC_nonbuffered_cor,CSC_random_cor,CSC_buffering_cor_endo, CSC_buffered_TB_Score2_cor_endo, CSC_nonbuffered_cor_endo,CSC_random_cor_endo)
head(CSC_cor_all)



cor(cf_buffering2_fre$Ratios, cf_buffering2_fre$X293T_endo, method= "pearson")
ggplot(cf_buffering_fre, aes(X293T_endo, Ratios))+geom_point()+geom_smooth(method = "lm", se = FALSE)+xlim(-0.25,0.25)+ylim(0,0.05)+labs(title = "Correlation between CSC and codon frequency (TB_score=1)",
                                                                                                                                         x = " CSC (HEK293T)",
                                                                                                                                         y = "Frequency of codon")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                         axis.title = element_text(size = 12),
                                                                                                                                                                         plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                         panel.spacing = unit(0.5, "lines"))+geom_text(aes(label = codon), position = position_nudge(y = 0.001),
                                                                                                                                                                                                                       size = 3,  
                                                                                                                                                                                                                       color = "blue",angle = 45,  # Rotate the text
                                                                                                                                                                                                                       hjust = 0.5,  # Horizontal justification
                                                                                                                                                                                                                       vjust = 0.5) # Vertical justificatio
ggplot(cf_buffering2_fre, aes(X293T_endo, Ratios))+geom_point()+geom_smooth(method = "lm", se = FALSE)+xlim(-0.25,0.25)+ylim(0,0.05)+labs(title = "Correlation between CSC and codon frequency (TB_score=2)",
                                                                                                                                          x = " CSC (HEK293T)",
                                                                                                                                          y = "Frequency of codon")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                          axis.title = element_text(size = 12),
                                                                                                                                                                          plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                          panel.spacing = unit(0.5, "lines"))+geom_text(aes(label = codon), position = position_nudge(y = 0.001),
                                                                                                                                                                                                                        size = 3,  
                                                                                                                                                                                                                        color = "blue",angle = 45,  # Rotate the text
                                                                                                                                                                                                                        hjust = 0.5,  # Horizontal justification
                                                                                                                                                                                                                        vjust = 0.5) # Vertical justificatio
ggplot(cf_nonbuffered_fre, aes(X293T_endo, Ratios))+geom_point()+geom_smooth(method = "lm", se = FALSE)+xlim(-0.25,0.25)+ylim(0,0.05)+labs(title = "Correlation between CSC and codon frequency (TB_score=3)",
                                                                                                                                           x = " CSC (HEK293T)",
                                                                                                                                           y = "Frequency of codon")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                           axis.title = element_text(size = 12),
                                                                                                                                                                           plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                           panel.spacing = unit(0.5, "lines"))+geom_text(aes(label = codon), position = position_nudge(y = 0.001),
                                                                                                                                                                                                                         size = 3,  
                                                                                                                                                                                                                         color = "blue",angle = 45,  # Rotate the text
                                                                                                                                                                                                                         hjust = 0.5,  # Horizontal justification
                                                                                                                                                                                                                         vjust = 0.5) # Vertical justificatio


cor(cf_nonbuffered_fre$Ratios, cf_nonbuffered_fre$X293T_ORFome, method= "pearson")

cor(cf_random_fre$Ratios, cf_random_fre$X293T_endo, method= "pearson")

# make  a table frequency of codons from each condition. Arrange from low CSC to high CSC. plot as histograms. 
cf_freq_all= cbind(cf_buffering_fre$codon, cf_buffering_fre$Ratios, cf_buffering2_fre$Ratios, cf_nonbuffered_fre$Ratios, cf_random_fre$Ratios)
colnames(cf_freq_all)=c("Codon", "TB_score1", "TB_score2", "TB_score3", "Random")
cf_freq_all=as.data.table(cf_freq_all)
cf_freq_all=cf_freq_all%>%
  mutate_at(vars(-all_of(c("Codon"))), as.numeric)
cf_freq_all=merge(cf_freq_all, CSC_table, by.x= "Codon", by.y= "codon")

# take the difference and between the frequency of buffered set and plot it as a function of CSC. 
head(cf_freq_all)
cf_freq_all$Buffered_nonBuffered_diff= (cf_freq_all$TB_score1- cf_freq_all$TB_score3)*100
ggplot(cf_freq_all, aes(X293T_endo, Buffered_nonBuffered_diff))+geom_point(size=1)+geom_smooth(method = "lm", se = FALSE)+xlim(-0.25,0.25)+ylim(-0.8,0.8)+labs(title = "Correlation between CSC and codon frequency difference TB score 1&3)",
                                                                                                                                                               x = " CSC (HEK293T)",
                                                                                                                                                               y = "Difference in codon frequency (%)")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                              axis.title = element_text(size = 12),
                                                                                                                                                                                                              plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                                                              panel.spacing = unit(0.5, "lines"))+geom_text_repel(aes(label = Name),size = 3,  
                                                                                                                                                                                                                                                                  color = "#3aa12b")
cf_freq_all$Buffered2_nonBuffered_diff= (cf_freq_all$TB_score2- cf_freq_all$TB_score3)*100
ggplot(cf_freq_all, aes(X293T_endo, Buffered2_nonBuffered_diff))+geom_point(size=1)+geom_smooth(method = "lm", se = FALSE)+xlim(-0.25,0.25)+ylim(-0.8,0.8)+labs(title = "Correlation between CSC and codon frequency difference TB score 2&3)",
                                                                                                                                                                x = " CSC (HEK293T)",
                                                                                                                                                                y = "Difference in codon frequency (%)")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                               axis.title = element_text(size = 12),
                                                                                                                                                                                                               plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                                                               panel.spacing = unit(0.5, "lines"))+geom_text_repel(aes(label = Name),size = 3,  
                                                                                                                                                                                                                                                                   color = "#3aa12b")
cf_freq_all$Random_nonBuffered_diff= (cf_freq_all$Random- cf_freq_all$TB_score3)*100
ggplot(cf_freq_all, aes(X293T_endo, Random_nonBuffered_diff))+geom_point(size=1)+geom_smooth(method = "lm", se = FALSE)+xlim(-0.25,0.25)+ylim(-0.8,0.8)+labs(title = "Correlation between CSC and codon frequency difference Random-Non Buffered)",
                                                                                                                                                             x = " CSC (HEK293T)",
                                                                                                                                                             y = "Difference in codon frequency (%)")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                            axis.title = element_text(size = 12),
                                                                                                                                                                                                            plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                                                            panel.spacing = unit(0.5, "lines"))+geom_text_repel(aes(label = Name),size = 3,  
                                                                                                                                                                                                                                                                color = "#3aa12b")
#plot the difference in frequency as fucntion of CSC 
cf_freq_all$CombinedCategory <- paste(cf_freq_all$Codon, cf_freq_all$Name, sep = "_")
cf_freq_all$CombinedCategory <- factor(cf_freq_all$CombinedCategory, levels = cf_freq_all$CombinedCategory[order(-cf_freq_all$Buffered_nonBuffered_diff)])
CFD_TBS1=ggplot(cf_freq_all, aes(x= CombinedCategory, y= Buffered_nonBuffered_diff, fill=X293T_endo))+geom_bar(stat = "identity")+labs(,
                                                                                                                                       x = " Codon ",
                                                                                                                                       y = "ΔCodon Frequency")+theme(axis.text = element_text(size = 12, angle = 90),
                                                                                                                                                                     axis.title = element_text(size = 12),
                                                                                                                                                                     plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"), panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                     panel.spacing = unit(0.5, "lines"))+scale_fill_gradient2(low = "purple", mid = "#e9cb2e", high = "green", midpoint = 0)+ylim(-1, 1.0)+theme(legend.position =c(0.8,0.8))

cor.test(cf_freq_all$Buffered_nonBuffered_diff, cf_freq_all$X293T_endo, method= "pearson")
#plot the difference in frequency as fucntion of CSC fr TB score 2 

cf_freq_all$CombinedCategory <- factor(cf_freq_all$CombinedCategory, levels = cf_freq_all$CombinedCategory[order(-cf_freq_all$Buffered2_nonBuffered_diff)])
CFD_TBS2=ggplot(cf_freq_all, aes(x= CombinedCategory, y= Buffered2_nonBuffered_diff, fill=X293T_endo))+geom_bar(stat = "identity")+labs(
  x = " Codon ",
  y = "Difference in codon frequency")+theme(axis.text = element_text(size = 12, angle = 90),
                                             axis.title = element_text(size = 12),
                                             plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"), panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                             panel.spacing = unit(0.5, "lines"))+ scale_fill_gradient2(low = "purple",mid="#e9cb2e", high = "green", midpoint= 0)+ylim(-0.5, 1.0)+theme(legend.position =c(0.8,0.8),legend.text = element_text(size =6),legend.title = element_text(size =12))


CFD_TBS2
cor.test(cf_freq_all$Buffered2_nonBuffered_diff, cf_freq_all$X293T_endo, method= "pearson")
#plot the difference in frequency as fucntion of CSC fr random

cf_freq_all$CombinedCategory <- factor(cf_freq_all$CombinedCategory, levels = cf_freq_all$CombinedCategory[order(-cf_freq_all$Random_nonBuffered_diff)])
CFD_Random=ggplot(cf_freq_all, aes(x= CombinedCategory, y= Random_nonBuffered_diff, fill=X293T_endo))+geom_bar(stat = "identity")+labs(title ="Codon freequency difference  (Random vs TB_score3)",
                                                                                                                                       x = " Codon ",
                                                                                                                                       y = "Difference in codon frequency")+theme(axis.text = element_text(size = 12, angle = 90),
                                                                                                                                                                                  axis.title = element_text(size = 12),
                                                                                                                                                                                  plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"), panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                                  panel.spacing = unit(0.5, "lines"))+ scale_fill_gradient2(low = "blue",mid="white", high = "red", midpoint= 0)+ylim(-0.6, 1)+theme(legend.position ="none")



cor.test(cf_freq_all$Random_nonBuffered_diff, cf_freq_all$X293T_endo, method= "pearson")




#MAD for proteomics

RNA_cancer_cell= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/rnaseq_fpkm_20220624.csv", header = TRUE)
dim(RNA_cancer_cell)
RNA_cancer_cell=as.data.table(RNA_cancer_cell)
RNA_cancer_cell[100:115,1:5]
Protein_cancer_cell= read_excel("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/1-s2.0-S1535610822002744-mmc3 (1).xlsx", sheet=1, skip=1)
dim(Protein_cancer_cell)
Protein_cancer_cell=as.data.table(Protein_cancer_cell)
Protein_cancer_cell[1:5,1:5]
#make a table with gene names as column and cell line as rows so as to extract our gene of interest.
#split the columns into  cell lines
Protein_cancer_cell[, c("Project_Identifier", "Cell_line_human") := tstrsplit(Project_Identifier, ";", type.convert = TRUE)]

Protein_cancer_cell=as.data.table(Protein_cancer_cell)
#extract uniprot ID. 

extracted_part <- gsub("^([^;]+);.*", "\\1", colnames(Protein_cancer_cell))
extracted_part


# from the uniport get all gene names
UniportID= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/idmapping_2024_03_25.csv", header = TRUE)

head(UniportID)
#considering the first is the gene name, extract along with Entry.name. replace SOX21_human with SOX2
UniportID= as.data.table(UniportID)

extracted_part <- gsub("^(\\S+).*", "\\1", UniportID$Gene.Names)
extracted_part
#make a data table of extracted part and Entry.Name.# 
UniportID_table= cbind(extracted_part, UniportID$Entry.Name)
UniportID_table=as.data.table(UniportID_table)
UniportID_table
colnames(UniportID_table)= c("UniportIDGene.Names", "UniportEntry.Name")
#remove the human.
UniportID_table$UniportEntry.Name= gsub("_HUMAN.*$", "",UniportID_table$UniportEntry.Name)

# extract the gene name and replace the column names
pattern <- ".*;(.*)"
extracted_names <- sub(pattern, "\\1", colnames(Protein_cancer_cell))
extracted_names
colnames(Protein_cancer_cell)= (extracted_names)
colnames(Protein_cancer_cell[,8000:8050])

#extract those that have "human"
cols_with_human <- grep("human", names(Protein_cancer_cell), ignore.case = TRUE,value = TRUE)
Protein_cancer_cell <- subset(Protein_cancer_cell, select= cols_with_human)
dim(Protein_cancer_cell)
# remove "human" from column names
colnames(Protein_cancer_cell) <- gsub("_HUMAN.*$", "", colnames(Protein_cancer_cell))
dim(Protein_cancer_cell)
Protein_cancer_cell[,8467:8468]
#make the last column the first column
Protein_cancer_cell= Protein_cancer_cell[, c(8468,1:8467)]
Protein_cancer_cell_numeric<-Protein_cancer_cell[, apply(.SD,2, as.numeric),.SDcols =-"Cell_line_human"]
# calculate the MAD of all genes ignoring NA 
mad_ignore_na <- function(x) {
  mad(x, na.rm = TRUE)
}
#each column is a gene and each row is a cell line. Calculate MAD value for each across cell lines
mad_values <- apply(Protein_cancer_cell_numeric, 2, mad_ignore_na)
# make a table of MAD values
mad_values= as.data.table(mad_values)
mad_values$genes= colnames(Protein_cancer_cell[,2:8468])
mad_values

# replace the gene names from that of uniport. 
# first merge the table with gene names 
mad_values_merge=merge(mad_values,UniportID_table,by.x = "genes", by.y = "UniportEntry.Name") 

# take the median of the protein abundance 
Protein_cancer_cell_median= Protein_cancer_cell[, apply(.SD,2, median,na.rm = TRUE),.SDcols =-"Cell_line_human"]
Protein_cancer_cell_median=as.data.table(Protein_cancer_cell_median)
dim(Protein_cancer_cell_median)

Protein_cancer_cell_median$genes= colnames(Protein_cancer_cell[,2:8468])
Protein_cancer_cell_median[1:2,1:2]
colnames(Protein_cancer_cell_median)= c("median", "genes")


# but these gene name is Uniprot ID 
# add gene names with uniprot ID
Protein_cancer_cell_median= merge(Protein_cancer_cell_median,UniportID_table,by.x = "genes", by.y = "UniportEntry.Name")

#Buffered or not
Protein_cancer_cell_median <- Protein_cancer_cell_median %>%
  mutate(Buffering = if_else(UniportIDGene.Names %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(UniportIDGene.Names %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))
Protein_cancer_cell_median

Protein_cancer_cell_median=Protein_cancer_cell_median[, treat := ifelse(Buffering == "TB_score=1", 1,
                                                                        ifelse(Buffering == "TB_score=2", 1, 0))]

Protein_cancer_cell_median <- na.omit(Protein_cancer_cell_median)







# also calculate on MAD of RNA 
# convert this into numeric
RNA_cancer_cell=RNA_cancer_cell[5:37606, 3:1433 := lapply(.SD, as.numeric), .SDcols = 3:1433]

dim(RNA_cancer_cell)

RNA_cancer_cell= RNA_cancer_cell[5:37606, MAD_RNA_cancer_cell := apply(.SD, 1, function(x) mad_ignore_na(as.numeric(x))), .SDcols = 3:1433]


#merge with RNA MAD with protein median=#make a separate table
Protein_cancer_cell_median_RNA <- merge(Protein_cancer_cell_median, RNA_cancer_cell[, .(X,  MAD_RNA_cancer_cell)], by.x = "UniportIDGene.Names", by.y = "X")
#how many na in MAD 


#
#
#

# find rows with NA in MAD RNA 
sum(is.na(Protein_cancer_cell_median_RNA$MAD_RNA_cancer_cell))
#12 genes with NA in MAD 

dim(Protein_cancer_cell_median_RNA)
# remove NA (thouse that do not have both RNA and protein)
Protein_cancer_cell_median_RNA=na.omit(Protein_cancer_cell_median_RNA)
# 12is removed as no correspoding RNA values 



# 
mad_RNA_cancer=pivot_longer(Protein_cancer_cell_median_RNA, cols = c("MAD_RNA_cancer_cell"), names_to = "Variable", values_to = "Value") 
ggplot(mad_RNA_cancer, aes(x=Variable, y=Value, color= Buffering))+geom_boxplot()+coord_cartesian(ylim=c(0,2000))
#merge the calculated mad value 
Protein_cancer_cell_median_mad_RNA= merge(Protein_cancer_cell_median_RNA,mad_values_merge, by.x=  "UniportIDGene.Names", by.y= "UniportIDGene.Names")
# without matching
Protein_cancer_cell_median_long=pivot_longer(Protein_cancer_cell_median_mad_RNA, cols = c("mad_values"), names_to = "Variable", values_to = "Value")     
mad_values_plot_unmatched= ggplot(Protein_cancer_cell_median_long, aes(x=  Variable, y= Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c( "#045275","#089099", "#B7E6A5"))+labs( x = "Condition", y = "Median Absolute Deviation")+theme(axis.text = element_text(size = 12),
                                                                                                                                              axis.title = element_text(size = 12),
                                                                                                                                              plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                              panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,1.5))+scale_x_discrete(expand = c(0.2, 0.2))

mad_values_plot_unmatched
install.packages("ggbreak")  # if not already installed
library(ggbreak)
p_main=ggplot(Protein_cancer_cell_median_long, aes(x = Value,fill= Buffering)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", 
                 alpha = 0.4, 
                 bins = 100)+scale_fill_manual(values = c( "#045275","#089099", "#B7E6A5"))+labs( y = "Density", x = "Protein abundance (MAD)")+theme(axis.text = element_text(size = 12),
                                                                                                                                                           axis.title = element_text(size = 12),
                                                                                                                                                           plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                           panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(xlim=c(0.2,1.5))

#
p_inset <- ggplot(Protein_cancer_cell_median_long, aes(x = Value, fill = Buffering)) +
geom_histogram(aes(y = after_stat(density)), bins = 100, alpha = 0.4, position = "identity") +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank())+scale_fill_manual(values = c( "#045275","#089099", "#B7E6A5"))
      
   

final_plot <- ggdraw() +
  draw_plot(p_main) +
  draw_plot(p_inset, x = 0.65, y = 0.65, width = 0.3, height = 0.3)
final_plot
significance_data <- Protein_cancer_cell_median_long %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, Protein_cancer_cell_median_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, Protein_cancer_cell_median_long= .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, Protein_cancer_cell_median_long = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )
significance_data


#
library(ggridges)
ggplot(Protein_cancer_cell_median_long, aes(x = Value, y = Buffering, fill = Buffering)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2, color = "black") +
  theme_minimal() 

# match with both 

m.out_abundnace <- matchit( treat ~  MAD_RNA_cancer_cell, Protein_cancer_cell_median_mad_RNA, method = "nearest")

summary(m.out_abundnace)
data_matched_Cancer_cell2020= match.data(m.out_abundnace)
data_matched_Cancer_cell2020


data_matched_Cancer_cell2020=pivot_longer(data_matched_Cancer_cell2020, cols = c("mad_values"), names_to = "Variable", values_to = "Value")     

mad_values_plot_matched= ggplot(data_matched_Cancer_cell2020, aes(x=  Variable, y= Value, fill= Buffering)) +
  geom_boxplot(outlier.shape = NA)+scale_fill_manual(values = c("TB_score=1"= "#045275","TB_score=2"="#089099", "TB_score=3"="#B7E6A5"))+labs(x = "Condition", y = "Median Absolute Deviation")+theme(axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
  plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
panel.spacing = unit(0.1, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+ scale_y_continuous(expand = c(0.0, 0))+coord_cartesian(ylim= c(0,1.5))


mad_values_plot_matched

p_main=ggplot(data_matched_Cancer_cell2020, aes(x = Value,fill= Buffering)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", 
                 alpha = 0.4, 
                 bins = 100)+scale_fill_manual(values = c( "#045275","#089099", "#B7E6A5"))+labs( y = "Density", x = "Protein abundance (MAD)")+theme(axis.text = element_text(size = 12),
                                                                                                                                                      axis.title = element_text(size = 12),
                                                                                                                                                      plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                      panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(xlim=c(0.2,1.5))

#
p_inset <- ggplot(data_matched_Cancer_cell2020, aes(x = Value, fill = Buffering)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100, alpha = 0.4, position = "identity") +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank())+
  scale_fill_manual(values = c( "#045275","#089099", "#B7E6A5"))


final_plot2 <- ggdraw() +
  draw_plot(p_main) +
  draw_plot(p_inset, x = 0.65, y = 0.65, width = 0.3, height = 0.3)
final_plot2
significance_data <- data_matched_Cancer_cell2020 %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data_matched_Cancer_cell2020 = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data_matched_Cancer_cell2020 = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data_matched_Cancer_cell2020 = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )
significance_data






# get this without matched (buffered genes compared to all other genes)

#add mad values





##cell_2020_tissue (for all genes: no matching)
Cell_2020_protein =  read_excel("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/mmc3.xlsx", sheet=7, skip=1)
Cell_2020_protein[1:5,1:5]
dim(Cell_2020_protein)
GeneName_cell2020= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/Gene_name_Cell_2020.csv")
GeneName_cell2020= as.data.table(GeneName_cell2020)
GeneName_cell2020
dim(GeneName_cell2020)
#add gene names using ensembl ID
Cell_2020_protein_MAD=merge(GeneName_cell2020, Cell_2020_protein , by.x=  "Gen_ID", by.y= "gene.id")
Cell_2020_protein_MAD
dim(Cell_2020_protein_MAD)

#calculate MAD across all rows. 
Cell_2020_protein_MAD_numeric= Cell_2020_protein_MAD[, apply(.SD,2, as.numeric),.SDcols =-"Gene_Name"]
Cell_2020_protein_MAD=Cell_2020_protein_MAD[, 2:34]
mad_values_Cell_2020 <- apply(Cell_2020_protein_MAD_numeric, 1, mad_ignore_na)
mad_values_Cell_2020
length(mad_values_Cell_2020) 
mad_values_Cell_2020=as.data.table(mad_values_Cell_2020)
mad_values_Cell_2020$genes= Cell_2020_protein_MAD$Gene_Name
dim(mad_values_Cell_2020)
mad_values_Cell_2020[1:5,]
# whether buffered or not
mad_values_Cell_2020 <- mad_values_Cell_2020 %>%
  mutate(Buffering = if_else(genes %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(genes %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))
mad_values_Cell_2020
mad_values_Cell_2020_long=pivot_longer(mad_values_Cell_2020, cols = c("mad_values_Cell_2020"), names_to = "Variable", values_to = "Value")     
mad_values_Cell_2020_plot= ggplot(mad_values_Cell_2020_long, aes(x=  Variable, y= Value, fill= Buffering)) +
  geom_boxplot(width=0.6)+scale_fill_manual(values = c("#045275","#089099", "#B7E6A5"))+labs(x = "Condition", y = "MAD")+theme(axis.text = element_text(size = 12),
                                                                                                                      axis.title = element_text(size = 12),                                                                                                                                                                                                     plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                      panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,1.5))+scale_x_discrete(expand = c(0.16, 0.16))


mad_values_Cell_2020_plot



significance_data <- mad_values_Cell_2020_long %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, mad_values_Cell_2020_long = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, mad_values_Cell_2020_long= .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, mad_values_Cell_2020_long = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )
significance_data







# can we match  with median abundance


Cell_2020_protein_MAD_median= apply(Cell_2020_protein_MAD_numeric, 1, median, na.rm = TRUE)
mad_values_Cell_2020$median=Cell_2020_protein_MAD_median
mad_values_Cell_2020
mad_values_Cell_2020 <- na.omit(mad_values_Cell_2020)
sum(is.na(mad_values_Cell_2020))
mad_values_Cell_2020=as.data.table(mad_values_Cell_2020)
# add MAD vaues from mRNA
Cell_2020_mRNA =  read_excel("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Protein_RNA_ratio/Cell_2020/mmc4.xlsx", sheet=3,skip=1)
Cell_2020_mRNA[1:5,1:5]

# calculate the MAD for tissue
Cell_2020_mRNA=as.data.table(Cell_2020_mRNA)
dim(Cell_2020_mRNA)
Cell_2020_mRNA= Cell_2020_mRNA[, MAD_RNA_tissue := apply(.SD, 1, mad, na.rm = TRUE), .SDcols = 2:33]
Cell_2020_mRNA$MAD_RNA_tissue
#add gene names using ensembl ID
Cell_2020_mRNA=merge(GeneName_cell2020, Cell_2020_mRNA , by.x=  "Gen_ID", by.y= "gene.id")
Cell_2020_mRNA
dim(Cell_2020_mRNA)



# add the mad column to the median protein values by merging uisng common gene names
mad_values_Cell_2020= mad_values_Cell_2020[Cell_2020_mRNA, MAD_RNA_tissue := i.MAD_RNA_tissue, on = .(genes == Gene_Name)]
mad_values_Cell_2020=na.omit(mad_values_Cell_2020)



# match it
# for Matchit
mad_values_Cell_2020[, treat := ifelse(Buffering == "TB_score=1", 1,
                                       ifelse(Buffering == "TB_score=2", 1, 0))]



### plot the mAD mRNA

mad_RNA_tissue=pivot_longer(mad_values_Cell_2020, cols = c("MAD_RNA_tissue"), names_to = "Variable", values_to = "Value") 
ggplot(mad_RNA_tissue, aes(x=Variable, y=Value, color= Buffering))+geom_boxplot()





# match for MAD of mRNA abundance
m.out_abundnace <- matchit( treat ~ MAD_RNA_tissue, mad_values_Cell_2020, method = "nearest")
summary(m.out_abundnace)
data_matched_Cell_2020= match.data(m.out_abundnace)
data_matched_Cell_2020
mad_values_Cell_2020=pivot_longer(data_matched_Cell_2020, cols = c("mad_values_Cell_2020"), names_to = "Variable", values_to = "Value")     
mad_values_Cell_2020_plot_matched= ggplot(mad_values_Cell_2020, aes(x=  Variable, y= Value, fill= Buffering)) +
  geom_boxplot(outlier.shape = NA)+scale_fill_manual(values = c("#045275","#089099", "#B7E6A5"))+labs(x = "Condition", y = "Median Absolute Deviation")+theme(axis.text = element_text(size = 12),
  axis.title = element_text(size = 12),
  plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
  panel.spacing = unit(0.1, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+ scale_y_continuous(expand = c(0.01, 0))+coord_cartesian(ylim= c(0,1.5))


mad_values_Cell_2020_plot_matched

significance_data <- data_matched_Cell_2020 %>%
  summarise(
    TB1_TB2 = wilcox.test(mad_values_Cell_2020~ Buffering, data_matched_Cell_2020 = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(mad_values_Cell_2020 ~ Buffering, data_matched_Cell_2020 = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(mad_values_Cell_2020 ~ Buffering, data_matched_Cell_2020 = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )
significance_data









# pLI score
# Happloinsuffciency
pLI= read.table("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/gnomad.v2.1.1.lof_metrics.by_gene.txt", header = TRUE, sep="\t")
head(pLI)
colnames(pLI)
pLI_genes=cbind(pLI$gene, pLI$pLI)
colnames(pLI_genes)= c("Gene", "pLI")
pLI_genes=as.data.table(pLI_genes)

# label as buffered vs non buffered
pLI_genes= pLI_genes %>%
  mutate(Buffering = if_else(Gene %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(Gene %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))

pLI_genes= na.omit(pLI_genes)
#plot a normalized distribution
# make another colomn with % of the total
pLI_genes= pLI_genes%>%
  mutate_at(vars(-all_of(c("Gene","Buffering"))), as.numeric)
#  Calculate total counts of Buffering 
# Calculate total score for each variable
pLI_genes=na.omit(pLI_genes)
# make pLI score bin
# Define score ranges
score_breaks <- seq(0, 1, by = 0.1)
score_labels <- paste0(score_breaks[-length(score_breaks)], "-", score_breaks[-1])

# Cut scores into bins and calculate percentage of each variable within each bin
summary_data <- pLI_genes %>%
  mutate(score_range = cut(pLI, breaks = score_breaks, labels = score_labels)) %>%
  group_by(Buffering, score_range) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  complete(Buffering, score_range, fill = list(count = 0)) %>%
  group_by(Buffering) %>%
  mutate(percentage = count / sum(count) * 100)

# Plot histogram
pLI_plot=ggplot(summary_data, aes(x = score_range, y = percentage, fill = Buffering)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "pLI score range",
    y = "Percentage of genes") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_fill_manual(values = c("#045275","#089099", "#B7E6A5"))

pLI_plot

#median 
median(subset(pLI_genes, pLI_genes$Buffering=="TB_score=1")$pLI)
median(subset(pLI_genes, pLI_genes$Buffering=="TB_score=2")$pLI)
median(subset(pLI_genes, pLI_genes$Buffering=="TB_score=3")$pLI)

# Check the phaplo and ptriplo from 2022 data set
Collins_2022= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/non deduplicated data/final_TE_8_2_23/1-s2.0-S0092867422007887-mmc7.csv", header= TRUE)
head(Collins_2022)
Collins_2022=as.data.table(Collins_2022)
Collins_2022= Collins_2022 %>%
  mutate(Buffering = if_else(Gene %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(Gene %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))
Collins_2022= na.omit(Collins_2022)


# Cut scores into bins and calculate percentage of each variable within each bin
summary_data <- Collins_2022 %>%
  mutate(score_range = cut(pTriplo, breaks = score_breaks, labels = score_labels)) %>%
  group_by(Buffering, score_range) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  complete(Buffering, score_range, fill = list(count = 0)) %>%
  group_by(Buffering) %>%
  mutate(percentage = count / sum(count) * 100)
summary_data=na.omit(summary_data)
# Plot histogram
pTriplo_plot= ggplot(summary_data, aes(x = score_range, y = percentage, fill = Buffering)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Triplo score range",
    y = "Percentage of genes") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_fill_manual(values = c("#045275","#089099", "#B7E6A5"))
pTriplo_plot

# Plot histogram (pHaplo)

# Cut scores into bins and calculate percentage of each variable within each bin
summary_data <- Collins_2022 %>%
  mutate(score_range = cut(pHaplo, breaks = score_breaks, labels = score_labels)) %>%
  group_by(Buffering, score_range) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  complete(Buffering, score_range, fill = list(count = 0)) %>%
  group_by(Buffering) %>%
  mutate(percentage = count / sum(count) * 100)
summary_data=na.omit(summary_data)


pHaplo_plot= ggplot(summary_data, aes(x = score_range, y = percentage, fill = Buffering)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "pHaplo score range",
    y = "Percentage of genes") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_fill_manual(values = c("#045275","#089099", "#B7E6A5"))
pHaplo_plot

#median 
median(subset(Collins_2022, Collins_2022$Buffering=="TB_score=1")$pTriplo)
median(subset(Collins_2022, Collins_2022$Buffering=="TB_score=2")$pTriplo)
median(subset(Collins_2022, Collins_2022$Buffering=="TB_score=3")$pTriplo)


# analyse  if variation within different cell types exhibits buffering
Cenik_RNA= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/GSE65912_RNASeq.csv.gz")
Cenik_Ribo= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/GSE65912_RibosomeProfiling.csv.gz")
dim(Cenik_RNA)
Cenik_RNA[1:5, ]
dim(Cenik_Ribo)
Cenik_Ribo[1:5, ]
# there are more ribo  than RNA. 
# extract those that a have matching RNA to Ribo
# each GM is an indiviual why we do not have RNA seq for all
# extract columns common between the two data table
length(intersect(colnames(Cenik_Ribo), colnames(Cenik_RNA)))
# is it the same from Ribo base
# more samples have RNA match 
# make a list of common columns
Cenik_common= c("GM12878", "GM12891", "GM12892", "GM19238", "GM19239", "GM19240","X")
# extract columns which have the columns 
Cenik_RNA= Cenik_RNA %>%
  dplyr::select(starts_with(Cenik_common))
Cenik_Ribo= Cenik_Ribo %>%
  dplyr::select(starts_with(Cenik_common))

# take the cpm of the counts. 
Cenik_RNA_cpm=calculate_cpm(Cenik_RNA[,1:18])
Cenik_Ribo_cpm=calculate_cpm(Cenik_Ribo[,1:13])
# add the transcript name
Cenik_RNA= cbind(Cenik_RNA_cpm, Cenik_RNA$X)
Cenik_Ribo= cbind(Cenik_Ribo_cpm, Cenik_Ribo$X)

# avresult <- your_data %>%
Cenik_RNA= Cenik_RNA %>%
  mutate(
    GM12878_mean = dplyr::select(., starts_with("GM12878")) %>% rowMeans(na.rm = TRUE),
    GM12891_mean = dplyr::select(., starts_with("GM12891")) %>% rowMeans(na.rm = TRUE),
    GM12892_mean = dplyr::select(., starts_with("GM12892")) %>% rowMeans(na.rm = TRUE), 
    GM19238_mean = dplyr::select(., starts_with("GM19238")) %>% rowMeans(na.rm = TRUE),
    GM19239_mean = dplyr::select(., starts_with("GM19239")) %>% rowMeans(na.rm = TRUE),
    GM19240_mean = dplyr::select(., starts_with("GM19240")) %>% rowMeans(na.rm = TRUE)
  )
#
Cenik_Ribo= Cenik_Ribo %>%
  mutate(
    GM12878_mean_Ribo = dplyr::select(., starts_with("GM12878")) %>% rowMeans(na.rm = TRUE),
    GM12891_mean_Ribo = dplyr::select(., starts_with("GM12891")) %>% rowMeans(na.rm = TRUE),
    GM12892_mean_Ribo = dplyr::select(., starts_with("GM12892")) %>% rowMeans(na.rm = TRUE), 
    GM19238_mean_Ribo = dplyr::select(., starts_with("GM19238")) %>% rowMeans(na.rm = TRUE),
    GM19239_mean_Ribo = dplyr::select(., starts_with("GM19239")) %>% rowMeans(na.rm = TRUE),
    GM19240_mean_Ribo = dplyr::select(., starts_with("GM19240")) %>% rowMeans(na.rm = TRUE)
  )

# calculate the CV 
Cenik_RNA= Cenik_RNA %>% 
  mutate(CV_Cenik_RNA = apply(dplyr::select(.,20,21,22,23,24,25), 1, function(x) sd(x) / mean(x) * 100))
# 
Cenik_Ribo= Cenik_Ribo %>% 
  mutate(CV_Cenik_Ribo = apply(dplyr::select(.,15,16,17,18,19,20), 1, function(x) sd(x) / mean(x) * 100))

# merge
Cenik_RNA_Ribo= merge(Cenik_Ribo, Cenik_RNA, by.x= "Cenik_Ribo$X", by.y= "Cenik_RNA$X")
# modify the names by extracting the chracyer before -
Cenik_RNA_Ribo <-Cenik_RNA_Ribo %>%
  mutate(genes = sub("-.*", "", `Cenik_Ribo$X`))
# add bufered or non buffered
Cenik_RNA_Ribo = Cenik_RNA_Ribo %>% 
  mutate(Buffering = if_else(genes %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(genes %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))
# take the ratio
Cenik_RNA_Ribo = Cenik_RNA_Ribo %>% 
  mutate(Cenik_RNA_Ribo_ratio= CV_Cenik_Ribo/CV_Cenik_RNA)
Cenik_RNA_Ribo=as.data.table(Cenik_RNA_Ribo)
# match it with vCG
Cenik_RNA_Ribo[, treat := ifelse(Buffering == "TB_score=1", 1,
                                 ifelse(Buffering == "TB_score=2", 1, 0))]
Cenik_RNA_Ribo=na.omit(Cenik_RNA_Ribo)
m.out <- matchit(treat ~  CV_Cenik_RNA, Cenik_RNA_Ribo, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched


# plot
Cenik_2015= ggplot(data_matched, aes(y= Cenik_RNA_Ribo_ratio, x= Buffering, fill= Buffering))+geom_point(position=position_jitterdodge(jitter.width=0.2), alpha=0.5)+ ylim(0,2.5)+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( x = "Condition", y = "CV Ribo/CV RNA")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                                                                                                                                    axis.title = element_text(size = 12),
                                                                                                                                                                                                                                                                                                                    plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")
Cenik_2015= ggplot(data_matched, aes(y= Cenik_RNA_Ribo_ratio, x= Buffering, color= Buffering))+geom_point(position=position_jitterdodge(jitter.width=0.2), alpha= 0.5)+scale_color_manual(values = c("#045275", "#089099", "#B7E6A5"))+coord_cartesian(ylim = c(0, 5))+stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black")+labs( x = "Condition", y = "CV Ribo/CV RNA")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                                                                                                                                                                                                          axis.title = element_text(size = 12),
                                                                                                                                                                                                                                                                                                                                                                                          plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")


Cenik_2015
significance_data <- data_matched %>%
  summarise(
    TB1_TB2 = wilcox.test(Cenik_RNA_Ribo_ratio ~ Buffering, data_matched  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Cenik_RNA_Ribo_ratio ~ Buffering, data_matched = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Cenik_RNA_Ribo_ratio ~ Buffering, data_matched = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )
significance_data



#siRNA





Ribosomeprofiling = Ribo_clr
RNAsequencing =RNA_clr
Translationefficiency= TE_GBM_clr 
Ribosomeprofiling=as.data.frame(Ribosomeprofiling)
RNAsequencing=as.data.frame(RNAsequencing)
Translationefficiency=as.data.frame(Translationefficiency)
#Makes the transcripts the row names and removes the transcript column and an unnecessary index column
row.names(Ribosomeprofiling)= Ribo_clr$transcript
row.names(RNAsequencing)= RNA_clr$transcript
row.names(Translationefficiency)= TE_GBM_clr$transcript
Ribosomeprofiling$transcript= NULL
RNAsequencing$transcript =NULL
Translationefficiency$transcript =NULL



Ribosomeprofiling[1:5,1:5]



row_headers_rpl5 <- "RPL5"
row_headers_rplp1 <- "RPLP1"
row_headers_rplp2 <- "RPLP2"
row_headers_rps19 <- "RPS19"
row_headers_eif1 <- "EIF1"
row_headers_hnrnpc <- "HNRNPC"
row_headers_eif3e <- "EIF3E"
row_headers_dhx36 <- "DHX36"
row_headers_dhx9 <- "DHX9"
row_headers_ythdf3 <- "YTHDF3"
row_headers_setd2 <- "SETD2"
row_headers_eif2b5 <- "EIF2B5"
row_headers_ythdf1 <- "YTHDF1"
row_headers_arntl <- "ARNTL"
row_headers_mettl3 <- "METTL3"
row_headers_kbtbd8 <- "KBTBD8"
row_headers_tcof1 <- "TCOF1"
row_headers_fbl <- "FBL"
row_headers_ythdf1_2 <- "YTHDF1"
row_headers_ythdf2 <- "YTHDF2"
row_headers_mettl14 <- "METTL14"
row_headers_fkbp10 <- "FKBP10"
row_headers_pdcd4 <- "PDCD4"
row_headers_slc3a2 <- "SLC3A2"
row_headers_eif2a <- "EIF2A"
row_headers_ddx3x <- "DDX3X"



col_headers_rpl5 <- c("GSM2360175", "GSM2360176", "GSM2360179", "GSM2360180")

control_rpl5 <- c("notcontrol", "notcontrol", "control", "control")

#--------------------------------------------------------------

col_headers_rplp1 <- c("GSM3900208", "GSM3900209", "GSM3900210", "GSM3900211")

control_rplp1 <- c("control", "control", "control", "notcontrol")

#--------------------------------------------------------------

col_headers_rplp2 <- c("GSM3900208", "GSM3900209", "GSM3900210", "GSM3900212", "GSM3900213")

control_rplp2 <- c("control", "control", "control", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_rps19 <- c("GSM2360175", "GSM2360176", "GSM2360177", "GSM2360178")

control_rps19 <- c("notcontrol", "notcontrol", "control", "control")

#--------------------------------------------------------------

col_headers_eif1 <- c("GSM2327826", "GSM2327828")

control_eif1 <- c("notcontrol", "control")

#--------------------------------------------------------------

col_headers_hnrnpc <- c("GSM2204389", "GSM2204390", "GSM2204391", "GSM2204392", "GSM2204393", "GSM2204394", "GSM2204395", "GSM2204396")

control_hnrnpc <- c("control", "control", "control", "control", "notcontrol", "notcontrol", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_eif3e <- c("GSM3762993", "GSM3762994", "GSM3762995", "GSM3762996", "GSM3762997", "GSM3762998")

control_eif3e <- c("notcontrol", "notcontrol", "notcontrol", "control", "control", "control")

#--------------------------------------------------------------

col_headers_dhx36 <- c("GSM2817679", "GSM2817680", "GSM2817681", "GSM2817682", "GSM2817683", "GSM2817684", "GSM2817685", "GSM2817686")

control_dhx36 <- c("control", "control", "control", "control", "control", "notcontrol", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_dhx9 <- c("GSM2817679", "GSM2817680", "GSM2817681", "GSM2817682", "GSM2817683", "GSM2817687", "GSM2817688", "GSM2817689")

control_dhx9 <- c("control", "control", "control", "control", "control", "notcontrol", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_ythdf3 <- c("GSM3944607", "GSM3944615", "GSM3944616", "GSM3944617")

control_ythdf3 <- c("control", "notcontrol", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_setd2 <- c("GSM3450419", "GSM3450420", "GSM3450424")

control_setd2 <- c("control", "control", "notcontrol")

#--------------------------------------------------------------

col_headers_eif2b5 <- c("GSM2883304", "GSM2883305", "GSM2883306", "GSM2883307", "GSM2883313", "GSM2883314", "GSM2883315", "GSM2883320", "GSM2883321", "GSM2883322", "GSM2883323")

control_eif2b5 <- c("control", "notcontrol", "control", "notcontrol", "notcontrol", "control", "notcontrol", "control", "notcontrol", "control", "notcontrol")

#--------------------------------------------------------------

col_headers_ythdf1 <- c("GSM4054749", "GSM4054750", "GSM4054751")

control_ythdf1 <- c("notcontrol", "notcontrol", "control")

#--------------------------------------------------------------

col_headers_arntl <- c("GSM1371443", "GSM1371444", "GSM1371445", "GSM1371446", "GSM1371447", "GSM1371448", "GSM1371449", "GSM1371450", "GSM1371451", "GSM1371452", "GSM1371453", "GSM1371455", "GSM1371456", "GSM1371457", "GSM1371458", "GSM1371459", "GSM1371460", "GSM1371461", "GSM1371462", "GSM1371463", "GSM1371464", "GSM1371465", "GSM1371466", "GSM1371467", "GSM1371468", "GSM1371469", "GSM1371470", "GSM1371471", "GSM1371472", "GSM1371473", "GSM1371474", "GSM1371475", "GSM1371476", "GSM1371477", "GSM1371478", "GSM1371479", "GSM1371480", "GSM1371481", "GSM1371482", "GSM1371483", "GSM1371484", "GSM1371485", "GSM1371486", "GSM1371487", "GSM1371488", "GSM1371489", "GSM1371490")

control_arntl <- c("control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "control", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_mettl3 <- c("GSM2602082", "GSM2602083", "GSM2602084", "GSM2602088", "GSM2602089", "GSM2602090")

control_mettl3 <- c("control", "control", "control", "notcontrol", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_kbtbd8 <- c("GSM1782874", "GSM1782876", "GSM1782879")

control_kbtbd8 <- c("control", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_tcof1 <- c("GSM1782874", "GSM1782877", "GSM1782880")

control_tcof1 <- c("control", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_fbl <- c("GSM2825128", "GSM2825129", "GSM2825131", "GSM2825132")

control_fbl <- c("control", "notcontrol", "notcontrol", "control")

#--------------------------------------------------------------

col_headers_ythdf1_2 <- c("GSM3944607", "GSM3944610", "GSM3944611")

control_ythdf1_2 <- c("control", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_ythdf2 <- c("GSM3944607", "GSM3944614")

control_ythdf2 <- c("control", "notcontrol")

#--------------------------------------------------------------

col_headers_mettl14 <- c("GSM3450419", "GSM3450420", "GSM3450421", "GSM3450422")

control_mettl14 <- c("control", "control", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_fkbp10 <- c("GSM3718424", "GSM3718425", "GSM3718426", "GSM3718427")

control_fkbp10 <- c("control", "control", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_pdcd4 <- c("GSM4110755", "GSM4110757", "GSM4110758", "GSM4110759", "GSM4110760")

control_pdcd4 <- c("notcontrol", "notcontrol", "control", "control", "control")

#--------------------------------------------------------------

col_headers_slc3a2 <- c("GSM1446854", "GSM1446856")

control_slc3a2 <- c("control", "notcontrol")

#--------------------------------------------------------------

col_headers_eif2a <- c("GSM5291928", "GSM5291929", "GSM5291930", "GSM5291931", "GSM5291932", "GSM5291933")

control_eif2a <- c("control", "control", "notcontrol", "notcontrol", "notcontrol", "notcontrol")

#--------------------------------------------------------------

col_headers_ddx3x <- c("GSM4258310", "GSM4258311", "GSM4258315", "GSM4258316", "GSM4258317", "GSM4258319")

control_ddx3x <- c("control", "control", "notcontrol", "notcontrol", "notcontrol", "notcontrol")





#Function for Processing Gene Data
create_gene_data <- function(gene_name, col_headers, row_headers, Ribosomeprofiling, RNAsequencing, Translationefficiency, control) {
  
  # Create empty vectors to store the data
  Riboprofdata <- vector("list", length(col_headers))
  RNAseqdata <- vector("list", length(col_headers))
  Translationdata <- vector("list", length(col_headers))
  
  # Loop over each column header
  for (i in seq_along(col_headers)) {
    # Check lengths of data before subsetting
    if (nrow(Ribosomeprofiling) != nrow(RNAsequencing) || nrow(RNAsequencing) != nrow(Translationefficiency)) {
      stop("Input data matrices (Ribosomeprofiling, RNAsequencing, Translationefficiency) have different numbers of rows.")
    }
    
    Riboprofdata[[i]] <- Ribosomeprofiling[row_headers, col_headers[i]]
    #print(Ribosomeprofiling[row_headers, col_headers[i]])
    RNAseqdata[[i]] <- RNAsequencing[row_headers, col_headers[i]]
    #print(RNAsequencing[row_headers, col_headers[i]])
    Translationdata[[i]] <- Translationefficiency[row_headers, col_headers[i]]
    #print(Translationdata[row_headers, col_headers[i]])
  }
  
  # Create named list for data frame columns
  df_list <- setNames(
    list(
      unlist(RNAseqdata),
      unlist(Riboprofdata),
      unlist(Translationdata),
      control
    ),
    c(
      paste0("RNAseq_", gene_name),
      paste0("Riboprof_", gene_name),
      paste0("Transeff_", gene_name),
      "Control"
    )
  )
  
  # Convert list to data frame
  data <- as.data.frame(df_list)
  
  return(data)
}



data_rpl5 <- create_gene_data("rpl5", col_headers_rpl5, row_headers_rpl5, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_rpl5)
#--------------------------------------------------------------
data_rplp1 <- create_gene_data("rplp1", col_headers_rplp1, row_headers_rplp1, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_rplp1)
#--------------------------------------------------------------
data_rplp2 <- create_gene_data("rplp2", col_headers_rplp2, row_headers_rplp2, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_rplp2)
#--------------------------------------------------------------
data_rps19 <- create_gene_data("rps19", col_headers_rps19, row_headers_rps19, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_rps19)
#--------------------------------------------------------------
data_eif1 <- create_gene_data("eif1", col_headers_eif1, row_headers_eif1, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_eif1)
#--------------------------------------------------------------
data_hnrnpc <- create_gene_data("hnrnpc", col_headers_hnrnpc, row_headers_hnrnpc, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_hnrnpc)
#--------------------------------------------------------------
data_eif3e <- create_gene_data("eif3e", col_headers_eif3e, row_headers_eif3e, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_eif3e)
#--------------------------------------------------------------
data_dhx36 <- create_gene_data("dhx36", col_headers_dhx36, row_headers_dhx36, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_dhx36)
#--------------------------------------------------------------
data_dhx9 <- create_gene_data("dhx9", col_headers_dhx9, row_headers_dhx9, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_dhx9)
#--------------------------------------------------------------
data_ythdf3 <- create_gene_data("ythdf3", col_headers_ythdf3, row_headers_ythdf3, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_ythdf3)
#--------------------------------------------------------------
data_setd2 <- create_gene_data("setd2", col_headers_setd2, row_headers_setd2, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_setd2)
#--------------------------------------------------------------
data_eif2b5 <- create_gene_data("eif2b5", col_headers_eif2b5, row_headers_eif2b5, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_eif2b5)
#--------------------------------------------------------------
data_ythdf1 <- create_gene_data("ythdf1", col_headers_ythdf1, row_headers_ythdf1, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_ythdf1)
#--------------------------------------------------------------
data_arntl <- create_gene_data("arntl", col_headers_arntl, row_headers_arntl, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_arntl)
#--------------------------------------------------------------
data_mettl3 <- create_gene_data("mettl3", col_headers_mettl3, row_headers_mettl3, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_mettl3)
#--------------------------------------------------------------
data_kbtbd8 <- create_gene_data("kbtbd8", col_headers_kbtbd8, row_headers_kbtbd8, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_kbtbd8)
#--------------------------------------------------------------
data_tcof1 <- create_gene_data("tcof1", col_headers_tcof1, row_headers_tcof1, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_tcof1)
#--------------------------------------------------------------
data_fbl <- create_gene_data("fbl", col_headers_fbl, row_headers_fbl, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_fbl)
#--------------------------------------------------------------
data_ythdf1_2 <- create_gene_data("ythdf1", col_headers_ythdf1_2, row_headers_ythdf1_2, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_ythdf1_2)
#--------------------------------------------------------------
data_ythdf2 <- create_gene_data("ythdf2", col_headers_ythdf2, row_headers_ythdf2, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_ythdf2)
#--------------------------------------------------------------
data_mettl14 <- create_gene_data("mettl14", col_headers_mettl14, row_headers_mettl14, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_mettl14)
#--------------------------------------------------------------
data_fkbp10 <- create_gene_data("fkbp10", col_headers_fkbp10, row_headers_fkbp10, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_fkbp10)
#--------------------------------------------------------------
data_pdcd4 <- create_gene_data("pdcd4", col_headers_pdcd4, row_headers_pdcd4, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_pdcd4)
#--------------------------------------------------------------
data_slc3a2 <- create_gene_data("slc3a2", col_headers_slc3a2, row_headers_slc3a2, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_slc3a2)
#--------------------------------------------------------------
data_eif2a <- create_gene_data("eif2a", col_headers_eif2a, row_headers_eif2a, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_eif2a)
#--------------------------------------------------------------
data_ddx3x <- create_gene_data("ddx3x", col_headers_ddx3x, row_headers_ddx3x, Ribosomeprofiling, RNAsequencing, Translationefficiency, control_ddx3x)




#Calculate the Mean of The Data and Center it at x = 0
# Create a function to calculate summary data
calculate_summary <- function(data, gene) {
  data %>%
    group_by(Control) %>%
    summarise(
      mean_RNAseq = mean(get(paste0("RNAseq_", gene))),
      mean_Riboprof = mean(get(paste0("Riboprof_", gene))),
      mean_Transeff = mean(get(paste0("Transeff_", gene)))
    ) %>%
    ungroup() %>%
    mutate(
      # Normalizes the data to the control grouup value, making it better to visualize the knockdown perturbaration
      mean_RNAseq = ifelse(Control == "notcontrol", mean_RNAseq - mean(mean_RNAseq[Control == "control"]), -1 * mean(mean_RNAseq[Control == "control"])),
      mean_Riboprof = ifelse(Control == "notcontrol", mean_Riboprof - mean(mean_Riboprof[Control == "control"]), -1 * mean(mean_Riboprof[Control == "control"])),
      mean_Transeff = ifelse(Control == "notcontrol", mean_Transeff - mean(mean_Transeff[Control == "control"]), -1 * mean(mean_Transeff[Control == "control"])),
      mean_RNAseq = ifelse(Control == "control", 0, mean_RNAseq),
      mean_Riboprof = ifelse(Control == "control", 0, mean_Riboprof),
      mean_Transeff = ifelse(Control == "control", 0, mean_Transeff)
    )
}

summary_data_rpl5 <- calculate_summary(data_rpl5, "rpl5")
summary_data_rplp1 <- calculate_summary(data_rplp1, "rplp1")
summary_data_rplp2 <- calculate_summary(data_rplp2, "rplp2")
summary_data_rps19 <- calculate_summary(data_rps19, "rps19")
summary_data_eif1 <- calculate_summary(data_eif1, "eif1")
summary_data_hnrnpc <- calculate_summary(data_hnrnpc, "hnrnpc")
summary_data_eif3e <- calculate_summary(data_eif3e, "eif3e")
summary_data_dhx36 <- calculate_summary(data_dhx36, "dhx36")
summary_data_dhx9 <- calculate_summary(data_dhx9, "dhx9")
summary_data_ythdf3 <- calculate_summary(data_ythdf3, "ythdf3")
summary_data_setd2 <- calculate_summary(data_setd2, "setd2")
summary_data_eif2b5 <- calculate_summary(data_eif2b5, "eif2b5")
summary_data_ythdf1 <- calculate_summary(data_ythdf1, "ythdf1")
summary_data_arntl <- calculate_summary(data_arntl, "arntl")
summary_data_mettl3 <- calculate_summary(data_mettl3, "mettl3")
summary_data_kbtbd8 <- calculate_summary(data_kbtbd8, "kbtbd8")
summary_data_tcof1 <- calculate_summary(data_tcof1, "tcof1")
summary_data_fbl <- calculate_summary(data_fbl, "fbl")
summary_data_ythdf1_2 <- calculate_summary(data_ythdf1_2, "ythdf1")
summary_data_ythdf2 <- calculate_summary(data_ythdf2, "ythdf2")
summary_data_mettl14 <- calculate_summary(data_mettl14, "mettl14")
summary_data_fkbp10 <- calculate_summary(data_fkbp10, "fkbp10")
summary_data_pdcd4 <- calculate_summary(data_pdcd4, "pdcd4")
summary_data_slc3a2 <- calculate_summary(data_slc3a2, "slc3a2")
summary_data_eif2a <- calculate_summary(data_eif2a, "eif2a")
summary_data_ddx3x <- calculate_summary(data_ddx3x, "ddx3x")

#Standardize the Columns Naming
colnames(summary_data_rpl5) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_rplp1) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_rplp2) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_rps19) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_eif1) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_hnrnpc) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_eif3e) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_dhx36) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_dhx9) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_ythdf3) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_setd2) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_eif2b5) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_ythdf1) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_arntl) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_mettl3) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_kbtbd8) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_tcof1) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_fbl) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_ythdf1_2) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_ythdf2) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_mettl14) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_fkbp10) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_pdcd4) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_slc3a2) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_eif2a) <- c("Control", "RNAseq", "Riboprof", "Transeff")
colnames(summary_data_ddx3x) <- c("Control", "RNAseq", "Riboprof", "Transeff")

#Combine All the Data to One Table
combined_data <- rbind(
  cbind(summary_data_rpl5, Gene = "RPL5", Bufferingscore = "1"),
  cbind(summary_data_rplp1, Gene = "RPLP1", Bufferingscore = "1"),
  cbind(summary_data_rplp2, Gene = "RPLP2", Bufferingscore = "1"),
  cbind(summary_data_rps19, Gene = "RPS19", Bufferingscore = "1"),
  cbind(summary_data_eif1, Gene = "EIF1", Bufferingscore = "1"),
  cbind(summary_data_hnrnpc, Gene = "HNRNPC", Bufferingscore = "1"),
  cbind(summary_data_eif3e, Gene = "EIF3E", Bufferingscore = "2"),
  cbind(summary_data_dhx36, Gene = "DHX36", Bufferingscore = "3"),
  cbind(summary_data_dhx9, Gene = "DHX9", Bufferingscore = "3"),
  cbind(summary_data_ythdf3, Gene = "YTHDF3", Bufferingscore = "3"),
  cbind(summary_data_setd2, Gene = "SETD2", Bufferingscore = "3"),
  cbind(summary_data_eif2b5, Gene = "EIF2B5", Bufferingscore = "3"),
  cbind(summary_data_ythdf1, Gene = "YTHDF1", Bufferingscore = "3"),
  cbind(summary_data_arntl, Gene = "ARNTL", Bufferingscore = "3"),
  cbind(summary_data_mettl3, Gene = "METTL3", Bufferingscore = "3"),
  cbind(summary_data_kbtbd8, Gene = "KBTBD8", Bufferingscore = "3"),
  cbind(summary_data_tcof1, Gene = "TCOF1", Bufferingscore = "3"),
  cbind(summary_data_fbl, Gene = "FBL", Bufferingscore = "3"),
  cbind(summary_data_ythdf1_2, Gene = "YTHDF1_2", Bufferingscore = "3"),
  cbind(summary_data_ythdf2, Gene = "YTHDF2", Bufferingscore = "3"),
  cbind(summary_data_mettl14, Gene = "METTL14", Bufferingscore = "3"),
  cbind(summary_data_fkbp10, Gene = "FKBP10", Bufferingscore = "3"),
  cbind(summary_data_pdcd4, Gene = "PDCD4", Bufferingscore = "3"),
  cbind(summary_data_slc3a2, Gene = "SLC3A2", Bufferingscore = "3"),
  cbind(summary_data_eif2a, Gene = "EIF2A", Bufferingscore = "3"),
  cbind(summary_data_ddx3x, Gene = "DDX3X", Bufferingscore = "3"))


#Filter combined_data based on conditions
filtered_combined_data <- combined_data %>%
  filter(!(Control == "notcontrol" & RNAseq > 0))

#If the not control is greater than 0, then the knockdown did not work, exclude it from the graph


# Translation Efficiency Average vs RNA Abundance Graph ~ Facet Wrappped

ggplot(filtered_combined_data, aes(x = RNAseq, y = Transeff, color = Bufferingscore)) +
  geom_point() +
  geom_smooth(aes(group = interaction(Bufferingscore, Gene)), method = 'lm', formula = y ~ x, se = FALSE) +
  scale_x_reverse(expand = c(0.1, 0.1)) +  # Reverse the x-axis
  scale_y_continuous(expand = c(0.1, 0.1)) +  # Zoom out on the y-axis
  labs(title = "Translation Efficiency Average vs RNA Abundance Average",
       x = "RNA Abundance",
       y = "Translation Efficiency") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "top"  # Position the legend at the top
  )  +
  facet_wrap(~ Bufferingscore, ncol = 3)

siRNA_RNA_Ribo= ggplot(filtered_combined_data, aes(x = RNAseq, y = Riboprof, color = Bufferingscore)) +
  geom_point() +
  geom_smooth(aes(group = interaction(Bufferingscore, Gene)), method = 'lm', formula = y ~ x, se = FALSE) +
  scale_x_reverse(expand = c(0.1, 0.1)) +  # Reverse the x-axis
  scale_y_continuous(expand = c(0.1, 0.1)) +  # Zoom out on the y-axis
  labs(
    x = "Δ mRNA abundance",
    y = "Δ Ribosome Occuppancy") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines")) +scale_color_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  facet_wrap(~ Bufferingscore, ncol = 3)+theme(legend.position ="none")
siRNA_RNA_Ribo





library(gghalves)

# Figure 5 
Array.ribo = Ribo("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/all.ribo")
rnaseq_CDS <- get_rnaseq(ribo.object = Array.ribo,
                         tidy        = TRUE,
                         region      = "CDS",
                         compact     = FALSE,
                         
)


Array_RNA= dcast(as.data.table(rnaseq_CDS), transcript ~ experiment)
head(Array_RNA)
dim(Array_RNA)
# take counts per million. 
# take cpm of counts
calculate_cpm <- function(expr_matrix) {
  # Calculate the column sums (total counts per sample)
  col_sums <- colSums(expr_matrix)
  
  # Calculate CPM
  cpm <- sweep(expr_matrix, 2, col_sums, FUN = "/") * 1e6
  
  return(cpm)
}

all_counts_cpm= calculate_cpm(Array_RNA[,2:9])
all_counts_cpm$transcript=Array_RNA$transcript
head(all_counts_cpm)
dim(all_counts_cpm)
all_counts_cpm[1:5,]
# extract the columns 
all_counts_cpm$gene_symbol <- sapply(str_split(all_counts_cpm$transcript, "\\|"), `[`, 6)


all_counts_cpm=all_counts_cpm[,c(1:8,10)]
all_counts_cpm=as.data.table(all_counts_cpm)

head(all_counts_cpm, 10)
dim(all_counts_cpm)
#19736 genes
 
sum(rowSums(all_counts_cpm[,1:8])<10) # 9015 genes

# remove genes fi cpm is less than 10 in more than 6 samples
# Keep genes where genes have more than 10 counts in atleast 6 samples
all_counts_cpm_filtered <- all_counts_cpm |> 
  filter(rowSums(all_counts_cpm[, 1:8] > 10) >= 6)
dim(all_counts_cpm_filtered)
# 8754  remaining

#all_counts_cpm_filtered <- all_counts_cpm_filtered %>%
 # rowwise() %>%
  #mutate(slope = coef(lm(c(AR_1_RF_POLY, AR_2_RF_POLY, AR_3_RF_POLY, AR_4_RF_POLY) ~ c(AR_1_RF, AR_2_RF, AR_3_RF, AR_4_RF)))[2]) %>% # [2] indicates slope: 2nd coefficient [1] being the intercept
 # ungroup()

# Also label it as TB score 1, 2,3 
all_counts_cpm_filtered = all_counts_cpm_filtered %>% 
  mutate(Buffering = if_else(gene_symbol %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(gene_symbol %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))

# # calculate CV
cv <- function(x) {
  (sd(x) / mean(x)) * 100
}


all_counts_cpm_filtered= all_counts_cpm_filtered %>%
  mutate(CV_RNA_Array= apply(all_counts_cpm_filtered[, c(1,3,5,7)], 1, cv))

head(all_counts_cpm_filtered)



# assesing the goodness of fit
m <- all_counts_cpm_filtered %>%
  rowwise() %>%
  mutate(
    model = list(lm( y ~ x,
                     data= data.frame(y= c(AR_1_RF_POLY, AR_2_RF_POLY, AR_3_RF_POLY, AR_4_RF_POLY),
                                      x= c(AR_1_RF, AR_2_RF, AR_3_RF, AR_4_RF)))),
    # Extract slope (coefficient of x)
    slope = coef(model)[2],  
    
    # Extract R-squared value
    r2 = summary(model)$r.squared  
  ) %>%
  ungroup() 

# distribution of r2 
hist(m$r2)

# assessing the goodness of fit chi square test

m <- m %>%
  rowwise() %>%
  mutate(expected_values = list(fitted(model))) %>%
  ungroup()


m_chi_sq <- m %>%
  rowwise() %>%
  mutate(
    # Extract observed values dynamically for each row
    observed_values = list(c(AR_1_RF_POLY, AR_2_RF_POLY, AR_3_RF_POLY, AR_4_RF_POLY)),
    
    # Compute residuals (Observed - Expected)
    residuals = list(observed_values - expected_values),  # Element-wise subtraction
    
    # Compute Chi-square statistic (for each row)
    chi_sq = sum((residuals^2) / expected_values, na.rm = TRUE),  # Avoid division by zero
    
    # Compute degrees of freedom
    df_residual = length(expected_values) - 2,  # df = number of observations - estimated parameters
    
    # Compute p-value (only if df_residual > 0)
    p_value = ifelse(df_residual > 0 & all(expected_values > 0),
                     pchisq(chi_sq, df = df_residual, lower.tail = FALSE),
                     NA_real_)  # Assign NA if invalid
  ) %>%
  ungroup() 

ggplot(m_chi_sq, aes(x = chi_sq, y = p_value, color = r2)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +  # Color shows R² values
  labs(title = "Chi-Square vs. p-value (Colored by R²)",
       x = "Chi-Square Statistic",
       y = "p-value",
       color = "R² Value") +
  theme_minimal()



# filetr samples with low r2 
summary(m_chi_sq$r2) #median is 0.27
all_counts_cpm_filtered_r2_cv= subset(m_chi_sq, m_chi_sq$r2 > median(m_chi_sq$r2) & m_chi_sq$CV_RNA_Array > median(m_chi_sq$CV_RNA_Array) )
dim(m_chi_sq) # 8754
dim(all_counts_cpm_filtered_r2_cv)
#2257 genes
ggplot(all_counts_cpm_filtered_r2_cv, aes(y= slope, x= Buffering, fill= Buffering))+ geom_boxplot()+geom_point()+ggtitle("filter: r2-_cv_median")+scale_y_continuous(breaks = seq(-1,2, by = 1))


# 
Translating_fraction=ggplot(all_counts_cpm_filtered_r2_cv, aes(y= slope, x= Buffering, fill= Buffering)) + 
geom_boxplot(width = .6, fill=NA, color=c("#045275", "#089099", "#B7E6A5"), linewidth=1.5) + 
  scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( y = "Slope")+theme(axis.text = element_text(size = 12),
                                                                                                       axis.title = element_text(size = 12),
                                                                                                       plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+
  coord_cartesian(ylim=c(-1.8,2.5))


Translating_fraction





Translating_fraction
# significance
all_counts_cpm_filtered_r2_cv %>%
  summarise(
    TB1_TB2 = wilcox.test(slope ~ Buffering, all_counts_cpm_filtered_r2_cv  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(slope ~ Buffering, all_counts_cpm_filtered_r2_cv  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(slope ~ Buffering, all_counts_cpm_filtered_r2_cv  = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )






# if matced with CV
all_counts_cpm_filtered_r2_cv=as.data.table(all_counts_cpm_filtered_r2_cv)
all_counts_cpm_filtered_r2_cv[, treat := ifelse(Buffering == "TB_score=1", 1,
                                                ifelse(Buffering == "TB_score=2", 1, 0))]
m.out <- matchit(treat ~  CV_RNA_Array, all_counts_cpm_filtered_r2_cv, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched
ggplot(data_matched, aes(y= slope, x= Buffering, fill= Buffering))+geom_boxplot()+geom_point()+ggtitle("filter: r2_cv_median,matched with CV")+scale_y_continuous(breaks = seq(-1,2, by = 1))
dim(data_matched)
# tere are 200 genes # 48 TB score =1 , 52 TB score =2
nrow(data_matched[data_matched$Buffering=="TB_score=1",])

Translating_fraction_matched=ggplot(data_matched, aes(y= slope, x= Buffering, fill= Buffering, color=Buffering)) + 
  geom_boxplot(width = .6, fill=NA, color=c("#045275", "#089099", "#B7E6A5"),linewidth=1.5) + 
  geom_jitter(width = .05, alpha = .3, shape= 21)+
  scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+scale_color_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( y = "Slope")+theme(axis.text = element_text(size = 12),
                                                                                          axis.title = element_text(size = 12),
                                                                                          plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+
  coord_cartesian(ylim=c(-2.5,2.5))




Translating_fraction_matched


data_matched %>%
  summarise(
    TB1_TB2 = wilcox.test(slope ~ Buffering, data_matched  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(slope ~ Buffering, data_matched  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(slope ~ Buffering, data_matched  = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )

# show this at the ribosome occupancy level
# # calculate the CV of ribo which is filtered for CV at mRNA
rc_CDS <- get_region_counts(ribo.object    = Array.ribo,
                            range.lower = 25,
                            range.upper = 32,
                            tidy       = TRUE,
                            transcript = FALSE,
                            region     = "CDS",
                            compact    = FALSE
)
head(rc_CDS)
Array_Ribo= dcast(as.data.table(rc_CDS), transcript ~ experiment)
head(Array_Ribo)
dim(Array_Ribo)

# the counts are duplicated so remove 
Array_Ribo=Array_Ribo[,c(1,3,5,7,9)]
# take counts per million. 
# take cpm of counts
dim(Array_Ribo)
all_counts_cpm_ribo= calculate_cpm(Array_Ribo[,2:5])
all_counts_cpm_ribo$transcript=Array_Ribo$transcript
head(all_counts_cpm_ribo)
dim(all_counts_cpm_ribo)
all_counts_cpm_ribo[1:5,]
# extract the columns 
all_counts_cpm_ribo$gene_symbol <- sapply(str_split(all_counts_cpm_ribo$transcript, "\\|"), `[`, 6)

dim(all_counts_cpm_ribo)
all_counts_cpm_ribo=all_counts_cpm_ribo[,c(1:4,6)]
all_counts_cpm_ribo=as.data.table(all_counts_cpm_ribo)
all_counts_cpm_ribo[1:5,]

# subset the gene set.

test_Ribo2= all_counts_cpm_ribo[gene_symbol %in% all_counts_cpm_filtered_r2_cv$gene_symbol]
test_Ribo2=as.data.table(test_Ribo2)
test_Ribo2 <- test_Ribo2 %>%
  mutate(CV_Ribo2 = apply(dplyr::select(., 1, 2, 3, 4), 1, function(x) sd(x) / mean(x) * 100))
# add buffereing score
test_Ribo2 <- test_Ribo2 %>%
  mutate(Buffering = if_else(gene_symbol %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(gene_symbol %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))

# plot the ggplot
test_Ribo2_longer= pivot_longer(test_Ribo2, cols = c("CV_Ribo2"), names_to = "Variable", values_to = "Value")

sum(is.na(test_Ribo2$CV_Ribo2))
#there six with NA
#why
test_Ribo2[c(which(is.na(test_Ribo2$CV_Ribo2))),]
# remove these genes
test_Ribo2= na.omit(test_Ribo2)

Ribo_MAD= ggplot(test_Ribo2_longer, aes(x= Buffering, y= Value, fill=Buffering)) + 
  geom_boxplot(width = .6, fill=NA, color=c("#045275", "#089099", "#B7E6A5"), linewidth=1.5) + 
scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( y = "CV of Ribosome occupancy")+theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 12),
                  plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                  panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,50))

Ribo_MAD

library(ggforce)
Ribo_MAD= ggplot(test_Ribo2_longer, aes(x= Buffering, y= Value, color=Buffering)) + 
  geom_sina(alpha = 0.5, maxwidth = 0.8) +
  scale_color_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( y = "CV of Ribosome occupancy")+stat_summary(fun = median, geom = "crossbar", 
                                                                                                                      width = 0.2, fill = "black", color = "black")+theme(axis.text = element_text(size = 12),
                                                                                                             axis.title = element_text(size = 12),
                                                                                                             plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                             panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,50))
  
test_Ribo2_longer %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, test_Ribo2_longer  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering,  test_Ribo2_longer  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering,  test_Ribo2_longer = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )
# do the same for matched set
test_Ribo2= all_counts_cpm_ribo[gene_symbol %in% data_matched$gene_symbol]
test_Ribo2=as.data.table(test_Ribo2)
dim(test_Ribo2)
test_Ribo2 <- test_Ribo2 %>%
  mutate(CV_Ribo2 = apply(dplyr::select(., 1, 2, 3, 4), 1, function(x) sd(x) / mean(x) * 100))
# add buffereing score
test_Ribo2 <- test_Ribo2 %>%
  mutate(Buffering = if_else(gene_symbol %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(gene_symbol %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))

# plot the ggplot
test_Ribo2_longer_matched= pivot_longer(test_Ribo2, cols = c("CV_Ribo2"), names_to = "Variable", values_to = "Value")
test_Ribo2_longer_matched$Buffering=factor(test_Ribo2_longer_matched$Buffering)
Ribo_MAD_matched= ggplot(test_Ribo2_longer_matched, aes(x= Buffering, y= Value, fill=Buffering)) + 
  geom_boxplot(width = .6, outlier.shape = NA, fill=NA, color=c("#045275", "#089099", "#B7E6A5"), linewidth=1.5) + 
  geom_jitter(width = .05, alpha = .3, shape= 21)+
  scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( y = "CV of Ribosome occupancy")+theme(axis.text = element_text(size = 12),
                                                                                                             axis.title = element_text(size = 12),
                                                                                                             plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                             panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,50))

Ribo_MAD_matched=ggplot(test_Ribo2_longer_matched, aes(x= Buffering, y= Value, color=Buffering)) + 
  geom_sina(alpha = 0.5, maxwidth = 0.8) +
  scale_color_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( y = "CV of Ribosome occupancy")+stat_summary(fun = median, geom = "crossbar", 
                                                                                                                     width = 0.2, fill = "black", color = "black")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                         axis.title = element_text(size = 12),
                                                                                                                                                                         plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                                                                         panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,50))






Ribo_MAD_matched

test_Ribo2_longer_matched %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, test_Ribo2_longer_matched  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering,  test_Ribo2_longer_matched  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering,  test_Ribo2_longer_matched = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )




# Is there any difference in fraction of genes entering the translational pool
# Differential expression of total and polysomal RNA
# make a box plot of the log FC according to the buffering score. 
# differential analysis

head(rnaseq_CDS)

CRISPR_RNA= dcast(as.data.table(rnaseq_CDS), transcript ~ experiment)

head(CRISPR_RNA)
dim(CRISPR_RNA)
# extract the columns 
CRISPR_RNA$gene_symbol <- sapply(str_split(CRISPR_RNA$transcript, "\\|"), `[`, 6)
CRISPR_RNA= CRISPR_RNA[,-1]
exptype <- factor(c("CRISPR_Total", "CRISPR_Polysome", "CRISPR_Total", "CRISPR_Polysome","CRISPR_Total", "CRISPR_Polysome","CRISPR_Total", "CRISPR_Polysome"))

CRISPR_RNA= na.omit(CRISPR_RNA)
y <- DGEList(counts=CRISPR_RNA[,-9],group=exptype, genes = CRISPR_RNA[,9])
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~0+exptype)
colnames(design) <- levels(exptype)
y <- estimateDisp(y,design)

plotBCV(y)
plotMDS(y)
fit <- glmQLFit(y,design)
my.contrasts <- makeContrasts(
  PolyvsTotal = (CRISPR_Polysome - CRISPR_Total),
  levels = design
)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"PolyvsTotal"])  
summary(decideTests(qlf, p.value = 0.05, adjust.method = "fdr"))





# match it with expression of RNA . Therefore add cpm value

polysomal_total_matched= cbind(qlf$table$logFC, qlf$table$PValue, qlf$genes, qlf$table$logCPM)
polysomal_total_matched= as.data.table(polysomal_total_matched)

polysomal_total_matched <- polysomal_total_matched %>%
  mutate(Buffering = if_else(gene_symbol %in% Top_buffered_human$transcript, "TB_score=1",
                             if_else(gene_symbol %in% Top2_buffered_human$transcript, "TB_score=2", "TB_score=3")))
#unmatched
polyvsTotalFC= ggplot(polysomal_total_matched, aes(y= `qlf$table$logFC`, x= Buffering, fill= Buffering))+geom_boxplot()+scale_fill_manual(values = c("TB_score=1"="#045275","TB_score=2"="#089099","TB_score=3"="#B7E6A5"))+labs(y = "Log2FC  translating RNA vs Total fraction")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                   axis.title = element_text(size = 12),
                                                                                                                                                                                                                                                   plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")



polyvsTotalFC= ggplot(polysomal_total_matched, aes(y= `qlf$table$logFC`, x= Buffering, fill= Buffering))+geom_violin(width=.6, alpha=.5)+stat_summary(fun = median, geom = "crossbar", 
                                                                                                                                                      width = 0.2, fill = "black", color = "black")+scale_fill_manual(values = c("TB_score=1"="#045275","TB_score=2"="#089099","TB_score=3"="#B7E6A5"))+labs(y = "Log2FC  translating RNA vs Total fraction")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                   axis.title = element_text(size = 12),
                                                                                                                                                                                                                                                   plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")


polyvsTotalFC
# The fraction of RNA that enters the polysomal pool is higher in buffered set compared to non buffered set.
significance_data <- polysomal_total_matched %>%
  summarise(
    TB1_TB2 = wilcox.test(`qlf$table$logFC` ~ Buffering, polysomal_total_matched  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(`qlf$table$logFC` ~ Buffering, polysomal_total_matched = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(`qlf$table$logFC` ~ Buffering, polysomal_total_matched = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )
polysomal_total_matched= as.data.table(polysomal_total_matched)
polysomal_total_matched[, treat := ifelse(Buffering == "TB_score=1", 1,
                                          ifelse(Buffering == "TB_score=2", 1, 0))]
m.out <- matchit(treat ~  `qlf$table$logCPM`, polysomal_total_matched, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched
polyvsTotalFC_matched= ggplot(data_matched, aes(y= `qlf$table$logFC`, x= Buffering, fill= Buffering))+geom_boxplot()+geom_point(position = position_jitter(width = 0.1), size = 1)+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs(y = "Log2FC  translating RNA vs Total fraction")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                axis.title = element_text(size = 12),
                                                                                                                                                                                                                                                plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")


polyvsTotalFC_matched= ggplot(data_matched, aes(y= `qlf$table$logFC`, x= Buffering, fill= Buffering))+geom_violin(width=.6, alpha=.5)+stat_summary(fun = median, geom = "crossbar", 
                                                                                                                                                   width = 0.2, fill = "black", color = "black")+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs(y = "Log2FC  translating RNA vs Total fraction")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                              axis.title = element_text(size = 12),
                                                                                                                                                                                                                                                                                                              plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")
polyvsTotalFC_matched
significance_data <- data_matched %>%
  summarise(
    TB1_TB2 = wilcox.test(`qlf$table$logFC` ~ Buffering, data_matched  = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(`qlf$table$logFC` ~ Buffering, data_matched = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(`qlf$table$logFC` ~ Buffering, data_matched= .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )



# half life of buffered genes across sublocations
Human_Data<- read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/Human_Data.csv")


#-------------------------------------------------------------------------------------------#----------------------------------------------------------------------------------------------------
#Renaming and Formatting

#So I can left join it with the same column name

Human_Data <- Human_Data |>
  dplyr::rename(
     "Chromatin"= Chromatin.Average,
    "Nucleoplasm"= Nucleoplasm.Average ,
    "Cytoplasm"= Cytoplasm.Average,
    "Untranslated Cytoplasm"=Polysome.Average,
    "Whole Cell"= Whole.Cell.Average
  )

Human_Data_matchit <-  Human_Data %>%
  mutate(Buffering = if_else(Symbol %in% Top_buffered_human$transcript, "TB1",
                             if_else(Symbol %in% Top2_buffered_human$transcript, "TB2", "TB3")))


#Expands the data so it can be graphed
data_long <- pivot_longer(Human_Data_matchit, 
                          cols = c("Chromatin", "Nucleoplasm", "Cytoplasm", "Untranslated Cytoplasm"), 
                          names_to = "Variable", 
                          values_to = "Value") |>
  dplyr::select(Gene, Symbol, PUND, Variable, Value,Buffering)# Convert Value column to numeric


data_long$Value <- as.numeric(data_long$Value)

# Identify the rows with NA values in the Value column
na_rows <- which(is.na(data_long$Value))
# View the rows with NA values
na_data <- data_long[na_rows, ]
# Remove rows where Value column has NA
data_long_clean <- data_long[!is.na(data_long$Value), ]


data_long_clean$Buffering <- factor(data_long_clean$Buffering, levels = c("TB1","TB2","TB3"))

half_life_subcellular= ggplot(data_long_clean, aes(Variable, Value, fill=factor(Buffering)))+geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( y = "Average minutes")+theme(axis.text = element_text(size = 12),
                                                                                                                                                                                                              
                                                                                                                                                                                                              axis.title = element_text(size = 12),
                                                                                                                                                                                                              plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,300))+scale_x_discrete(limits = c("Chromatin", "Nucleoplasm", "Cytoplasm", "Untranslated Cytoplasm"))
half_life_subcellular


# calculate the median 
data_long_clean_TB3= data_long_clean[data_long_clean$Buffering=="TB3",]
median(subset(data_long_clean_TB3, data_long_clean_TB3$Variable=="Cytoplasm")$Value)


# average half life of mRNA. 
Human_Data_matchit_Wholecell <- Human_Data_matchit |>
  dplyr::select(Symbol, Buffering, `Whole Cell`)
Human_Data_matchit_Wholecell$`Whole Cell` <- as.numeric(Human_Data_matchit_Wholecell$`Whole Cell`)
# Identify the rows with NA values in the Value column
na_rows <- which(is.na(Human_Data_matchit$Cell))
# View the rows with NA values
na_data <- Human_Data_matchit[na_rows, ] 

Human_Data_matchit_Wholecell$Buffering <- factor(Human_Data_matchit_Wholecell$Buffering, levels = c("TB1","TB2","TB3"))

# the plot represents data from 12311 genes and others were not calculated. 

half_life_av=ggplot(Human_Data_matchit_Wholecell, aes(x = factor(Buffering), y = `Whole Cell`, fill = factor(Buffering))) +
  geom_boxplot(notch = TRUE)+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+labs( y = "Average minutes")+theme(axis.text = element_text(size = 12),
                                                                                                                   
                                                                                                                   axis.title = element_text(size = 12),
                                                                                                                   plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")



half_life_av
half_life_av=ggplot(Human_Data_matchit_Wholecell, aes(x =`Whole Cell`, color = factor(Buffering))) +
  stat_ecdf(size=1)+scale_color_manual(values = c("TB1"= "#045275","TB2"= "#089099", "TB3"="#B7E6A5"))+labs( y = "Cumulative probability", x="Half-Life (min)")+theme(axis.text = element_text(size = 12),
                                                                                                                               
                                                                                                                               axis.title = element_text(size = 12),
                                                                                                                               plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+
  coord_cartesian(xlim=c(0,2000))



half_life_av

# calculate themedian 
Human_Data_matchit_Wholecell_TB3= Human_Data_matchit_Wholecell[Human_Data_matchit_Wholecell$Buffering=="TB3",]

median(Human_Data_matchit_Wholecell_TB3$`Whole Cell`, na.rm="TRUE")
Human_Data_matchit_Wholecell_TB2= Human_Data_matchit_Wholecell[Human_Data_matchit_Wholecell$Buffering=="TB2",]

median(Human_Data_matchit_Wholecell_TB2$`Whole Cell`, na.rm="TRUE")

Human_Data_matchit_Wholecell_TB1= Human_Data_matchit_Wholecell[Human_Data_matchit_Wholecell$Buffering=="TB1",]

median(Human_Data_matchit_Wholecell_TB1$`Whole Cell`, na.rm="TRUE")

Human_Data_matchit_Wholecell %>%
  summarise(
    TB1_TB2 = wilcox.test(`Whole Cell` ~ Buffering, Human_Data_matchit_Wholecell  = .,
                          subset = Buffering %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(`Whole Cell` ~ Buffering, Human_Data_matchit_Wholecell = .,
                          subset = Buffering %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(`Whole Cell` ~ Buffering, Human_Data_matchit_Wholecell= .,
                          subset = Buffering %in% c("TB2", "TB3"))$p.value
  )


# if we control for the total half life.

#-----------------------------------------------------------------------------------------------------
# Matches the data based on Whole Cell Average
Human_Data_matchit= as.data.table(Human_Data_matchit)
# for matching convert to numeric
cols_to_convert <- c("Chromatin", "Nucleoplasm", "Cytoplasm", "Untranslated Cytoplasm","Whole Cell","NotPund.NucExp","PUND.NucExp")
Human_Data_matchit[, (cols_to_convert) := lapply(.SD, as.numeric), .SDcols = cols_to_convert]
Human_Data_matchit

Human_Data_matchit[, treat := ifelse(Buffering == "TB1", 1,
                                     ifelse(Buffering == "TB2", 1, 0))]



#remove NA# this is removing many genes???
is.na(Human_Data_matchit$`Whole Cell`)
Human_Data_matchit <- Human_Data_matchit[!is.na(Human_Data_matchit$`Whole Cell`)]
sum(is.na(Human_Data_matchit))

#set seed
set.seed(123)

m.out <- matchit(treat ~ `Whole Cell`, Human_Data_matchit, method = "nearest")
summary(m.out)
data_matched= match.data(m.out)
data_matched



data_matched_long <- pivot_longer(data_matched, cols = c("Chromatin", "Nucleoplasm", "Cytoplasm", "Untranslated Cytoplasm"), names_to = "Variable", values_to = "Value")                   
data_matched_long                                                         
Subcellular_matched_plot=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#045275", "#089099", "#B7E6A5"))+
  labs( y = "Time in minutes")+theme(axis.text = element_text(size = 12),
                                     axis.title = element_text(size = 12),
                                     plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,350))+scale_x_discrete(limits = c("Chromatin", "Nucleoplasm", "Cytoplasm", "Untranslated Cytoplasm"))
Subcellular_matched_plot
# calculate the median 
data_matched_long_TB1= data_matched_long[data_matched_long$Buffering=="TB1",]
median(subset(data_matched_long_TB1, data_matched_long_TB1$Variable=="Nucleoplasm")$Value)


#parallel plot

library(GGally)
ggparcoord(data_matched, 
           columns = 3:7,    # The 5 gene-related variables
           groupColumn = "Buffering",  # Column to color by (7th column: "Group")
           scale = "globalminmax",    # Standardizes values for better comparison
           alphaLines = 0.6  # Transparency of lines
) +
  theme_minimal() +
  labs(title = "Parallel Coordinate Plot of Gene Data",
       x = "Gene-related Variables",
       y = "Standardized Values") +
  scale_color_manual(values = c("TB1" = "blue", "TB2" = "red", "TB3" = "green"))+ylim(0,300)+
  facet_wrap(~ Buffering)


#MOUSE
RNA_clr_mouse=read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/RNA_clr_mouse.csv")
Ribo_clr_mouse=read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/Ribo_clr_mouse.csv")
TE_GBM_clr_mouse=read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/TE_GBM_clr_mouse.csv")
dim(RNA_clr_mouse)
dim(Ribo_clr_mouse)
dim(TE_GBM_clr_mouse)


RNA_clr_mouse[1:5,1:5]
RNA_clr_mouse=RNA_clr_mouse[,-1]
Ribo_clr_mouse=Ribo_clr_mouse[,-1]
TE_GBM_clr_mouse= TE_GBM_clr_mouse[,-1]

RNA_clr_mouse[1:5,1:5]
Ribo_clr_mouse[1:5,1:5]
TE_GBM_clr_mouse[1:5,1:5]
RNA_clr_mouse=as.data.table(RNA_clr_mouse)
Ribo_clr_mouse=as.data.table(Ribo_clr_mouse)
TE_GBM_clr_mouse=as.data.table(TE_GBM_clr_mouse)


# it has 8135 gene and 550 samples 
#each row is a gene and each column is a sample
RNA_mouse_nond= RNA_clr_mouse
mouse_TE_nond=  TE_GBM_clr_mouse
Ribo_mouse_nond= Ribo_clr_mouse
transcript=RNA_clr_mouse$transcript

numeric_RNA_nond_mouse<-RNA_mouse_nond[, apply(.SD,1, as.numeric),.SDcols =-"transcript"]

numeric_TE_nond_mouse <-mouse_TE_nond[, apply(.SD,1, as.numeric),.SDcols =-"transcript"]
# in the below code each column is corelated to each other column such that the corleation between each geneis stored as diagonal
TE_RNA_cor_nond_mouse= cor(numeric_TE_nond_mouse, numeric_RNA_nond_mouse, method = "spearman")
TE_RNA_cor_value_nond_mouse= diag(TE_RNA_cor_nond_mouse)
TE_RNA_sample_corl_nond_mouse= data.table(transcript, TE_RNA_cor_value_nond_mouse)
median_TE_RNA_sample_corl_nond_mouse = median(TE_RNA_sample_corl_nond_mouse$TE_RNA_cor_value_nond_mouse)


TE_RNA_nond_mouse=ggplot(TE_RNA_sample_corl_nond_mouse, aes(x = TE_RNA_cor_value_nond_mouse)) +
  geom_histogram(binwidth = 0.02, fill = "#FCE1A4", alpha=0.8)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines")) +
  geom_vline(xintercept = median_TE_RNA_sample_corl_nond_mouse, linetype = "dashed", color = c("blue"), linewidth = 1)+
  scale_x_continuous(name = "Spearman Correlation  coeffcient", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(name = "Freequency", breaks = scales::pretty_breaks(n = 10),expand = c(0, 0))+coord_cartesian(ylim = c(0, 400), xlim=c(-0.8
                                                                                                                                                               , 0.6))
TE_RNA_nond_mouse

numeric_Ribo_nond_mouse <-Ribo_mouse_nond[, apply(.SD,1, as.numeric),.SDcols =-"transcript"]
Ribo_RNA_cor_value_nond_mouse= cor(numeric_Ribo_nond_mouse, numeric_RNA_nond_mouse, method = "spearman")
Ribo_RNA_cor_nond_mouse= diag(Ribo_RNA_cor_value_nond_mouse)
length(Ribo_RNA_cor_nond_mouse)
Ribo_RNA_cor_nond_mouse=as.data.table(Ribo_RNA_cor_nond_mouse)
Ribo_RNA_cor_nond_mouse$transcript= transcript
median_Ribo_RNA_cor_nond_mouse = median(Ribo_RNA_cor_nond_mouse$Ribo_RNA_cor_nond_mouse)

Ribo_RNA_nond_plot_mouse=ggplot(Ribo_RNA_cor_nond_mouse, aes(x = Ribo_RNA_cor_nond_mouse)) +
  geom_histogram(binwidth = 0.02, fill = "#FCE1A4", alpha=0.8)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines")) +
  geom_vline(xintercept = median_Ribo_RNA_cor_nond_mouse, linetype = "dashed", color = c("blue"), linewidth = 1)+
  scale_x_continuous(name = "Spearman Correlation  coeffcient", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(name = "Frequency", breaks = scales::pretty_breaks(n = 10),expand = c(0, 0))+coord_cartesian(ylim = c(0, 450), xlim=c(-0.1
                                                                                                                                                               , 1))

Ribo_RNA_nond_plot_mouse


# Repeat this for mouse. 
#mouse MAD calculation. 
#add if buffered or not
MAD_RNA1_mouse= apply(numeric_RNA_nond_mouse, 2, mad)
MAD_Ribo1_mouse= apply(numeric_Ribo_nond_mouse, 2, mad)
length(MAD_RNA1_mouse)
length(MAD_Ribo1_mouse)
MAD_RNA_Ribo1_mouse= cbind(MAD_RNA1_mouse, MAD_Ribo1_mouse, transcript)
MAD_RNA_Ribo1_mouse=as.data.table(MAD_RNA_Ribo1_mouse)
colnames(MAD_RNA_Ribo1_mouse)= c("MAD_RNA_mouse","MAD_Ribo_mouse","transcript")

#remove difference in the case
MAD_RNA_Ribo1_mouse= MAD_RNA_Ribo1_mouse[,transcript := toupper(transcript)]
MAD_RNA_Ribo1_mouse
#ratui
MAD_RNA_Ribo1_mouse_ratio= mutate(MAD_RNA_Ribo1_mouse, MAD_ratio_mouse= as.numeric(MAD_Ribo_mouse)/as.numeric(MAD_RNA_mouse))
TE_RNA_sample_corl_nond_mouse= TE_RNA_sample_corl_nond_mouse[,transcript := toupper(transcript)]
Buffering_mouse=merge(TE_RNA_sample_corl_nond_mouse, MAD_RNA_Ribo1_mouse_ratio, by.x= "transcript", by.y= "transcript")


Buffering_mouse=Buffering_mouse%>%
  mutate_at(vars(-all_of(c("transcript"))), as.numeric)

dim(Buffering_mouse)

head(Buffering_mouse)
Ribo_RNA_cor_nond_mouse= Ribo_RNA_cor_nond_mouse[,transcript := toupper(transcript)]
# add RNA_Ribo corr
Buffering_mouse=merge(Buffering_mouse, Ribo_RNA_cor_nond_mouse,  by.x= "transcript", by.y= "transcript")
tail(Buffering_mouse)

# Give weight to each criteria. Multiple the rank with the weight and sum it up . Rank according to score again (lowest to highest)
#
# rank according to the parameter $# Multiply each parameter by 2,1, 0. 5 and  0.5 then take the sum
Buffering_mouse_rank= Buffering_mouse %>% mutate(
  TE_RNA_rank_mouse = rank(TE_RNA_cor_value_nond_mouse, ties.method = "first"),
  MAD_rank_mouse = rank(MAD_ratio_mouse, ties.method = "first")) %>%
  mutate(
    TE_RNA_rank_adj_mouse = TE_RNA_rank_mouse * 1,
    MAD_rank_adj_mouse = MAD_rank_mouse * 1) %>%
  rowwise() %>%
  mutate(Row_Sum = sum(c(TE_RNA_rank_adj_mouse,  MAD_rank_adj_mouse))) %>%
  ungroup() %>%
  mutate(TB_rank_mouse=rank(Row_Sum, ties.method="first")) %>%
  arrange(TB_rank_mouse)


as.list(Buffering_mouse_rank[1:500,1])


# Top 250 as the TB score 1, 251-5-- TB score 2 and rest as TB score3 

Top2_buffered_mouse= subset(Buffering_mouse_rank, Buffering_mouse_rank$TB_rank_mouse > 250 & Buffering_mouse_rank$TB_rank_mouse <= 500) %>% arrange(TB_rank_mouse)
Top_buffered_mouse= subset(Buffering_mouse_rank, Buffering_mouse_rank$TB_rank_mouse > 0 & Buffering_mouse_rank$TB_rank_mouse <= 250) %>% arrange(TB_rank_mouse)

# Give TB score
Buffering_mouse_rank <- Buffering_mouse_rank %>%
  mutate(Buffering = if_else(transcript %in% Top_buffered_mouse$transcript, "TB_score=1",
                             if_else(transcript %in% Top2_buffered_mouse$transcript, "TB_score=2", "TB_score=3")))



#mouse
MAD_ratio_plot_mouse= ggplot(Buffering_mouse, aes(x= MAD_ratio_mouse))+geom_histogram(binwidth = 0.05, fill = "#FCE1A4",alpha= 0.8)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines")) +
  scale_x_continuous(name = "MAD Ribo/MAD RNA ratio", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(name = "Frequency", breaks = scales::pretty_breaks(n = 10),expand = c(0, 0))+coord_cartesian(ylim = c(0, 700), xlim=c(0
                                                                                                                                                              , 3))


MAD_ratio_plot_mouse


# Make the MAD plot with Buffering score
MAD_table_plot_mouse =ggplot(Buffering_mouse_rank, aes(x = MAD_RNA_mouse, y = MAD_Ribo_mouse, shape = Buffering, color = as.factor(Buffering))) +
  # Outer layer for border effect with a larger size
  geom_point(color = NA, stroke = 1.5) +
  # Inner layer for filled points
  geom_point(size = 1, aes(fill = Buffering, alpha = ifelse(Buffering== "TB_score=1", 0.4, ifelse(Buffering == "TB_score=2", 0.4, 0.1)))) +
  geom_abline(intercept = 0, slope = 1, color = "#003147", linetype= "dashed")+
  scale_shape_manual(values = c(24, 23, 21)) +
  scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4")) +  # Set fill colors
  scale_color_manual(values = c("#861350", "#d62929", "#fbd783"))+# Use shapes with borders
  labs(shape = "Condition", fill = "Fill Color", color = "Border Color")+guides(alpha = "none", color="none", fill="none")+theme(axis.text = element_text(size = 12),
                                                                                                                                 axis.title = element_text(size = 12),
                                                                                                                                 plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
                                                                                                                                 panel.spacing = unit(0.5, "lines"))+
  scale_x_continuous(name = "Median Absolute Deviation (mRNA abundance)") +scale_y_continuous(name = "Median Absolute Deviation (Ribosome Occupancy)")+theme(
    axis.ticks = element_line(color = "black"),  # Customize tick marks
    axis.ticks.length = unit(0.2, "cm"),  # Adjust length of tick marks
    axis.text.x = element_text(hjust = 0.5),  # Rotate x-axis text
    axis.text.y = element_text(hjust = 0.5)  # Customize y-axis text size
  )+labs(shape= "Buffering_score")+theme(legend.position ="none")+coord_cartesian(ylim = c(0, 3), xlim=c(0
                                                                                                            , 3))
MAD_table_plot_mouse

#example mouse # Fus as it is the top candidate
# for mouse
RNA_Ribo_plot_mouse=  function(gene) {
  gene_RNA_nond <- RNA_mouse_nond[transcript ==gene] 
  gene_RNA_nond[,1:5]
  gene_TE_nond<- mouse_TE_nond[transcript ==gene]
  gene_TE_nond[,1:5]
  gene_Ribo_nond = Ribo_mouse_nond[transcript==gene]
  gene_cor_nond_Ribo= cor(as.numeric(gene_Ribo_nond[,2:550]), as.numeric(gene_RNA_nond[,2:550]), method = "spearman")
  gene_cor_nond_Ribo
  gene_RNA_Ribo_nond = data.table(as.numeric(gene_Ribo_nond[,2:550]), as.numeric(gene_RNA_nond[,2:550]))
  colnames(gene_RNA_Ribo_nond) = c("gene_Ribo_nond", "gene_RNA_nond")
  gene_RNA_TE_nond_plot= ggplot(gene_RNA_Ribo_nond, aes(x=gene_RNA_nond, y=gene_Ribo_nond))+geom_point(alpha=0.5,size= 1, shape=24, color="#861350", fill="#AB1866")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))+
    scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Ribosome occupancy", breaks = scales::pretty_breaks(n = 5))+
    theme(legend.position = "none")
  gene_RNA_TE_nond_plot
}
Fus_RNA_Ribo_plot= RNA_Ribo_plot_mouse(gene="Fus")
Fus_RNA_Ribo_plot
eIF1_RNA_Ribo_plot= RNA_Ribo_plot_mouse(gene="Eif1")
eIF1_RNA_Ribo_plot

RNA_TE_plot_mouse=  function(gene) {
  gene_RNA_nond <- RNA_mouse_nond[transcript ==gene] 
  gene_RNA_nond[,1:5]
  gene_TE_nond<- mouse_TE_nond[transcript ==gene]
  gene_TE_nond[,1:5]
  gene_Ribo_nond = Ribo_mouse_nond[transcript==gene]
  gene_cor_nond_TE= cor(as.numeric(gene_TE_nond[,2:550]), as.numeric(gene_RNA_nond[,2:550]), method = "spearman")
  gene_cor_nond_TE
  gene_RNA_TE_nond = data.table(as.numeric(gene_TE_nond[,2:550]), as.numeric(gene_RNA_nond[,2:550]))
  colnames(gene_RNA_TE_nond) = c("gene_TE_nond", "gene_RNA_nond")
  gene_RNA_TE_nond_plot= ggplot(gene_RNA_TE_nond, aes(x=gene_RNA_nond, y=gene_TE_nond))+geom_point(alpha=0.5,size= 1, shape=24, color="#861350", fill="#AB1866")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,face="plain",size= 8),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))+
    scale_x_continuous(name = "mRNA abundance", breaks = scales::pretty_breaks(n = 5)) +scale_y_continuous(name = "Translation efficiency", breaks = scales::pretty_breaks(n = 5))+
    theme(legend.position = "none")
  gene_RNA_TE_nond_plot
}
Fus_RNA_TE_plot=RNA_TE_plot_mouse(gene="Fus")
Fus_RNA_TE_plot


eIF1_RNA_TE_plot=RNA_TE_plot_mouse(gene="Eif1")
eIF1_RNA_TE_plot

Mouse_human=merge(Buffering_human_rank, Buffering_mouse_rank, by.x= "transcript", by.y = "transcript") 

dim(Mouse_human)
Mouse_human$overlap <- with(Mouse_human, ifelse(Buffering.x == "TB_score=1" & Buffering.y == "TB_score=1", 
                                                "Common",
                                                ifelse(Buffering.x == "TB_score=1", 
                                                       "Human_Top", 
                                                       ifelse(Buffering.y == "TB_score=1", 
                                                              "Mouse_Top", 
                                                              "Others"))))


Human_mouse_overlap_plot= ggplot(Mouse_human, aes(x= TE_RNA_cor_value_nond, y=TE_RNA_cor_value_nond_mouse,color=overlap,fill=overlap))+geom_point(size=1, alpha=0.5, shape=21)+labs(
  x = "mRNA~TE (Human)",
  y = "mRNA~TE (Mouse)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none", legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"))+
  scale_color_manual(values= c("Others"= "#def5d6" , "Common"="#FC600A", "Human_Top"= "#045275", "Mouse_Top"= "#AB1866"))+
  scale_fill_manual(values= c("Others"= "#def5d6" , "Common"="#FC600A", "Human_Top"= "#045275", "Mouse_Top"= "#AB1866"))
Human_mouse_overlap_plot
Mouse_human_cor=cor.test(Mouse_human$TE_RNA_cor_value_nond,Mouse_human$TE_RNA_cor_value_nond_mouse, method="spearman")

#howmany are human with TB1 # 229
dim(Mouse_human[Mouse_human$Buffering.x=="TB_score=1",])
#howmany are mouse with TB1 # 225
dim(Mouse_human[Mouse_human$Buffering.y=="TB_score=1",])
dim(Mouse_human)#7113


# fisher test (how many of genes in mouse and human are common)
dim(Mouse_human)
fisher.test(matrix(c(126,103,99, 6785),nrow = 2, byrow = TRUE))

#2B.
# gene ontology of the overlap
Buffering_score =merge(Buffering_human_rank, Buffering_mouse_rank, by.x= "transcript", by.y = "transcript") 
Buffering_score= cbind(Buffering_score$transcript, Buffering_score$Buffering.x, Buffering_score$Buffering.y)
colnames(Buffering_score)= c("transcript", "Human Buffering score", "Mouse Buffering score")
Buffering_score=as.data.table(Buffering_score)





common_buffering126= Buffering_score[Buffering_score$`Human Buffering score`=="TB_score=1" & Buffering_score$`Mouse Buffering score`=="TB_score=1"]

#need to manually assign function. 
#write.csv(common_buffering126, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/common_buffering126.csv")
#assign genefunction



#common_buffering53_function= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/common_buffering53_function.csv", header = FALSE)
#common_buffering53_function
#merge the columns 
#merged_126= merge(common_buffering53_function, common_buffering126, by.x= "V1", by.y= "transcript", all =TRUE )
#write.csv(merged_126, "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/merged_126.csv")
merged_126= read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/merged_126.csv")
dim(merged_126)
# Group them and take a percentage of each term and make a pie chart. 
function_percent <- merged_126 %>%
  group_by(V2) %>%
  summarise(Count = n()) %>%
  mutate(Percent = Count / sum(Count) * 100)

# Create a tree map
library(viridis)
treemap_table=treemap(function_percent,
                      index = "V2",
                      vSize = "Percent",
                      type = "index",
                      fontsize.title = 14,
                      fontsize.labels = 12,
                      fontcolor.labels = "white",
                      bg.labels = "#444444",
                      align.labels = c("center", "center"),
                      border.col = "#FFFFFF",
                      border.lwds = 3,
                      palette = viridis(10, option = "virdis")
)

treemap_table
library(treemapify)

treemap_plot= ggplot(function_percent, aes(area = Percent, fill = V2, label = V2 )) +
  geom_treemap() +                             # Draw the treemap
  geom_treemap_text(colour = "white",          # Add text labels
                    place = "center", 
                    grow = TRUE) +
  theme_minimal()+theme(legend.position ="none")

gene_counts <- merged_126 |> 
  group_by(V2) |> 
  summarise(Gene_Count = n()) |> 
  arrange(desc(Gene_Count))
gene_counts

# Create a pie chart
pie_chart=
  ggplot(function_percent, aes(x = "", y = Percent, fill = V2)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Percentage of Genes by Function") +
  theme_void() + # Remove background, gridlines, and axes
  geom_text(aes(label = paste0(round(Percent, 1), "%")), 
            position = position_stack(vjust = 0.5))


pie_chart




#Figure 2mouse


set.seed(123)
Mouse_stats_path = "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/Buffering MS/Mouse_statsR (3).csv"
Mouse_stats <- read.csv(Mouse_stats_path)

# Capitalize the 'transcript' column
Mouse_stats <- Mouse_stats %>%
  mutate(transcript = toupper(transcript))

# Change 0's to 1's in specific columns (e.g., Column1, Column2, Column3) so log(1) = 0
Mouse_stats <- Mouse_stats %>%
  mutate(UTR3_GC = ifelse(UTR3_GC == 0, 1, UTR3_GC),
         UTR3_len = ifelse(UTR3_len == 0, 1, UTR3_len),
         UTR5_GC = ifelse(UTR5_GC == 0, 1, UTR5_GC),
         UTR5_len = ifelse(UTR5_len == 0, 1, UTR5_len),
         CDS_GC = ifelse(CDS_GC == 0, 1, CDS_GC),
         CDS_len = ifelse(CDS_len == 0, 1, CDS_len)
  )

Mouse_stats <- Mouse_stats %>%
  mutate(UTR3_GC = UTR3_GC_Perc/100,
         UTR3_len = log10(UTR3_len),
         UTR5_GC = UTR5_GC_Perc/100,
         UTR5_len = log10(UTR5_len),
         CDS_GC = CDS_GC_Perc/100,
         CDS_len = log10(CDS_len)
  )




# CDS length mouse
Mouse_stats_Match_Prep <- Mouse_stats |>
  mutate(Buffering = if_else(transcript %in% Top_buffered_mouse$transcript, "TB1",
                             if_else(transcript %in% Top2_buffered_mouse$transcript, "TB2", "TB3"))) |>
  dplyr::select(Buffering, UTR5_len, UTR3_len, CDS_len, transcript)


# Binary for Matchit TB1 & TB2 = 1 | TB3 = 0
Mouse_stats_Match_Prep <- as.data.table(Mouse_stats_Match_Prep)
Mouse_stats_Match_Prep[, treat := ifelse(Buffering %in% c("TB1", "TB2"), 1, 0)]

#remove NA
Mouse_stats_Match_Prep <- na.omit(Mouse_stats_Match_Prep)

# Matches based on UTR3_GC_Perc + UTR5_GC_Perc
m.out <- matchit(treat ~ UTR3_len + UTR5_len, Mouse_stats_Match_Prep, method = "nearest")
#summary(m.out)

data_matched <- match.data(m.out)

data_matched_long <- pivot_longer(data_matched, cols = c(CDS_len), names_to = "Variable", values_to = "Value") |>
  dplyr::select(Value, Variable, Buffering, transcript)

#-----------------------------------------------------------

# Analysis

#-----------------------------------------------------------
# Define the levels for comparison
data_matched_long$Buffering <- factor(data_matched_long$Buffering, levels = c("TB1","TB2","TB3"))
data_matched_long_CDS_length_mouse= data_matched_long

CDS_length_plot_mouse=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+
  labs( y = "log 10 CDS length")+theme(axis.text = element_text(size = 8),
                                       axis.title = element_text(size = 8),
                                       plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0,4.5))


CDS_length_plot_mouse

my_comparison = list(c("TB1", "TB2"), c("TB1", "TB3"), c("TB2", "TB3"))

significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB2", "TB3"))$p.value
  )

significance_data



# 5UTR Len with (CDS Len Matchit)

# Matches based on UTR3_GC_Perc + UTR5_GC_Perc
m.out <- matchit(treat ~ CDS_len, Mouse_stats_Match_Prep, method = "nearest")
#summary(m.out)

data_matched <- match.data(m.out)

data_matched_long <- pivot_longer(data_matched, cols = c(UTR5_len), names_to = "Variable", values_to = "Value") |>
  dplyr::select(Value, Variable, Buffering, transcript)



# Define the levels for comparison
data_matched_long$Buffering <- factor(data_matched_long$Buffering, levels = c("TB1","TB2","TB3"))
data_matched_long_5UTR_length_mouse= data_matched_long

UTR5_length_plot_mouse=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+
  labs(y = "log 10 5UTR length")+theme(axis.text = element_text(size = 8),
                                       axis.title = element_text(size = 8),
                                       plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_x_discrete(expand = c(0.2, 0.2))+coord_cartesian(ylim = c(0,4.5))
UTR5_length_plot_mouse
my_comparison = list(c("TB1", "TB2"), c("TB1", "TB3"), c("TB2", "TB3"))

significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB2", "TB3"))$p.value
  )

significance_data



# UTR 3 

m.out <- matchit(treat ~ CDS_len, Mouse_stats_Match_Prep, method = "nearest")
#summary(m.out)

data_matched <- match.data(m.out)

data_matched_long <- pivot_longer(data_matched, cols = c(UTR3_len), names_to = "Variable", values_to = "Value") |>
  dplyr::select(Value, Variable, Buffering, transcript)



# Define the levels for comparison
data_matched_long$Buffering <- factor(data_matched_long$Buffering, levels = c("TB1","TB2","TB3"))
data_matched_long_3UTR_length_mouse= data_matched_long
UTR3_length_plot_mouse=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+
  labs( y = "log 10 3UTR length")+theme(axis.text = element_text(size = 8),
                                        axis.title = element_text(size = 8),
                                        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim = c(0,4.5))+scale_x_discrete(expand = c(0.2, 0.2))
UTR3_length_plot_mouse
my_comparison = list(c("TB1", "TB2"), c("TB1", "TB3"), c("TB2", "TB3"))

significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB2", "TB3"))$p.value
  )

significance_data

# mouse combined figure
length_new_mouse= rbind(data_matched_long_CDS_length_mouse,data_matched_long_5UTR_length_mouse, data_matched_long_3UTR_length_mouse)
length_new_mouse$Variable <- factor(length_new_mouse$Variable, levels = c("CDS_len", "UTR5_len", "UTR3_len"))

Combined_length_mouse= ggplot(length_new_mouse, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+
  labs( y = "length (log10)")+theme(axis.text = element_text(size = 12),
                                    axis.title = element_text(size = 12),
                                    plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim = c(0.5,4.3))+scale_x_discrete(expand = c(0.2, 0.2))



Combined_length_mouse



# supplementary 2
# 2B CDS GC percentage
# Matches based on UTR3_GC_Perc + UTR5_GC_Perc
Mouse_stats_Match_Prep <- Mouse_stats |>
  mutate(Buffering = if_else(transcript %in% Top_buffered_mouse$transcript, "TB1",
                             if_else(transcript %in% Top2_buffered_mouse$transcript, "TB2", "TB3"))) |>
  dplyr::select(Buffering, UTR3_GC, UTR5_GC, CDS_GC, transcript)


# Binary for Matchit TB1 & TB2 = 1 | TB3 = 0
Mouse_stats_Match_Prep <- as.data.table(Mouse_stats_Match_Prep)
Mouse_stats_Match_Prep[, treat := ifelse(Buffering %in% c("TB1", "TB2"), 1, 0)]

#remove NA
Mouse_stats_Match_Prep <- na.omit(Mouse_stats_Match_Prep)

# Matches based on UTR3_GC_Perc + UTR5_GC_Perc
m.out <- matchit(treat ~ UTR3_GC + UTR5_GC, Mouse_stats_Match_Prep, method = "nearest")
#summary(m.out)

data_matched <- match.data(m.out)

data_matched_long <- pivot_longer(data_matched, cols = c(CDS_GC), names_to = "Variable", values_to = "Value") |>
  dplyr::select(Value, Variable, Buffering, transcript)

#-----------------------------------------------------------

# Analysis

#-----------------------------------------------------------
# Define the levels for comparison

data_matched_long$Buffering <- factor(data_matched_long$Buffering, levels = c("TB1","TB2","TB3"))
data_matched_long_CDS_GC_mouse= data_matched_long
CDS_GC_plot_mouse=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+
  labs( y = "CDS GC content")+theme(axis.text = element_text(size = 8),
                                    axis.title = element_text(size =8),
                                    plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim = c(0,1))+scale_x_discrete(expand = c(0.2, 0.2))



CDS_GC_plot_mouse
significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB2", "TB3"))$p.value
  )

significance_data


# 2B mouse CDS GC
# 2F mouse 3 UTR GC
m.out <- matchit(treat ~ CDS_GC, Mouse_stats_Match_Prep, method = "nearest")
#summary(m.out)

data_matched <- match.data(m.out)

data_matched_long <- pivot_longer(data_matched, cols = c(UTR3_GC), names_to = "Variable", values_to = "Value") |>
  dplyr::select(Value, Variable, Buffering, transcript)

#-----------------------------------------------------------

# Analysis

#-----------------------------------------------------------
# Define the levels for comparison
data_matched_long$Buffering <- factor(data_matched_long$Buffering, levels = c("TB1","TB2","TB3"))
data_matched_long_UTR3_GC_mouse= data_matched_long
UTR3_GC_plot_mouse=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+
  labs(y = "UTR3 GC content")+theme(axis.text = element_text(size = 8),
                                    axis.title = element_text(size = 8),
                                    plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim = c(0,1))+scale_x_discrete(expand = c(0.2, 0.2))


UTR3_GC_plot_mouse
significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB2", "TB3"))$p.value
  )

significance_data

#2D mouse 5 UTR GC

m.out <- matchit(treat ~ CDS_GC, Mouse_stats_Match_Prep, method = "nearest")
#summary(m.out)

data_matched <- match.data(m.out)

data_matched_long <- pivot_longer(data_matched, cols = c(UTR5_GC), names_to = "Variable", values_to = "Value") |>
  dplyr::select(Value, Variable, Buffering, transcript)

#-----------------------------------------------------------

# Analysis

#-----------------------------------------------------------
# Define the levels for comparison
data_matched_long$Buffering <- factor(data_matched_long$Buffering, levels = c("TB1","TB2","TB3"))
data_matched_long_5UTR_GC_mouse= data_matched_long
UTR5_GC_plot_mouse=ggplot(data_matched_long, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+
  labs(y = "UTR5 GC content")+theme(axis.text = element_text(size = 8),
                                    axis.title = element_text(size = 8),
                                    plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim = c(0,1))+scale_x_discrete(expand = c(0.2, 0.2))


UTR5_GC_plot_mouse
significance_data <- data_matched_long %>%
  group_by(Variable) %>%
  summarise(
    TB1_TB2 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(Value ~ Buffering, data = ., subset = Buffering %in% c("TB2", "TB3"))$p.value
  )

significance_data




#combine
# mouse combined figure
length_new_mouse= rbind(data_matched_long_CDS_GC_mouse,data_matched_long_5UTR_GC_mouse, data_matched_long_UTR3_GC_mouse)
length_new_mouse$Variable <- factor(length_new_mouse$Variable, levels = c("CDS_GC", "UTR5_GC", "UTR3_GC"))


Combined_GC_mouse= ggplot(length_new_mouse, aes(x = Variable, y = Value, fill= Buffering)) +
  geom_boxplot()+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+
  labs( y = "GC content")+theme(axis.text = element_text(size = 12),
                                    axis.title = element_text(size = 12),
                                    plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim = c(0.2,0.9))+scale_x_discrete(expand = c(0.2, 0.2))



Combined_GC_mouse




# for mouse
numeric_RNA_mouse=RNA_mouse_nond%>%
  mutate_at(vars(-all_of(c("transcript"))), as.numeric)
# Add the medians as a new column to the data table
medians_mouse <- numeric_RNA_mouse[, apply(.SD,1, median),.SDcols =-"transcript"]
# 
RNA_median_mouse<- numeric_RNA_mouse %>%
  mutate(median = medians_mouse)
#change the case
RNA_median_mouse= RNA_median_mouse[,transcript := toupper(transcript)]
# also add if buffering or non buffering
RNA_median_mouse <-RNA_median_mouse%>%
  mutate(Buffering = if_else(transcript %in% Top_buffered_mouse$transcript, "TB_score=1",
                             if_else(transcript %in% Top2_buffered_mouse$transcript, "TB_score=2", "TB_score=3")))



# Effect of expression levels on buffering.
#with all genes in x axis and y xis mRNA median levels
RNA_median_mouse= cbind(RNA_median_mouse$transcript, RNA_median_mouse$median, RNA_median_mouse$Buffering)
colnames(RNA_median_mouse)= c("transcript","median", "Buffering")
RNA_median_mouse= as.data.table(RNA_median_mouse)
RNA_median_mouse=RNA_median_mouse%>%
  mutate_at(vars(-all_of(c("transcript", "Buffering"))), as.numeric)
#plot the median 
Median_RNA_expression_mouse= ggplot(RNA_median_mouse, aes(x= Buffering, y= median, fill= Buffering))+geom_boxplot()+theme(axis.text = element_text(size = 12),
                                                                                                                          axis.title = element_text(size = 12),
                                                                                                                          plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                          panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+labs(x = "Buffering",y = "Median RNA expression")+ylim(-1,5)


Median_RNA_expression_mouse
Median_RNA_expression_mouse= ggplot(RNA_median_mouse, aes(x= median, color= Buffering))+geom_density(size=1)+theme(axis.text = element_text(size = 12),
                                                                                                                   axis.title = element_text(size = 12),
                                                                                                                   plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                   panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_color_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+labs(x = "Median RNA expression",y = "Frequency")+ylim(0,0.5)+xlim(-5.5,5)


Median_RNA_expression_mouse

#calculate the significance


significance_data <- RNA_median_mouse %>%
  summarise(
    TB1_TB2 = wilcox.test(median ~ Buffering, RNA_median_mouse = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=2"))$p.value,
    TB1_TB3 = wilcox.test(median ~ Buffering, RNA_median_mouse = .,
                          subset = Buffering %in% c("TB_score=1", "TB_score=3"))$p.value,
    TB2_TB3 = wilcox.test(median ~ Buffering, RNA_median_mouse = .,
                          subset = Buffering %in% c("TB_score=2", "TB_score=3"))$p.value
  )
significance_data



# For mouse
Buffering_mouse <-Buffering_mouse %>%
  mutate(Buffering = if_else(transcript %in% Top_buffered_mouse$transcript, "TB_score=1",
                             if_else(transcript %in% Top2_buffered_mouse$transcript, "TB_score=2", "TB_score=3")))
RNA_MAD_mouse= ggplot(Buffering_mouse, aes(x= Buffering, y= MAD_RNA_mouse, fill= Buffering))+geom_violin()+stat_summary(fun = median, geom = "crossbar", 
                                                                                                                         width = 0.2, fill = "black", color = "black")+theme(axis.text = element_text(size = 12),
                                                                                                                 axis.title = element_text(size = 12),
                                                                                                                 plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA),
                                                                                                                 panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+scale_fill_manual(values = c( "#AB1866", "#E05C5C", "#FCE1A4"))+labs(
                                                                                                                   x = "Buffering",y = "MAD RNA expression")+coord_cartesian(ylim=c(0,2.6))
RNA_MAD_mouse






Figure1a= plot_grid(TE_RNA_nond, TE_RNA_nond_mouse+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y= element_blank() ), MAD_ratio_plot_human, MAD_ratio_plot_mouse+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y= element_blank() ),MAD_table_plot_human, MAD_table_plot_mouse+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y= element_blank() ),
                    ncol = 2,
                    rel_widths = c(1.2,1,1.2,1,1.2,1),  # Adjust widths
                    rel_heights = c(1, 1,1,1),
                    labels = c('A','B', 'C','D','E','F'))

Figure1a
Figure1b= plot_grid(G6PD_RNA_Ribo_nond_plot,G6PD_RNA_TE_nond_plot,FUS_RNA_Ribo_nond_plot, FUS_RNA_TE_nond_plot,Fus_RNA_Ribo_plot, Fus_RNA_TE_plot,
                    ncol = 2,
                    rel_widths = c(1, 1,1,1),  # Adjust widths
                    rel_heights = c(1, 1,1,1),
                    labels = c('G','H','I','J','K','L'))


Figure1b

Figure1 <- plot_grid(Figure1a, Figure1b, axis = 'l', align= 'v')
Figure1

Sup_Figure1= plot_grid(Ribo_RNA_nond_plot, Ribo_RNA_nond_plot_mouse, JUNB_RNA_Ribo_nond_plot, JUNB_RNA_TE_nond_plot,PPIA_RNA_Ribo_nond_plot,  PPIA_RNA_TE_nond_plot, eIF1_RNA_Ribo_plot, eIF1_RNA_TE_plot,
                       ncol = 2,
                       rel_widths = c(1, 1,1,1),  # Adjust widths
                       rel_heights = c(1, 1,1,1),
                       labels = c('A','B', 'C','F','D','G','E','H')# Adjust heights
)

Sup_Figure1


top_row_2 <- plot_grid(Human_mouse_overlap_plot, treemap_plot, Combined_length, Combined_length_mouse+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y= element_blank()), labels = c('A', 'B','C', 'D'), label_size =12, rel_widths = c(1, 1), ncol = 2, # Adjust widths
                       rel_heights = c(1, 1))

#bottom_row_2 <- plot_grid(CDS_length_plot,UTR5_length_plot, UTR3_length_plot,  CDS_length_plot_mouse  ,  UTR5_length_plot_mouse,  UTR3_length_plot_mouse, labels = c('B', 'D','F','C','E','G'), label_size =12, rel_widths = c(1, 0.95),  # Adjust widths
#rel_heights = c(1, 1))
bottom_row_2 = plot_grid(Median_RNA_expression_human, Median_RNA_expression_mouse, labels = c('E', 'F'), label_size =12, rel_widths = c(1, 1),  # Adjust widths
                         rel_heights = c(1, 1), ncol=1)

bottom_row_3= plot_grid( bottom_row_2, CFD_TBS1, label_size = 12, ncol = 2 , rel_widths=c(0.4,1), rel_heights = c(1,1))

Figure2=plot_grid(top_row_2,  bottom_row_3, nrow=2, rel_heights =  c(1.5, 1))
Figure2







Sup_Figure2 <- plot_grid(
  plot_grid(Combined_length_GC, Combined_GC_mouse+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y= element_blank()), RNA_MAD_human, RNA_MAD_mouse+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y= element_blank()), ncol=2, labels = c('A', 'B','C','D'), label_size =12, rel_widths = c(1, 1),  # Adjust widths
            rel_heights = c(1, 1)) 
  # Adjust heights
)



Sup_Figure2

Sup_Figure3 <- plot_grid(
  cai_buf_ran_long_plot, CFD_TBS2, labels = c('A', 'B'), label_size = 12, ncol = 2 ,rel_widths = c(0.4, 1),  # Adjust widths
  rel_heights = c(1, 1)
  # Adjust heights
)
Sup_Figure3




Figure3 <- plot_grid(NULL, final_plot, mad_values_Cell_2020_plot, mad_mouse,
  ncol = 2,
  labels = c('A','B', 'C','E'), label_size =12 
  # Adjust heights
)



Sup_Figure4 = plot_grid(final_plot2, mad_values_Cell_2020_plot_matched+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y= element_blank()) ,mad_mouse_matched_proteome +theme(axis.title.y = element_blank()),
                        ncol = 3,
                        rel_widths = c(1.15, 0.9,1),  # Adjust widths
                    
                        labels = c('A','B', 'C'), label_size =12 
                        # Adjust heights
)

Sup_Figure4


Figure4 <- plot_grid(
  pLI_plot, pTriplo_plot, siRNA_RNA_Ribo,Cenik_2015,
  ncol = 2,
  rel_widths = c(1, 1,1,1),  # Adjust widths
  rel_heights = c(1, 1,1,1),
  labels = c('B', 'C','D','E')
  # Adjust heights
)
Figure4













Sup_Figure5= plot_grid(pHaplo_plot,
                       ncol = 2,
                       rel_widths = c(1, 1,1,1),  # Adjust widths
                       rel_heights = c(1, 1,1,1) 
                       # Adjust heights
)
Sup_Figure5






Figure5_top= plot_grid(
  Translating_fraction, Ribo_MAD,polyvsTotalFC,
  ncol = 3,
  labels = c('B', 'C','D')
  # Adjust heights
)
Figure5_bottom= plot_grid(
  half_life_av, half_life_subcellular,
  ncol = 2,
  rel_widths = c(0.5,1),
  labels = c( 'E','F')
  # Adjust heights
)
Figure5= plot_grid(Figure5_top, Figure5_bottom,
                   nrow=2)
Figure5
Sup_Figure6= plot_grid(
  Translating_fraction_matched, Ribo_MAD_matched, polyvsTotalFC_matched,Subcellular_matched_plot, 
  ncol = 2
  # Adjust height
)
Sup_Figure6



path = "C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Imputation methods/"
ggsave("Figure1.pdf", plot = Figure1, path = path, width = 10, height = 8,dpi=300)

ggsave("Sup_Figure1.pdf", plot =Sup_Figure1, path = path, width = 10, height = 10,dpi=300)
ggsave("Figure2.pdf", plot = Figure2, path = path, width = 10, height = 10,dpi=300)
ggsave("Sup_Figure2.pdf", plot =Sup_Figure2, path = path, width = 10, height = 10,dpi=300)
ggsave("Figure3.pdf", plot = Figure3, path = path, width = 10, height = 10,dpi=300)
ggsave("Sup_Figure3.pdf", plot =Sup_Figure3, path = path, width = 15, height = 6,dpi=300)
ggsave("Sup_Figure4.pdf", plot =Sup_Figure4, path = path, width = 10, height = 5,dpi=300)
ggsave("Sup_Figure5.pdf", plot =Sup_Figure5, path = path, width = 10, height = 10,dpi=300)
ggsave("Figure4.pdf", plot = Figure4, path = path, width = 14, height = 10,dpi=300)
ggsave("Figure5.pdf", plot = Figure5, path = path, width = 8, height = 5,dpi=300)
ggsave("Sup_Figure6.pdf", plot = Sup_Figure6, path = path, width = 10, height = 10,dpi=300)











# Mouse MAD
# Imports the protein abundance data

Protein_abundance_table <- read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Aden/Mad Mouse Score ~ New/Sequencing_data/Protein_abundance/Protein_abundance.csv")

# Imports the longest appris transcripts data
Appris_data <-read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Aden/Mad Mouse Score ~ New/Longest_transcript/longest_appris_transcripts.csv")




# Imports the data needed to convert protein ID's to transcript and gene codes
Protein_to_transcript <- read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Aden/Mad Mouse Score ~ New/G_profiler/gProfiler_mmusculus_protein_transcript.csv")

#Convert Protein ID's to Transcripts

# Formats the G_profiler File
Protein_to_transcript <- Protein_to_transcript |>
  rename(initial_alias= "PROTEIN_ID",
        converted_alias= "Transcript_code")

Protein_to_transcript=  Protein_to_transcript[, c("PROTEIN_ID","Transcript_code", "name")]


# Left joins the two tables
Protein_abundance_transcript <- left_join(Protein_to_transcript, Protein_abundance_table, by = "PROTEIN_ID")


#Select only Principal Longest APPRIS Isoforms
# Retrieve only Principal:1 Isoforms

# Formats the Appris data
Appris_data <- Appris_data |>
  rename(Transcript_Code = "Transcript_code")

Protein_abundance_table <- left_join(Appris_data, Protein_abundance_transcript, by = "Transcript_code")

# Remove any NA's during the transfer process and only takes the Gene names from the Appris data
Protein_abundance_table <- filter(Protein_abundance_table, !is.na(PROTEIN_ID), PROTEIN_ID != "N/A") |>
  dplyr::select(!c("Length", "GENE_NAME")) |>
  rename(name="GENE_NAME")

# Checks for duplicate genes
Duplicate_genes <- Protein_abundance_table |>
  group_by(GENE_NAME) |>
  summarize(count = n()) |>
  filter(count >= 2) |>
  dplyr::select(GENE_NAME)

# Remove gene names that exist in Duplicate_genes. i.e one protein ID from two transcripts 
Protein_abundance_table <- Protein_abundance_table |>
  anti_join(Duplicate_genes, by = "GENE_NAME")

#-----------------------------------------------------------------
# Remove Excess Data
remove(Duplicate_genes)
remove(Appris_data)
remove(Protein_abundance_transcript)


# Calculates the NA count per sample across each gene
Na_count <- apply(Protein_abundance_table, 1, function(x) sum(is.na(x)))

# Concatenates this NA_count to each of the genes
Protein_abundance_table <- cbind(Protein_abundance_table, Na_count)

# Groups the NA Counts to visualize the Distribution of Invalid Samples
Protein_abundance_table_Na <- Protein_abundance_table |>   
  group_by(Na_count) |>   
  summarize(count = n()) 

# Checks for duplicate genes
# Protein_abundance_table |>
#   group_by(GENE_NAME) |>
#   filter(n() != 1) |>
#   arrange(GENE_NAME)

#-------------------------------------------------------------------
# Visulize the distributions of Invalid Samples Per Gene

ggplot(Protein_abundance_table_Na, aes(x = factor(Na_count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  
  labs(
    title = "Gene Count by NA Count",  # Updated title to match the data
    x = "NA Count",
    y = "Gene Count"
  ) +
  
  scale_x_discrete(breaks = 0:41) +
  # Fixes the bars to the x-axis while giving space above
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  
  
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    legend.position = "top"  # Position the legend at the top
  )

#-----------------------------------------------------------------

# Filters for genes that have atleast 27% valid samples
Protein_abundance_table <- filter(Protein_abundance_table, Na_count <= 30)

# Checks for duplicate genes
# Protein_abundance_table |>
#    group_by(GENE_NAME) |>
#    filter(n() != 1) |>
#    arrange(GENE_NAME)


#-----------------------------------------------------------------
# Remove Excess Data
remove(Na_count)
remove(Protein_abundance_table_Na)

#Perform Median Absolute Deviation Calculation
# Slices the Protein ID and Gene Name
Protein_transcript_id <- Protein_abundance_table[,1:4]

# Remove naming columns so MAD can be performed
Protein_abundance_table <- Protein_abundance_table |>
  dplyr::select(!c(Transcript_code, Gene_Code, PROTEIN_ID, GENE_NAME, Na_count))

# Calculates the NA count per sample across each gene
Protein_NA <- apply(Protein_abundance_table, 1, function(x) sum(is.na(x)))

# Caclulates the Median RNA Abundance
Median_protein <- apply(Protein_abundance_table, 1 , median, na.rm = TRUE)

# Caclulates the MAD of the Protein
Mad_protein <- apply(Protein_abundance_table, constant = 1.4826, 1, mad, na.rm = TRUE)

# Restore the Identification
Mad_protein_table <- data.table(Protein_transcript_id, Median_protein, Mad_protein, Protein_NA) |>
  mutate(GENE_NAME = toupper(GENE_NAME))

#-------------------------------------------------------------------
# Remove Excess Data
remove(Protein_NA)
remove(Median_protein)
remove(Mad_protein)
remove(Protein_abundance_table)
remove(Protein_transcript_id)

#Import Mouse Buffered Gene List

Buffered_list <- Buffering_mouse_rank


# Gives the order 1 for the first 250 rows then 2 for the next 250 rows
Buffered_list <- Buffered_list |>
  rename(
     transcript="GENE_NAME"
  ) |>
  mutate(Order = case_when(row_number() <= 250 ~ 1,
                        row_number() <= 500 ~ 2,
                        TRUE~ 3)) |>
  dplyr::select(GENE_NAME, Order)


#Map MAD Genes to Buffering Score
# Joins the table and if Order is empty, give the rest of the genes the score 3
Mad_buffered_genes <- merge(Mad_protein_table, Buffered_list, by = "GENE_NAME")
  
dim(Mad_buffered_genes)
#6501 genes are common between th eprotein expression and our Ribobase list
#Format Data
# Renaming the Buffering Score

Mad_buffered_genes <- Mad_buffered_genes |>
  mutate(Order = case_when(
    Order == 1 ~ "TB1",
    Order == 2 ~ "TB2",
    Order == 3 ~ "TB3"
  ))


#Visualization without Match-it
# Define the levels for comparison
Mad_buffered_genes$Order <- factor(Mad_buffered_genes$Order, levels = c("TB1","TB2","TB3"))

mad_mouse=ggplot(data = Mad_buffered_genes, aes(x = factor(Order), y = Mad_protein, fill = factor(Order))) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+labs(x = "Condition", y = "Median absolute Deviation")+theme(axis.text = element_text(size = 12),
                                                                                                                               axis.title = element_text(size = 12),
                                                                                                                               plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA), panel.spacing = unit(0.5, "lines"))+theme(legend.position ="none")+coord_cartesian(ylim=c(0,0.6))+scale_x_discrete(expand = c(0.2, 0.2))



mad_mouse

significance_data <- Mad_buffered_genes %>%
  summarise(
    TB1_TB2 = wilcox.test(Median_protein ~ Order, Mad_buffered_genes = .,
                          subset = Order %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(Median_protein ~ Order, Mad_buffered_genes = .,
                          subset = Order %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(Median_protein~ Order, Mad_buffered_genes= .,
                          subset = Order %in% c("TB2", "TB3"))$p.value
  )
significance_data




#-------------------------------------------------------------------

#Import RNA Sequencing TPM
Rna_abundance_table <- read.csv("C:/Users/sjr2797/Box/Cenik lab_Shilpa/FUS/FUS_2023/Aden/Mad Mouse Score ~ New/Sequencing_data/RNA_abundance.csv")
 
#Quality Control Data
# Calculates the NA count per sample across each gene
Na_count <- apply(Rna_abundance_table, 1, function(x) sum(is.na(x)))

# Concatenates this NA_count to each of the genes
Rna_abundance_table <- cbind(Rna_abundance_table, Na_count)

# Groups the NA Counts to visualize the Distribution of Invalid Samples
Rna_abundance_table_Na <- Rna_abundance_table |>   
  group_by(Na_count) |>   
  summarize(count = n()) 

#-------------------------------------------------------------------
# Visulize the distributions of Invalid Samples Per Gene

ggplot(Rna_abundance_table_Na, aes(x = factor(Na_count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  
  labs(
    title = "Gene Count by NA Count",  # Updated title to match the data
    x = "NA Count",
    y = "Gene Count"
  ) +
  
  scale_x_discrete(breaks = 0:41) +
  # Fixes the bars to the x-axis while giving space above
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  
  
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    legend.position = "top"  # Position the legend at the top
  )

#-------------------------------------------------------------------
# Filters for genes that have atleast 70% valid samples

Rna_abundance_table <- filter(Rna_abundance_table, Na_count <= 20)



#Perform Median Absolute Deviation Calculations
# Slices the ENSG Code
Rna_transcript_id <- Rna_abundance_table[,1]

# Remove naming columns so MAD can be performed
Rna_abundance_table <- Rna_abundance_table |>
  dplyr::select(!c(Gene.ID))

# Calculates the NA count per sample across each gene
Rna_NA <- apply(Rna_abundance_table, 1, function(x) sum(is.na(x)))

# Caclulates the Median RNA Abundance
Median_rna <- apply(Rna_abundance_table, 1 , median, na.rm = TRUE)

# Caclulates the MAD of the Protein
Mad_rna <- apply(Rna_abundance_table, constant = 1.4826, 1, mad, na.rm = TRUE)

# Restore the Identification
Mad_rna_table <- data.table(Rna_transcript_id, Median_rna, Mad_rna, Rna_NA)



#Join the RNA and Protein Abundance Tables
# Format the name to match the current Ribo table for matching
Mad_rna_table <- Mad_rna_table |>
  rename(Rna_transcript_id="Gene_Code")

Mad_buffered_genes_joined <- left_join(Mad_buffered_genes, Mad_rna_table, by = "Gene_Code")

# Checks for duplicates
Mad_buffered_genes_joined |>
  dplyr::count(GENE_NAME) |>
  filter(n > 1)

# Remove unsuccessful joining
Mad_buffered_genes_joined <- filter(Mad_buffered_genes_joined, !is.na(Median_rna), Median_rna != "N/A")

#--
# Define the levels for comparison
Mad_buffered_genes_joined$Order <- factor(Mad_buffered_genes_joined$Order, levels = c("TB1","TB2","TB3"))

#match it with MAD
Mad_buffered_genes_matched <- Mad_buffered_genes_joined |>
  mutate(match_it_order = case_when(
    Order == "TB1" ~ 1,
    Order == "TB2" ~ 1,                 # Assign 0 for Order value 3
    Order == "TB3" ~ 0,                
  ))

m.out <- matchit(match_it_order ~ Mad_rna, Mad_buffered_genes_matched, method = "nearest")
Mad_buffered_genes_matched = match.data(m.out)
#Visualize Match-It MAD Protein (MAD RNA)

## Define the levels for comparison
Mad_buffered_genes_matched$Order <- factor(Mad_buffered_genes_matched$Order, levels = c("TB1","TB2","TB3"))

mad_mouse_matched_proteome= ggplot(data = Mad_buffered_genes_matched, aes(x = Order, y = Mad_protein, fill = Order)) +
  geom_boxplot(outlier.shape = NA)+
scale_fill_manual(values = c("#AB1866", "#E05C5C", "#FCE1A4"))+labs(x = "Condition", y = "Median absolute Deviation")+theme(axis.text = element_text(size = 12),
                                                                                                                              axis.title = element_text(size = 12),
                                                                                                                              plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 1, fill = NA), panel.spacing = unit(0.1, "lines"))+
  theme(legend.position ="none")+scale_x_discrete(expand = c(0.2,0.2))+ scale_y_continuous(expand = c(0.01, 0))+coord_cartesian(ylim =c(0,0.8))



mad_mouse_matched_proteome

significance_data <- Mad_buffered_genes_matched %>%
  summarise(
    TB1_TB2 = wilcox.test(Median_protein ~ Order, Mad_buffered_genes_matched = .,
                          subset = Order %in% c("TB1", "TB2"))$p.value,
    TB1_TB3 = wilcox.test(Median_protein ~ Order, Mad_buffered_genes_matched = .,
                          subset = Order %in% c("TB1", "TB3"))$p.value,
    TB2_TB3 = wilcox.test(Median_protein~ Order, Mad_buffered_genes_matched= .,
                          subset = Order %in% c("TB2", "TB3"))$p.value
  )
significance_data


