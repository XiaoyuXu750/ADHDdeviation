# This script generates the barplot of S-A ranks for 17 networks, S-A connectional axis ranks, 
# and save out the RGB for the networks.

library(ggplot2)
library(reshape)
library(pals)
library(scales)
library(grDevices)
library(tidyverse)
library(openxlsx)
rm(list=ls())


wd <- getwd()
homepath <- str_split_i(wd, "Normative_model", 1)
demopath <- paste0(homepath, '/Normative_model/demography')
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_ABCD')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_ABCD")
FigureFolder <- paste0(homepath, "/Normative_model/Figures_ABCD")

Yeo17_SAvertex <- read.csv(paste0(resultFolder, '/Yeo17S400_SArank_vertex.csv'))
Yeo17_SAtab <- read.csv(paste0(resultFolder, '/Yeo17S400_SArank.csv'))
Du15_SAvertex <- read.csv(paste0(resultFolder, '/Du15_SArank_vertex.csv'))
Du15_SAtab <- read.csv(paste0(resultFolder, '/Du15_SArank.csv'))
Du15_RGB <- read.xlsx("D:/xuxiaoyu/Atlas/DU15NET/DU15NET/RGBindex.xlsx")
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))

Yeo17Schaefer400_SA <- schaefer400_index_SA %>% group_by(network_label_17network) %>%
  summarise(MedianSA = median(finalrank.wholebrain))
source(paste0(functionFolder, "/plotmatrix.R"))

# 1. S-A connectional rank, 15*15
#############################
axeslabels=c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.C", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B") 
Matrix.tmp <- matrix(NA, nrow = 15, ncol=15)
for (x in 1:15){
  for (y in 1:15){
    Matrix.tmp[x,y] <- (x+y)^2+(x-y)^2
  }
}
SCrank <- rank(Matrix.tmp[lower.tri(Matrix.tmp, diag = T)], ties.method = "average")
SCrank.df <- data.frame(SCrank=SCrank)
PaletteSet <- list(Name="RdBu", drirection=-1, lmmin = 1, lmmax = max(SCrank), anglex=0, angley=0,
                   hjustx = 1, hjusty = 1, vjustx = 1, vjusty=0.5)
axeslabels=c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.C", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B") 

FigSCrank <- plotmatrix(dataname="SCrank.df", "SCrank", ds.resolution=15, 
                          Pvar=NA, NAcol="white", lmthr=NA, axeslabels, axeslabelsGap=F, 
                          linerange_frame=NA, PaletteSet=PaletteSet)

filename<-paste0(FigureFolder, "/Connectionrank_squareDistance_15_Yeo17.tiff")
ggsave(filename, height = 14, width = 16, units = "cm")

#####################################

# Yeo network mean SA axis rank
######################################
Yeo17_SAtab$index <- c(0:17)
limbicindex <- Yeo17_SAtab$index[str_detect(Yeo17_SAtab$NetworkLabel, "Limbic")]
Yeo17_SAtab <- Yeo17_SAtab %>% filter(! index %in% limbicindex)
Yeo17_SAtab$RGB <- c("white", "#771183", "#FF0200", "#4680B3", "#2ACCA5", "#499C3D", "#007512", "#C334F4", 
                            "#F997D8", "#778CB3", "#E89321", "#8A314B", "#0D28FD", "#000083", "#F9FF04", "#D13F52")
                            
Yeo17_SAtab <- Yeo17_SAtab[order(Yeo17_SAtab$MedianSA), ]

Yeo17_SAvertex <- Yeo17_SAvertex %>% filter(! Yeo17S400_Group_Label_NetLabel %in% limbicindex)
Yeo17_SAvertex$Yeo17_Group_Label <- factor(Yeo17_SAvertex$Yeo17S400_Group_Label_NetLabel, levels = Yeo17_SAtab$index, 
                                           labels=Yeo17_SAtab$NetworkLabel)

Yeo17_SAvertex <- Yeo17_SAvertex %>% filter(Yeo17_Group_Label != "None") %>% droplevels()
Yeo17_SAtab <- Yeo17_SAtab %>% filter(NetworkLabel != "None") %>% droplevels()
axeslabels=c("VS.P", "SM.A", "VS.C", "SM.B", "DA.A", "DA.B", "FP.C", "DM.C", "VA.A", "TP", "FP.A", "VA.B", "DM.A", "FP.B", "DM.B") 
ggplot(data=Yeo17_SAvertex)+
  geom_boxplot(aes(x=Yeo17_Group_Label, y=SAaxis_Label, fill=Yeo17_Group_Label))+
  #scale_fill_manual(values = c("#8E3896", "#679DD0", "#298230", "#D35CFF", "#F6AE4C", "#DA6073"))+
  scale_fill_manual(values = Yeo17_SAtab$RGB)+
  scale_x_discrete(labels= axeslabels)+
  scale_y_continuous(breaks=c(5, 40000), labels=c("S pole", "A pole"), position = "left")+
  labs(x=NULL, y=NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=21, color="black"), 
        axis.text.x = element_text(angle = 90, hjust=1, vjust = 1),
        axis.text.y = element_text(angle = 90, hjust=0),
        axis.ticks.y = element_blank(),
        axis.title =element_text(size=21),aspect.ratio = 0.55,
        plot.title = element_text(size=21, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
filename<-paste0(FigureFolder, "/Yeo17S400network_SAaxis_vertex.tiff")
ggsave(filename, width=11.5, height =10, units = "cm")

## Yeo17 bar
Yeo17_SAtab$SArank <- as.numeric(c(1:15))
Yeo17_SAtab$SArank_lab <- factor(Yeo17_SAtab$SArank, levels=c(1:15))
Yeo17_SAtab$long <- 1
cols <- rev(COL2(diverging = "RdBu", n=17))
cols <- cols[2:16]
#[1] "#950F26" "#BB2933" "#D15648" "#E58267" "#F5AB8A" "#FACDB5" "#FDE8DB" "#FFFFFF" "#E2EEF5" "#C1DDEB" "#99C9E0"
#[13] "#6AACD0" "#3E8DC0" "#2971B1" "#16518F"
fig <- ggplot(data=Yeo17_SAtab)+
  geom_bar(aes(x=long, y=long, fill = SArank_lab, group=SArank_lab), color="black", 
           alpha=1, stat = "identity", position = "stack", width=0.1)+
  scale_fill_manual(values=cols)+
  scale_y_continuous(breaks=seq(from=0.5, to=14.5, by=1), labels = rev(axeslabels), expand = c(0, 0))+
  scale_x_continuous(breaks=NULL, expand = c(0, 0), position = "bottom")+
  labs(x=NULL, y=NULL, color=NULL, fill=NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=16, color='black'),
        axis.text.y = element_text(angle = 0, vjust=0.5, hjust = 1),
        axis.title = element_text(size = 16),aspect.ratio =14.5,
        axis.line = element_blank(), 
        plot.title=element_text(size=16, color='black', hjust=0.5),
        legend.position = "none")
fig
  
ggsave(paste0(FigureFolder, '/Yeo17_colorLabel.tiff'),fig, width=2.5, height =14, units = "cm")



hex_to_rgb <- function(hex_color) {
  if (!is.character(hex_color) || length(hex_color) != 1) {
    stop("输入必须是一个字符向量")
  }
  
  # 移除可能存在的#号
  hex_color_clean <- gsub("^#", "", hex_color)
  
  # 检查十六进制格式是否有效（3位或6位）
  if (!grepl("^[0-9A-Fa-f]{3}$|^[0-9A-Fa-f]{6}$", hex_color_clean)) {
    stop("无效的十六进制颜色格式，请使用3位或6位格式（如FFF或FF0000）")
  }
  
  # 使用R内置函数转换
  rgb_matrix <- col2rgb(hex_color)
  
  # 转换为命名向量，更直观
  rgb_vector <- c(r = rgb_matrix[1, 1], 
                  g = rgb_matrix[2, 1], 
                  b = rgb_matrix[3, 1])
  
  return(rgb_vector)
}

for (i in 1:15){
  index <- Yeo17_SAtab$SArank[i]
  rbg <- hex_to_rgb(cols[index])
  rbg <- rbg / 256
  
  Yeo17_SAtab$r[i] <- rbg[1]
  Yeo17_SAtab$b[i] <- rbg[3]
  Yeo17_SAtab$g[i] <- rbg[2]
}



write.csv(Yeo17_SAtab, paste0(resultFolder, "/Yeo17_RdBu15.csv"))
######################################


