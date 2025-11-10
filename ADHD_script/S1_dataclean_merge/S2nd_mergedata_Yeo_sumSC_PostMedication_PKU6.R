## This script is to generate a dataframe for the scans after medication treatment, 
## in which each column is the strength for an edge in large-scale network.

library(R.matlab)
library(ggplot2)
library(tidyverse)
library(parallel)
library(openxlsx)
library(rjson)
library(reshape)
library(psych)
rm(list = ls())
wdpath <- getwd()
homepath <- str_split_i(wdpath, "Normative_model", 1)
# input Yeo resolution of 7 or 17.
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
elementnum <- Yeoresolution.delLM*(Yeoresolution.delLM+1) /2

SC_path_PKU6 <-'D:/xuxiaoyu/open_dataset_information/PKU6/SCmatpost'
Volume_path_PKU6 <-'D:/xuxiaoyu/open_dataset_information/PKU6/schaefer400_nodevolumepost'
qc_path <- 'D:/xuxiaoyu/open_dataset_information/PKU6/qc_jsonpost'

demopath <- paste0(homepath, '/Normative_model/demography')
interfileFolder <- paste0(homepath, '/Normative_model/interfileFolder_EFNYnoCCNP')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder <- paste0(homepath, "/Normative_model/results_EFNYnoCCNP")
FigureFolder <- paste0(homepath, "/Normative_model/Figures_EFNYnoCCNP")

Behavior <- read.csv(paste0(demopath, '/ALL_SITES_basic_demo.csv'))
Behavior <- Behavior %>% distinct(ID, .keep_all = TRUE)
Behavior <- Behavior %>% filter(age >= 6.5, age <= 15.5)

#### import schaefer400 index
schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.csv'))
## qsiprep output matrix is in Yeo 7 order, so reorder schaefer400 index to Yeo 7 order
schaefer400_index<-schaefer400_index_SA[order(schaefer400_index_SA$index),]
limbicindex <- which(str_detect(schaefer400_index$label_17network, "Limbic"))
schaefer400_index <- schaefer400_index[-limbicindex, ]
schaefer376_delLM <- schaefer400_index$index
## Rearrange left and right regions.
schaefer400_index$index_7network_LRmixed <- schaefer400_index$index
schaefer400_index$index_7network_LRmixed[str_detect(schaefer400_index$label, "RH_")] <- schaefer400_index$index_7network_LRmixed[str_detect(schaefer400_index$label, "RH_")] - 200
orderYeo_7<-order(schaefer400_index$index_7network_LRmixed)

schaefer400_index$index_17network_LRmixed <- schaefer400_index$index_17network
schaefer400_index$index_17network_LRmixed[str_detect(schaefer400_index$label_17network, "RH_")] <- schaefer400_index$index_17network_LRmixed[str_detect(schaefer400_index$label_17network, "RH_")] - 200
orderYeo_17<-order(schaefer400_index$index_17network_LRmixed)

# filter index of P75th of CV.
deleteindex75 <- readRDS(paste0(interfileFolder, '/CV75_deleteindex.Yeo', Yeoresolution,'.delLM.rds'))

# assign each region to Yeo.resolution*Yeo.resolution network.
schaefer400_index.Yeo7 <- schaefer400_index[order(schaefer400_index$index_7network_LRmixed),]
schaefer400_index.Yeo17 <- schaefer400_index[order(schaefer400_index$index_17network_LRmixed),]
schaefer400_index.Yeo7 <- schaefer400_index.Yeo7 %>% mutate(Yeo.resolutionnode = recode_factor(network_label,
                                                                                               "Vis" = 1,
                                                                                               "SomMot" = 2,
                                                                                               "DorsAttn" = 3,
                                                                                               "SalVentAttn" = 4,
                                                                                               "Cont" = 5,
                                                                                               "Default" = 6,
                                                                                               "Limbic" = 0))

summary(schaefer400_index.Yeo7$Yeo.resolutionnode)
schaefer400_index.Yeo17 <- schaefer400_index.Yeo17 %>% 
  mutate(Yeo.resolutionnode = recode_factor(network_label_17network,
                                            "VisCent" = 3,
                                            "VisPeri" = 1,
                                            "SomMotA" =2,
                                            "SomMotB" = 4,
                                            "DorsAttnA" =5,
                                            "DorsAttnB" = 6,
                                            "SalVentAttnA" =9,
                                            "SalVentAttnB" =12,
                                            "ContA" =11,
                                            "ContB" =14,
                                            "ContC" =7,
                                            "DefaultA" = 13,
                                            "DefaultB" = 15,
                                            "DefaultC" = 8,
                                            "TempPar" = 10))

summary(schaefer400_index.Yeo17$Yeo.resolutionnode)

if (Yeoresolution == 7){
  Yeo.resolutionnode <- schaefer400_index.Yeo7$Yeo.resolutionnode
  Yeo.resolutionnode <- factor(Yeo.resolutionnode, levels = c(1:Yeoresolution.delLM))
}else if (Yeoresolution == 17){
  Yeo.resolutionnode <- schaefer400_index.Yeo17$Yeo.resolutionnode
  Yeo.resolutionnode <- factor(Yeo.resolutionnode, levels = c(1:Yeoresolution.delLM))
}else{
  print("Invalid Yeoresolution!")
}

# SC 376*376 --> Yeoresolution.delLM*Yeoresolution.delLM
matrixYeo.resolution <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
matrixYeo.resolution[lower.tri(matrixYeo.resolution, diag = T)] <- c(1:elementnum)
matrixYeo.resolution[upper.tri(matrixYeo.resolution)] <- t(matrixYeo.resolution)[upper.tri(matrixYeo.resolution)]
matrix_SCYeo.resolution <- matrix(NA, 376, 376)
for (x in 1:Yeoresolution.delLM){
  for (y in 1:Yeoresolution.delLM){
    xindex <- which(Yeo.resolutionnode==x)
    yindex <- which(Yeo.resolutionnode==y)
    matrix_SCYeo.resolution[xindex, yindex] <- matrixYeo.resolution[x,y]
  }
}
# an index telling how 376*376 map to 12*12
Yeo.resolution.index <- matrix_SCYeo.resolution[lower.tri(matrix_SCYeo.resolution, diag = T)]
#################################################

#### import SC data
#### Yeoresolution.delLM regions, (Yeoresolution.delLM+1)*Yeoresolution.delLM/2=elementnum SCs
#### extract a dataframe containing elementnum columns, each represents an edge.
#################################################
colname <- character(length = elementnum)
for (i in 1:elementnum){
  colname[i] <- paste0('SC.', as.character(i))
}

Behavior <- Behavior %>% filter(site=="PKU6")
SCdata.sum <- list()
for (i in 1:nrow(Behavior)){
  ID <- Behavior$ID[i]
  site <- Behavior$site[i]
  
  SCname <- paste0(ID, '.mat')
  SC_file_path <- paste0(SC_path_PKU6, '/', SCname)
  volumefile <- paste0(Volume_path_PKU6, '/', ID, '_Volume7.txt')
  
  
  if (file.exists(SC_file_path)){
    qcfile <- fromJSON(file=paste0(qc_path, '/', ID, '.json'))
    mean_fd <- qcfile$subjects[[1]]$mean_fd
    Behavior$mean_fd[i] <- mean_fd
    
    SCmat <- readMat(SC_file_path)
    # load steamline counts matrix & fiber length matrix
    SCmat_raw <- SCmat$schaefer400.sift.radius2.count.connectivity[schaefer376_delLM, schaefer376_delLM]
    length_raw <- SCmat$schaefer400.radius2.meanlength.connectivity[schaefer376_delLM, schaefer376_delLM]
    if (Yeoresolution == 7){
      SCmat_raw <- SCmat_raw[orderYeo_7, orderYeo_7] # 376*376 nodes sorted by Yeo index
      length_raw <- length_raw[orderYeo_7, orderYeo_7] # 376*376 nodes sorted by Yeo index
    }else if (Yeoresolution == 17){
      SCmat_raw <- SCmat_raw[orderYeo_17, orderYeo_17]
      length_raw <- length_raw[orderYeo_17, orderYeo_17]
    }else{
      print("Invalid Yeoresolution!")
    }
    
    totallength_raw <- length_raw * SCmat_raw
    indexup <- upper.tri(SCmat_raw)
    indexsave <- !indexup
    SCmat_raw <- SCmat_raw[indexsave] # 1*70876 each element represents streamline counts
    SCmat_raw75 <- SCmat_raw
    SCmat_raw75[deleteindex75]<-0 # remove top 1/4 inconsistent connetions
    totallength_raw <- totallength_raw[indexsave]
    totallength_raw75 <- totallength_raw
    totallength_raw75[deleteindex75]<-0
    df <- data.frame(
      group = Yeo.resolution.index,
      value75 = SCmat_raw75,
      length75 = totallength_raw75
    )
    # compute the sum of streamline counts / length for each fraction, in total of elementnum.
    result <- df %>%
      group_by(group) %>%
      summarise(sum_value75 = sum(value75), sum_length75=sum(length75))
    mean_length75 <- (result$sum_length75 / result$sum_value75)[1:elementnum]
    sumSC.raw75 <- result$sum_value75[1:elementnum]
    ## node volume
    if (file.exists(volumefile)){
      nodevolume <- read_table(volumefile, col_names=F)
      if (nrow(nodevolume)==453){
        nodevolume <- as.numeric(nodevolume$X2[schaefer376_delLM]) # delete limbic regions
        if (Yeoresolution == 7){
          nodevolume <- nodevolume[orderYeo_7] # sorted by Yeo index
        }else if (Yeoresolution == 17){
          nodevolume <- nodevolume[orderYeo_17]
        }else{
          print("Invalid Yeoresolution!")
        }
        
        df2 <- data.frame(
          group = Yeo.resolutionnode,
          value = nodevolume
        )
        result2 <- df2 %>% arrange(group) %>%
          group_by(group) %>% 
          summarise(sum_value = sum(value))
        nodevolume_sum <- result2$sum_value[1:Yeoresolution.delLM] # sum of nodes' volume for each node fraction (Yeo.resolution).
        
        ### Yeo.resolution*Yeo.resolution
        volumeSC <- matrix(NA, Yeoresolution.delLM, Yeoresolution.delLM)
        for (x in 1:Yeoresolution.delLM){
          for (y in 1:Yeoresolution.delLM){
            volumeSC[x,y] <- (nodevolume_sum[x]+nodevolume_sum[y])/2
          }
        }
        volumeSC <- volumeSC[lower.tri(volumeSC, diag = T)] # the scale values of node volume for each edge.
        sumSC.invnode75 <- sumSC.raw75 / volumeSC
      }else{
        sumSC.invnode75 <- rep(NA, elementnum)
      }}else{
        sumSC.invnode75 <- rep(NA, elementnum)
      }
    ###keep lower triangle and diagonal
    SCdat75 <- as.data.frame(sumSC.invnode75)
    SCdat75 <- as.data.frame(t(SCdat75), row.names = NULL)
    names(SCdat75) <- colname
    row.names(SCdat75) <- NULL
    SCdat75$ID[1] <- ID
    #SCdata.sum75<-rbind(SCdata.sum75, SCdat75)
    
    # not inverse node volume
    SCdat75_noInvNode <- as.data.frame(sumSC.raw75)
    SCdat75_noInvNode <- as.data.frame(t(SCdat75_noInvNode), row.names = NULL)
    names(SCdat75_noInvNode) <- colname
    row.names(SCdat75_noInvNode) <- NULL
    SCdat75_noInvNode$ID[1] <- ID
    
    print(i)
    SCdata.sum[[i]] <- data.frame(SCdat75, SCdat75_noInvNode)
  }
  
}

SCdata.sum2 <- SCdata.sum[!sapply(SCdata.sum, is.null)]

print("The parallel Mat reading is over!")
saveRDS(SCdata.sum2, paste0(interfileFolder, "/SCdata.sum.list_Yeo", Yeoresolution, "_postmedication.rds"))

SCdata.sum75 <- do.call(rbind, lapply(SCdata.sum2, function(x) as.data.frame(x[,c(1:(elementnum+1))])))
SCdata.sum75_noInvNode <- do.call(rbind, lapply(SCdata.sum2, function(x) as.data.frame(x[,c((elementnum+2):(elementnum*2+2))])))
names(SCdata.sum75_noInvNode) <- c(colname, "ID")

SCdata.sum75.merge <- merge(SCdata.sum75, Behavior, by="ID") %>% drop_na(SC.1)
SCdata.sum75_noInvNode.merge <- merge(SCdata.sum75_noInvNode, Behavior, by="ID")

# convert variables' types
SCdata.sum75.merge$ID <- as.factor(SCdata.sum75.merge$ID) ; SCdata.sum75.merge$site <- as.factor(SCdata.sum75.merge$site)
SCdata.sum75_noInvNode.merge$ID <- as.factor(SCdata.sum75_noInvNode.merge$ID) ; SCdata.sum75_noInvNode.merge$site <- as.factor(SCdata.sum75_noInvNode.merge$site)

saveRDS(SCdata.sum75.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.sum.msmtcsd.merge_postmedication.rds'))
saveRDS(SCdata.sum75_noInvNode.merge, paste0(interfileFolder, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSC.sum.msmtcsd.merge_postmedication.rds'))


