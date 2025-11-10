# This script generates age distribution bar plots and demographic tables.

rm(list=ls())
library(ggplot2)
library(tidyverse)
library(tableone)

# set resolution
Yeoresolution <- 17
if (Yeoresolution == 7){
  Yeoresolution.delLM = 6
}else if (Yeoresolution == 17){
  Yeoresolution.delLM = 15
}
element_num <- Yeoresolution.delLM*(Yeoresolution.delLM+1)/2

# input directory
wd <- getwd()
homepath <- str_split_i(wd, "Normative_model", 1)
interfileFolder_EFNY <- paste0(homepath, '/Normative_model/interfileFolder_EFNYnoCCNP')
interfileFolder_ABCD <- paste0(homepath, '/Normative_model/interfileFolder_ABCD')
functionFolder <- paste0(homepath, "/Normative_model/functions")
resultFolder_EFNY <- paste0(homepath, "/Normative_model/results_EFNYnoCCNP")
resultFolder_ABCD <- paste0(homepath, "/Normative_model/results_ABCD")
FigureFolder_EFNY <- paste0(homepath, "/Normative_model/Figures_EFNYnoCCNP")
FigureFolder_ABCD <- paste0(homepath, "/Normative_model/Figures_ABCD")
demopath <- paste0(homepath, '/Normative_model/demography')

deviation.EFNY <- readRDS(paste0(interfileFolder_EFNY, '/SCdata_Yeo',Yeoresolution, '_CV75.deviations_ADHD_All.rds'))
deviation.EFNY <- deviation.EFNY %>% filter(visit=="0")
SCdata.TD.trainset <- readRDS(paste0(interfileFolder_EFNY, "/SCdata.allTD_EFNYnoCCNP.SCYeo", element_num, ".rds"))
data.EFNY <- rbind(deviation.EFNY[,c("age", "mean_fd", "sex", "diagnosis", "site", "ID")], 
                   SCdata.TD.trainset[,c("age", "mean_fd", "sex", "diagnosis", "site", "ID")])


deviation.ABCD <- readRDS(paste0(interfileFolder_ABCD, '/SCdata_Yeo',Yeoresolution, '_CV75_sumSCinvnode.deviations_TDtest_ADHD.rds'))
deviation.ABCD <- deviation.ABCD %>% mutate(
  Setlabel = case_when(if_TD=="TD" ~ "TDtest",
                       .default = "ADHD")
)
SCdata.TD.trainset.ABCD <- readRDS(paste0(interfileFolder_ABCD, "/SCdata.TD.trainset_SCYeo", element_num, ".rds"))
SCdata.TD.trainset.ABCD$if_TD <- "TD"
SCdata.TD.trainset.ABCD$Setlabel <- "TDtrain"
SCdata.ABCD <- rbind(deviation.ABCD[,names(SCdata.TD.trainset.ABCD)], SCdata.TD.trainset.ABCD)

## Figure Age distribution
############################

## Chinese Cohort
data.EFNY$diagnosis <- factor(data.EFNY$diagnosis, levels=c("TD", "ADHD"))
fillcolor = brewer.pal(3, "Paired")

ggplot(data = data.EFNY, aes(age, y = ..count.., fill = site)) +
  geom_histogram(binwidth = 1, color = "black", position = "stack",linewidth=0.5) +
  facet_wrap(~diagnosis, ncol = 2, nrow = 1, scales = "free", 
             labeller = labeller(diagnosis=c("TD"="TD, N=292", "ADHD"="ADHD, N=468"))) +
  labs(x = "Age (years)", y = "Frequency", fill="Study") +
  #scale_x_continuous(limits = c(9, 16), breaks = c(9,10,11,12,13,14,15,16)) +
  #scale_y_continuous(limits = c(0, 50), breaks = c(0,10, 20,30,40, 50)) +
  scale_fill_manual(values = rep(fillcolor,2))+
  #geom_hline(aes(yintercept = 10), colour = "red", linetype="dashed")+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.9,
        plot.title = element_text(color = "black", size = 15, hjust = 0.5),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "transparent",colour = NA),
        axis.title = element_text(color = "black", size = 15),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 15), 
        legend.text = element_text(color = "black", size = 15),
        legend.title = element_text(color = "black", size = 15))

ggsave(paste(FigureFolder_EFNY, '/FigS_Age_distribution_count.tiff', sep = ''),  width = 18.5, height = 12, units = "cm")

## ABCD
SCdata.ABCD$if_TD <- factor(SCdata.ABCD$if_TD, levels=c("TD", "ADHD"))
fillcolor = brewer.pal(3, "Paired")
SCdata.ABCD$eventname2 <- factor(SCdata.ABCD$eventname2, levels=c("baselineYear1Arm1", "2YearFollowUpYArm1", "4YearFollowUpYArm1"))

print(paste(length(which(table(SCdata.ABCD$subID)==3)), "participants have 3 visits,", length(which(table(SCdata.ABCD$subID)==2)), "participants have 2 visits."))
print(paste(length(unique(SCdata.ABCD$subID)), "participants in total."))

fillcolor = brewer.pal(3, "Paired")[c(2, 1, 3)]

ggplot(data = SCdata.ABCD, aes(age, y = ..count.., fill = eventname2)) +
  geom_histogram(binwidth = 0.5, color = "black", position = "stack",linewidth=0.5) +
  facet_wrap(~if_TD, ncol = 2, nrow = 1, scales = "free", 
             labeller = labeller(if_TD=c("TD"="TD, N=6,687", "ADHD"="ADHD, N=1,114"))) +
  labs(x = "Age (years)", y = "Frequency", fill="Visit") +
  #scale_x_continuous(limits = c(9, 16), breaks = c(9,10,11,12,13,14,15,16)) +
  #scale_y_continuous(limits = c(0, 50), breaks = c(0,10, 20,30,40, 50)) +
  scale_fill_manual(values = rep(fillcolor,2), labels=c("Baseline", "2-Year follow-up", "4-Year follow-up"))+
  #geom_hline(aes(yintercept = 10), colour = "red", linetype="dashed")+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.9,
        plot.title = element_text(color = "black", size = 14.8, hjust = 0.5),
        strip.text = element_text(size = 14.8, color = "black"),
        strip.background = element_rect(fill = "transparent",colour = NA),
        axis.title = element_text(color = "black", size = 14.8),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14.8), 
        legend.text = element_text(color = "black", size = 14.8),
        legend.title = element_text(color = "black", size = 14.8))

ggsave(paste(FigureFolder_ABCD, '/Fig1_Age_distribution_count_ABCD_ADHDall.tiff', sep = ''),  width = 21, height = 12, units = "cm")


## Describe demographic information
####################
## ABCD
SCdata.ABCD$sex <- factor(SCdata.ABCD$sex, levels=c(1, 2), labels = c("M", "F"))
Interest.vars <- c("cbcl_scr_syn_external_t", 
                   "cbcl_scr_syn_attention_t", "cbcl_scr_dsm5_adhd_t", "sex", "age", "meanFD", "race_ethnicity", "siteID", "eventname2")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "Setlabel", data = SCdata.ABCD, factorVars=c("sex", "race_ethnicity", "siteID", "eventname2"),test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder_ABCD, "/demoinfo_testtrainsets.csv"), row.names = T)

# ADHD VS TD test
SCdata.ABCD.test <- SCdata.ABCD %>% filter(Setlabel %in% c("ADHD", "TDtest"))
Interest.vars <- c("cbcl_scr_syn_external_t", 
                   "cbcl_scr_syn_attention_t", "cbcl_scr_dsm5_adhd_t", "sex", "age", "meanFD", "race_ethnicity", "siteID", "eventname2")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "Setlabel", data = SCdata.ABCD.test, factorVars=c("sex", "race_ethnicity", "siteID", "eventname2"), test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder_ABCD, "/demoinfo_testsets.csv"), row.names = T)

# TD train VS TD test
SCdata.ABCD.TD <- SCdata.ABCD %>% filter(Setlabel %in% c("TDtrain", "TDtest"))
Interest.vars <- c("cbcl_scr_syn_external_t", 
                   "cbcl_scr_syn_attention_t", "cbcl_scr_dsm5_adhd_t", "sex", "age", "meanFD", "race_ethnicity", "siteID", "eventname2")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "Setlabel", data = SCdata.ABCD.TD, factorVars=c("sex", "race_ethnicity", "siteID", "eventname2"), test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder_ABCD, "/demoinfo_TDsets.csv"), row.names = T)

## Chinese Cohort
Behavior <- read.csv(paste0(demopath, "/basic_demo_PKU6.csv"))
Behavior$ID <- paste0("sub-", Behavior$ID)
data.EFNY <- data.EFNY %>% left_join(dplyr::select(Behavior, ID, medication, IA_B, HI_B, TO_B, IA_F, HI_F, TO_F, ICV), by="ID")
data.EFNY <- data.EFNY %>% mutate(
  IA_B = ifelse(is.na(IA_F), NA, IA_B),
  HI_B = ifelse(is.na(HI_F), NA, HI_B)
)

Interest.vars <- c("sex", "age", "mean_fd", "IA_B", "HI_B", "TO_B", "IA_F", "HI_F", "TO_F", "site", "medication")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "diagnosis", data = data.EFNY,factorVars=c("sex","site"), test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder_EFNY, "/demoinfo_diagnosissets.csv"), row.names = T)

data.EFNY.medication <- data.EFNY %>% drop_na(IA_B)
Interest.vars <- c("sex", "age", "mean_fd", "IA_B", "HI_B", "TO_B", "IA_F", "HI_F", "TO_F", "site")
tableone.df <- CreateTableOne(vars=Interest.vars, strata = "medication", data = data.EFNY.medication,
                              factorVars=c("sex","site"), test = TRUE)
tableone.df <- print(tableone.df)
write.csv(tableone.df, paste0(resultFolder_EFNY, "/demoinfo_diagnosissets_medication.csv"), row.names = T)








