


# --------------------------- Ferrara et al. 2024 - R script N. 3 --------------------------- #
#                                                                                             #
#                                                                                             #
#  Effect of preconditioning treatments on photosynthetic efficiency and                      #
#  tissue color at time points "day-30"                                                       #
#                                                                                             #
#  In this script we calculate the effect of preconditioning on recovery and survival rate    #
#  to estimate the increase in long term resilience after the heat assay                      #
#                                                                                             #
#  The effects are calculated on both metrics (photosynthetic efficiency and tissue color)    #
#  and on survival score                                                                      #
#-------------------------------------------------------------------------------------------- #


# Libraries ----


library(here)
library(lubridate)
library(scales)
library(ggpubr)
library(bestNormalize)
library(rstatix)
library(multcomp)
library(easystats)
library(gridExtra)
library(lme4)
library(ggbeeswarm)
library(ggtext)
library(dplyr)   
library(tidyverse)
library(emmeans)
library(car)
library(lmtest)
library(ggnewscale)
library(survminer)
library(patchwork)

# Load data ----

# Load the photosynthetic efficiency df (which include also the survival data) and the coral tissue color df.
# These df contain the data from the "post-preconditioning", "post-heat", "day-15", and "day-30" time points. 
# DATA ANALYSIS of recovery phase is shown here.
# At each time point, each fragment effective quantum yield was measured in three different spots to obtained a comprehensive overview 
# of coral fragments health status. Before any analysis, the mean of these three measurements was calculated. 

#Statistical analysis of effective quantum yield and tissue color data are performed in a separate script (Ferrara_effect_size_analysis)


#Load the efficiency data (effective quantum yield and surival data)
H_rcv <- read_csv2(here("Data", "Ferrara_PAM_master_df.csv"))


# To perform the tissue color analyses, we calculated the complementary values (255 - x) to show the 0 as white and 255 as bleak.
Clr  <- read_csv2(here("Data", "Ferrara_tissue_color_master_df.csv"))







# RECOVERY------------------------------------------------------------------------------------------

## PAM Data subest ----
Recovery_data <- H_rcv %>%
  
  dplyr::filter(
    Species_ID == "Pru"& Day != 1| 
      Species_ID == "Gfa"& !(Day %in% c("1","2"))| 
      Species_ID == "Mdi"& !(Day %in% c("1","2"))|
      Species_ID == "Pve"|
      Species_ID == "Spi"|
      Species_ID == "Amu") %>%
  
  mutate(Species_short = case_when(str_detect(Species_ID, "Pru") ~ "P. rus",
                                  str_detect(Species_ID, "Gfa") ~ "G. fascicularis",
                                  str_detect(Species_ID, "Mdi") ~ "M. digitata",
                                  str_detect(Species_ID, "Spi") ~ "S. pistillata",
                                  str_detect(Species_ID, "Pve") ~ "P. verrucosa",
                                  str_detect(Species_ID, "Amu") ~ "A. muricata"))%>%
  
  mutate(across(c(Species_ID, Prec, Colony, Frg_n, Treatment, Species_short),
                as.factor),
         across(c(YII, F, Fm), as.numeric))%>%
  
  mutate (Species =  factor(Species, level = c("Galaxea fascicularis", "Porites rus", "Acropora muricata", "Montipora digitata", 
                                               "Pocillopora verrucosa", "Stylophora pistillata"
  )))%>%
  # view()
  
  group_by(Species_ID, Day, Colony, Treatment, Prec, Frg_n, Species_short )%>%
  dplyr::summarise(F_mean = mean(F),F_sd = sd(F),
            Fm_mean = mean(Fm),Fm_sd = sd(Fm),
            YII_mean = mean(YII),YII_sd = sd(YII)
  )%>%
  
  mutate_if(is.numeric, round, 3) %>%
  ungroup() 

## Theme model----
### Theme PAM ----
theme_recovery <- theme_classic()+
  theme (
         axis.text = element_text(size=20, face="bold"),
         axis.text.x=element_blank(),
         axis.title = element_text(size = 17, face="bold"),
         axis.title.y = element_text( margin = margin(t = 0, r = 10, b = 0, l = 0)),
         legend.position= "none", 
         legend.justification = c(0, 0),
         legend.direction = "horizontal",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_rect(colour = "black"),
         strip.background = element_blank(),
         strip.text.x = element_blank(),
         strip.text.y = element_blank(), 
         plot.margin=unit(c(-0.4,0.1,-0.4,-0.9), "cm"))

  
  

### Theme Tissue color ----

theme_recovery_IMG <- theme_classic()+
  
  theme (
         axis.text = element_text(size=20, face="bold"),
         axis.title = element_text(size = 17, face="bold"),
         axis.title.y = element_text( margin = margin(t = 0, r = 10, b = 0, l = 0)),
         legend.position= "none", 
         legend.justification = c(0, 0),
         legend.direction = "horizontal",
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         strip.text.y = element_blank(),
         strip.text.x = element_text(size = 20, face = "bold.italic"),
         strip.placement = "outside",
         plot.margin=unit(c(-0.4,0.1,0,-0.9), "cm") ,
         panel.border = element_rect(color = "black", fill = NA)  
         
         )
         
  
### Theme survival ----
theme_recovery_surv <- theme_classic()+
  theme (
    axis.text = element_text(size=20, face="bold"),
    axis.text.x=element_blank(),
    axis.title = element_text(size = 20, face="bold"),
    axis.title.y = element_text( margin = margin(t = 0, r = 10, b = 0, l = 0)),
    legend.position= "none", 
    legend.justification = c(0, 0),
    legend.direction = "horizontal",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(), 
    plot.margin=unit(c(-0.4,0.1,-0.4,-0.9), "cm"))      







##  PLOTS - Effective quantum yield (PAM) -------------------------------------------
### PAM - Porites ----

(Pru_PAM_recovery <- 
    Recovery_data %>%
    filter(Species_ID == "Pru"
           )%>%
    mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                           str_ends(Day, "2") ~ "0",
                           str_detect(Day, "15") ~ "10",
                           str_starts(Day, "30") ~ "20"))%>%
    
    group_by(Day, Prec, Treatment) %>% 
    dplyr::mutate(YII_mean_all = mean(YII_mean)) %>%
    
    mutate (Days =  as.numeric(Days))%>%
    
    
    ggplot(mapping = aes(x = Days))+
    
    geom_point(aes(y = YII_mean, color = Prec, shape = Treatment), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
               alpha = 0.6, fill = NA, set.seed(1) ) +
    geom_line(aes(y = YII_mean_all, color = Prec, linetype = Treatment), linewidth = 1) +
    geom_point(aes(y = YII_mean_all, color = Prec), size = 3)+
    scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
    scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
    
    annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "#1AD45E")+
    annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "#F75840")+
    
    ggeasy::easy_center_title() +
    ylim(0, 0.72)+
    
    theme_recovery+
    theme(
      axis.text.y=element_blank()
    )+
    
    
    scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                        label = c("-3", "0", "15", "30"))+
    
    labs(x = "", y = "")+ #y = "ΔF/Fm'"
    facet_wrap(Species_short ~., strip.position ="bottom"))


### PAM - Galaxea ----
(Gfa_PAM_recovery <- 
  Recovery_data %>%
  filter(Species_ID == "Gfa"
    
         )%>%
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "3") ~ "0",
                          str_detect(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
  
   group_by(Day, Prec, Treatment) %>% 
   dplyr::mutate(YII_mean_all = mean(YII_mean)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = YII_mean, color = Prec, shape = Treatment), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   geom_line(aes(y = YII_mean_all, color = Prec, linetype = Treatment), linewidth = 1) +
   geom_point(aes(y = YII_mean_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   ggeasy::easy_center_title() +
   ylim(0, 0.72)+
   
   theme_recovery+

   # theme(
   #   plot.margin=unit(c(-0.4,0.1,-0.4,0), "cm")
   # )+
   # 
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-3", "0", "15", "30"))+
   
  labs(x = "", y = "ΔF/Fm'")+ #y = "ΔF/Fm'"
   # labs(x = "", y = "")+
   # theme(
   #   axis.text.y=element_blank(),
   # )+
   facet_wrap(Species_short ~., strip.position ="bottom"))






### PAM - Montipora ----
(Mdi_PAM_recovery <- 
    Recovery_data %>%
    filter(Species_ID == "Mdi")%>%
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "3") ~ "0",
                          str_detect(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
    
    group_by(Day, Prec, Treatment) %>% 
    dplyr::mutate(YII_mean_all = mean(YII_mean)) %>%
    
    mutate (Days =  as.numeric(Days))%>%
    
    
    ggplot(mapping = aes(x = Days))+
    
    geom_point(aes(y = YII_mean, color = Prec, shape = Treatment), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
               alpha = 0.6, fill = NA, set.seed(1) ) +
    
    # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
    #               ) +
    geom_line(aes(y = YII_mean_all, color = Prec, linetype = Treatment), linewidth = 1) +
    geom_point(aes(y = YII_mean_all, color = Prec), size = 3)+
    scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
    scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
    
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   ylim(0, 0.72)+
   
   ggeasy::easy_center_title() +
   
   theme_recovery+
   
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-3", "0", "15", "30"))+
   
   labs(x = "", y = "ΔF/Fm'")+ #y = "ΔF/Fm'"
   labs(x = "", y = "")+
   theme(
     axis.text.y=element_blank(),
     # legend.position= "bottom",
   )+
   facet_wrap(Species_short ~., strip.position ="bottom"))




### PAM - Pocillopora ----
(Pve_PAM_recovery <- 
   Recovery_data %>%
   filter(Species_ID == "Pve", Treatment == "Heat"& Day <= 1 | Treatment == "Control")%>%
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "1") ~ "0",
                          str_ends(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
   
   group_by(Day, Prec, Treatment) %>% 
   dplyr::mutate(YII_mean_all = mean(YII_mean)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = YII_mean, color = Prec, shape = Treatment), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = YII_mean_all, color = Prec, linetype = Treatment), size = 1) +
   geom_point(aes(y = YII_mean_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   ggeasy::easy_center_title() +
   ylim(0, 0.72)+
   
   theme_recovery+
   theme(
     axis.text.y=element_blank()
   )+
   
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-1", "0", "15", "30"))+
   
   labs(x = "", y = "")+ #y = "ΔF/Fm'"
   facet_wrap(Species_short ~., strip.position ="bottom"))



### PAM - Stylophora ----
(Spi_PAM_recovery <- 
   Recovery_data %>%
   filter(Species_ID == "Spi", Treatment == "Heat"& Day <= 15 | Treatment == "Control")%>%
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "1") ~ "0",
                          str_ends(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
   
   group_by(Day, Prec, Treatment) %>% 
   dplyr::mutate(YII_mean_all = mean(YII_mean)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = YII_mean, color = Prec, shape = Treatment), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = YII_mean_all, color = Prec, linetype = Treatment), size = 1) +
   geom_point(aes(y = YII_mean_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   ggeasy::easy_center_title() +
   
   ylim(0, 0.72)+
   theme_recovery+
   theme(
     axis.text.y=element_blank()
   )+
   
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-1", "0", "15", "30"))+
   
   labs(x = "", y = "")+ #y = "ΔF/Fm'"
   facet_wrap(Species_short ~., strip.position ="bottom"))




### PAM - Acropora ----
(Amu_PAM_recovery <- 
   Recovery_data %>%
   filter(Species_ID == "Amu", Treatment == "Control"& Day <= 1 | Treatment == "Heat"& Day <= 1 )%>%
   
 
   
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "1") ~ "0",
                          str_ends(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
   
   group_by(Day, Prec, Treatment) %>% 
   dplyr::mutate(YII_mean_all = mean(YII_mean)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = YII_mean, color = Prec, shape = Treatment), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = YII_mean_all, color = Prec, linetype = Treatment), size = 1) +
   geom_point(aes(y = YII_mean_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   ggeasy::easy_center_title() +
   ylim(0, 0.72)+
   theme_recovery+
   theme(
     axis.text.y=element_blank()
   )+
   
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-1", "0", "15", "30"), limits = c(-4.5, 20.55))+
   
   labs(x = "", y = "")+ #y = "ΔF/Fm'"
   facet_wrap(Species_short ~., strip.position ="bottom"))

   

## PLOTS - Tissue color (Clr) -----------------------------------------------------------------------





### Clr - Galaxea ----
(Gfa_Clr_recovery <- 
   Clr %>%
   filter(ID == "Gfa")%>%
   
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "3") ~ "0",
                          str_ends(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
   
   group_by(Day, Prec, Trtm) %>% 
   dplyr::mutate(grays_sum_all = mean(gray_sum)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = gray_sum, color = Prec, shape = Trtm), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = grays_sum_all, color = Prec, linetype = Trtm), size = 1) +
   geom_point(aes(y = grays_sum_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   ylim(30, 205)+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   theme_recovery_IMG+
   # theme(
   #   plot.margin=unit(c(-0.4,0.1,0,0), "cm")
   # )+

   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-5", "0", "15", "30"))+
   
   labs(x = "", y = "Tissue color intensity")+
   # labs(x = "", y = "")+
   # labs(x = "", y = "")+
   # theme(
   #   axis.text.y=element_blank(),
   # )+
   facet_wrap(Species_short ~., strip.position ="bottom"))







### Clr - Montipora ----
(Mdi_Clr_recovery <- 
   Clr %>%
   filter(ID == "Mdi")%>%
   
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "3") ~ "0",
                          str_ends(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
   
   group_by(Day, Prec, Trtm) %>% 
   dplyr::mutate(grays_sum_all = mean(gray_sum)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = gray_sum, color = Prec, shape = Trtm), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = grays_sum_all, color = Prec, linetype = Trtm), size = 1) +
   geom_point(aes(y = grays_sum_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   ylim(30, 205)+
   theme_recovery_IMG+
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-5", "0", "15", "30"))+
   
   labs(x = "", y = "Tissue color intensity")+
   labs(x = "", y = "")+
   theme(
     axis.text.y=element_blank(),
 
   )+
   facet_wrap(Species_short ~., strip.position ="bottom"))








### Clr - Pocillopora ----
(Pve_Clr_recovery <- 
   Clr %>%
  filter(ID == "Pve")%>%
  
  mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                         str_ends(Day, "1") ~ "0",
                         str_ends(Day, "15") ~ "10",
                         str_starts(Day, "30") ~ "20"))%>%
  
   mutate (Days =  as.numeric(Days))%>%
   
   filter(ID == "Pve" & Trtm == "Heat"& Days <= 10 | ID == "Pve" & Trtm == "Control" & Days <= 20) %>%
   
   
   
   group_by(Day, Prec, Trtm) %>% 
   dplyr::mutate(grays_sum_all = mean(gray_sum)) %>%
   
   mutate (Day =  as.numeric(Day))%>%
   mutate (Days =  as.numeric(Days))%>%
   # 
   # filter(Trtm == "Heat"& Day <= 15| Trtm == "Control"& Day <= 30)%>%
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = gray_sum, color = Prec, shape = Trtm), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = grays_sum_all, color = Prec, linetype = Trtm), size = 1) +
   geom_point(aes(y = grays_sum_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   ylim(30, 205)+
   theme_recovery_IMG+
   theme(
     axis.text.y=element_blank(),
   )+
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-1", "0", "15", "30"))+
   
   labs(x = "", y = "")+
   facet_wrap(Species_short ~., strip.position ="bottom"))








### Clr - Stylopora ----
(Spi_Clr_recovery <- 
   Clr %>%
   filter(ID == "Spi")%>%
   
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "1") ~ "0",
                          str_ends(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
   
   group_by(Day, Prec, Trtm) %>% 
   dplyr::mutate(grays_sum_all = mean(gray_sum)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   filter(ID == "Spi" & Trtm == "Heat"& Days <= 10 | ID == "Spi" & Trtm == "Control" & Days <= 20) %>%
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = gray_sum, color = Prec, shape = Trtm), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = grays_sum_all, color = Prec, linetype = Trtm), size = 1) +
   geom_point(aes(y = grays_sum_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   ylim(30, 206)+
   theme_recovery_IMG+
   theme(
     axis.text.y=element_blank()
   )+
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-1", "0", "15", "30"))+
   
   labs(x = "", y = "")+
   facet_wrap(Species_short ~., strip.position ="bottom"))






### Clr - Porites ----
(Pru_Clr_recovery <- 
   Clr%>%
   # Clr2_fixed %>%
   filter(ID == "Pru")%>%
   
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "2") ~ "0",
                          str_ends(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
   
   group_by(Day, Prec, Trtm) %>% 
   dplyr::mutate(grays_sum_all = mean(gray_sum)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = gray_sum, color = Prec, shape = Trtm), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = grays_sum_all, color = Prec, linetype = Trtm), size = 1) +
   geom_point(aes(y = grays_sum_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+

   ylim(30, 205)+
   theme_recovery_IMG+
   theme(
     axis.text.y=element_blank()
   )+
   
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-3", "0", "15", "30"))+

   # labs(x = "Days after heat stress", y = "Tissue color intensity")+
   labs(x = "", y = "")+
   facet_wrap(Species_short ~., strip.position ="bottom"))





### Clr - Acropora ----
(Amu_Clr_recovery <- 
   Clr %>%
   filter(ID == "Amu"& Day<=1)%>%
   
   mutate(Days =case_when(str_starts(Day, "0") ~ "-3",
                          str_ends(Day, "1") ~ "0",
                          str_ends(Day, "15") ~ "10",
                          str_starts(Day, "30") ~ "20"))%>%
   
   group_by(Day, Prec, Trtm) %>% 
   dplyr::mutate(grays_sum_all = mean(gray_sum)) %>%
   
   mutate (Days =  as.numeric(Days))%>%
   
   filter(ID == "Amu" & Trtm == "Heat"& Days <= 10 | ID == "Amu" & Trtm == "Control" & Days <= 20) %>%
   
   
   ggplot(mapping = aes(x = Days))+
   
   geom_point(aes(y = gray_sum, color = Prec, shape = Trtm), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
              alpha = 0.6, fill = NA, set.seed(1) ) +
   
   # geom_errorbar(test_pru, mapping = aes(ymin=YII_mean-ci, ymax=YII_mean+ci,color = Prec ), width= 0.8, stat = "identity"
   #               ) +
   geom_line(aes(y = grays_sum_all, color = Prec, linetype = Trtm), size = 1) +
   geom_point(aes(y = grays_sum_all, color = Prec), size = 3)+
   scale_linetype_manual(breaks=c("Control","Heat"), values=c(2,1))+
   scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
   
   annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#1AD45E")+
   annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
            alpha = .1,fill = "#F75840")+
   
   ylim(27, 205)+
   theme_recovery_IMG+
   theme(
     axis.text.y=element_blank()
   )+
   scale_x_continuous (breaks = c(-3, 0 , 10, 20),
                       label = c("-1", "0", "15", "30"), limits = c(-4.5, 20.55))+
   
   # labs(x = "Days after heat stress", y = "Tissue color intensity")+
   labs(x = "", y = "")+
   
   facet_wrap(Species_short ~., strip.position ="bottom"))






# SURVIVAL RATE ----------------------------------------------------------

# Cleaning and table preparation. Fragments 3 and 4 at "day-30 time point are selected together with the survival information

## Df subset ----
survival_table <- H_rcv %>%
  
  dplyr::filter(Frg_n %in% c("3", "4") & Day == "30"& obs=="1")%>%
  select(c(Species_ID, Colony, Prec, Frg_n, Survived_day, Status))%>%
  mutate(across(c(Species_ID, Prec, Colony, Frg_n),
                as.factor),
         across(c(Survived_day, Status), as.numeric))

## Porites rus (Pru)----

#subset the df to include only porites rus
Pru_survival <- survival_table %>%
  dplyr::filter(Species_ID == "Pru", Frg_n =="4")

#fit data to build the survival plot (see below)
Pru.sfit <- survfit(Surv(Survived_day, Status)~Prec, data=Pru_survival)

### Cox model ----
# test statistical difference between survival rates
Pru.cfit <- coxph(Surv(Survived_day, Status)~ Prec, data= Pru_survival)
Pru.cfit

### Post-hoc test----
Pru.cfit_pairwise <-pairwise_survdiff(Surv(Survived_day, Status) ~ Prec, data = Pru_survival)
Pru.cfit_pairwise

#plot 
(Surv_Pru <- ggsurvplot(sfit, conf.int=F, pval=F, risk.table=F, 
                       legend.labs=c("Ambient", "ST", "VT"), legend.title="Preconditioning", 
                       ggtheme = theme_classic(base_size=20),
                       # xlim = c(0, 30), break.x.by = 10,
                       palette = c("#4682B4", "#B4B446", "#D4711C"),
                       #surv.median.line = "hv",
                       xlab = "Days after heat-stress", ylab = c("Survival rate (%)"),

                       #title = expression(paste("Kaplan-Meier Curve for coral Survival rate, ", italic("Porites rus"))), 
                       risk.table.height=.15))



  #Extracting the table with e values from the survival model
  Surv_Pru_df<-Surv_Pru[["plot"]][["data"]]


  #Delete the not necessary columns
  Surv_Pru_df <- Surv_Pru_df[, -c(2:4, 6:9)]
  
  #Add the Pre-heat points
  Surv_Pru_df[nrow(Surv_Pru_df) + 1,] = c(-3,1.000, "Ambient")
  Surv_Pru_df[nrow(Surv_Pru_df) + 1,] = c(-3,1.000, "ST")
  Surv_Pru_df[nrow(Surv_Pru_df) + 1,] = c(-3,1.000, "VT")

  ### GGPLOT ----
  
  (Pru_surv <- Surv_Pru_df %>%
     
     mutate(across( c( time, surv), as.numeric))%>%
     
     mutate(ID = "P. rus")%>%
     
     group_by(Prec, time)%>%
     
     mutate(Days = ifelse(time > 0, (time * 20)/30, time ))%>%
     
     ungroup()%>%
     
     mutate(across(c(time, surv, Days), as.numeric))%>%
     
     ggplot(mapping = aes(x = Days, y = surv))+
     
     geom_point(aes(y =surv, color = Prec), alpha = 1, fill = NA , set.seed(1),size = 3)+ 
     
     geom_step(aes(y = surv, color = Prec), linewidth = 1)+
     
     scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
     ylim(0,1)+
     
     new_scale("color")+

     geom_point(aes(y =surv, x= Days, color = Prec),
                alpha = 0, fill = NA, set.seed(1) )+
     
     annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#1AD45E")+
     annotate("rect", xmin = -1.5, xmax = 1.35, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#F75840")+
     
     scale_x_continuous(breaks = c(-3,0,10,20),
                        labels = c(-2,0,15,30))+
     labs(x = "", y = "")+
     theme_recovery +
     theme (
       plot.margin=unit(c(0.5,0.1,-0.4,-0.9), "cm"))+
     theme(
       axis.text.y=element_blank()
     )+
     
     facet_wrap(ID ~., strip.position ="bottom")
   
  )



### ggsave ----

#ggsave(here ("Output", "Pru_survival.png"), dpi = 400, units = "cm", width = 15, height = 15)


## Acropora muricata (Amu) ----
  
#subset the df to include only Acropora muricata
Amu_survival <- survival_table %>%
  dplyr::filter(  Species_ID == "Amu", Frg_n =="4") 

Amu.sfit <- survfit(Surv(Survived_day, Status)~Prec, data= Amu_survival)

### Cox model ----

Amu.cfit <- coxph(Surv(Survived_day, Status)~ Prec, data=Amu_survival)
Amu.cfit

### Post-hoc test----
Amu.cfit_pairwise <-pairwise_survdiff(Surv(Survived_day, Status) ~ Prec, data = Amu_survival)
Amu.cfit_pairwise

#plot 
(Surv_Amu <- 
  ggsurvplot(sfit, conf.int=F, pval=F, risk.table=F, 
                       legend.labs=c("Ambient", "ST", "VT"), legend.title="Preconditioning", 
                       ggtheme = theme_classic2(base_size=20),
            xlim = c(0, 30), break.x.by = 10,
            #surv.median.line = "hv",
            
                       palette = c("#4682B4", "#B4B446", "#D4711C"),
                       
                       
                       #title = expression(paste("Kaplan-Meier Curve for coral Survival rate, ", italic("Porites rus"))), 
                       risk.table.height=.15))


#Extracting the table with e values from the survival model 
Surv_Amu_df<-Surv_Amu[["plot"]][["data"]]

#Delete the not necessary columns
Surv_Amu_df <- Surv_Amu_df[, -c(2:4, 6:9)]

#Add the Pre-heat points
Surv_Amu_df[nrow(Surv_Amu_df) + 1,] = c(-3,1.000, "Ambient")
Surv_Amu_df[nrow(Surv_Amu_df) + 1,] = c(-3,1.000, "ST")
Surv_Amu_df[nrow(Surv_Amu_df) + 1,] = c(-3,1.000, "VT")


### GGPLOT ----

(Amu_surv <- Surv_Amu_df %>%
 
    mutate(across( c( time, surv), as.numeric))%>%
  
    mutate(ID = "A. muricata")%>%
    
    group_by(Prec, time)%>%
   
    mutate(Days = ifelse(time > 0, (time * 25)/30, time ))%>% # attempt to create plots with the same size
    
    ungroup()%>%
    
    mutate(across(c(time, surv, Days), as.numeric))%>%
    
    ggplot(mapping = aes(x = Days, y = surv))+
    
    geom_point(aes(y =surv, color = Prec), alpha = 1, fill = NA , set.seed(1),size = 3)+ 
    
    geom_step(aes(y = surv, color = Prec), size = 1)+
   
    scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
    
    new_scale("color")+
    
    geom_point(aes(y =surv, x= Days, color = Prec), position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.15),
               alpha = 0, fill = NA, set.seed(1) )+
    
   ylim(0, 1)+
    annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "#1AD45E")+
    annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
             alpha = .1,fill = "#F75840")+
    
    scale_x_continuous(breaks = c(-3,0,10,20),
                       labels = c(-2,0,15,30), limits = c(-4.5, 20.55))+
    labs(x = "", y = "")+
    theme_recovery +
   theme(
     axis.text.y=element_blank()
     ,
     plot.margin=unit(c(0.5,0.1,-0.4,-0.9), "cm")
     )+
    
    facet_wrap(ID ~., strip.position ="bottom")
    
    )

### ggsave ----

#ggsave(here ("Output", "Amu_survival.png"), dpi = 400, units = "cm", width = 30, height = 15)


## Pocillopora verrucosa (Pve) ----

#subset the df to include only Pocillopora verrucosa
  Pve_survival <- survival_table %>%
    dplyr::filter(  Species_ID == "Pve", Frg_n =="4") 
  
  Pve.sfit <- survfit(Surv(Survived_day, Status)~Prec, data= Pve_survival)
  
  ### Cox model ----
  
  Pve.cfit <- coxph(Surv(Survived_day, Status)~ Prec, data=Pve_survival)
  Pve.cfit
  
  ### Post-hoc test----
  Pve.cfit_pairwise <-pairwise_survdiff(Surv(Survived_day, Status) ~ Prec, data = Pve_survival)
  Pve.cfit_pairwise 
  
  #plot 
  (Surv_Pve <- 
  ggsurvplot(sfit, conf.int=F, pval=F, risk.table=F, 
             legend.labs=c("Ambient", "ST", "VT"), legend.title="Preconditioning", 
             ggtheme = theme_classic2(base_size=20),
             xlim = c(0, 30), break.x.by = 10,
             palette = c("#4682B4", "#B4B446", "#D4711C"),
             risk.table.height=.15))
 
  
  #Extracting the table with e values from the survival model 
  Surv_Pve_df<-Surv_Pve[["plot"]][["data"]]
  
  #Delete the not necessary columns 
  Surv_Pve_df <- Surv_Pve_df[, -c(2:4, 6:9)]

  
  #Add the Pre-heat points
  Surv_Pve_df[nrow(Surv_Pve_df) + 1,] = c(-3,1.000, "Ambient")
  Surv_Pve_df[nrow(Surv_Pve_df) + 1,] = c(-3,1.000, "ST")
  Surv_Pve_df[nrow(Surv_Pve_df) + 1,] = c(-3,1.000, "VT")
  
  
  ### GGPLOT ----
  
  (Pve_surv <- Surv_Pve_df %>%
     
     mutate(across( c( time, surv), as.numeric))%>%
     
     mutate(ID = "P. verrucosa")%>%
     
     group_by(Prec, time)%>%
     
     mutate(Days = ifelse(time > 0, (time * 20)/30, time ))%>%
     
     ungroup()%>%
     
     mutate(across(c(time, surv, Days), as.numeric))%>%
     
     ggplot(mapping = aes(x = Days, y = surv))+
     
     geom_point(aes(y =surv, color = Prec), alpha = 1, fill = NA , set.seed(1),size = 3)+ 
     
     geom_step(aes(y = surv, color = Prec), size = 1)+
     
     scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
     
     new_scale("color")+
     ylim(0, 1)+
     
     geom_point(aes(y =surv, x= Days, color = Prec), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
                alpha = 0, fill = NA, set.seed(1) )+
     
     annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#1AD45E")+
     annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#F75840")+
     
     scale_x_continuous(breaks = c(-3,0,10,20),
                        labels = c(-2,0,15,30))+
     labs(x = "", y = "")+
     theme_recovery +
     theme(
       axis.text.y=element_blank()
     )+
     theme (
       plot.margin=unit(c(0.5,0.1,-0.4,-0.9), "cm"))+
     
     facet_wrap(ID ~., strip.position ="bottom")
   
  )
  
  ### ggsave ----
  
  #ggsave(here ("Output", "Pve_mortality.png"), dpi = 400, units = "cm", width = 20, height = 15)
  
  
   
## Montopora digitata (Mdi) ----
  
  #subset the df to include only Montopora digitata
  Mdi_survival <- survival_table %>%
    dplyr::filter(  Species_ID == "Mdi", Frg_n =="4") 
  
  Mdi.sfit <- survfit(Surv(Survived_day, Status)~Prec, data= Mdi_survival)
  
  ### Cox model ----
  
  Mdi.cfit <- coxph(Surv(Survived_day, Status)~ Prec, data=Mdi_survival) # problems with samples size and number of events (dead vs survived)
  Mdi.cfit
  
  ### Post-hoc test----
  Mdi.cfit_pairwise <-pairwise_survdiff(Surv(Survived_day, Status) ~ Prec, data = Mdi_survival)
  Mdi.cfit_pairwise 
  
  #plot 
  (Surv_Mdi <- 
    ggsurvplot(sfit, conf.int=F, pval=F, risk.table=F, 
               legend.labs=c("Ambient", "ST", "VT"), legend.title="Preconditioning", 
               ggtheme = theme_classic2(base_size=20),
               xlim = c(0, 30), break.x.by = 10,
               palette = c("#4682B4", "#B4B446", "#D4711C"),
               risk.table.height=.15))
  
  
  #Extracting the table with e values from the survival model 
  Surv_Mdi_df<-Surv_Mdi[["plot"]][["data"]]
  
  #Delete the not necessary columns 
  Surv_Mdi_df <- Surv_Mdi_df[, -c(2:4, 6:9)]

  
  #Add the Pre-heat points
  Surv_Mdi_df[nrow(Surv_Mdi_df) + 1,] = c(-3,1.000, "Ambient")
  Surv_Mdi_df[nrow(Surv_Mdi_df) + 1,] = c(-3,1.000, "ST")
  Surv_Mdi_df[nrow(Surv_Mdi_df) + 1,] = c(-3,1.000, "VT")
  
  
  ### GGPLOT ----
  
  (Mdi_surv <- Surv_Mdi_df %>%
     
     mutate(across( c( time, surv), as.numeric))%>%
     
     mutate(ID = "M. digitata")%>%
     
     group_by(Prec, time)%>%
     
     mutate(Days = ifelse(time > 0, (time * 20)/30, time ))%>%
     
     ungroup()%>%
     
     mutate(across(c(time, surv, Days), as.numeric))%>%
     
     ggplot(mapping = aes(x = Days, y = surv))+
     
     geom_point(aes(y =surv, color = Prec), alpha = 1, fill = NA , set.seed(1),size = 3)+ 
     
     geom_step(aes(y = surv, color = Prec), size = 1)+
     
     scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
     ylim(0,1)+
     
     new_scale("color")+
     
     
     geom_point(aes(y =surv, x= Days, color = Prec), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
                alpha = 0, fill = NA, set.seed(1) )+
     
     annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#1AD45E")+
     annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#F75840")+
     
     scale_x_continuous(breaks = c(-3,0,10,20),
                        labels = c(-2,0,15,30))+
     labs(x = "", y = "Survival rate (%)")+
     theme_recovery +
     theme (
       plot.margin=unit(c(0.5,0.1,-0.4,-0.9), "cm"))+
     labs(x = "", y = "")+
     theme(
       axis.text.y=element_blank(),
     )+
     facet_wrap(ID ~., strip.position ="bottom")
   
  )
  ### ggsave ----
  
  #ggsave(here ("Output", "Mdi_mortality.png"), dpi = 400, units = "cm", width = 20, height = 15)
  

  
## Galaxea fascicularis (Gfa) ----
  
  #subset the df to include only Galaxea fascicularis 
  Gfa_survival <- survival_table %>%
    dplyr::filter(  Species_ID == "Gfa", Frg_n =="4") 
  
  Gfa.sfit <- survfit(Surv(Survived_day, Status)~Prec, data= Gfa_survival)
  
  ### Cox model ----
  
  Gfa.cfit <- coxph(Surv(Survived_day, Status)~ Prec, data=Gfa_survival)# problems with samples size and number of events (dead vs survived)
  Gfa.cfit
  
  ### Post-hoc test----
  Gfa.cfit_pairwise <-pairwise_survdiff(Surv(Survived_day, Status) ~ Prec, data = Gfa_survival) 
  Gfa.cfit_pairwise

  
  #plot 
  (Surv_Gfa <- 
    ggsurvplot(sfit, conf.int=F, pval=F, risk.table=F, 
               legend.labs=c("Ambient", "ST", "VT"), legend.title="Preconditioning", 
               ggtheme = theme_classic2(base_size=20),
               xlim = c(0, 30), break.x.by = 10,
               #surv.median.line = "hv",
               
               palette = c("#4682B4", "#B4B446", "#D4711C"),
               
               
               #title = expression(paste("Kaplan-Meier Curve for coral Survival rate, ", italic("Porites rus"))), 
               risk.table.height=.15))
  
  #Extracting the table with e values from the survival model 
  Surv_Gfa_df<-Surv_Gfa[["plot"]][["data"]]
  
  #Delete the not necessary columns 
  Surv_Gfa_df <- Surv_Gfa_df[, -c(2:4, 6:9)]

  
  #Add the Pre-heat points
  Surv_Gfa_df[nrow(Surv_Gfa_df) + 1,] = c(-3,1.000, "Ambient")
  Surv_Gfa_df[nrow(Surv_Gfa_df) + 1,] = c(-3,1.000, "ST")
  Surv_Gfa_df[nrow(Surv_Gfa_df) + 1,] = c(-3,1.000, "VT")
  
  
  ### GGPLOT ----
  
  (Gfa_surv <- Surv_Gfa_df %>%
     
     mutate(across( c( time, surv), as.numeric))%>%
     
     mutate(ID = "M. digitata")%>%
     
     group_by(Prec, time)%>%
     
     mutate(Days = ifelse(time > 0, (time * 20)/30, time ))%>%
     
     ungroup()%>%
     
     mutate(across(c(time, surv, Days), as.numeric))%>%
     
     ggplot(mapping = aes(x = Days, y = surv))+
     
     geom_point(aes(y =surv, color = Prec), alpha = 1, fill = NA , set.seed(1),size = 3)+ 
     
     geom_step(aes(y = surv, color = Prec), size = 1)+
     
     scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
     ylim(0,1)+
     
     new_scale("color")+
     
     geom_point(aes(y =surv, x= Days, color = Prec), position = position_jitterdodge(dodge.width = 1, jitter.width = 0.15),
                alpha = 0, fill = NA, set.seed(1) )+
     
     annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#1AD45E")+
     annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#F75840")+
     
     scale_x_continuous(breaks = c(-3,0,10,20),
                        labels = c(-2,0,15,30))+
     labs(x = "", y = "Survival rate (%)")+
     theme_recovery +
     theme (
       plot.margin=unit(c(0.5,0.1,-0.4,0), "cm"))+
     facet_wrap(ID ~., strip.position ="bottom")
   
  )
  
  ### ggsave ----
  
  #ggsave(here ("Output", "Gfa_mortality.png"), dpi = 400, units = "cm", width = 20, height = 15)

  
  
  ## Stylophora pistillata (Spi) ----
 
  #subset the df to include only Stylophora pistillata  
  Spi_survival <- survival_table %>%
    dplyr::filter(  Species_ID == "Spi", Frg_n =="4") 
  
  Spi.sfit <- survfit(Surv(Survived_day, Status)~Prec, data= Spi_survival)
  
  ### Cox model ----
  
  Spi.cfit <- coxph(Surv(Survived_day, Status)~ Prec, data=Spi_survival)
  Spi.cfit
  
  ### Post-hoc test----
  Spi.cfit_pairwise <-pairwise_survdiff(Surv(Survived_day, Status) ~ Prec, data = Spi_survival) 
  Spi.cfit_pairwise
  
  #plot 
  (Surv_Spi <- 
      ggsurvplot(sfit, conf.int=F, pval=F, risk.table=F, 
                 legend.labs=c("Ambient", "ST", "VT"), legend.title="Preconditioning", 
                 ggtheme = theme_classic2(base_size=20),
                 xlim = c(0, 30), break.x.by = 10,
                 palette = c("#4682B4", "#B4B446", "#D4711C"),
                 risk.table.height=.15))
  
  
  #Extracting the table with e values from the survival model 
  Surv_Spi_df<-Surv_Spi[["plot"]][["data"]]
  
  #Delete the not necessary columns 
  Surv_Spi_df <- Surv_Spi_df[, -c(2:4, 6:9)]

  
  #Add the Pre-heat points
  Surv_Spi_df[nrow(Surv_Spi_df) + 1,] = c(-3,1.000, "Ambient")
  Surv_Spi_df[nrow(Surv_Spi_df) + 1,] = c(-3,1.000, "ST")
  Surv_Spi_df[nrow(Surv_Spi_df) + 1,] = c(-3,1.000, "VT")
  
  
  ### GGPLOT ----
  
  (Spi_surv <- Surv_Spi_df %>%
     
     mutate(across( c( time, surv), as.numeric))%>%
     
     mutate(ID = "M. digitata")%>%
     
     group_by(Prec, time)%>%
     
     mutate(Days = ifelse(time > 0, (time * 25)/30, time ))%>%
     
     ungroup()%>%
     
     mutate(across(c(time, surv, Days), as.numeric))%>%
     
     ggplot(mapping = aes(x = Days, y = surv))+
     
     geom_point(aes(y =surv, color = Prec), alpha = 1, fill = NA , set.seed(1),size = 3)+ 
     
     geom_step(aes(y = surv, color = Prec), size = 1)+
     
     scale_color_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
     ylim(0,1)+
     
     new_scale("color")+
     
     geom_point(aes(y =surv, x= Days, color = Prec), position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.15),
                alpha = 0, fill = NA, set.seed(1) )+
     
     annotate("rect", xmin = -4.5, xmax = -1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#1AD45E")+
     annotate("rect", xmin = -1.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
              alpha = .1,fill = "#F75840")+
     
     scale_x_continuous(breaks = c(-3,0,10,20),
                        labels = c(-2,0,15,30), limits = c(-4.5, 20.55))+
     labs(x = "", y = "")+
     theme_recovery +
     theme(
       axis.text.y=element_blank(),
     )+
     theme (
       plot.margin=unit(c(0.5,0.1,-0.4,-0.9), "cm"))+
     
     facet_wrap(ID ~., strip.position ="bottom")
   
  )
  
  ### ggsave ----
  
  #ggsave(here ("Output", "Pve_mortality.png"), dpi = 400, units = "cm", width = 20, height = 15)
  
  
  
  
  
  
  
# Plot Assembly  --------------------------------------------
  
  # In this section all recovery plots (Survival, Photosynthetic efficiency, and tissue color) are combined first by species and then all together

  # Pru

  
  (Pru <- ggarrange(Pru_surv,Pru_PAM_recovery,Pru_Clr_recovery, 
                           align = "v", 
                           heights= c(1, 0.85, 1.2),
                           ncol=1))
  

  #ggsave(here ("Output/Recovery", "Pru_patch_new_tissue.png"), dpi = 400, units = "cm", width = 15, height = 25)

  #Gfa
  
    (Gfa <- ggarrange(Gfa_surv,Gfa_PAM_recovery,Gfa_Clr_recovery, 
                  align = "v", 
                  heights= c(1, 0.85, 1.2),
                  ncol=1))
  

  #ggsave(here ("Output/Recovery", "Gfa_patch_new_tissue.png"), dpi = 400, units = "cm", width = 15, height = 25)
  
  
  #Pve
  
  (Pve <- ggarrange(Pve_surv,Pve_PAM_recovery,Pve_Clr_recovery, 
                  align = "v",  
                  heights= c(1, 0.85, 1.2),
                  ncol=1))
  
  
  #ggsave(here ("Output/Recovery", "Pve_patch_new_tissue.png"), dpi = 400, units = "cm", width = 15, height = 25)
  
  
  #Amu

  (Amu <- ggarrange(Amu_surv,Amu_PAM_recovery,Amu_Clr_recovery, 
                  align = "v", 
                  heights= c(1, 0.85, 1.2),
                  ncol=1))
  
  
  #ggsave(here ("Output/Recovery", "Amu_patch_new_tissue.png"), dpi = 400, units = "cm", width = 15, height = 25)
  
  #Spi
  
  (Spi <- ggarrange(Spi_surv,Spi_PAM_recovery,Spi_Clr_recovery, 
                  align = "v", 
                  heights= c(1, 0.85, 1.2),
                  ncol=1))
  
  # ggsave(here ("Output/Recovery", "Spi_patch_new_tissue.png"), dpi = 400, units = "cm", width = 15, height = 25)
  
  #Mdi

  (Mdi <- ggarrange(Mdi_surv,Mdi_PAM_recovery,Mdi_Clr_recovery, 
                  align = "v",
                  heights= c(1, 0.85, 1.2),
                  # legend = "bottom",
                  ncol=1))
  
  
  # ggsave(here ("Output/Recovery", "Mdi_patch_new_tissue.png"), dpi = 400, units = "cm", width = 9, height = 28)
  


 ## All species ----
  
  ggarrange(Gfa, Pru, Amu, Mdi, Pve, Spi, 
            align = "h",  
            widths= c(1.33, 1,1,1,1,1),
            legend = "bottom",
            ncol=6)
  
  
 
  
  #ggsave(here ("Output", "Recovery_phase.png"), dpi = 400, units = "cm", width = 45, height = 25)
  
   
