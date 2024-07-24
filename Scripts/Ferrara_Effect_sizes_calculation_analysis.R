 

# --------------------------- Ferrara et al. 2024 - R script N. 4 --------------------------- #
#                                                                                             #
#                                                                                             #
#  Effect size calculation using dabestR package                                              #
#                                                                                             #
#  In this script we calculate the effect of preconditioning on baseline physioligy,          #
#  heat stress response, and recovery capacity improvements as Hedges'g effect size           #
#                                                                                             #
#  The effects are calculated on both metrics (photosynthetic efficiency and tissue color)    #
#                                                                                             #
#-------------------------------------------------------------------------------------------- #

#0 Library ----

library(magrittr)
library(dabestr)
library(tidyverse)
library(here)
options(scipen = 999)# avoid scientific annotations


## Load data ----

# Load effective quantum yield data
H_rcv <- read_csv(here("Data", "Ferrara_PAM_master_df.csv"))

# Load tissue color data
Clr  <- read_csv(here("Data", "Ferrara_tissue_color_master_df.csv"))


## Data preparation ----
# To calculate the effect of preconditioning as effect size, we first have to calculate the mean value for each measurement as explained in 
# previous scripts

H_rcv_mean <-H_rcv %>%
  
  mutate(across(c(Species_ID, Prec, Tank_n, Colony, Frg_n, Region, Treatment), as.factor),
         across(c(YII, F, Fm), as.numeric))%>%
  
  mutate (Species =  factor(Species, level = c("Galaxea fascicularis", "Porites rus", "Acropora muricata", "Montipora digitata", 
                                               "Pocillopora verrucosa", "Stylophora pistillata"
  )))%>%
  # mutate(Species= ifelse(Species_ID == "Mdi", "Montipora digitata", Species))%>%
  
  group_by(Species_ID,Prec, Colony, Treatment, Frg_n, Day,)%>%
  
  summarise(F_mean = mean(F),F_sd = sd(F),
            Fm_mean = mean(Fm),Fm_sd = sd(Fm),
            YII_mean = mean(YII),YII_sd = sd(YII)
  )%>%
  ungroup()

# 1.1 Effect of PRECONDITIONING on coral physiology -------------------------------------------------------------------------------------------

## Photosynthetic efficiency (PAM) ----
 

## Grouping ----
Prec.effect.H_rcv_mean<- H_rcv_mean%>%
  mutate(YII_mean = pmax(YII_mean + rnorm(length(YII_mean), mean = 0, sd = 1e-6), 0))%>% # to avoid the presence of ties, a small positive values is added to each measurement
  

  mutate(t= NULL, Date= NULL, Time= NULL, Survived_day=NULL, Status=NULL)%>%
    # to ran the analysis with the dabestR package is necessary to create a new variable to group samples (the function doesn't have a grouping option)
  mutate(Day_Trtm_Prec = 
           case_when(

       
              (Species_ID == "Amu" & Day == "0" & Prec == "Ambient")~ "Amu_AT", # Species ID - Preconditioning regime (AT = Ambient temperature)
              (Species_ID == "Amu" & Day == "0" & Prec == "ST")~ "Amu_ST",
              (Species_ID == "Amu" & Day == "0" & Prec == "VT")~ "Amu_VT",
              
              
              (Species_ID == "Gfa" & Day == "0" & Prec == "Ambient")~ "Gfa_AT",
              (Species_ID == "Gfa" & Day == "0" & Prec == "ST")~ "Gfa_ST",
              (Species_ID == "Gfa" & Day == "0" & Prec == "VT")~ "Gfa_VT",
              
              
              (Species_ID == "Mdi" & Day == "0" & Prec == "Ambient")~ "Mdi_AT",
              (Species_ID == "Mdi" & Day == "0" & Prec == "ST")~ "Mdi_ST",
              (Species_ID == "Mdi" & Day == "0" & Prec == "VT")~ "Mdi_VT",
              
              
              (Species_ID == "Pve" & Day == "0" & Prec == "Ambient")~ "Pve_AT",
              (Species_ID == "Pve" & Day == "0" & Prec == "ST")~ "Pve_ST",
              (Species_ID == "Pve" & Day == "0" & Prec == "VT")~ "Pve_VT",
              
              
              (Species_ID == "Pru" & Day == "0" & Prec == "Ambient")~ "Pru_AT",
              (Species_ID == "Pru" & Day == "0" & Prec == "ST")~ "Pru_ST",
              (Species_ID == "Pru" & Day == "0" & Prec == "VT")~ "Pru_VT",
              
              
              (Species_ID == "Spi" & Day == "0" & Prec == "Ambient")~ "Spi_AT",
              (Species_ID == "Spi" & Day == "0" & Prec == "ST")~ "Spi_ST",
              (Species_ID == "Spi" & Day == "0" & Prec == "VT")~ "Spi_VT"))
              



## Loading model ----
# Create a list of contrasts

groups<- list(
  c("Gfa_AT","Gfa_ST", "Gfa_VT"),
  c("Pru_AT","Pru_ST", "Pru_VT"),
  c("Amu_AT","Amu_ST", "Amu_VT"),
  c("Mdi_AT","Mdi_ST", "Mdi_VT"),
  c("Pve_AT","Pve_ST", "Pve_VT"),
  c("Spi_AT","Spi_ST", "Spi_VT")
)

# Load the data into the dabstR function
Load.Prec_Pam <- 
  Prec.effect.H_rcv_mean%>%
  filter(Day==0)%>%
  load( Day_Trtm_Prec, YII_mean, 
        idx= groups,
        colour = Prec)

## Calculate the effect size ----




Effect.Prec_Pam <- hedges_g(Load.Prec_Pam)

## Plot ----
dabest_plot(Effect.Prec_Pam,
            raw_marker_spread=0.5,
            swarm_label = "", contrast_label = "Hedges g")


# Save the plot
ggsave(here ("Output", "Baseline_Pam.png"), dpi = 400, units = "cm", width = 45, height = 20)


## Results extraction ----

# The exact effect size values are listed under "difference"
Prec_PAM_eff_g<- Effect.Prec_Pam[["boot_result"]]

write_excel_csv(Prec_PAM_eff_g, here("Output","Statistics/Prec_PAM_eff_g.csv"))





# 1.2 Effect of PRECONDITIONING on coral physiology -------------------------------------------------------------------------------------------

## Tissue color (Clr) ----


## Grouping ----

Prec.effect.Clr<- Clr%>%
 
  # to ran the analysis with the dabestR package is necessary to create a new variable to group samples (the function doesn't have a grouping option)
  mutate(Day_Trtm_Prec = 
           case_when(
             
             
             (ID == "Amu" & Day == "0" & Prec == "Ambient")~ "Amu_AT", # ID - Preconditioning regime (AT = Ambient temperature)
             (ID == "Amu" & Day == "0" & Prec == "ST")~ "Amu_ST",
             (ID == "Amu" & Day == "0" & Prec == "VT")~ "Amu_VT",
             
             
             (ID == "Gfa" & Day == "0" & Prec == "Ambient")~ "Gfa_AT",
             (ID == "Gfa" & Day == "0" & Prec == "ST")~ "Gfa_ST",
             (ID == "Gfa" & Day == "0" & Prec == "VT")~ "Gfa_VT",
             
             
             (ID == "Mdi" & Day == "0" & Prec == "Ambient")~ "Mdi_AT",
             (ID == "Mdi" & Day == "0" & Prec == "ST")~ "Mdi_ST",
             (ID == "Mdi" & Day == "0" & Prec == "VT")~ "Mdi_VT",
             
             
             (ID == "Pve" & Day == "0" & Prec == "Ambient")~ "Pve_AT",
             (ID == "Pve" & Day == "0" & Prec == "ST")~ "Pve_ST",
             (ID == "Pve" & Day == "0" & Prec == "VT")~ "Pve_VT",
             
             
             (ID == "Pru" & Day == "0" & Prec == "Ambient")~ "Pru_AT",
             (ID == "Pru" & Day == "0" & Prec == "ST")~ "Pru_ST",
             (ID == "Pru" & Day == "0" & Prec == "VT")~ "Pru_VT",
             
             
             (ID == "Spi" & Day == "0" & Prec == "Ambient")~ "Spi_AT",
             (ID == "Spi" & Day == "0" & Prec == "ST")~ "Spi_ST",
             (ID == "Spi" & Day == "0" & Prec == "VT")~ "Spi_VT"))




## Loading model ----
# Create a list of contrasts

groups<- list(
  c("Gfa_AT","Gfa_ST", "Gfa_VT"),
  c("Pru_AT","Pru_ST", "Pru_VT"),
  c("Amu_AT","Amu_ST", "Amu_VT"),
  c("Mdi_AT","Mdi_ST", "Mdi_VT"),
  c("Pve_AT","Pve_ST", "Pve_VT"),
  c("Spi_AT","Spi_ST", "Spi_VT")
)

# Load the data into the dabstR function
Load.Prec_Clr <- 
  Prec.effect.Clr%>%
  filter(Day==0)%>%
  group_by(Species)%>%
  load( Day_Trtm_Prec, gray_sum, 
        idx= groups,
        colour = Prec)

## Calculate the effect size ----




Effect.Prec_Clr <- hedges_g(Load.Prec_Clr)

## Plot ----
dabest_plot(Effect.Prec_Clr,
            raw_marker_spread=0.5,
            swarm_label = "", contrast_label = "Hedges g")


# Save the plot
ggsave(here ("Output", "Baseline_IMG.png"), dpi = 400, units = "cm", width = 45, height = 20)


## Results extraction ----

# The exact effect size values are listed under "difference"
Prec_Clr_eff_g<- Effect.Prec_Clr[["boot_result"]]

write_excel_csv(Prec_Clr_eff_g, here("Output","Statistics/Prec_Clr_eff_g.csv"))










# 2.1 Effect of HEAT on thermal tolerance -------------------------------------------------------------------------------------------

## Photosynthetic efficiency (PAM) ----

## Grouping ----
Heat.effect.H_rcv_mean<- H_rcv_mean%>%
  
  # to avoid the presence of ties, a small positive values is added to each measurement
  mutate(YII_mean = as.numeric(YII_mean+ round(runif(length(YII_mean), 0.00001, 0.00003),7)))%>% 
  
  mutate(t= NULL, Date= NULL, Time= NULL, Survived_day=NULL, Status=NULL)%>%
  
  # to ran the analysis with the dabestR package is necessary to create a new variable to group samples (the function doesn't have a grouping option)
  mutate(Day_Trtm_Prec = 
           case_when(
             (Species_ID == "Amu" & Day == "0" & Treatment == "Heat" & Prec == "Ambient")~ "Amu_AT_C", # H and C at the end of the ID indicate the "post-preconditioning" (C) and "Post-heat" (H) time points 
             (Species_ID == "Amu" & Day == "0" & Treatment == "Heat" & Prec == "ST")~ "Amu_ST_C",
             (Species_ID == "Amu" & Day == "0" & Treatment == "Heat" & Prec == "VT")~ "Amu_VT_C",
             (Species_ID == "Amu" & Day == "1" & Treatment == "Heat" & Prec == "Ambient")~ "Amu_AT_H",
             (Species_ID == "Amu" & Day == "1" & Treatment == "Heat" & Prec == "ST")~ "Amu_ST_H",
             (Species_ID == "Amu" & Day == "1" & Treatment == "Heat" & Prec == "VT")~ "Amu_VT_H",
             
         
             (Species_ID == "Gfa" & Day == "0" & Treatment == "Heat" & Prec == "Ambient")~ "Gfa_AT_C",
             (Species_ID == "Gfa" & Day == "0" & Treatment == "Heat" & Prec == "ST")~ "Gfa_ST_C",
             (Species_ID == "Gfa" & Day == "0" & Treatment == "Heat" & Prec == "VT")~ "Gfa_VT_C",
             (Species_ID == "Gfa" & Day == "3" & Treatment == "Heat" & Prec == "Ambient")~ "Gfa_AT_H",
             (Species_ID == "Gfa" & Day == "3" & Treatment == "Heat" & Prec == "ST")~ "Gfa_ST_H",
             (Species_ID == "Gfa" & Day == "3" & Treatment == "Heat" & Prec == "VT")~ "Gfa_VT_H",
             
      
             (Species_ID == "Mdi" & Day == "0" & Treatment == "Heat" & Prec == "Ambient")~ "Mdi_AT_C",
             (Species_ID == "Mdi" & Day == "0" & Treatment == "Heat" & Prec == "ST")~ "Mdi_ST_C",
             (Species_ID == "Mdi" & Day == "0" & Treatment == "Heat" & Prec == "VT")~ "Mdi_VT_C",
             (Species_ID == "Mdi" & Day == "3" & Treatment == "Heat" & Prec == "Ambient")~ "Mdi_AT_H",
             (Species_ID == "Mdi" & Day == "3" & Treatment == "Heat" & Prec == "ST")~ "Mdi_ST_H",
             (Species_ID == "Mdi" & Day == "3" & Treatment == "Heat" & Prec == "VT")~ "Mdi_VT_H",
             
    
             (Species_ID == "Pve" & Day == "0" & Treatment == "Heat" & Prec == "Ambient")~ "Pve_AT_C",
             (Species_ID == "Pve" & Day == "0" & Treatment == "Heat" & Prec == "ST")~ "Pve_ST_C",
             (Species_ID == "Pve" & Day == "0" & Treatment == "Heat" & Prec == "VT")~ "Pve_VT_C",
             (Species_ID == "Pve" & Day == "1" & Treatment == "Heat" & Prec == "Ambient")~ "Pve_AT_H",
             (Species_ID == "Pve" & Day == "1" & Treatment == "Heat" & Prec == "ST")~ "Pve_ST_H",
             (Species_ID == "Pve" & Day == "1" & Treatment == "Heat" & Prec == "VT")~ "Pve_VT_H",
             

             (Species_ID == "Pru" & Day == "0" & Treatment == "Heat" & Prec == "Ambient")~ "Pru_AT_C",
             (Species_ID == "Pru" & Day == "0" & Treatment == "Heat" & Prec == "ST")~ "Pru_ST_C",
             (Species_ID == "Pru" & Day == "0" & Treatment == "Heat" & Prec == "VT")~ "Pru_VT_C",
             (Species_ID == "Pru" & Day == "2" & Treatment == "Heat" & Prec == "Ambient")~ "Pru_AT_H",
             (Species_ID == "Pru" & Day == "2" & Treatment == "Heat" & Prec == "ST")~ "Pru_ST_H",
             (Species_ID == "Pru" & Day == "2" & Treatment == "Heat" & Prec == "VT")~ "Pru_VT_H",
             

             (Species_ID == "Spi" & Day == "0" & Treatment == "Heat" & Prec == "Ambient")~ "Spi_AT_C",
             (Species_ID == "Spi" & Day == "0" & Treatment == "Heat" & Prec == "ST")~ "Spi_ST_C",
             (Species_ID == "Spi" & Day == "0" & Treatment == "Heat" & Prec == "VT")~ "Spi_VT_C",
             (Species_ID == "Spi" & Day == "1" & Treatment == "Heat" & Prec == "Ambient")~ "Spi_AT_H",
             (Species_ID == "Spi" & Day == "1" & Treatment == "Heat" & Prec == "ST")~ "Spi_ST_H",
             (Species_ID == "Spi" & Day == "1" & Treatment == "Heat" & Prec == "VT")~ "Spi_VT_H"))





## Loading model ----
# Create a list of contrasts
# we contrast the "Post-preconditioning"  with the "Post-heat" time points within each preconditioning regime, across each species
Heat.groups <- list(
  c("Gfa_AT_C","Gfa_AT_H"),
  c("Gfa_ST_C","Gfa_ST_H"),
  c("Gfa_VT_C","Gfa_VT_H")
  ,
  c("Pru_AT_C","Pru_AT_H"),
  c("Pru_ST_C","Pru_ST_H"),
  c("Pru_VT_C","Pru_VT_H")
  ,
  
  c("Amu_AT_C","Amu_AT_H"),
  c("Amu_ST_C","Amu_ST_H"),
  c("Amu_VT_C","Amu_VT_H")
  ,
  
  c("Mdi_AT_C","Mdi_AT_H"),
  c("Mdi_ST_C","Mdi_ST_H"),
  c("Mdi_VT_C","Mdi_VT_H")
  ,
  
  c("Pve_AT_C","Pve_AT_H"),
  c("Pve_ST_C","Pve_ST_H"),
  c("Pve_VT_C","Pve_VT_H")
  ,
  c("Spi_AT_C","Spi_AT_H"),
  c("Spi_ST_C","Spi_ST_H"),
  c("Spi_VT_C","Spi_VT_H")
)


Load.Heat_Pam <- load(Heat.effect.H_rcv_mean, Day_Trtm_Prec, YII_mean, 
                                  idx= Heat.groups, 
                                  colour = Prec)  



## Calculate the effect size ----


Effect.Heat_Pam <- hedges_g(Load.Heat_Pam)

## Plot ----
dabest_plot(Effect.Heat_Pam,
            raw_marker_spread=0.5,
            swarm_label = "ΔF/Fm’", contrast_label = "Hedges g",
            contrast_ylim =c(-6,1))


# Save the plot
ggsave(here ("Output", "Heat_respone_PAM.png"), dpi = 400, units = "cm", width = 45, height = 20)


## Results extraction ----

# The exact effect size values are listed under "difference"
Heat_Pam_eff_g<- Effect.Heat_Pam[["boot_result"]]

write_excel_csv(Heat_Pam_eff_g, here("Output","Statistics/Heat_Pam_eff_g.csv"))






# 2.2 Effect  of HEAT on  thermal tolerance -------------------------------------------------------------------------------------------

## Tissue color (Clr) ----

## Grouping ----

Heat.effect.Clr<- Clr%>%
  
  # to ran the analysis with the dabestR package is necessary to create a new variable to group samples (the function doesn't have a grouping option)
  mutate(Day_Trtm_Prec = 
           case_when(
             (ID == "Amu" & Day == "0" & Trtm == "Heat" & Prec == "Ambient")~ "Amu_AT_C", # H and C at the end of the ID indicate the "post-preconditioning" (C) and "Post-heat" (H) time points 
             (ID == "Amu" & Day == "0" & Trtm == "Heat" & Prec == "ST")~ "Amu_ST_C",
             (ID == "Amu" & Day == "0" & Trtm == "Heat" & Prec == "VT")~ "Amu_VT_C",
             (ID == "Amu" & Day == "1" & Trtm == "Heat" & Prec == "Ambient")~ "Amu_AT_H",
             (ID == "Amu" & Day == "1" & Trtm == "Heat" & Prec == "ST")~ "Amu_ST_H",
             (ID == "Amu" & Day == "1" & Trtm == "Heat" & Prec == "VT")~ "Amu_VT_H",
             
             
             (ID == "Gfa" & Day == "0" & Trtm == "Heat" & Prec == "Ambient")~ "Gfa_AT_C",
             (ID == "Gfa" & Day == "0" & Trtm == "Heat" & Prec == "ST")~ "Gfa_ST_C",
             (ID == "Gfa" & Day == "0" & Trtm == "Heat" & Prec == "VT")~ "Gfa_VT_C",
             (ID == "Gfa" & Day == "3" & Trtm == "Heat" & Prec == "Ambient")~ "Gfa_AT_H",
             (ID == "Gfa" & Day == "3" & Trtm == "Heat" & Prec == "ST")~ "Gfa_ST_H",
             (ID == "Gfa" & Day == "3" & Trtm == "Heat" & Prec == "VT")~ "Gfa_VT_H",
             
             
             (ID == "Mdi" & Day == "0" & Trtm == "Heat" & Prec == "Ambient")~ "Mdi_AT_C",
             (ID == "Mdi" & Day == "0" & Trtm == "Heat" & Prec == "ST")~ "Mdi_ST_C",
             (ID == "Mdi" & Day == "0" & Trtm == "Heat" & Prec == "VT")~ "Mdi_VT_C",
             (ID == "Mdi" & Day == "3" & Trtm == "Heat" & Prec == "Ambient")~ "Mdi_AT_H",
             (ID == "Mdi" & Day == "3" & Trtm == "Heat" & Prec == "ST")~ "Mdi_ST_H",
             (ID == "Mdi" & Day == "3" & Trtm == "Heat" & Prec == "VT")~ "Mdi_VT_H",
             
             
             (ID == "Pve" & Day == "0" & Trtm == "Heat" & Prec == "Ambient")~ "Pve_AT_C",
             (ID == "Pve" & Day == "0" & Trtm == "Heat" & Prec == "ST")~ "Pve_ST_C",
             (ID == "Pve" & Day == "0" & Trtm == "Heat" & Prec == "VT")~ "Pve_VT_C",
             (ID == "Pve" & Day == "1" & Trtm == "Heat" & Prec == "Ambient")~ "Pve_AT_H",
             (ID == "Pve" & Day == "1" & Trtm == "Heat" & Prec == "ST")~ "Pve_ST_H",
             (ID == "Pve" & Day == "1" & Trtm == "Heat" & Prec == "VT")~ "Pve_VT_H",
             
             
             (ID == "Pru" & Day == "0" & Trtm == "Heat" & Prec == "Ambient")~ "Pru_AT_C",
             (ID == "Pru" & Day == "0" & Trtm == "Heat" & Prec == "ST")~ "Pru_ST_C",
             (ID == "Pru" & Day == "0" & Trtm == "Heat" & Prec == "VT")~ "Pru_VT_C",
             (ID == "Pru" & Day == "2" & Trtm == "Heat" & Prec == "Ambient")~ "Pru_AT_H",
             (ID == "Pru" & Day == "2" & Trtm == "Heat" & Prec == "ST")~ "Pru_ST_H",
             (ID == "Pru" & Day == "2" & Trtm == "Heat" & Prec == "VT")~ "Pru_VT_H",
             
             
             (ID == "Spi" & Day == "0" & Trtm == "Heat" & Prec == "Ambient")~ "Spi_AT_C",
             (ID == "Spi" & Day == "0" & Trtm == "Heat" & Prec == "ST")~ "Spi_ST_C",
             (ID == "Spi" & Day == "0" & Trtm == "Heat" & Prec == "VT")~ "Spi_VT_C",
             (ID == "Spi" & Day == "1" & Trtm == "Heat" & Prec == "Ambient")~ "Spi_AT_H",
             (ID == "Spi" & Day == "1" & Trtm == "Heat" & Prec == "ST")~ "Spi_ST_H",
             (ID == "Spi" & Day == "1" & Trtm == "Heat" & Prec == "VT")~ "Spi_VT_H"))




## Loading model ----
# Create a list of contrasts
# we contrast the "Post-preconditioning"  with the "Post-heat" time points within each preconditioning regime, across each species
Heat.groups <- list(
  c("Gfa_AT_C","Gfa_AT_H"),
  c("Gfa_ST_C","Gfa_ST_H"),
  c("Gfa_VT_C","Gfa_VT_H")
  ,
  c("Pru_AT_C","Pru_AT_H"),
  c("Pru_ST_C","Pru_ST_H"),
  c("Pru_VT_C","Pru_VT_H")
  ,
  
  c("Amu_AT_C","Amu_AT_H"),
  c("Amu_ST_C","Amu_ST_H"),
  c("Amu_VT_C","Amu_VT_H")
  ,
  
  c("Mdi_AT_C","Mdi_AT_H"),
  c("Mdi_ST_C","Mdi_ST_H"),
  c("Mdi_VT_C","Mdi_VT_H")
  ,
  
  c("Pve_AT_C","Pve_AT_H"),
  c("Pve_ST_C","Pve_ST_H"),
  c("Pve_VT_C","Pve_VT_H")
  ,
  c("Spi_AT_C","Spi_AT_H"),
  c("Spi_ST_C","Spi_ST_H"),
  c("Spi_VT_C","Spi_VT_H")
)


Load.Heat_Clr <- load(Heat.effect.Clr, Day_Trtm_Prec, gray_sum, 
                      idx= Heat.groups, 
                      colour = Prec)  



## Calculate the effect size ----


Effect.Heat_Clr <- hedges_g(Load.Heat_Clr)

## Plot ----
dabest_plot(Effect.Heat_Clr,
            raw_marker_spread=0.5,
            swarm_label = "ΔF/Fm’", contrast_label = "Hedges g",
            contrast_ylim =c(-14,1)
            )


# Save the plot
ggsave(here ("Output", "Heat_respone_Clr.png"), dpi = 400, units = "cm", width = 45, height = 20)


## Results extraction ----

# The exact effect size values are listed under "difference"
Heat_Clr_eff_g<- Effect.Heat_Clr[["boot_result"]]

write_excel_csv(Heat_Clr_eff_g, here("Output","Statistics/Heat_Clr_eff_g.csv"))







# 3.1 Effect on Recovery -------------------------------------------------------------------------------------------



## Photosynthetic efficiency (PAM) -----

## Grouping ----

Recovery.effect.H_rcv_mean <- H_rcv_mean%>%
  filter(Day== 30)%>%
  mutate(YII_mean = as.numeric(YII_mean+ round(runif(length(YII_mean), 0.00001, 0.00003),7)))%>% 
  mutate(Day_Treatment_Prec = 
           case_when(
             
             (Species_ID == "Amu" & Day == "30" & Treatment == "Control" & Prec == "Ambient")~ "Amu_CT_C",
             (Species_ID == "Amu" & Day == "30" & Treatment == "Control" & Prec == "ST")~ "Amu_ST_C",
             (Species_ID == "Amu" & Day == "30" & Treatment == "Control" & Prec == "VT")~ "Amu_VT_C",
             (Species_ID == "Amu" & Day == "30" & Treatment == "Heat" & Prec == "Ambient")~ "Amu_CT_H",
             (Species_ID == "Amu" & Day == "30" & Treatment == "Heat" & Prec == "ST")~ "Amu_ST_H",
             (Species_ID == "Amu" & Day == "30" & Treatment == "Heat" & Prec == "VT")~ "Amu_VT_H",
             
             
             (Species_ID == "Gfa" & Day == "30" & Treatment == "Control" & Prec == "Ambient")~ "Gfa_CT_C",
             (Species_ID == "Gfa" & Day == "30" & Treatment == "Control" & Prec == "ST")~ "Gfa_ST_C",
             (Species_ID == "Gfa" & Day == "30" & Treatment == "Control" & Prec == "VT")~ "Gfa_VT_C",
             (Species_ID == "Gfa" & Day == "30" & Treatment == "Heat" & Prec == "Ambient")~ "Gfa_CT_H",
             (Species_ID == "Gfa" & Day == "30" & Treatment == "Heat" & Prec == "ST")~ "Gfa_ST_H",
             (Species_ID == "Gfa" & Day == "30" & Treatment == "Heat" & Prec == "VT")~ "Gfa_VT_H",
             
             (Species_ID == "Mdi" & Day == "30" & Treatment == "Control" & Prec == "Ambient")~ "Mdi_CT_C",
             (Species_ID == "Mdi" & Day == "30" & Treatment == "Control" & Prec == "ST")~ "Mdi_ST_C",
             (Species_ID == "Mdi" & Day == "30" & Treatment == "Control" & Prec == "VT")~ "Mdi_VT_C",
             (Species_ID == "Mdi" & Day == "30" & Treatment == "Heat" & Prec == "Ambient")~ "Mdi_CT_H",
             (Species_ID == "Mdi" & Day == "30" & Treatment == "Heat" & Prec == "ST")~ "Mdi_ST_H",
             (Species_ID == "Mdi" & Day == "30" & Treatment == "Heat" & Prec == "VT")~ "Mdi_VT_H",
             
             
             (Species_ID == "Pve" & Day == "30" & Treatment == "Control" & Prec == "Ambient")~ "Pve_CT_C",
             (Species_ID == "Pve" & Day == "30" & Treatment == "Control" & Prec == "ST")~ "Pve_ST_C",
             (Species_ID == "Pve" & Day == "30" & Treatment == "Control" & Prec == "VT")~ "Pve_VT_C",
             (Species_ID == "Pve" & Day == "30" & Treatment == "Heat" & Prec == "Ambient")~ "Pve_CT_H",
             (Species_ID == "Pve" & Day == "30" & Treatment == "Heat" & Prec == "ST")~ "Pve_ST_H",
             (Species_ID == "Pve" & Day == "30" & Treatment == "Heat" & Prec == "VT")~ "Pve_VT_H",
             
             
             (Species_ID == "Pru" & Day == "30" & Treatment == "Control" & Prec == "Ambient")~ "Pru_CT_C",
             (Species_ID == "Pru" & Day == "30" & Treatment == "Control" & Prec == "ST")~ "Pru_ST_C",
             (Species_ID == "Pru" & Day == "30" & Treatment == "Control" & Prec == "VT")~ "Pru_VT_C",
             (Species_ID == "Pru" & Day == "30" & Treatment == "Heat" & Prec == "Ambient")~ "Pru_CT_H",
             (Species_ID == "Pru" & Day == "30" & Treatment == "Heat" & Prec == "ST")~ "Pru_ST_H",
             (Species_ID == "Pru" & Day == "30" & Treatment == "Heat" & Prec == "VT")~ "Pru_VT_H",
             
             
             (Species_ID == "Spi" & Day == "30" & Treatment == "Control" & Prec == "Ambient")~ "Spi_CT_C",
             (Species_ID == "Spi" & Day == "30" & Treatment == "Control" & Prec == "ST")~ "Spi_ST_C",
             (Species_ID == "Spi" & Day == "30" & Treatment == "Control" & Prec == "VT")~ "Spi_VT_C",
             (Species_ID == "Spi" & Day == "30" & Treatment == "Heat" & Prec == "Ambient")~ "Spi_CT_H",
             (Species_ID == "Spi" & Day == "30" & Treatment == "Heat" & Prec == "ST")~ "Spi_ST_H",
             (Species_ID == "Spi" & Day == "30" & Treatment == "Heat" & Prec == "VT")~ "Spi_VT_H"))




## Loading model ----
# Create a list of contrasts
# we contrast the "Heat" with the "Control" treatments at "day-30" time point within each preconditioning regime, across each species


Recovery.groups <- list(c("Gfa_CT_C","Gfa_CT_H"),
                        c("Gfa_ST_C","Gfa_ST_H"),
                        c("Gfa_VT_C","Gfa_VT_H")
                        ,
                        c("Pru_CT_C","Pru_CT_H"),
                        c("Pru_ST_C","Pru_ST_H"),
                        c("Pru_VT_C","Pru_VT_H")
                        ,
                        c("Mdi_CT_C","Mdi_CT_H"),
                        c("Mdi_ST_C","Mdi_ST_H"),
                        c("Mdi_VT_C","Mdi_VT_H"))
                        
                        
  



Load.recovery_Pam <- load(Recovery.effect.H_rcv_mean, Day_Treatment_Prec, YII_mean, 
                                   idx= Recovery.groups,
                                   colour = Prec )

## Calculate the effect size ----


Effect.recovery_Pam <- hedges_g(Load.recovery_Pam)

## Plot ----
dabest_plot(Effect.recovery_Pam,
            raw_marker_spread=0.5,
            swarm_label = "ΔF/Fm’", contrast_label = "Hedges g",
            contrast_ylim =c(-14,1)
)


# Save the plot
ggsave(here ("Output", "Heat_respone_Clr.png"), dpi = 400, units = "cm", width = 45, height = 20)



## Results extraction ----

# The exact effect size values are listed under "difference"
# The p value showed in table S6 are obtained from the following df

Recovery_Pam_g<- Effect.recovery_Pam[["boot_result"]]

write_excel_csv(Recovery_Pam_g, here("Output","Statistics/Recovery_Pam_g.csv"))

## Results extraction Table S6----

# The p value showed in table S6 are obtained from the following df

Recovery_Pam_pvalues_g <- Effect.recovery_Pam[["permtest_pvals"]]



# 3.2 Effect on Recovery -------------------------------------------------------------------------------------------

## Tissue color (Clr) -----

## Grouping ----

Recovery.effect.Clr <- Clr%>%
  filter(Day== 30)%>%
  mutate(Day_Treatment_Prec = 
           case_when(
             
             (ID == "Amu" & Day == "30" & Trtm == "Control" & Prec == "Ambient")~ "Amu_CT_C",
             (ID == "Amu" & Day == "30" & Trtm == "Control" & Prec == "ST")~ "Amu_ST_C",
             (ID == "Amu" & Day == "30" & Trtm == "Control" & Prec == "VT")~ "Amu_VT_C",
             (ID == "Amu" & Day == "30" & Trtm == "Heat" & Prec == "Ambient")~ "Amu_CT_H",
             (ID == "Amu" & Day == "30" & Trtm == "Heat" & Prec == "ST")~ "Amu_ST_H",
             (ID == "Amu" & Day == "30" & Trtm == "Heat" & Prec == "VT")~ "Amu_VT_H",
             
             
             (ID == "Gfa" & Day == "30" & Trtm == "Control" & Prec == "Ambient")~ "Gfa_CT_C",
             (ID == "Gfa" & Day == "30" & Trtm == "Control" & Prec == "ST")~ "Gfa_ST_C",
             (ID == "Gfa" & Day == "30" & Trtm == "Control" & Prec == "VT")~ "Gfa_VT_C",
             (ID == "Gfa" & Day == "30" & Trtm == "Heat" & Prec == "Ambient")~ "Gfa_CT_H",
             (ID == "Gfa" & Day == "30" & Trtm == "Heat" & Prec == "ST")~ "Gfa_ST_H",
             (ID == "Gfa" & Day == "30" & Trtm == "Heat" & Prec == "VT")~ "Gfa_VT_H",
             
             (ID == "Mdi" & Day == "30" & Trtm == "Control" & Prec == "Ambient")~ "Mdi_CT_C",
             (ID == "Mdi" & Day == "30" & Trtm == "Control" & Prec == "ST")~ "Mdi_ST_C",
             (ID == "Mdi" & Day == "30" & Trtm == "Control" & Prec == "VT")~ "Mdi_VT_C",
             (ID == "Mdi" & Day == "30" & Trtm == "Heat" & Prec == "Ambient")~ "Mdi_CT_H",
             (ID == "Mdi" & Day == "30" & Trtm == "Heat" & Prec == "ST")~ "Mdi_ST_H",
             (ID == "Mdi" & Day == "30" & Trtm == "Heat" & Prec == "VT")~ "Mdi_VT_H",
             
             
             (ID == "Pve" & Day == "30" & Trtm == "Control" & Prec == "Ambient")~ "Pve_CT_C",
             (ID == "Pve" & Day == "30" & Trtm == "Control" & Prec == "ST")~ "Pve_ST_C",
             (ID == "Pve" & Day == "30" & Trtm == "Control" & Prec == "VT")~ "Pve_VT_C",
             (ID == "Pve" & Day == "30" & Trtm == "Heat" & Prec == "Ambient")~ "Pve_CT_H",
             (ID == "Pve" & Day == "30" & Trtm == "Heat" & Prec == "ST")~ "Pve_ST_H",
             (ID == "Pve" & Day == "30" & Trtm == "Heat" & Prec == "VT")~ "Pve_VT_H",
             
             
             (ID == "Pru" & Day == "30" & Trtm == "Control" & Prec == "Ambient")~ "Pru_CT_C",
             (ID == "Pru" & Day == "30" & Trtm == "Control" & Prec == "ST")~ "Pru_ST_C",
             (ID == "Pru" & Day == "30" & Trtm == "Control" & Prec == "VT")~ "Pru_VT_C",
             (ID == "Pru" & Day == "30" & Trtm == "Heat" & Prec == "Ambient")~ "Pru_CT_H",
             (ID == "Pru" & Day == "30" & Trtm == "Heat" & Prec == "ST")~ "Pru_ST_H",
             (ID == "Pru" & Day == "30" & Trtm == "Heat" & Prec == "VT")~ "Pru_VT_H",
             
             
             (ID == "Spi" & Day == "30" & Trtm == "Control" & Prec == "Ambient")~ "Spi_CT_C",
             (ID == "Spi" & Day == "30" & Trtm == "Control" & Prec == "ST")~ "Spi_ST_C",
             (ID == "Spi" & Day == "30" & Trtm == "Control" & Prec == "VT")~ "Spi_VT_C",
             (ID == "Spi" & Day == "30" & Trtm == "Heat" & Prec == "Ambient")~ "Spi_CT_H",
             (ID == "Spi" & Day == "30" & Trtm == "Heat" & Prec == "ST")~ "Spi_ST_H",
             (ID == "Spi" & Day == "30" & Trtm == "Heat" & Prec == "VT")~ "Spi_VT_H"))




## Loading model ----
# Create a list of contrasts
# we contrast the "Heat" with the "Control" treatments at "day-30" time point within each preconditioning regime, across each species


Recovery.groups <- list(c("Gfa_CT_C","Gfa_CT_H"),
                        c("Gfa_ST_C","Gfa_ST_H"),
                        c("Gfa_VT_C","Gfa_VT_H")
                        ,
                        c("Pru_CT_C","Pru_CT_H"),
                        c("Pru_ST_C","Pru_ST_H"),
                        c("Pru_VT_C","Pru_VT_H")
                        ,
                        c("Mdi_CT_C","Mdi_CT_H"),
                        c("Mdi_ST_C","Mdi_ST_H"),
                        c("Mdi_VT_C","Mdi_VT_H"))






Load.recovery_Clr <- load(Recovery.effect.Clr, Day_Treatment_Prec, gray_sum, 
                          idx= Recovery.groups,
                          colour = Prec )

## Calculate the effect size ----


Effect.recovery_Clr <- hedges_g(Load.recovery_Clr)

## Plot ----
dabest_plot(Effect.recovery_Clr,
            raw_marker_spread=0.5,
            swarm_label = "ΔF/Fm’", contrast_label = "Hedges g",
            contrast_ylim =c(-14,1)
)


# Save the plot
ggsave(here ("Output", "Recovery_Clr.png"), dpi = 400, units = "cm", width = 45, height = 20)



## Results extraction ----

# The exact effect size values are listed under "difference"
Recovery_Clr_g<- Effect.recovery_Clr[["boot_result"]]

write_excel_csv(Recovery_Clr_g, here("Output","Statistics/Recovery_Clr_g.csv"))


## Results extraction Table S6----

# The p value showed in table S6 are obtained from the following df

Recovery_Clr_pvalues_g <- Effect.recovery_Clr[["permtest_pvals"]]


