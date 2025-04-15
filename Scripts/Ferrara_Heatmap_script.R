
# --------------------------- Ferrara et al. 2024 - R script N. 5 --------------------------- #
#                                                                                             #
#                                                                                             #
#                Heatmap the resume the overall coral species response                        #
#                                                                                             #
#  In this script we create a heatmap using the effect size of post-heat coral responses,     #
#  the recovery rate (measured as the difference, showed as effect size, between the Heat     #
#  treatment and the respective Control treatment within each preconditioning regime), and    #
#  survival rate                                                                              #
#                                                                                             #
#                                                                                             #
#-------------------------------------------------------------------------------------------- #

# Libraries ----

library(here)
library(ggpubr)
library(gridExtra)
library(ggtext)
library(dplyr)   
library(tidyverse)
library(ggnewscale)
library(RColorBrewer)

library(here)
library(lubridate)
library(scales)
library(ggpubr)

library(rstatix)
library(patchwork)
library(multcomp)

library(easystats)
library(gridExtra)
library(lme4)

library(ggbeeswarm)

library(ggtext)

library(dplyr)   
library(tidyverse)
library(emmeans)
library(ggnewscale)
library(survival)
library(survminer)
library(car)





# 1 (PAM) Post-heat Photosynthetic efficiency effect sizes -----------------------------------------------------------------------------------------------


# Heat_Pam_eff_g the effect size table obtained in script 4 about the coral response (photosynthetic efficiency) to heat is modified to create the heatmap

Heat_Pam_eff_heatmap <-Heat_Pam_eff_g%>%
  
  separate(col=test_group, into=c('Species_ID', 'Prec',"Treatment"), sep='_', remove = FALSE)%>%
  mutate(Prec = if_else(Prec == "AT", "Ambient", Prec),
         Treatment =if_else(Treatment == "H", "Heat", Treatment))%>%
  
  mutate(Analysis = "PAM", Value = "PAM")%>%
  
  mutate(Species =case_when(str_starts(Species_ID, "Pru") ~ "P. rus",
                            str_starts(Species_ID, "Gfa") ~ "G. fascicularis",
                            str_starts(Species_ID, "Mdi") ~ "M. digitata",
                            str_starts(Species_ID, "Spi") ~ "S. pistillata",
                            str_starts(Species_ID, "Pve") ~ "P. verrucosa",
                            str_starts(Species_ID, "Amu") ~ "A. muricata"),
         
         Species =  factor(Species, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata", 
                                              "P. verrucosa", "S. pistillata")))%>%

  dplyr::select(Species_ID, Species, Prec, Treatment, Analysis, difference, Value)%>%
  rename(eff = difference)


# 2 (Clr) Post-heat tissue color effect sizes --------------------------------------------------------------------------------------------------------------

# Heat_Clr_eff_g the effect size table obtained in script 4 about the coral response (tissue color) to heat is modified to create the heatmap

Heat_Clr_eff_heatmap <-Heat_Clr_eff_g %>%
  
  separate(col=test_group, into=c('Species_ID', 'Prec',"Treatment"), sep='_', remove = FALSE)%>%
  mutate(Prec = if_else(Prec == "AT", "Ambient", Prec),
         
         Treatment =if_else(Treatment == "H", "Heat", Treatment))%>%
  
  mutate(Analysis = "Clr", Value = "Clr")%>%
  
  mutate(Species =case_when(str_starts(Species_ID, "Pru") ~ "P. rus",
                            str_starts(Species_ID, "Gfa") ~ "G. fascicularis",
                            str_starts(Species_ID, "Mdi") ~ "M. digitata",
                            str_starts(Species_ID, "Spi") ~ "S. pistillata",
                            str_starts(Species_ID, "Pve") ~ "P. verrucosa",
                            str_starts(Species_ID, "Amu") ~ "A. muricata"),
         
         Species =  factor(Species, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata", 
                                              "P. verrucosa", "S. pistillata")))%>%
  
  dplyr::select(Species_ID, Species, Prec, Treatment, Analysis,  difference, Value)%>%

  rename(eff = difference)%>%
  
  view()


# 3 (PAM-Recovery) Photosynthetic efficiency effect sizes during recovery phase  ---------------------------------------------------------------------------

#  Recovery_Pam_g the effect size table obtained in script 4 about the coral recovery rate (photosynthetic efficiency) is modified to create the heatmap

Recovery_Pam_eff_heatmap <- Recovery_Pam_g %>% 
  
  separate(col=test_group, into=c('Species_ID', 'Prec',"Treatment"), sep='_', remove = FALSE)%>%
  mutate(Prec = if_else(Prec == "CT", "Ambient", Prec),
         
         Treatment =if_else(Treatment == "H", "Heat", Treatment))%>%
  
  mutate(Analysis = "PAM_Recovery", Value = "PAM_Recovery")%>%
  
 
  
  mutate(Species =case_when(str_starts(Species_ID, "Pru") ~ "P. rus",
                            str_starts(Species_ID, "Gfa") ~ "G. fascicularis",
                            str_starts(Species_ID, "Mdi") ~ "M. digitata",
                            str_starts(Species_ID, "Spi") ~ "S. pistillata",
                            str_starts(Species_ID, "Pve") ~ "P. verrucosa",
                            str_starts(Species_ID, "Amu") ~ "A. muricata"),
         
         Species =  factor(Species, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata", 
                                              "P. verrucosa", "S. pistillata"))
         )%>%
    
  rename(eff = difference) %>%
  
  dplyr::select(Species_ID, Species, Prec, Treatment, Analysis, eff, Value)

## Missing values ----

# Here we add manually arbitrary values for the species that did not survived till day 30, to have the same table structure.

PAM_miss <- data.frame(Species_ID  = c("Amu","Amu", "Amu", "Spi","Spi","Spi","Pve","Pve","Pve"),
                           Species     = c("A. muricata","A. muricata","A. muricata", "S. pistillata", "S. pistillata", "S. pistillata","P. verrucosa","P. verrucosa","P. verrucosa"),
                           Prec        = c("Ambient", "ST", "VT","Ambient", "ST", "VT","Ambient", "ST", "VT"),
                           Treatment   = c("Heat","Heat","Heat","Heat","Heat","Heat","Heat","Heat","Heat"),
                           Analysis    = c("PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery"),
                           eff         = c(-25,-25,-25,-25,-25,-25,-25,-25,-25),
                           Value       = c("PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery","PAM_Recovery"))

# Merge the two df
Recovery_Pam_eff_heatmap <- bind_rows(Recovery_Pam_eff_heatmap,PAM_miss)






# 4 (Clr-Recovery) Tissue color effect sizes during recovery phase  ----


# Recovery_Clr_g the effect size table obtained in script 4 about the coral recovery rate (tissue color) is modified to create the heatmap

Recovery_Clr_eff_heatmap <-Recovery_Clr_g %>% 
  
  separate(col=test_group, into=c('Species_ID', 'Prec',"Treatment"), sep='_', remove = FALSE)%>%
  mutate(Prec = if_else(Prec == "CT", "Ambient", Prec),
         
         Treatment =if_else(Treatment == "H", "Heat", Treatment))%>%
  mutate(Analysis = "Clr_Recovery", Value = "Clr_Recovery")%>%
  
  mutate(Species =case_when(str_starts(Species_ID, "Pru") ~ "P. rus",
                            str_starts(Species_ID, "Gfa") ~ "G. fascicularis",
                            str_starts(Species_ID, "Mdi") ~ "M. digitata",
                            str_starts(Species_ID, "Spi") ~ "S. pistillata",
                            str_starts(Species_ID, "Pve") ~ "P. verrucosa",
                            str_starts(Species_ID, "Amu") ~ "A. muricata"),
         Species =  factor(Species, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata", 
                                              "P. verrucosa", "S. pistillata")))%>%
  
  dplyr::select(Species_ID, Species, Prec, Treatment, Analysis, difference, Value)%>%
  
  rename(eff = difference)


## Missing values ----

# Here we add manually arbitrary values for the species that did not survived till day 30, to have the same table structure.

Clr_miss <- data.frame(Species_ID  = c("Amu","Amu", "Amu", "Spi","Spi","Spi","Pve","Pve","Pve"),
                           Species     = c("A. muricata","A. muricata","A. muricata", "S. pistillata", "S. pistillata", "S. pistillata","P. verrucosa","P. verrucosa","P. verrucosa"),
                           Prec        = c("Ambient", "ST", "VT","Ambient", "ST", "VT","Ambient", "ST", "VT"),
                           Treatment   = c("Heat","Heat","Heat","Heat","Heat","Heat","Heat","Heat","Heat"),
                           Analysis    = c("Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery"),
                           eff         = c(-25,-25,-25,-25,-25,-25,-25,-25,-25),
                           Value       = c("Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery","Clr_Recovery"))

Recovery_Clr_eff_heatmap <- bind_rows(Recovery_Clr_eff_heatmap,Clr_miss) 







# 5 Survival rate  ----------------------------------------------------------------------------------------------------

Surv_Pru_df <- Surv_Pru_df%>%
  mutate(ID = "Pru", Species ="P.rus")

Surv_Amu_df <- Surv_Amu_df%>%
  mutate(ID = "Amu", Species ="A.muricata")

Surv_Mdi_df <- Surv_Mdi_df%>%
  mutate(ID = "Mdi", Species ="M. digitata")

Surv_Spi_df <- Surv_Spi_df%>%
  mutate(ID = "Spi", Species ="S.pistillata")

Surv_Gfa_df <- Surv_Gfa_df%>%
  mutate(ID = "Gfa", Species ="G. fascicularis")

Surv_Pve_df <- Surv_Pve_df%>%
  mutate(ID = "Pve", Species ="P. verrucosa")


Survival_heatmap.df <- bind_rows(Surv_Pru_df,
                           Surv_Amu_df,
                           Surv_Pve_df,
                           Surv_Gfa_df,
                           Surv_Spi_df,
                           Surv_Mdi_df) 


Survival_heatmap.df.clean <-  Survival_heatmap.df %>%
  
  rename( Species_ID = ID)%>%
  
  group_by(Species_ID, Prec, )%>%
  
  slice_min(surv, n = 1, with_ties = FALSE)%>%
 
  mutate(Analysis = "Survive rate", time =NULL, Treatment = "Heat", Value = "Survive rate")%>%
  
  mutate(Species =case_when(str_starts(Species_ID, "Pru") ~ "P. rus",
                            str_starts(Species_ID, "Gfa") ~ "G. fascicularis",
                            str_starts(Species_ID, "Mdi") ~ "M. digitata",
                            str_starts(Species_ID, "Spi") ~ "S. pistillata",
                            str_starts(Species_ID, "Pve") ~ "P. verrucosa",
                            str_starts(Species_ID, "Amu") ~ "A. muricata"),
         
         Species =  factor(Species, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata", 
                                              "P. verrucosa", "S. pistillata")))%>%
  
  rename(eff = surv)%>%
  mutate(across (c(eff), as.numeric))%>%
  mutate_if(is.numeric, round, 3) %>%
  view()


# 6 Merging df for heatmap ---------------------------------------------------------------------------------

Heatmap.df <- bind_rows(Heat_Pam_eff_heatmap, Heat_Clr_eff_heatmap, Recovery_Pam_eff_heatmap, Recovery_Clr_eff_heatmap, Survival_heatmap.df.clean)

Heatmap.df <-Heatmap.df %>%
  mutate_if(is.numeric, round,3)%>%
  mutate (Species =  factor(Species, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata", 
                                               "P. verrucosa", "S. pistillata")),
          Analysis =  factor(Analysis, level = c("PAM", "Clr", "PAM_Recovery","Clr_Recovery", "Survive rate")))








# 7 HEATMAP-------------------------------------------------------------


# Color palette for PAM and tissue color
colors <- c("white", "#FFEE96FE", "#f3c999", "#da8d56", "#c36f44", "#9b5636", "#7b3e2a", "#5e2a06", "#4b2e16", "black")
colors1 <- c("white", "#d2eaa3", "#b7dc74", "#74a94e", "#507d32", "#2A6337FC", "#145C0E","#064011FE")

# Breakpoints
values <- scales::rescale(c(-10,-9.8, -7, -4.7, -3, -2.5, -2, -1.5, -1,-0.5, -0.3, 0), to = c(0, 1))
values1 <- scales::rescale(c(-15, -8, -7, -6, -5, -4, -3, -2, -1, -0.5, 0), to = c(0, 1))



# The heatmap is build on 5 different layers, one for each of the different set of data

Heatmap.df %>%


  mutate(Value =  factor(Value, level = c("Survive rate","Clr_Recovery", "Clr","PAM_Recovery","PAM" )
                         ))%>%
  # PAM
  ggplot(mapping = aes( x = Prec , y = Value)) +
  geom_tile(data = ~ filter(.x,Heatmap.df$Value == "PAM") %>% rename(`Changes in ΔF/Fm'` = eff), mapping = aes(fill = `Changes in ΔF/Fm'`), color = "lightgray") +

  scale_fill_gradientn(colors = colors1,
                       values = values1,
                       guide = guide_colorbar(order = 1,title.position = "top",barwidth = 20),
                       na.value = "#F6FFE3FE", # To ensure the rest of the values remain readable, the limit was set to -10. The effect size for S. pistillata from the Ambient regime, which was -17, was colored using a light shade color that was not white.
                       limits = c(-10,0)
  ) +
  
  # Tissue color
  new_scale("fill") +
  
  geom_tile(data = ~ filter(.x,Heatmap.df$Value == "Clr") %>% rename(`Changes in tissue color` = eff), mapping = aes(fill = `Changes in tissue color`), color = "lightgray") +
  scale_fill_gradientn(colors = colors,
                       values = values,
                       guide = guide_colorbar(order = 3,title.position = "top",barwidth = 20),
                       limits = c(-10,0),oob = scales::squish) +

  # Pam recovery
  new_scale("fill") +
  
  geom_tile(data = ~ filter(.x,Heatmap.df$Value == "PAM_Recovery") %>% rename(`Change in ΔF/Fm'(Recovery)` = eff),
            
            mapping = aes(fill = `Change in ΔF/Fm'(Recovery)`), color = "lightgray") +
  
  scale_fill_gradientn(colors = colors1,
                       values = values1,
                       guide =  FALSE, 
                       limits = c(-10,0),oob = scales::squish) +
  
  # Tissue color recovery
  new_scale("fill") +
  
  geom_tile(data = ~ filter(.x,Heatmap.df$Value == "Clr_Recovery") %>% rename(`Changes in tissue color (Recovery)` = eff),
            mapping = aes(fill = `Changes in tissue color (Recovery)`), color = "lightgray") +

  scale_fill_gradientn(colors = colors,
                       values = values,
                       guide =  FALSE,
                       limits = c(-10,0),oob = scales::squish) +
  
  new_scale("fill") +
  
  geom_tile(data = ~ filter(.x,Heatmap.df$Value == "Survive rate") %>% rename(`Survive rate (%)` = eff)%>% mutate(`Survive rate (%)`=(100*`Survive rate (%)`)), mapping = aes(fill = `Survive rate (%)`), color = "lightgray") +
  scale_fill_gradientn(colors = rev(c("#0D00FF",  "white")), guide = guide_colorbar(order = 5,title.position = "top",barwidth = 20)) +
  
  scale_y_discrete(label = c("Clr_Recovery" = "Tissue color<br>(recovery)",
                             "PAM_Recovery" = "ΔF/Fm'<br>(recovery)",
                             "Survive rate" = "Surive rate (%)",
                             "Clr" = "Tissue color",
                             "PAM" = "ΔF/Fm'"
                             ),
                   drop = FALSE # this command allow to maintain the level information despite this heatmap is build on separate layers
                   )+
  
  scale_x_discrete(position = "top")+
  
  
  guides(
    
    fill = guide_colourbar(title.position = "top",barwidth = 15))+
  
  theme_recovery_IMG+
  theme(
    plot.margin=unit(c(0.5,1,0,0), "cm"),
    legend.position = "bottom", 
    axis.text.y = element_markdown(angle = 45, vjust = 0),
    legend.justification = c(0.5, 0.5),
    legend.direction = "horizontal",
    legend.title = element_text(size = 22, angle = 0, hjust = 0.5,), 
    legend.text = element_text(size = 18),
    legend.spacing.x = unit(1, "cm")
  )+

  labs(x="",y= "")+
  
  coord_fixed() +
  facet_wrap(~ Species, strip.position ="bottom", ncol = 6)


###ggsave ----

ggsave(here ("Output", "Heatmap_hedges_manuscript.png"), dpi = 400, units = "cm", width = 53, height = 20)












