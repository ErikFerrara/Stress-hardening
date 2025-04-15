

# --------------------------- Ferrara et al. 2024 - R script N. 1 --------------------------- #
#                                                                                             #
#                                                                                             #
#  Effect of preconditioning treatments on photosynthetic efficiency at time points           #
#  "Post-preconditioning" and "Post-heat"                                                     #
#                                                                                             #
#  In this script we calculate the effect of preconditioning on baseline physiology           #
#  and the increase in thermal tolerance after the heat assay                                 #
#                                                                                             #
#  The effects are calculated on photosynthetic efficiency metrics only                       #
#                                                                                             #
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


# Load data ----

# Load the photosynthetic efficiency data (effective quantum yield). This df contain the effective quantum yield values from the "post-preconditioning",
# "post-heat", "day-15", and "day-30" time points. However, DATA ANALYSIS OF "POST-PRECONDITIONING" AND "POST-HEAT" are shown here.
# At each time point, each fragment effective quantum yield was measured in three different spots to obtained a comprehensive overview 
# of coral fragments health status. Before any analysis, the mean of these three measurements was calculated. 


#Load the photosynthetic efficiency data (effective quantum yield and surival data)

H_rcv <- read_csv2(here("Data", "Ferrara_etal_PAM.csv"))


#                            1 Effect of Prec. on baseline physiology -----
############################################################################################################# #


#In this first chapter, the effect of thermal preconditioning regimes on coral baseline effective quantum yield
#is tested using LMER models for each species separately. 


## Preconditioning df subset----

#Subset the main df to include only the "post-preconditioning" values.  Mean values (each fragment
#was measured in three different spots) are calculated. 

PAM_Preconditioning <- H_rcv %>%
  
  mutate(across(c(Species_ID, Prec, Treatment, Tank_n, Colony, Day, Frg_n, Region, Treatment), as.factor),
         across(c(YII, F, Fm), as.numeric))%>%
  
  dplyr::filter(Day == 0)%>%
  
  mutate (Species =  factor(Species, level = c("Galaxea fascicularis", "Porites rus", "Acropora muricata", "Montipora digitata", 
                                               "Pocillopora verrucosa", "Stylophora pistillata"
  )))%>%
  
  group_by(Species_ID, Species ,Prec, Colony, Region, Treatment, Tank_n, Frg_n, Day,)%>%
  
  summarise(F_mean = mean(F),F_sd = sd(F),
            Fm_mean = mean(Fm),Fm_sd = sd(Fm),
            YII_mean = mean(YII),YII_sd = sd(YII)
  )%>%
  mutate_if(is.numeric, round, 3) %>%
  ungroup()

### Shapiro test  with mean data----

#Test the normal distribution within each preconditioning treatment across species
PAM_Preconditioning %>%
  ungroup() %>%
  group_by(Species_ID, Prec) %>%
  rstatix::shapiro_test(YII_mean) %>%
  mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))

### Levene test ----

PAM_Preconditioning %>%
    ungroup() %>%
  group_by(Species_ID) %>%
  rstatix::levene_test(YII_mean ~ Prec) %>%
  mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))









#                          Post-Preconditioning Linear Mixed-Effects Models ----
################################################################################################################### #
# The effect of the two preconditioning treatments is tested for each species separately



## Post-Prec - MONTIPORA DIGITATA (Mdi) ----------------------------------
#### Subset ----


Mdi.PAM_Preconditioning <-
  PAM_Preconditioning %>%
  dplyr::filter(Species_ID == "Mdi")



#### Models ----


Mdi.lmer.Prec <- lmerTest::lmer(YII_mean ~ Prec + (1|Colony)+ (1|Tank_n), data= Mdi.PAM_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Mdi.lmer.Prec))

# Residual normality
check_normality(Mdi.lmer.Prec)

# Check homogeneity
check_homogeneity(Mdi.lmer.Prec)

#Check heteroscedasticity
## Residual plot
plot(residuals(Mdi.lmer.Prec) ~ fitted(Mdi.lmer.Prec), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Mdi.lmer.Prec), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Mdi.lmer.Prec)
fitted_lmer <- fitted(Mdi.lmer.Prec)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Mdi_Prec_PAM_anova<-as_tibble(anova(Mdi.lmer.Prec))
Mdi_Prec_PAM_anova

write.table(Mdi_Prec_PAM_anova, sep = ";", "Statistics/Mdi_Prec_PAM_anova.csv")



#### Post-hoc analysis ---------

Mdi.prec.result <- tidy(glht(Mdi.lmer.Prec, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>%  add_significance("adj.p.value")
Mdi.prec.result

write.table(Mdi.prec.result, sep = ";", "Statistics/Mdi_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Mdi.lmer.Prec)
report(Mdi.lmer.Prec)











## Post-Prec - PORITES RUS (Pru) -------------------------------------------------------------
#### Subset ----


Pru.PAM_Preconditioning <-
  PAM_Preconditioning %>%
  dplyr::filter(Species_ID == "Pru")



### Transformation ----

bestN <- bestNormalize(Pru.PAM_Preconditioning$YII_mean, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN # visualize the best trasformation methods

norm_res <-bestNormalize::boxcox(Pru.PAM_Preconditioning$YII_mean)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(YII_norm_bcx = value)

Pru.PAM_Preconditioning <- bind_cols(Pru.PAM_Preconditioning, norm_dat)



#### Models ----

#Data were transformed using the box-cox transformation method
Pru.lmer.Prec <- lmerTest::lmer(YII_norm_bcx ~ Prec + (1|Colony), data= Pru.PAM_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Pru.lmer.Prec))

# Residual normality
check_normality(Pru.lmer.Prec)

# Check homogeneity
check_homogeneity(Pru.lmer.Prec)

#Check heteroscedasticity
check_heteroscedasticity(Pru.lmer.Prec)

## Residual plot
plot(residuals(Pru.lmer.Prec) ~ fitted(Pru.lmer.Prec), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Pru.lmer.Prec), main = "Histogram of residuals")


#### F statistic (ANOVA)----
Pru_Prec_PAM_anova<-as_tibble(anova(Pru.lmer.Prec))
Pru_Prec_PAM_anova

# write.table(Pru_Prec_PAM_anova, sep = ";", "Statistics/Pru_Prec_PAM_anova.csv")




#### Post-hoc analysis -------

Pru.prec.result <- tidy(glht(Pru.lmer.Prec, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>%  add_significance("adj.p.value")
Pru.prec.result

# write.table(Pru.prec.result, sep = ";", "Statistics/Pru_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Pru.lmer.Prec)
report(Pru.lmer.Prec)







## Post-Prec - GALAXEA FASCICULARIS (Gfa) ------------------------------------------------------------------------


#### Subset ----


Gfa.PAM_Preconditioning <-
  PAM_Preconditioning %>%
  dplyr::filter(Species_ID == "Gfa")




#### Models ----


Gfa.lmer.Prec <- lme4::lmer(YII_mean ~ Prec + (1|Colony)+ (1|Tank_n), data= Gfa.PAM_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Gfa.lmer.Prec))

# Residual normality
check_normality(Gfa.lmer.Prec)

# Check homogeneity
check_homogeneity(Gfa.lmer.Prec)

#Check heteroscedasticity

## Residual plot
plot(residuals(Gfa.lmer.Prec) ~ fitted(Gfa.lmer.Prec), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Gfa.lmer.Prec), main = "Histogram of residuals")


#### F statistic (ANOVA)----
Gfa_Prec_PAM_anova<-as_tibble(anova(Gfa.lmer.Prec))
Gfa_Prec_PAM_anova

# write.table(Gfa_Prec_PAM_anova, sep = ";", "Statistics/Gfa_Prec_PAM_anova.csv")



#### Post-hoc analysis -------------

Gfa.prec.result <- tidy(glht(Gfa.lmer.Prec, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>%  add_significance("adj.p.value")
Gfa.prec.result

# write.table(Gfa.prec.result, sep = ";", "Statistics/Gfa_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Gfa.lmer.Prec)
report(Gfa.lmer.Prec)










## Post-Prec - ACROPORA MURICATA (Amu) ------------------------------------------------------------------

#### Subset ----

Amu.PAM_Preconditioning <-
  PAM_Preconditioning %>%
  dplyr::filter(Species_ID == "Amu")


### Transformation ----

# This function compare different transformation methods to find the best one
bestN <- bestNormalize(Amu.PAM_Preconditioning$YII_mean, loo = T, warn = T, allow_orderNorm = T, 
                      allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)

bestN


norm_res <-orderNorm(Amu.PAM_Preconditioning$YII_mean)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(YII_norm = value)

Amu.PAM_Preconditioning <- bind_cols(Amu.PAM_Preconditioning, norm_dat)



#### Models ----

# Data are transformed using the Ordered Quantile Normalizing Transformation
Amu.lmer.Prec <- lmerTest::lmer(YII_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Amu.PAM_Preconditioning, REML = T)


# Residual plot
qqPlot(residuals(Amu.lmer.Prec))

# Residual normality
check_normality(Amu.lmer.Prec)

# Check homogeneity
check_homogeneity(Amu.lmer.Prec)

#Check heteroscedasticity
check_heteroscedasticity(Amu.lmer.Prec)

## Residual plot
plot(residuals(Amu.lmer.Prec) ~ fitted(Amu.lmer.Prec), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Amu.lmer.Prec), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Amu.lmer.Prec)
fitted_lmer <- fitted(Amu.lmer.Prec)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Amu_Prec_PAM_anova<-as_tibble(anova(Amu.lmer.Prec))
Amu_Prec_PAM_anova

# write.table(Amu_Prec_PAM_anova, sep = ";", "Statistics/Amu_Prec_PAM_anova.csv")




#### Post-hoc analysis --------

Amu.prec.result <- tidy(glht(Amu.lmer.Prec, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>%
  add_significance("adj.p.value", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                   symbols = c("***", "**", "*", "ns"))
Amu.prec.result

 #write.table(Amu.prec.result, sep = ";", "Amu_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Amu.lmer.Prec)
report(Amu.lmer.Prec)











## Post-Prec - POCILLOPORA VERRUCOSA (Pve) --------------------------------------------------------

#### Subset ----

# To test the effect of preconditioning treatments on baseline effective quantum yield, 
# the main df is subset to include data of one species at the time.
Pve.PAM_Preconditioning <-
  PAM_Preconditioning %>%
  dplyr::filter(Species_ID == "Pve")



### Transformation ----

bestN <- bestNormalize(Pve.PAM_Preconditioning$YII_mean, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <-boxcox(Pve.PAM_Preconditioning$YII_mean)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(YII_norm = value)

Pve.PAM_Preconditioning <- bind_cols(Pve.PAM_Preconditioning, norm_dat)




#### Models ----

# Data are transformed using the box-cox Transformation
Pve.lmer.Prec <- lmerTest::lmer(YII_norm ~ Prec + (1|Colony), data= Pve.PAM_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Pve.lmer.Prec))

# Residual normality
check_normality(Pve.lmer.Prec)

# Check homogeneity
check_homogeneity(Pve.lmer.Prec)


#Check heteroscedasticity
check_heteroscedasticity(Pve.lmer.Prec)

## Residual plot
plot(residuals(Pve.lmer.Prec) ~ fitted(Pve.lmer.Prec), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Pve.lmer.Prec), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Pve.lmer.Prec)
fitted_lmer <- fitted(Pve.lmer.Prec)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Pve_Prec_PAM_anova<-as_tibble(anova(Pve.lmer.Prec))
Pve_Prec_PAM_anova

#write.table(Pve_Prec_PAM_anova, sep = ";", "Statistics/Pve_Prec_PAM_anova.csv")





#### Post-hoc analysis ----

Pve.prec.result <- tidy(glht(Pve.lmer.Prec, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>% 
  add_significance("adj.p.value", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                   symbols = c("***", "**", "*", "ns"))%>% tibble()
Pve.prec.result

write.table(Pve.prec.result, sep = ";", "Pve_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Pve.lmer.Prec)
report(Pve.lmer.Prec)







## Post-Prec - STYLOPHORA PISTILLATA (Spi) --------------------------------------------------------


#### Subset ----

# To test the effect of preconditioning treatments on baseline effective quantum yield, 
# the main df is subset to include data of one species at the time.
Spi.PAM_Preconditioning <-
  PAM_Preconditioning %>%
  dplyr::filter(Species_ID == "Spi")


### Transformation ----

bestN <- bestNormalize(Spi.PAM_Preconditioning$YII_mean, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <-boxcox(Spi.PAM_Preconditioning$YII_mean)
norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(YII_norm = value)
Spi.PAM_Preconditioning <- bind_cols(Spi.PAM_Preconditioning, norm_dat)




#### Models ------------


# Data are transformed using the box-cox Transformation
Spi.lmer.Prec <- lmerTest::lmer(YII_norm ~ Prec + (1|Colony) , data= Spi.PAM_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Spi.lmer.Prec))

# Residual normality
check_normality(Spi.lmer.Prec)

# Check homogeneity
check_homogeneity(Spi.lmer.Prec)


#Check heteroscedasticity
check_heteroscedasticity(Spi.lmer.Prec)

## Residual plot
plot(residuals(Spi.lmer.Prec) ~ fitted(Spi.lmer.Prec), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Spi.lmer.Prec), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Spi.lmer.Prec)
fitted_lmer <- fitted(Spi.lmer.Prec)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Spi_Prec_PAM_anova<-as_tibble(anova(Spi.lmer.Prec))
Spi_Prec_PAM_anova

#write.table(Spi_Prec_PAM_anova, sep = ";", "Statistics/Spi_Prec_PAM_anova.csv")




#### Post-hoc analysis --------------

Spi.prec.result <- tidy(glht(Spi.lmer.Prec, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>% 
  add_significance("adj.p.value", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                   symbols = c("***", "**", "*", "ns"))
Spi.prec.result

write.table(Spi.prec.result, sep = ";", "Spi_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Spi.lmer.Prec)
report(Spi.lmer.Prec)




## Significance plot position  ----------------------------

#To add the lmer p values the results for each model was joint in one single table and then, the position was
#obtained from another test (e.g Kruscal-Wallis) replasing the p values with those from lmer

Spi.prec.result1 <- Spi.prec.result %>%
 separate(contrast, sep = " - ",
          into = c("group2", "group1"))%>%
  add_column(Species = "Stylophora pistillata")

Amu.prec.result1 <- Amu.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Acropora muricata")

Pru.prec.result1 <- Pru.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Porites rus")

Mdi.prec.result1 <- Mdi.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Montipora digitata")

Gfa.prec.result1 <- Gfa.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Galaxea fascicularis")

Pve.prec.result1 <- Pve.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Pocillopora verrucosa")


# bind all rows together
Prec_significance <- bind_rows(Spi.prec.result1, Amu.prec.result1, Pru.prec.result1, Mdi.prec.result1, Gfa.prec.result1, Pve.prec.result1)%>%
  rename( p.adj.signif =adj.p.value.signif )%>%
  mutate(term =NULL)

#run a test to that include all species together, to obtain the coordinates where to add the annotations
Pwc_mean <- PAM_Preconditioning %>%
  mutate(across(c(Species_ID, Prec, Treatment, Tank_n, Colony, Day, Frg_n, Region, Treatment), as.factor))%>%
  mutate (Species =  factor(Species, level = c("Galaxea fascicularis", "Acropora muricata", "Porites rus", "Pocillopora verrucosa",
                                               "Montipora digitata", "Stylophora pistillata"
  )))%>%
  
  ungroup() %>%
  group_by(Species) %>%
  rstatix::wilcox_test(YII_mean ~ Prec, p.adjust.method = "BH") %>% #, ref.group = "Ambient"
  #mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))%>%
  rstatix::add_xy_position(x = "Species", dodge = 0.8)

Pwc1 <- Pwc_mean %>%
  mutate(p= NULL , p.adj= NULL, p.adj.signif = NULL)

### Annotations coordinates ----
#Join the two df to get the right annotation coordinates and significance
Prec_significance <- left_join(Pwc1, Prec_significance, by = c("Species", "group1", "group2"))%>% mutate_if(is.numeric, round, 3)


#write_csv2(Prec_significance, here("Prec_annotations.csv"))











## PLOT - PAM Post-preconditioning -----------------------------------------------------------------------------

# This plot shows the effect of preconditioning on baseline effective quantum yield
# Connecting lines and significance symbols are derived from the results of post-hoc tests performed on the mixed-effects model.
# This plot is the one showed in Fig.2 panel A of the manuscript

  PAM_Preconditioning %>%
    mutate (Species =  factor(Species, level = c("Galaxea fascicularis", "Acropora muricata",  "Porites rus", 
                                                 "Pocillopora verrucosa",  "Montipora digitata",   "Stylophora pistillata"
    )))%>%
  ggplot(mapping = aes(x = Species, y = YII_mean, )) +
 
 
  geom_quasirandom(aes(fill = Prec), pch = 21, col = "black", dodge.width = .75, cex = 2.7, width=0.1) +
  geom_boxplot(lwd = 1, aes(fill = Prec), alpha = 0.7) + # outlier.shape = NA
  scale_fill_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
    

  stat_pvalue_manual(
    Prec_significance,
    label = "p.adj.signif", label.size = 14, tip.length = 0.01, hide.ns = TRUE
  ) +
  ylim(0.33,0.72)+
   
   scale_x_discrete(
     label = c("Porites rus" = "*Porites<br>rus*",
               "Galaxea fascicularis" = "*Galaxea<br>fascicularis*",
               "Montipora digitata"= "*Montipora<br>digitata*",
               "Stylophora pistillata"= "*Stylophora<br>pistillata*",
               "Pocillopora verrucosa" ="*Pocillopora<br>verrucosa*",
               "Acropora muricata" = "*Acropora<br>muricata*"))+
  

  
  theme_classic() +
  
  ylab("ΔF/Fm'") +
  xlab("") +

  ggeasy::easy_center_title() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size = 25, face = "bold"),
    axis.text.x = element_markdown(),
    axis.title = element_markdown(size = 25, face = "bold"),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(),
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = 35, face = "bold"),
    plot.subtitle = element_markdown(color = "black", size = 25),
    legend.text = element_text(size = 20),
    plot.caption = element_markdown(size = 15),
    panel.grid.major.y = element_line(colour="lightgray"),
    legend.position = "bottom",
    strip.text = element_text(size = 20)
    )

### ggsave ----


ggsave(here ("output", "PAM_post-preconditioning.png"), dpi = 400, units = "cm", width = 40, height = 25)
















#                               2 Effect of prec. on THERMAL TOLERANCE ----
################################################################################################################################### #
#in the second chapter, we analyzed the effect of precondition on coral thermal tolerance.
# First, we tested the effect of heat exposure on corals and then we calculated the delta between
# post-heat and post-preconditioning effective quantum yield within each preconditioning regime and
# across species to estimate the change in thermal tolerance


## DF subset (Paired data) ------------------

# subset the master df to include only the "post-heat" and "post-preconditioning" values to calculated the paired difference (delta).
# As before, the mean values of effective quantum yield is calculated (each fragment was measured in three different points)


Pam_heat_paired  <- H_rcv %>%
  mutate(across(c(F, Fm, YII), as.double))%>%
  
  dplyr::filter(
    Species_ID == "Pru"& Day == 0| 
      Species_ID == "Pru"& Day == 2|
      Species_ID == "Gfa"& Day == 0| 
      Species_ID == "Gfa"& Day == 3|
      Species_ID == "Mdi"& Day == 0|
      Species_ID == "Mdi"& Day == 3|
      Species_ID == "Pve"& Day == 0|
      Species_ID == "Pve"& Day == 1|
      Species_ID == "Spi"& Day == 0|
      Species_ID == "Spi"& Day == 1|
      Species_ID == "Amu"& Day == 0|
      Species_ID == "Amu"& Day == 1) %>%
  
  #Shorten the species names
  mutate(Species =case_when(str_detect(Species, "Porites rus") ~ "P. rus",
                                  str_detect(Species, "Galaxea fascicularis") ~ "G. fascicularis",
                                  str_detect(Species, "Montipora digitata") ~ "M. digitata",
                                  str_detect(Species, "Stylophora pistillata") ~ "S. pistillata",
                                  str_detect(Species, "Pocillopora verrucosa") ~ "P. verrucosa",
                                  str_detect(Species, "Acropora muricata") ~ "A. muricata"))%>%
  
  mutate (Species =  factor(Species, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata", 
                                               "P. verrucosa", "S. pistillata"
  )))%>%
  
  #add a new variable to calculate the delta
  mutate(Heat_stress =ifelse(Day =="0", "Post-preconditioning","Post-heat"))%>%
  mutate(across(c(Species_ID, Prec, Treatment, Tank_n, Colony, Day, Frg_n, Region, Treatment, Heat_stress),
                as.factor))%>%
  mutate (Heat_stress =  factor(Heat_stress, level = c("Post-preconditioning", "Post-heat"
  )))%>%
  
  #calculate the mean values of each fragment
  group_by(Species_ID, Species ,Prec, Colony, Region, Tank_n, Frg_n, Day, Heat_stress,Treatment)%>%
  summarise(F_mean = mean(F),F_sd = sd(F),
            Fm_mean = mean(Fm),Fm_sd = sd(Fm),
            YII_mean = mean(YII),YII_sd = sd(YII)
  )%>%
  mutate_if(is.numeric, round, 3) %>%
  
  # create a new variable to allow the color pattern in the following ggplot
  mutate(Prec_Day = factor(Prec:Heat_stress))%>%

  ungroup() 








## Statistical test - effect of heat exposure ------------------





###  Kruskal-Wallis test ----------------------



# Kruskal-Wallis test was used to test whether the heat had an effect on the effective quantum yield within each species



# "Heat" treatment
Heat_PAM_kruskal <- Pam_heat_paired %>%
  filter(Treatment=="Heat")%>%
  
  group_by(Species) %>%
  rstatix::kruskal_test(YII_mean ~ Heat_stress) %>%
  mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                            symbols = c("***", "**", "*", "ns"))



# "Control" treatment
Heat_PAM_kruskal_control <- Pam_heat_paired %>%
  filter(Treatment=="Control")%>%
  
  group_by(Species) %>%
  rstatix::kruskal_test(YII_mean ~ Heat_stress) %>%
  mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                            symbols = c("***", "**", "*", "ns"))


### Wilcox test -----------------------------------

# post-hoc test - Heat 


Heat_PAM_Wilcox <- Pam_heat_paired %>%
  mutate (Species =  factor(Species, level = c("G. fascicularis", "A. muricata", "P. rus",  
                                               "P. verrucosa",  "M. digitata", "S. pistillata"
  )))%>%
  filter(Treatment=="Heat")%>%
  mutate(Prec_Day = factor(Prec:Heat_stress))%>%
  group_by(Species, Prec) %>%
  rstatix::pairwise_wilcox_test(YII_mean ~ Prec_Day, p.adjust.method = "BH") %>%
 
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))%>%
  mutate_if (is.numeric, round, 6)%>%
  
  rstatix::add_xy_position(x = "Prec_Day", dodge = 0.8)

#save the results
write_csv2(Heat_PAM_Wilcox,here("Output/Statistics","Wilcox_Heat_PAM.result.csv"))

here()


# post-hoc test - Control 


Heat_PAM_Wilcox_control <- Pam_heat_paired %>%
  mutate (Species =  factor(Species, level = c("G. fascicularis", "A. muricata", "P. rus",  
                                               "P. verrucosa",  "M. digitata", "S. pistillata"
  )))%>%
  filter(Treatment=="Control")%>%
  mutate(Prec_Day = factor(Prec:Heat_stress))%>%
  group_by(Species, Prec) %>%
  rstatix::pairwise_wilcox_test(YII_mean ~ Prec_Day, p.adjust.method = "BH") %>%
  
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))%>%
  mutate_if (is.numeric, round, 6)%>%
  
  rstatix::add_xy_position(x = "Prec_Day", dodge = 0.8)






##                            PLOT - POST-HEAT paired values ----
####################################################################################################################### #
#plot showing the coral response to heat (post-heat) in comparison to the paired values at the "post-preconditioning" time point





#### Plot - Paired "Heat" --------------------

#paired data of the "Heat" treatment


Pam_heat_paired%>%
  mutate(Prec_Day = factor(Prec:Heat_stress))%>%
  mutate (Species =  factor(Species, level = c("G. fascicularis",  "P. rus",  "A. muricata",
                                               "M. digitata", "P. verrucosa",  "S. pistillata"
  )))%>%
  
  dplyr::filter(Treatment=="Heat")%>%
  ggplot(mapping = aes(x = Prec_Day, y = YII_mean))+
  
  geom_quasirandom (aes( fill = Prec_Day, shape = Heat_stress ), 
                    col = "black",dodge.width= .9, size = 2 )+
  geom_boxplot(aes(fill = Prec_Day, shape = Heat_stress),lwd = 0.7)+ 

  scale_fill_manual(values=c("#4682B4B3", "#63B8FF40",  "#B4B446B3", "#EEEE0040",  "#D4711CB3", "#FF7F2440" ))+
  
  scale_shape_manual(values=c(21,24))+
  stat_pvalue_manual(Heat_PAM_Wilcox,
                     label = "p.adj.signif", label.size = 14, tip.length = 0.01, hide.ns = T
  ) +
  
  theme_classic() +
  
  ylab("ΔF/Fm'") +
  xlab("") +
  ggeasy::easy_center_title() +
  
  theme(
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(size = 20),
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = 35, face = "bold"),
    plot.subtitle = element_text(color = "black", size = 25),
    legend.text = element_text(size = 20),
    strip.text.x = element_text(size = 25, face = "bold.italic"),
    strip.text.y = element_text(size = 20, face = "bold"),
    strip.placement = "outside",
    panel.grid.major.y = element_line(colour="lightgray"), 
    strip.background = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position = "none") +
  
  facet_wrap(.~ Species, ncol = 6, strip.position ="bottom")

ggsave(here ("Plot", "Output/PAM_post-heat_paired-heat.png"), dpi = 400, units = "cm", width = 53.4, height = 11.7)




#### Plot - Paired "Control" ------------------------------------------


#paired data of the "Control" treatment


Pam_heat_paired%>%
  mutate(Prec_Day = factor(Prec:Heat_stress))%>%
  mutate (Species =  factor(Species, level = c("G. fascicularis", "A. muricata", "P. rus",  
                                               "P. verrucosa",  "M. digitata", "S. pistillata"
  )))%>%
  
  dplyr::filter(Treatment=="Control")%>%
  ggplot(mapping = aes(x = Prec_Day, y = YII_mean))+
  
  geom_quasirandom (aes( fill = Prec_Day, shape = Heat_stress ), 
                    col = "black",dodge.width= .9, size = 2 )+
  geom_boxplot(aes(fill = Prec_Day, shape = Heat_stress),lwd = 0.7)+ 
  
  scale_fill_manual(values=c("#4682B4B3", "#63B8FF40",  "#B4B446B3", "#EEEE0040",  "#D4711CB3", "#FF7F2440" ))+
  
  scale_shape_manual(values=c(21,24))+
  stat_pvalue_manual(Heat_PAM_Wilcox_control,
                     label = "p.adj.signif", label.size = 14, tip.length = 0.01, hide.ns = T
  ) +
  
  theme_classic() +
  
  ylab("ΔF/Fm'") +
  xlab("") +
  ggeasy::easy_center_title() +
  
  theme(
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(size = 20),
    legend.title = element_blank(),
    plot.title = element_text(color = "black", size = 35, face = "bold"),
    plot.subtitle = element_text(color = "black", size = 25),
    legend.text = element_text(size = 20),
    strip.text.x = element_text(size = 25, face = "bold.italic"),
    strip.text.y = element_text(size = 20, face = "bold"),
    strip.placement = "outside",
    panel.grid.major.y = element_line(colour="lightgray"), 
    strip.background = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.position = "bottom") +
  
  facet_wrap(.~ Species, ncol = 6, strip.position ="bottom")

ggsave(here ("Plot", "PAM_post-heat_paired-control.png"), dpi = 400, units = "cm", width = 53.5, height = 14)










## PAIRED DELTA calculation -----------------------------------------------------------------

#Delta values are calculated as Post-heat minus Post-preconditioning paired ΔF/Fm' values
PAM_heat_delta <- Pam_heat_paired %>%
  filter(Treatment  == "Heat") %>%
  group_by(Species_ID, Species ,Prec, Colony, Region, Treatment, Tank_n, Frg_n)%>%
  mutate(Delta_pair = lead(YII_mean) - YII_mean )


##                          THERMAL TOLERANCE CHANGES (Linear Mixed-Effects Models) ----
################################################################################################################## #


# The effect of preconditioning was evaluated comparing the paired delta difference 
# (Post-heat minus Post-preconditioning) of each preconditioning regime, separately for each species





## Heat - MONTIPORA DIGITATA (Mdi) -----------------------

#### Subset ----

Mdi.PAM_heat_delta <- PAM_heat_delta %>%
  dplyr::filter(Species_ID == "Mdi")%>%
  mutate(Delta_pair_1 = (Delta_pair+1))# add 1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Mdi.PAM_heat_delta$Delta_pair, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <- orderNorm(Mdi.PAM_heat_delta$Delta_pair)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Mdi.PAM_heat_delta <- bind_cols(Mdi.PAM_heat_delta, norm_dat)



#### Models ----

# Data are transformed using the Ordered Quantile Normalizing Transformation 
Mdi.lmer.heat <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony), data= Mdi.PAM_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Mdi.lmer.heat)) 

# Residual normality
check_normality(Mdi.lmer.heat)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Mdi.lmer.heat)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Mdi.lmer.heat)

## Residual plot
plot(residuals(Mdi.lmer.heat) ~ fitted(Mdi.lmer.heat), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Mdi.lmer.heat), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Mdi.lmer.heat)
fitted_lmer <- fitted(Mdi.lmer.heat)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals



#### F statistic (ANOVA)----
Mdi.lmer.heat_anova<-as_tibble(anova(Mdi.lmer.heat))
Mdi.lmer.heat_anova # ANOVA test is not significant

#write.table(Mdi.lmer.heat_anova, sep = ";", "Statistics/Mdi_Prec_PAM_anova.csv")

# summary and report for better interpretation
summary(Mdi.lmer.heat)
report(Mdi.lmer.heat)








## Heat - GALAXEA FASCICULARIS (Gfa) ----------------------------------------

#### Subset ----

Gfa.PAM_heat_delta <- PAM_heat_delta %>%
  dplyr::filter(Species_ID == "Gfa")%>%
  mutate(Delta_pair_1 = (Delta_pair+1))%>%
  drop_na(Delta_pair_1)# add 1 to make all values positive


### Transformation ----

bestN <- bestNormalize(Gfa.PAM_heat_delta$Delta_pair_1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <- orderNorm(Gfa.PAM_heat_delta$Delta_pair)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Gfa.PAM_heat_delta <- bind_cols(Gfa.PAM_heat_delta, norm_dat)



#### Models ----

# Data are transformed using the Ordered Quantile Normalizing Transformation

Gfa.lmer.heat <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony), data= Gfa.PAM_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Gfa.lmer.heat)) 

# Residual normality
check_normality(Gfa.lmer.heat)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Gfa.lmer.heat)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Gfa.lmer.heat)

## Residual plot
plot(residuals(Gfa.lmer.heat) ~ fitted(Gfa.lmer.heat), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Gfa.lmer.heat), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Gfa.lmer.heat)
fitted_lmer <- fitted(Gfa.lmer.heat)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals



#### F statistic (ANOVA)----
Gfa.lmer.heat_anova<-as_tibble(anova(Gfa.lmer.heat))
Gfa.lmer.heat_anova

#write.table(Gfa.lmer.heat_anova, sep = ";", "Statistics/Gfa_Prec_PAM_anova.csv")




#### Post-hoc analysis ----

Gfa.lmer.heat.ems <- emmeans::emmeans(Gfa.lmer.heat, list(pairwise ~ Prec), adjust = "bonferroni")

(Gfa.lmer.heat.ems.results <- Gfa.lmer.heat.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())


#write.table(Gfa.lmer.heat.result, sep = ";", "Statistics/Gfa_Result_Prec_Pam.csv")


# summary and report for better interpretation
summary(Gfa.lmer.heat)
report(Gfa.lmer.heat)






## Heat - ACROPORA MURICATA (Amu) ----------------------------------------------

#### Subset ----

Amu.PAM_heat_delta <- PAM_heat_delta %>%
  dplyr::filter(Species_ID == "Amu")%>%
  mutate(Delta_pair_1 = (Delta_pair+1))%>%
  drop_na(Delta_pair_1)# add 1 to make all values positive


### Transformation ----

bestN <- bestNormalize(Amu.PAM_heat_delta$Delta_pair, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

# Data are transformed using the Box-cox Transformation
norm_res <- bestNormalize::boxcox(Amu.PAM_heat_delta$Delta_pair_1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Amu.PAM_heat_delta <- bind_cols(Amu.PAM_heat_delta, norm_dat)

#### Models ----

# Data are transformed using the Box-cox Transformation

#two models were tested and the with random Intercept but fixed slope works better (Amu.lmer.heat) 
Amu.lmer.heat <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Amu.PAM_heat_delta, REML = T)

Amu.lmer.heat.slope <- lmerTest::lmer(Delta_pair_norm ~ Prec + (Prec|Colony) + (1|Tank_n), data= Amu.PAM_heat_delta, REML = T)

#Models comparison
AIC(Amu.lmer.heat, Amu.lmer.heat.slope)
anova(Amu.lmer.heat, Amu.lmer.heat.slope)

# Residual plot
qqPlot(residuals(Amu.lmer.heat)) 

# Residual normality
check_normality(Amu.lmer.heat)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Amu.lmer.heat)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Amu.lmer.heat)

## Residual plot
plot(residuals(Amu.lmer.heat) ~ fitted(Amu.lmer.heat), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Amu.lmer.heat), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Amu.lmer.heat)
fitted_lmer <- fitted(Amu.lmer.heat)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals




#### F statistic (ANOVA)----
Amu.lmer.heat_anova<-as_tibble(anova(Amu.lmer.heat))
Amu.lmer.heat_anova





#### Post-hoc analysis ----

Amu.lmer.heat.ems <- emmeans::emmeans(Amu.lmer.heat, list(pairwise ~ Prec), adjust = "bonferroni")

(Amu.lmer.heat.ems.results <- Amu.lmer.heat.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())



# summary and report for better interpretation
summary(Amu.lmer.heat)
report(Amu.lmer.heat)







##Heat -PORITES RUS (Pru) -----------------------------------------------------------

#### Subset ----

Pru.PAM_heat_delta <- PAM_heat_delta %>%
  dplyr::filter(Species_ID == "Pru")%>%
  mutate(Delta_pair_1 = (Delta_pair+1))%>%
  drop_na(Delta_pair_1)# add 1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Pru.PAM_heat_delta$Delta_pair_1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN


norm_res <- bestNormalize::boxcox(Pru.PAM_heat_delta$Delta_pair_1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Pru.PAM_heat_delta <- bind_cols(Pru.PAM_heat_delta, norm_dat)




#### Models ----

# Data are transformed using the Box-cox Transformation
Pru.lmer.heat <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Pru.PAM_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Pru.lmer.heat)) 

# Residual normality
check_normality(Pru.lmer.heat)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Pru.lmer.heat)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Pru.lmer.heat)

## Residual plot
plot(residuals(Pru.lmer.heat) ~ fitted(Pru.lmer.heat), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Pru.lmer.heat), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Pru.lmer.heat)
fitted_lmer <- fitted(Pru.lmer.heat)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals




#### F statistic (ANOVA)----
Pru.lmer.heat_anova<-as_tibble(anova(Pru.lmer.heat))
Pru.lmer.heat_anova

#write.table(Pru.lmer.heat_anova, sep = ";", "Statistics/Pru_Prec_PAM_anova.csv")





#### Post-hoc analysis ----

Pru.lmer.heat.ems <- emmeans::emmeans(Pru.lmer.heat, list(pairwise ~ Prec), adjust = "bonferroni")

(Pru.lmer.heat.ems.results <- Pru.lmer.heat.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())


#write.table(Pru.lmer.heat.result, sep = ";", "Statistics/Pru_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Pru.lmer.heat)
report(Pru.lmer.heat)








## Heat - POCILLOPORA VERRUCOSA (Pve) ---------------------------------------------------------

#### Subset ----

Pve.PAM_heat_delta <- PAM_heat_delta %>%
  dplyr::filter(Species_ID == "Pve")%>%
  mutate(Delta_pair_1 = (Delta_pair+1))%>%
  drop_na(Delta_pair_1)# add 1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Pve.PAM_heat_delta$Delta_pair_1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN


norm_res <- bestNormalize::boxcox(Pve.PAM_heat_delta$Delta_pair_1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Pve.PAM_heat_delta <- bind_cols(Pve.PAM_heat_delta, norm_dat)





#### Models ----

# Data are transformed using the Box-cox Transformation
Pve.lmer.heat <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Pve.PAM_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Pve.lmer.heat)) 

# Residual normality
check_normality(Pve.lmer.heat)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Pve.lmer.heat)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Pve.lmer.heat)

## Residual plot
plot(residuals(Pve.lmer.heat) ~ fitted(Pve.lmer.heat), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Pve.lmer.heat), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Pve.lmer.heat)
fitted_lmer <- fitted(Pve.lmer.heat)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals





#### F statistic (ANOVA)----
Pve.lmer.heat_anova<-as_tibble(anova(Pve.lmer.heat))
Pve.lmer.heat_anova

#write.table(Pve.lmer.heat_anova, sep = ";", "Statistics/Pve_Prec_PAM_anova.csv")





#### Post-hoc analysis ----

Pve.lmer.heat.ems <- emmeans::emmeans(Pve.lmer.heat, list(pairwise ~ Prec), adjust = "bonferroni")

(Pve.lmer.heat.ems.results <- Pve.lmer.heat.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())


#write.table(Pve.lmer.heat.result, sep = ";", "Statistics/Pve_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Pve.lmer.heat)
report(Pve.lmer.heat)







## Heat - STYLOPHORA PISTILLATA (Spi) ---------------------------

#### Subset ----

Spi.PAM_heat_delta <- PAM_heat_delta %>%
  dplyr::filter(Species_ID == "Spi")%>%
  mutate(Delta_pair_1 = (Delta_pair+1))%>%
  drop_na(Delta_pair_1)# add 1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Spi.PAM_heat_delta$Delta_pair_1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN


norm_res <- bestNormalize::boxcox(Spi.PAM_heat_delta$Delta_pair_1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Spi.PAM_heat_delta <- bind_cols(Spi.PAM_heat_delta, norm_dat)




#### Models ----

# Data are transformed using the Box-cox Transformation
Spi.lmer.heat <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Spi.PAM_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Spi.lmer.heat)) 

# Residual normality
check_normality(Spi.lmer.heat)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Spi.lmer.heat)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Spi.lmer.heat)

## Residual plot
plot(residuals(Spi.lmer.heat) ~ fitted(Spi.lmer.heat), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Spi.lmer.heat), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Spi.lmer.heat)
fitted_lmer <- fitted(Spi.lmer.heat)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals





#### F statistic (ANOVA)----
Spi.lmer.heat_anova<-as_tibble(anova(Spi.lmer.heat))
Spi.lmer.heat_anova

#write.table(Spi.lmer.heat_anova, sep = ";", "Statistics/Spi_Prec_PAM_anova.csv")





#### Post-hoc analysis ----

Spi.lmer.heat.ems <- emmeans::emmeans(Spi.lmer.heat, list(pairwise ~ Prec), adjust = "bonferroni")

(Spi.lmer.heat.ems.results <- Spi.lmer.heat.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())


#write.table(Spi.lmer.heat.result, sep = ";", "Statistics/Spi_Result_Prec_Pam.csv")

# summary and report for better interpretation
summary(Spi.lmer.heat)
report(Spi.lmer.heat)




