

# --------------------------- Ferrara et al. 2024 - R script N. 2 --------------------------- #
#                                                                                             #
#                                                                                             #
#  Effect of preconditioning treatments on tissue color at time points                        #
#  "Post-preconditioning" and "Post-heat"                                                     #
#                                                                                             #
#  In this script we calculate the effect of preconditioning on baseline physiology           #
#  and the increase in thermal tolerance after the heat assay                                 #
#                                                                                             #
#  The effects are calculated on tissue color metrics only                                    #
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

# Load the data of coral tissue color, used as proxy for bleaching. This df contain the tissue color values from the "post-preconditioning",
# "post-heat", "day-15", and "day-30" time points. However, Only DATA ANALYSIS OF "POST-PRECONDITIONING" AND "POST-HEAT" is shown here.

Clr  <- read_csv2(here("Data", "Ferrara_tissue_color_master_df.csv"))

# To perform the tissue color analyses, we calculated the complementary values (255 - x) to show the 0 as white and 255 as bleak.
Clr <- Clr %>%
  mutate(gray_sum = (255 - gray_sum))


#  1 Effect of Prec. on baseline physiology (tissue color) ---- 

#In this first chapter, the effect of thermal preconditioning regimes on coral baseline tissue color
#is tested using LMER models for each species. 

## Preconditioning df subset----

#subset the master df to include only the "post-preconditioning" tissue color values.

Clr_Preconditioning <- Clr %>%
  dplyr::filter(Day == 0) %>% 
  mutate(across(c(Species, ID, Prec, Trtm, Tank_n, Colony, Day, Frg, Region), as.factor),
         across(c(gray_sum), as.numeric))%>%
  
  mutate (Species =  factor(Species, level = c("Galaxea fascicularis", "Porites rus", "Acropora muricata", "Montipora digitata", 
                                               "Pocillopora verrucosa", "Stylophora pistillata"
  )))%>%
  mutate_if(is.numeric, round, 3) %>%
  ungroup()
  


### Shapiro test ----

#Test the normal distribution within each preconditioning treatment across species
Clr_Preconditioning %>%
  ungroup() %>%
  group_by(ID, Prec) %>%
  rstatix::shapiro_test(gray_sum) %>%
  mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))

### Levene test ----

Clr_Preconditioning %>%
  ungroup() %>%
  group_by(ID) %>%
  rstatix::levene_test(gray_sum ~ Prec) %>%
  mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))




# The effect of the two preconditioning treatments on coral tissue color is tested for each species separately
# Post-Preconditioning Linear Mixed-Effects Models ----





## Post-Prec - MONTIPORA DIGITATA (Mdi) ----------------------------------------------------------------------

#### Subset ----

Mdi.Clr_Preconditioning <-
  Clr_Preconditioning %>%
  dplyr::filter(ID == "Mdi")

### Transformation ----

bestN <- bestNormalize(Mdi.Clr_Preconditioning$gray_sum, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN # visualize the best trasformation methods

norm_res <-bestNormalize::boxcox(Mdi.Clr_Preconditioning$gray_sum)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(gray_sum_bcx = value)

Mdi.Clr_Preconditioning <- bind_cols(Mdi.Clr_Preconditioning, norm_dat)


#### Models ----

#Data were transformed using the box-cox transformation method
Mdi.Clr.Prec.lmer <- lmerTest::lmer(gray_sum_bcx ~ Prec + (1|Colony), data= Mdi.Clr_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Mdi.Clr.Prec.lmer))

# Residual normality
check_normality(Mdi.Clr.Prec.lmer)

# Check homogeneity
check_homogeneity(Mdi.Clr.Prec.lmer)

#Check heteroscedasticity
check_heteroscedasticity(Mdi.Clr.Prec.lmer)

## Residual plot
plot(residuals(Mdi.Clr.Prec.lmer) ~ fitted(Mdi.Clr.Prec.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Mdi.Clr.Prec.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Mdi.Clr.Prec.lmer)
fitted_lmer <- fitted(Mdi.Clr.Prec.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Mdi_Prec_Clr_anova<-as_tibble(anova(Mdi.Clr.Prec.lmer))
Mdi_Prec_Clr_anova


#### Post-hoc analysis ----
Mdi.Clr.prec.result <- tidy(glht(Mdi.Clr.Prec.lmer, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>% 
  add_significance("adj.p.value", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                   symbols = c("***", "**", "*", "ns"))
Mdi.Clr.prec.result

# summary and report for better interpretation
summary(Mdi.lmer.Prec)
report(Mdi.lmer.Prec)







## Post-Prec - PORITES RUS (Pru) --------------------------------------------------------------------------------

#### Subset ----

Pru.Clr_Preconditioning <-
  Clr_Preconditioning %>%
  dplyr::filter(ID == "Pru")

### Transformation ----

bestN <- bestNormalize(Pru.Clr_Preconditioning$gray_sum, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN # visualize the best trasformation methods

norm_res <-bestNormalize::boxcox(Pru.Clr_Preconditioning$gray_sum)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(gray_sum_bcx = value)

Pru.Clr_Preconditioning <- bind_cols(Pru.Clr_Preconditioning, norm_dat)


#### Models ----

#Data were transformed using the box-cox transformation method
Pru.Clr.Prec.lmer <- lmerTest::lmer(gray_sum_bcx ~ Prec + (1|Colony), data= Pru.Clr_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Pru.Clr.Prec.lmer))

# Residual normality
check_normality(Pru.Clr.Prec.lmer)

# Check homogeneity
check_homogeneity(Pru.Clr.Prec.lmer)

#Check heteroscedasticity
check_heteroscedasticity(Pru.Clr.Prec.lmer)

## Residual plot
plot(residuals(Pru.Clr.Prec.lmer) ~ fitted(Pru.Clr.Prec.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Pru.Clr.Prec.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Pru.Clr.Prec.lmer)
fitted_lmer <- fitted(Pru.Clr.Prec.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Pru_Prec_Clr_anova<-as_tibble(anova(Pru.Clr.Prec.lmer))
Pru_Prec_Clr_anova


# summary and report for better interpretation
summary(Pru.lmer.Prec)
report(Pru.lmer.Prec)







## Post-Prec - GALAXEA FASCICULARIS (Gfa) --------------------------------------------------------------------------------

#### Subset ----

Gfa.Clr_Preconditioning <-
  Clr_Preconditioning %>%
  dplyr::filter(ID == "Gfa")

### Transformation ----

bestN <- bestNormalize(Gfa.Clr_Preconditioning$gray_sum, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN # visualize the best trasformation methods

norm_res <-bestNormalize::boxcox(Gfa.Clr_Preconditioning$gray_sum)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(gray_sum_bcx = value)

Gfa.Clr_Preconditioning <- bind_cols(Gfa.Clr_Preconditioning, norm_dat)


#### Models ----

#Data were transformed using the box-cox transformation method
Gfa.Clr.Prec.lmer <- lmerTest::lmer(gray_sum_bcx ~ Prec + (1|Colony), data= Gfa.Clr_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Gfa.Clr.Prec.lmer))

# Residual normality
check_normality(Gfa.Clr.Prec.lmer)

# Check homogeneity
check_homogeneity(Gfa.Clr.Prec.lmer)

#Check heteroscedasticity
check_heteroscedasticity(Gfa.Clr.Prec.lmer)

## Residual plot
plot(residuals(Gfa.Clr.Prec.lmer) ~ fitted(Gfa.Clr.Prec.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Gfa.Clr.Prec.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Gfa.Clr.Prec.lmer)
fitted_lmer <- fitted(Gfa.Clr.Prec.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Gfa_Prec_Clr_anova<-as_tibble(anova(Gfa.Clr.Prec.lmer))
Gfa_Prec_Clr_anova

# summary and report for better interpretation
summary(Gfa.lmer.Prec)
report(Gfa.lmer.Prec)






## Post-Prec - ACROPORA MURICATA (Amu) -----------------------------------------------------------------------------------

#### Subset ----

Amu.Clr_Preconditioning <-
  Clr_Preconditioning %>%
  dplyr::filter(ID == "Amu")

### Transformation ----

bestN <- bestNormalize(Amu.Clr_Preconditioning$gray_sum, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN # visualize the best trasformation methods

norm_res <-bestNormalize::boxcox(Amu.Clr_Preconditioning$gray_sum)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(gray_sum_bcx = value)

Amu.Clr_Preconditioning <- bind_cols(Amu.Clr_Preconditioning, norm_dat)


#### Models ----

#Data were transformed using the box-cox transformation method
Amu.Clr.Prec.lmer <- lmerTest::lmer(gray_sum_bcx ~ Prec + (1|Colony) + (1|Tank_n), data= Amu.Clr_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Amu.Clr.Prec.lmer))

# Residual normality
check_normality(Amu.Clr.Prec.lmer)

# Check homogeneity
check_homogeneity(Amu.Clr.Prec.lmer)

#Check heteroscedasticity
check_heteroscedasticity(Amu.Clr.Prec.lmer)

## Residual plot
plot(residuals(Amu.Clr.Prec.lmer) ~ fitted(Amu.Clr.Prec.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Amu.Clr.Prec.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Amu.Clr.Prec.lmer)
fitted_lmer <- fitted(Amu.Clr.Prec.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Amu_Prec_Clr_anova<-as_tibble(anova(Amu.Clr.Prec.lmer))
Amu_Prec_Clr_anova


#### Post-hoc analysis ----
Amu.Clr.prec.result <- tidy(glht(Amu.Clr.Prec.lmer, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>% 
  add_significance("adj.p.value", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                   symbols = c("***", "**", "*", "ns"))
Amu.Clr.prec.result

# summary and report for better interpretation
summary(Amu.lmer.Prec)
report(Amu.lmer.Prec)







## Post-Prec - POCILLOPORA VERRUCOSA (Pve) ----------------------------------------------

#### Subset ----

Pve.Clr_Preconditioning <-
  Clr_Preconditioning %>%
  dplyr::filter(ID == "Pve")

### Transformation ----

bestN <- bestNormalize(Pve.Clr_Preconditioning$gray_sum, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN # visualize the best trasformation methods

norm_res <-bestNormalize::boxcox(Pve.Clr_Preconditioning$gray_sum)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(gray_sum_bcx = value)

Pve.Clr_Preconditioning <- bind_cols(Pve.Clr_Preconditioning, norm_dat)


#### Models ----

#Data were transformed using the box-cox transformation method
Pve.Clr.Prec.lmer <- lmerTest::lmer(gray_sum_bcx ~ Prec + (1|Colony) + (1|Tank_n), data= Pve.Clr_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Pve.Clr.Prec.lmer))

# Residual normality
check_normality(Pve.Clr.Prec.lmer)

# Check homogeneity
check_homogeneity(Pve.Clr.Prec.lmer)

#Check heteroscedasticity
check_heteroscedasticity(Pve.Clr.Prec.lmer)

## Residual plot
plot(residuals(Pve.Clr.Prec.lmer) ~ fitted(Pve.Clr.Prec.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Pve.Clr.Prec.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Pve.Clr.Prec.lmer)
fitted_lmer <- fitted(Pve.Clr.Prec.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Pve_Prec_Clr_anova<-as_tibble(anova(Pve.Clr.Prec.lmer))
Pve_Prec_Clr_anova


#### Post-hoc analysis ----
Pve.Clr.prec.result <- tidy(glht(Pve.Clr.Prec.lmer, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>% 
  add_significance("adj.p.value", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                   symbols = c("***", "**", "*", "ns"))
Pve.Clr.prec.result

# summary and report for better interpretation
summary(Pve.lmer.Prec)
report(Pve.lmer.Prec)






## Post-Prec - STYLOPHORA PISTILLATA (Spi) ----------------------------------------------------

#### Subset ----

Spi.Clr_Preconditioning <-
  Clr_Preconditioning %>%
  dplyr::filter(ID == "Spi")

### Transformation ----

bestN <- bestNormalize(Spi.Clr_Preconditioning$gray_sum, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN # visualize the best trasformation methods

norm_res <-bestNormalize::boxcox(Spi.Clr_Preconditioning$gray_sum)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(gray_sum_bcx = value)

Spi.Clr_Preconditioning <- bind_cols(Spi.Clr_Preconditioning, norm_dat)


#### Models ----

#Data were transformed using the box-cox transformation method
Spi.Clr.Prec.lmer <- lmerTest::lmer(gray_sum_bcx ~ Prec + (1|Colony) + (1|Tank_n), data= Spi.Clr_Preconditioning, REML = T)

# Residual plot
qqPlot(residuals(Spi.Clr.Prec.lmer))

# Residual normality
check_normality(Spi.Clr.Prec.lmer)

# Check homogeneity
check_homogeneity(Spi.Clr.Prec.lmer)

#Check heteroscedasticity
check_heteroscedasticity(Spi.Clr.Prec.lmer)

## Residual plot
plot(residuals(Spi.Clr.Prec.lmer) ~ fitted(Spi.Clr.Prec.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Spi.Clr.Prec.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Spi.Clr.Prec.lmer)
fitted_lmer <- fitted(Spi.Clr.Prec.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test)

#### F statistic (ANOVA)----
Spi_Prec_Clr_anova<-as_tibble(anova(Spi.Clr.Prec.lmer))
Spi_Prec_Clr_anova


#### Post-hoc analysis ----
Spi.Clr.prec.result <- tidy(glht(Spi.Clr.Prec.lmer, linfct = mcp(Prec = "Tukey"), test = adjusted("holm"))) %>% 
  add_significance("adj.p.value", cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                   symbols = c("***", "**", "*", "ns"))
Spi.Clr.prec.result

# summary and report for better interpretation
summary(Spi.lmer.Prec)
report(Spi.lmer.Prec)



##### Spi_IMG. Emmeans test ----

Spi_IMG.ems <- emmeans::emmeans(Spi_IMG.lmer4, list(pairwise ~ Prec), adjust = "holm")

(Spi_IMG.pairwise <- Spi_IMG.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 3)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())







#### Significance annotations position ---------------------------------------------------------------------------------

#To add the lmer p values the results for each model was joint in one single table and then, the position was
#obtained from another test (e.g Kruscal-Wallis) replasing the p values with those from lmer

Spi.Clr.prec.result1 <- Spi.Clr.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Stylophora pistillata")

Amu.Clr.prec.result1 <- Amu.Clr.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Acropora muricata")

Pru.Clr.prec.result1 <- Pru.Clr.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Porites rus") %>% 
  mutate(adj.p.value.signif = "ns")# the anova on the lmer model was non significant therefore we don't have to calculate the post-hoc test

Mdi.Clr.prec.result1 <- Mdi.Clr.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Montipora digitata")

Gfa.Clr.prec.result1 <- Gfa.Clr.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Galaxea fascicularis")

Pve.Clr.prec.result1 <- Pve.Clr.prec.result %>%
  separate(contrast, sep = " - ",
           into = c("group2", "group1"))%>%
  add_column(Species = "Pocillopora verrucosa")



# bind all rows together
Clr.Prec_significance <- bind_rows(Spi.Clr.prec.result1, Amu.Clr.prec.result1, Pru.Clr.prec.result1, Mdi.Clr.prec.result1, Gfa.Clr.prec.result1, Pve.Clr.prec.result1)%>%
  rename( p.adj.signif =adj.p.value.signif )%>%
  mutate(term =NULL)

#run a test to that include all species together, to obtain the coordinates where to add the annotations
Clr.Pos.bp <- Clr_Preconditioning %>%
  group_by(Species) %>%
  rstatix::wilcox_test(gray_sum ~ Prec, p.adjust.method = "BH") %>% #, ref.group = "Ambient"
  #mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))%>%
  rstatix::add_xy_position(x = "Species", dodge = 0.8)

Clr.Pos.bp <- Clr.Pos.bp %>%
  mutate(p= NULL , p.adj= NULL, p.adj.signif = NULL)

### Annotations coordinates ----
#Join the two df to get the right annotation coordinates and significance
Prec_Clr_significance <- left_join(Clr.Pos.bp, Clr.Prec_significance, by = c("Species", "group1", "group2"))%>% mutate_if(is.numeric, round, 3)


  








## PLOT - Post-preconditioning ----------------------------------------------------------------------------------------------------------------- 

# This plot shows the effect of preconditioning on baseline tissue color
# Connecting lines and significance symbols are derived from the results of post-hoc tests performed on the mixed-effects model.
# This plot is the one showed in Fig.2 panel C of the manuscript
 
Clr_Preconditioning %>%
 
  ggplot(mapping = aes(x = Species, y = gray_sum))+

  geom_quasirandom (aes( fill = Prec ),width=0.1, pch= 21, col = "black", dodge.width= .75, cex = 2.7 )+
  geom_boxplot(lwd=1,aes(fill = Prec), alpha = 0.7) +
  scale_fill_manual(values=c("#4682B4", "#B4B446", "#D4711C"))+
  
  stat_pvalue_manual(
    Prec_Clr_significance,  label = "p.adj.signif", label.size = 14, tip.length = 0.01, hide.ns = TRUE
  ) +
  theme_classic() +
  
  ylab("Tissue color intensity")+
  xlab("")+
  ggeasy::easy_center_title()+
  
  scale_x_discrete(
    label = c("Porites rus" = "P. rus",
              "Galaxea fascicularis" = "G. fascicularis",
              "Montipora digitata"= "M. digitata",
              "Stylophora pistillata"= "S. pistillata",
              "Pocillopora verrucosa" ="P. verrucosa",
              "Acropora muricata" = "A. muricata"))+
  
  theme(
    panel.grid.minor.y = element_blank(),
    
    
    axis.text = element_text(size=20, face="bold"),
    axis.title = element_text(size = 25, face="bold"),
    axis.title.y = element_text( margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(),
    legend.title = element_blank(),
    plot.title = element_markdown(color = "black", size = 35, face = "bold"),
    #plot.subtitle = element_text(color = "black", size = 25),
    legend.text = element_text(size = 20),
    plot.subtitle = element_markdown(color = "black", size = 25),
    panel.grid.major.y = element_line(colour="lightgray"),
    plot.caption = element_markdown(size = 15),
    legend.position="bottom",
    axis.text.x = element_text(size = 25,face = "bold.italic")
  ) +
  theme(strip.text = element_text(
    size = 20))

## ggsave ----

ggsave(here ("Plot", "CLR_Post-Preconditioning.png"), dpi = 400, units = "cm", width = 45, height = 18)


# 2 Effect of prec. on THERMAL TOLERANCE ----

# In the second chapter, we analyzed the effect of precondition on coral thermal tolerance.
# First, we tested the effect of heat exposure on corals and then we calculated the delta between
# post-heat and post-preconditioning tissue color values within each preconditioning regime and
# across species to estimate the change in thermal tolerance after preconditioning


## DF subset (Paired data calculation) ----

# Subset the master df to include only the "post-heat" and "post-preconditioning" values to calculated the paired difference (delta).
Clr_heat_paired <- Clr %>%

dplyr::filter(
  ID == "Pru"& Day == 0| 
    ID == "Pru"& Day == 2|
    ID == "Gfa"& Day == 0| 
    ID == "Gfa"& Day == 3|
    ID == "Mdi"& Day == 0|
    ID == "Mdi"& Day == 3|
    ID == "Pve"& Day == 0|
    ID == "Pve"& Day == 1|
    ID == "Spi"& Day == 0|
    ID == "Spi"& Day == 1|
    ID == "Amu"& Day == 0|
    ID == "Amu"& Day == 1) %>%
  ungroup()%>%
  
  mutate(Species_short =case_when(str_detect(Species, "Porites rus") ~ "P. rus",
                            str_detect(Species, "Galaxea fascicularis") ~ "G. fascicularis",
                            str_detect(Species, "Montipora digitata") ~ "M. digitata",
                            str_detect(Species, "Stylophora pistillata") ~ "S. pistillata",
                            str_detect(Species, "Pocillopora verrucosa") ~ "P. verrucosa",
                            str_detect(Species, "Acropora muricata") ~ "A. muricata"))%>%
  
  mutate (Species_short =  factor(Species_short, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata",
                                               "P. verrucosa", "S. pistillata"
  )))%>%
  
  mutate(Heat_stress =ifelse(Day =="0", "Pre-heat","Post-heat"))%>%
  mutate(across(c(ID, Prec, Trtm, Tank_n, Colony, Day, Frg, Region, Trtm, Heat_stress),
                as.factor),
         Heat_stress =  factor(Heat_stress, level = c("Pre-heat", "Post-heat"
         )),
         Prec_Day = factor(Prec:Heat_stress))%>%
  mutate_if(is.numeric, round, 3) %>%

  group_by(ID, Frg, Species ,Prec, Colony)%>%
  mutate(Delta_pair = lead(gray_sum) - gray_sum )%>%
  ungroup()

## Statistical test - effect of heat exposure ----

###  Kruskal-Wallis test ----
# Kruskal-Wallis test was used to test coral response in each preconditioning treatment to heat exposure within each species. 
# The test was performed on paired data from the post-preconditioning and post-heat time points

#### "Heat" treatment ----
Heat_Clr_kruskal <- Clr_heat_paired %>%
  filter(Trtm=="Heat")%>%
  
  group_by(Species) %>%
  rstatix::kruskal_test(gray_sum ~ Heat_stress) %>%
  mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                            symbols = c("***", "**", "*", "ns"))

#### "Control" treatment ----
Heat_Clr_kruskal_control <- Clr_heat_paired %>%
  filter(Trtm=="Control")%>%
  
  group_by(Species) %>%
  rstatix::kruskal_test(gray_sum ~ Heat_stress) %>%
  mutate_if(is.numeric, round, 3) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                            symbols = c("***", "**", "*", "ns"))


### Wilcox test ----

#### "Heat" post-hoc test ----
Heat_Clr_Wilcox <- Clr_heat_paired %>%
  
  filter(Trtm=="Heat")%>%
  group_by(Species_short, Prec) %>%
  rstatix::pairwise_wilcox_test(gray_sum ~ Prec_Day, p.adjust.method = "BH") %>%
  
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))%>%
  mutate_if (is.numeric, round, 6)%>%
  
  rstatix::add_xy_position(x = "Prec_Day", dodge = 0.8)

#save the results
#write_csv2(Heat_PAM_Wilcox,here("Statistics","Wilcox_Heat_PAM.result.csv"))


#### "Control" post-hoc test ---- 
Heat_Clr_Wilcox_control <- Clr_heat_paired %>%
  
  filter(Trtm=="Control")%>%
  mutate(Prec_Day = factor(Prec:Heat_stress))%>%
  group_by(Species_short, Prec) %>%
  rstatix::pairwise_wilcox_test(gray_sum ~ Prec_Day, p.adjust.method = "BH") %>%
  
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", "ns"))%>%
  mutate_if (is.numeric, round, 6)%>%
  
  rstatix::add_xy_position(x = "Prec_Day", dodge = 0.8)








##                            PLOT - POST-HEAT paired values ----
####################################################################################################################### #
#plot showing the coral response to heat (post-heat) in comparison to the paired values at the "post-preconditioning" time point

#### Plot - Paired "Heat" ----

#paired data of the "Heat" treatment


Clr_heat_paired%>%
  mutate (Species_short =  factor(Species_short, level = c("G. fascicularis", "P. rus", "A. muricata", "M. digitata", "P. verrucosa", "S. pistillata"
  )))%>%
  dplyr::filter(Trtm=="Heat")%>%
  ggplot(mapping = aes(x = Prec_Day, y = gray_sum))+
  
  geom_quasirandom (aes( fill = Prec_Day, shape = Heat_stress ), 
                    col = "black",dodge.width= .9, size = 2 )+
  geom_boxplot(aes(fill = Prec_Day, shape = Heat_stress), alpha = 0.7,lwd = 0.7)+ 
  
  scale_fill_manual(values=c("#4682B4", "#63B8FF",  "#B4B446", "#EEEE00",  "#D4711C", "#FF7F24" ))+
  
  scale_shape_manual(values=c(21,24))+
  stat_pvalue_manual(Heat_Clr_Wilcox,
                     label = "p.adj.signif", label.size = 14, tip.length = 0.01, hide.ns = T
  ) +

  theme_classic() +
  
  ylab("Tissue color intensity") +
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
  
  facet_wrap(.~ Species_short, ncol = 6, strip.position ="bottom", )

ggsave(here ("Plot", "Clr_post-heat_paired-heat.png"), dpi = 400, units = "cm", width = 53.5, height = 14)



#### Plot - Paired "Control" ----
#paired data of the "Control" treatment


Clr_heat_paired%>%

  
  dplyr::filter(Trtm=="Control")%>%
  ggplot(mapping = aes(x = Prec_Day, y = gray_sum))+
  
  geom_quasirandom (aes( fill = Prec_Day, shape = Heat_stress ), 
                    col = "black",dodge.width= .9, size = 2 )+
  geom_boxplot(aes(fill = Prec_Day, shape = Heat_stress), alpha = 0.7,lwd = 0.7)+ 
  
  scale_fill_manual(values=c("#4682B4", "#63B8FF",  "#B4B446", "#EEEE00",  "#D4711C", "#FF7F24" ))+
  
  scale_shape_manual(values=c(21,24))+
  stat_pvalue_manual(Heat_Clr_Wilcox_control,
                     label = "p.adj.signif", label.size = 14, tip.length = 0.01, hide.ns = T
  ) +
  
  theme_classic() +
  
  ylab("Tissue color intensity") +
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
  
  facet_wrap(.~ Species_short, ncol = 6, strip.position ="bottom")

ggsave(here ("Plot", "Clr_post-heat_paired-control.png"), dpi = 400, units = "cm", width = 53.5, height = 14)








##                          THERMAL TOLERANCE CHANGES (Linear Mixed-Effects Models) ----
################################################################################################################## #


# The effect of preconditioning was evaluated comparing the paired delta difference 
# (Post-heat minus Post-preconditioning) of each preconditioning regime, separately for each species





## Heat - MONTIPORA DIGITATA (Mdi) ----

#### Subset ----

Mdi.Clr_heat_delta <- Clr_heat_paired %>%
  dplyr::filter(ID == "Mdi", Trtm =="Heat")%>%
  drop_na(Delta_pair)%>%
  mutate(Delta_pair1 = (-1*(Delta_pair)))# multiply per -1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Mdi.Clr_heat_delta$Delta_pair1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <- bestNormalize::boxcox(Mdi.Clr_heat_delta$Delta_pair1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Mdi.Clr_heat_delta <- bind_cols(Mdi.Clr_heat_delta, norm_dat)

#### Models ----

# Data are transformed using the Box-Cox Transformation 
Mdi.Clr.Heat.lmer <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony), data= Mdi.Clr_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Mdi.Clr.Heat.lmer)) 

# Residual normality
check_normality(Mdi.Clr.Heat.lmer)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Mdi.Clr.Heat.lmer)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Mdi.Clr.Heat.lmer)

## Residual plot
plot(residuals(Mdi.Clr.Heat.lmer) ~ fitted(Mdi.Clr.Heat.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Mdi.Clr.Heat.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Mdi.Clr.Heat.lmer)
fitted_lmer <- fitted(Mdi.Clr.Heat.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals

#### F statistic (ANOVA)----
Mdi.Clr.Heat.lmer_anova<-as_tibble(anova(Mdi.Clr.Heat.lmer))
Mdi.Clr.Heat.lmer_anova # ANOVA test is not significant


#### Post-hoc analysis ----

Mdi.Clr.Heat.lmer.ems <- emmeans::emmeans(Mdi.Clr.Heat.lmer, list(pairwise ~ Prec), adjust = "bonferroni")

(Mdi.Clr.Heat.lmer.ems.results <- Mdi.Clr.Heat.lmer.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())

# summary and report for better interpretation
summary(Mdi.Clr.Heat.lmer)
report(Mdi.Clr.Heat.lmer)


## Heat - GALAXEA FASCICULARIS (Gfa) ----

#### Subset ----

Gfa.Clr_heat_delta <- Clr_heat_paired %>%
  dplyr::filter(ID == "Gfa", Trtm =="Heat")%>%
  drop_na(Delta_pair)%>%
  mutate(Delta_pair1 = (-1*(Delta_pair)))# multiply per -1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Gfa.Clr_heat_delta$Delta_pair1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <- bestNormalize::boxcox(Gfa.Clr_heat_delta$Delta_pair1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Gfa.Clr_heat_delta <- bind_cols(Gfa.Clr_heat_delta, norm_dat)

#### Models ----

# Data are transformed using the Box-Cox Transformation 
Gfa.Clr.Heat.lmer <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Gfa.Clr_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Gfa.Clr.Heat.lmer)) 

# Residual normality
check_normality(Gfa.Clr.Heat.lmer)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Gfa.Clr.Heat.lmer)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Gfa.Clr.Heat.lmer)

## Residual plot
plot(residuals(Gfa.Clr.Heat.lmer) ~ fitted(Gfa.Clr.Heat.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Gfa.Clr.Heat.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Gfa.Clr.Heat.lmer)
fitted_lmer <- fitted(Gfa.Clr.Heat.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals

#### F statistic (ANOVA)----
Gfa.Clr.Heat.lmer_anova<-as_tibble(anova(Gfa.Clr.Heat.lmer))
Gfa.Clr.Heat.lmer_anova # ANOVA test is not significant


#### Post-hoc analysis ----

Gfa.Clr.Heat.lmer.ems <- emmeans::emmeans(Gfa.Clr.Heat.lmer, list(pairwise ~ Prec), adjust = "bonferroni")

(Gfa.Clr.Heat.lmer.ems.results <- Gfa.Clr.Heat.lmer.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())

# summary and report for better interpretation
summary(Gfa.Clr.Heat.lmer)
report(Gfa.Clr.Heat.lmer)


## Heat - PORITES RUS (Pru) ----

#### Subset ----

Pru.Clr_heat_delta <- Clr_heat_paired %>%
  dplyr::filter(ID == "Pru", Trtm =="Heat")%>%
  #drop_na(Delta_pair)%>%
  mutate(Delta_pair1 = (-1*(Delta_pair-5)))# multiply per -1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Pru.Clr_heat_delta$Delta_pair1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <- bestNormalize::boxcox(Pru.Clr_heat_delta$Delta_pair1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Pru.Clr_heat_delta <- bind_cols(Pru.Clr_heat_delta, norm_dat)

#### Models ----

# Data are transformed using the Box-Cox Transformation 
Pru.Clr.Heat.lmer <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Pru.Clr_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Pru.Clr.Heat.lmer)) 

# Residual normality
check_normality(Pru.Clr.Heat.lmer)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Pru.Clr.Heat.lmer)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Pru.Clr.Heat.lmer)

## Residual plot
plot(residuals(Pru.Clr.Heat.lmer) ~ fitted(Pru.Clr.Heat.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Pru.Clr.Heat.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Pru.Clr.Heat.lmer)
fitted_lmer <- fitted(Pru.Clr.Heat.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals

#### F statistic (ANOVA)----
Pru.Clr.Heat.lmer_anova<-as_tibble(anova(Pru.Clr.Heat.lmer))
Pru.Clr.Heat.lmer_anova # ANOVA test is not significant


#### Post-hoc analysis ----

Pru.Clr.Heat.lmer.ems <- emmeans::emmeans(Pru.Clr.Heat.lmer, list(pairwise ~ Prec), adjust = "bonferroni")

(Pru.Clr.Heat.lmer.ems.results <- Pru.Clr.Heat.lmer.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())

# summary and report for better interpretation
summary(Pru.Clr.Heat.lmer)
report(Pru.Clr.Heat.lmer)



## Heat - ACROPORA MURICATA (Amu) ----

#### Subset ----

Amu.Clr_heat_delta <- Clr_heat_paired %>%
  dplyr::filter(ID == "Amu", Trtm =="Heat")%>%
  #drop_na(Delta_pair)%>%
  mutate(Delta_pair1 = (-1*(Delta_pair)))# multiply per -1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Amu.Clr_heat_delta$Delta_pair1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <- bestNormalize::boxcox(Amu.Clr_heat_delta$Delta_pair1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Amu.Clr_heat_delta <- bind_cols(Amu.Clr_heat_delta, norm_dat)

#### Models ----

# Data are transformed using the Box-Cox Transformation 
Amu.Clr.Heat.lmer <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Amu.Clr_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Amu.Clr.Heat.lmer)) 

# Residual normality
check_normality(Amu.Clr.Heat.lmer)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Amu.Clr.Heat.lmer)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Amu.Clr.Heat.lmer)

## Residual plot
plot(residuals(Amu.Clr.Heat.lmer) ~ fitted(Amu.Clr.Heat.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Amu.Clr.Heat.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Amu.Clr.Heat.lmer)
fitted_lmer <- fitted(Amu.Clr.Heat.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals

#### F statistic (ANOVA)----
Amu.Clr.Heat.lmer_anova<-as_tibble(anova(Amu.Clr.Heat.lmer))
Amu.Clr.Heat.lmer_anova # ANOVA test is not significant


#### Post-hoc analysis ----

Amu.Clr.Heat.lmer.ems <- emmeans::emmeans(Amu.Clr.Heat.lmer, list(pairwise ~ Prec), adjust = "bonferroni")

(Amu.Clr.Heat.lmer.ems.results <- Amu.Clr.Heat.lmer.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())

# summary and report for better interpretation
summary(Amu.Clr.Heat.lmer)
report(Amu.Clr.Heat.lmer)




## Heat - POCILLOPORA VERRUCOSA (Pve) ----

#### Subset ----

Pve.Clr_heat_delta <- Clr_heat_paired %>%
  dplyr::filter(ID == "Pve", Trtm =="Heat")%>%
  #drop_na(Delta_pair)%>%
  mutate(Delta_pair1 = (-1*(Delta_pair)))# multiply per -1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Pve.Clr_heat_delta$Delta_pair1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <- bestNormalize::boxcox(Pve.Clr_heat_delta$Delta_pair1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Pve.Clr_heat_delta <- bind_cols(Pve.Clr_heat_delta, norm_dat)

#### Models ----

# Data are transformed using the Box-Cox Transformation 
Pve.Clr.Heat.lmer <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Pve.Clr_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Pve.Clr.Heat.lmer)) 

# Residual normality
check_normality(Pve.Clr.Heat.lmer)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Pve.Clr.Heat.lmer)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Pve.Clr.Heat.lmer)

## Residual plot
plot(residuals(Pve.Clr.Heat.lmer) ~ fitted(Pve.Clr.Heat.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Pve.Clr.Heat.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Pve.Clr.Heat.lmer)
fitted_lmer <- fitted(Pve.Clr.Heat.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals

#### F statistic (ANOVA)----
Pve.Clr.Heat.lmer_anova<-as_tibble(anova(Pve.Clr.Heat.lmer))
Pve.Clr.Heat.lmer_anova # ANOVA test is not significant


#### Post-hoc analysis ----

Pve.Clr.Heat.lmer.ems <- emmeans::emmeans(Pve.Clr.Heat.lmer, list(pairwise ~ Prec), adjust = "bonferroni")

(Pve.Clr.Heat.lmer.ems.results <- Pve.Clr.Heat.lmer.ems$`pairwise differences of Prec` %>% as_tibble() %>% mutate_if (is.numeric, round, 5)%>%
    add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>% tibble())

# summary and report for better interpretation
summary(Pve.Clr.Heat.lmer)
report(Pve.Clr.Heat.lmer)



## Heat - STYLOPHORA PISTILLATA (Spi) ----

#### Subset ----

Spi.Clr_heat_delta <- Clr_heat_paired %>%
  dplyr::filter(ID == "Spi", Trtm =="Heat")%>%
  #drop_na(Delta_pair)%>%
  mutate(Delta_pair1 = (-1*(Delta_pair)))# multiply per -1 to make all values positive

### Transformation ----

bestN <- bestNormalize(Spi.Clr_heat_delta$Delta_pair1, loo = T, warn = T, allow_orderNorm = T, 
                       allow_lambert_s = T, allow_lambert_h = T, out_of_sample = T)
bestN

norm_res <- bestNormalize::boxcox(Spi.Clr_heat_delta$Delta_pair1)

norm_dat <- as_tibble(norm_res$x.t) %>%
  rename(Delta_pair_norm = value)

Spi.Clr_heat_delta <- bind_cols(Spi.Clr_heat_delta, norm_dat)

#### Models ----

# Data are transformed using the Box-Cox Transformation 
Spi.Clr.Heat.lmer <- lmerTest::lmer(Delta_pair_norm ~ Prec + (1|Colony) + (1|Tank_n), data= Spi.Clr_heat_delta, REML = T)

# Residual plot
qqPlot(residuals(Spi.Clr.Heat.lmer)) 

# Residual normality
check_normality(Spi.Clr.Heat.lmer)# the significant p value is due to the two dead corals

# Check homogeneity
check_homogeneity(Spi.Clr.Heat.lmer)# the significant p value is due to the two dead corals

#Check heteroscedasticity
check_heteroscedasticity(Spi.Clr.Heat.lmer)

## Residual plot
plot(residuals(Spi.Clr.Heat.lmer) ~ fitted(Spi.Clr.Heat.lmer), main = "Residuals vs Fitted")
abline(h = 0, col = "red")

## Histogram of residuals
hist(resid(Spi.Clr.Heat.lmer), main = "Histogram of residuals")

## Breusch-Pagan test
residuals_lmer <- residuals(Spi.Clr.Heat.lmer)
fitted_lmer <- fitted(Spi.Clr.Heat.lmer)
lm_test <- lm(residuals_lmer ~ fitted_lmer)
bptest(lm_test) # the significant p value is due to the two dead corals

#### F statistic (ANOVA)----
Spi.Clr.Heat.lmer_anova<-as_tibble(anova(Spi.Clr.Heat.lmer))
Spi.Clr.Heat.lmer_anova # ANOVA test is not significant


# summary and report for better interpretation
summary(Spi.Clr.Heat.lmer)
report(Spi.Clr.Heat.lmer)
