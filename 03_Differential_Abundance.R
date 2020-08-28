# Differential abundance metrics

# packages ####
library(tidyverse)
library(patchwork)
library(phyloseq)
library(corncob)
library(vegan)
library(broom)
library(purrr)

#functions
source("./plot_bar2.R")
source("./bbdml_helper.R")

pal.discrete <- c("#c1593c","#688e52","#643d91","#894e7d","#477887","#12aa91","#705f36","#8997b2","#c4a113",
                  "#753c2b","#3c3e44","#b3bf2d","#82b2a4","#820616","#a17fc1","#262a8e","#abb5b5","#000000",
                  "#493829","#816C5B","#A9A18C","#613318","#855723","#B99C6B","#8F3B1B","#D57500","#DBCA69",
                  "#404F24","#668D3C","#BDD09F","#4E6172","#83929F","#A3ADB8")


# Load phyloseq data ####
ps_FF <- readRDS("./output/final_phyloseq_object.RDS")

names(sample_data(ps_FF))

# glom taxa at various levels
ps_order <- tax_glom(ps_FF,"Order")
ps_family <- tax_glom(ps_FF,"Family")
ps_genus <- tax_glom(ps_FF,"Genus")



######################## Burned vs Unburned #########################


# Order-level ##

set.seed(123)
da_analysis <- differentialTest(formula = ~ FireTreatment, #abundance
                                phi.formula = ~ FireTreatment, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_order,
                                fdr_cutoff = 0.05)

fire_bbdml <- multi_bbdml(da_analysis,
                          ps_object = ps_family,
                          mu_predictor = "FireTreatment",
                          phi_predictor = "FireTreatment",
                          taxlevels = 2:7)
length(fire_bbdml)

plot_multi_bbdml(bbdml_list = fire_bbdml,
                 color = "FireTreatment",
                 pointsize = 3)

# Save figure of DA taxa
ls(pattern = "bbdml_plot_")
bbdml_plot_1 / bbdml_plot_2 / bbdml_plot_3 / bbdml_plot_4
ggsave("./output/DA_Plots_FireTreatment_Order.png",height = 8,width = 14,dpi=300)

# Save model summaries for individual taxa
sink("./output/bbdml_summaries_FireTreatment_Order.txt")
fire_bbdml
sink(NULL)

rm(list = ls(pattern = "bbdml_plot_"))


# family-level ##

set.seed(123)
da_analysis <- differentialTest(formula = ~ FireTreatment, #abundance
                                phi.formula = ~ FireTreatment, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_family,
                                fdr_cutoff = 0.05)

fire_bbdml <- multi_bbdml(da_analysis,
                          ps_object = ps_family,
                          mu_predictor = "FireTreatment",
                          phi_predictor = "FireTreatment",
                          taxlevels = 2:7)

length(fire_bbdml)

plot_multi_bbdml(bbdml_list = fire_bbdml,
                 color = "FireTreatment",
                 pointsize = 3)

# Save figure of DA taxa
ls(pattern = "bbdml_plot_")
bbdml_plot_1 / bbdml_plot_2 / bbdml_plot_3 / bbdml_plot_4
ggsave("./output/DA_Plots_FireTreatment_Family.png",height = 8,width = 14,dpi=300)

# Save model summaries for individual taxa
sink("./output/bbdml_summaries_FireTreatment_Family.txt")
fire_bbdml
sink(NULL)

rm(list = ls(pattern = "bbdml_plot_"))

# Genus-level ##

set.seed(123)
da_analysis <- differentialTest(formula = ~ FireTreatment, #abundance
                                phi.formula = ~ FireTreatment, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

fire_bbdml <- multi_bbdml(da_analysis,
                          ps_object = ps_family,
                          mu_predictor = "FireTreatment",
                          phi_predictor = "FireTreatment",
                          taxlevels = 2:7)

length(fire_bbdml)

plot_multi_bbdml(bbdml_list = fire_bbdml,
                 color = "FireTreatment",
                 pointsize = 3)

# Save figure of DA taxa
ls(pattern = "bbdml_plot_")
bbdml_plot_1
ggsave("./output/DA_Plots_FireTreatment_Genus.png",height = 8,width = 14,dpi=300)

# Save model summaries for individual taxa
sink("./output/bbdml_summaries_FireTreatment_Genus.txt")
fire_bbdml
sink(NULL)

rm(list = ls(pattern = "bbdml_plot_"))


######################## Temperature Groups #########################


# Genus-level ##

set.seed(123)
da_analysis <- differentialTest(formula = ~ TempGroup, #abundance
                                phi.formula = ~ TempGroup, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_order,
                                fdr_cutoff = 0.05)
plot(da_analysis)
da_analysis$significant_models


fire_bbdml <- multi_bbdml(da_analysis,
                          ps_object = ps_family,
                          mu_predictor = "TempGroup",
                          phi_predictor = "TempGroup",
                          taxlevels = 2:7)

plot_multi_bbdml(bbdml_list = fire_bbdml,
                 color = "TempGroup",
                 pointsize = 3)

# Save figure of DA taxa
ls(pattern = "bbdml_plot_")
bbdml_plot_1 / bbdml_plot_2 / bbdml_plot_3 / bbdml_plot_4
# ggsave("./output/DA_Plots_TempGroup_Order.png",height = 8,width = 14,dpi=300)

# Save model summaries for individual taxa
sink("./output/bbdml_summaries_TempGroup_Order.txt")
fire_bbdml
sink(NULL)

rm(list = ls(pattern = "bbdml_plot_"))

# family-level ##


set.seed(123)
da_analysis <- differentialTest(formula = ~ TempGroup, #abundance
                                phi.formula = ~ TempGroup, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_family,
                                fdr_cutoff = 0.05)
plot(da_analysis) + theme(axis.text.y = element_text(face="bold"))


fire_bbdml <- multi_bbdml(da_analysis,
                          ps_object = ps_family,
                          mu_predictor = "TempGroup",
                          phi_predictor = "TempGroup",
                          taxlevels = 2:7)

plot_multi_bbdml(bbdml_list = fire_bbdml,
                 color = "TempGroup",
                 pointsize = 3)

# Save figure of DA taxa
ls(pattern = "bbdml_plot_")
(bbdml_plot_1 + bbdml_plot_2 + bbdml_plot_3 + bbdml_plot_4) / (bbdml_plot_5 + bbdml_plot_6 + bbdml_plot_7 + bbdml_plot_8) / 
  (bbdml_plot_9 + bbdml_plot_10 + bbdml_plot_11 + bbdml_plot_12) / (bbdml_plot_13 + bbdml_plot_14)
ggsave("./output/DA_Plots_TempGroup_Family.png",height = 18,width = 24,dpi=300)

# Save model summaries for individual taxa
sink("./output/bbdml_summaries_TempGroup_Family.txt")
fire_bbdml
sink(NULL)

rm(list = ls(pattern = "bbdml_plot_"))

# Genus-level ##


set.seed(123)
da_analysis <- differentialTest(formula = ~ TempGroup, #abundance
                                phi.formula = ~ TempGroup, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

fire_bbdml <- multi_bbdml(da_analysis,
                          ps_object = ps_family,
                          mu_predictor = "TempGroup",
                          phi_predictor = "TempGroup",
                          taxlevels = 2:7)

plot_multi_bbdml(bbdml_list = fire_bbdml,
                 color = "TempGroup",
                 pointsize = 3)

# Save figure of DA taxa
ls(pattern = "bbdml_plot_")
bbdml_plot_1 / bbdml_plot_2 / bbdml_plot_3 / bbdml_plot_4
ggsave("./output/DA_Plots_TempGroup_Genus.png",height = 8,width = 14,dpi=300)

# Save model summaries for individual taxa
sink("./output/bbdml_summaries_TempGroup_Genus.txt")
fire_bbdml
sink(NULL)

rm(list = ls(pattern = "bbdml_plot_"))

# DA analysis for each burnyear group, based on 
set.seed(123)
da_analysis <- differentialTest(formula = ~ TempGroup, #abundance
                                phi.formula = ~ TempGroup, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)

fire_bbdml <- multi_bbdml(da_analysis,
                          ps_object = ps_family,
                          mu_predictor = "TempGroup",
                          phi_predictor = "TempGroup",
                          taxlevels = 2:7)

plot_multi_bbdml(bbdml_list = fire_bbdml,
                 color = "TempGroup",
                 pointsize = 3)




#########################  Stacked barcharts ##########################

# Merge samples by fire treatment and temp group
ps_FF@sam_data$TempGroupName <- 
ps_FF@sam_data$TempGroup %>% as.character() %>% 
  str_replace(pattern = "A: -7.9:-6.4", replacement = "Low") %>%
  str_replace(pattern = "B: -6.4:-4.9", replacement = "Med") %>%
  str_replace(pattern = "C: -4.9:-3.4", replacement = "High")
ps_FF@sam_data$TempGroupName <- factor(ps_FF@sam_data$TempGroupName,levels=c("Low","Med","High"))
NewPastedVar <- paste(ps_FF@sam_data$FireTreatment,ps_FF@sam_data$TempGroupName,sep = "_")
ps_FF@sam_data$FireAndTempGroup <- NewPastedVar
ps_FF_merged <- ps_FF %>% merge_samples("FireAndTempGroup")

# repair metadata
ps_FF_merged@sam_data$FireTreatment <- unlist(purrr::map(str_split(sample_names(ps_FF_merged),"_"),1))
ps_FF_merged@sam_data$TempGroupName <- unlist(purrr::map(str_split(sample_names(ps_FF_merged),"_"),2))
ps_FF_merged@sam_data$TempGroupName <- factor(ps_FF_merged@sam_data$TempGroupName,levels=c("Low","Med","High"))

sample_data(ps_FF_merged)

p <- ps_FF_merged %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(x="TempGroupName",fill = "Phylum") +
  scale_fill_manual(values=pal.discrete) + 
  labs(x="Annual mean temperature group",y="Relative abundance") 

p + facet_wrap(~FireTreatment) + theme(strip.background = element_blank(),
                                       strip.text = element_text(size=12,face="bold"),
                                       axis.title = element_text(size=12,face="bold"),
                                       axis.text.x = element_text(face="bold"),
                                       legend.title = element_text(size=12,face="bold"))
ggsave("./output/Barplot_Phylum_by_Burn_and_tempgroup.png",dpi=300,height = 8,width = 8)


# merge samples by fire treatment and burnyear

YSB <- 2020 - ps_FF@sam_data$BurnYear
ps_FF@sam_data$YearsSinceBurn <- YSB

NewPastedVar <- paste(ps_FF@sam_data$FireTreatment,ps_FF@sam_data$YearsSinceBurn,sep = "_")
ps_FF@sam_data$FireAndYearsSinceBurn <- NewPastedVar
ps_FF_merged <- ps_FF %>% merge_samples("FireAndYearsSinceBurn")

# repair metadata and convert new variable to factor
ps_FF_merged@sam_data$FireTreatment <- unlist(purrr::map(str_split(sample_names(ps_FF_merged),"_"),1))
ps_FF_merged@sam_data$YearsSinceBurn <- unlist(purrr::map(str_split(sample_names(ps_FF_merged),"_"),2))
ps_FF_merged@sam_data$YearsSinceBurn <- factor(ps_FF_merged@sam_data$YearsSinceBurn,levels = c("5","8","13","14","17","20"))

p <- ps_FF_merged %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(x="YearsSinceBurn",fill = "Phylum") +
  scale_fill_manual(values=pal.discrete) + 
  labs(x="Years Since Burn",y="Relative abundance") 

p + facet_wrap(~FireTreatment) + theme(strip.background = element_blank(),
                                       strip.text = element_text(size=12,face="bold"),
                                       axis.title = element_text(size=12,face="bold"),
                                       axis.text.x = element_text(face="bold"),
                                       legend.title = element_text(size=12,face="bold"))
ggsave("./output/Barplot_Phylum_by_Burn_and_YearsSinceBurn.png",dpi=300,height = 8,width = 8)


# is temperature correlated with burnyear?
ggplot(meta, aes(x=ANNUAL_MEAN_TEMP,y=BurnYear)) + geom_smooth() + geom_point()
cor(exp(meta$ANNUAL_MEAN_TEMP),meta$BurnYear)


