# ITS Forest Soil Fungi Burn Recovery Exploratory analyses

# Load packages ####
library(phyloseq) 
library(tidyverse)
library(vegan)
library(ade4)
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthdata)
# library(maptools)
library(tools)
library(ggbiplot)
library(microbiome)
library(RColorBrewer)
library(dada2)
library(ecodist)
<<<<<<< HEAD
# library(ggmap)
# library(maps)
library(colorblindr)
library(lme4)
library(corncob)
library(broom)
library(patchwork)
library(microbiome)
library(ggpubr)
library(lme4)
library(lmerTest)

source("./palettes.R")
source("./plot_bar2.R")
palette_plot(pal.discrete)

=======
library(modelr)
>>>>>>> 951dbbb763cdfb68db506575459e058ae34d8b10
#Set ggplot theme
theme_set(theme_bw())

# import phyloseq object and metadata ####
ps_FF = readRDS("./output/phyloseq_object_ITS.RDS")
ps_ITS <- ps_FF

sample_data(ps_FF)

meta <- sample_data(ps_FF)
sample_names(ps_FF)

meta_df <- meta(ps_FF)


# examine non-fungi ####
taxa = ps_ITS@tax_table
kingdoms = c(unique(taxa[,1]))

# plot kingdom abundance ####

ggplot(as.data.frame(taxa), aes(x=Kingdom)) +
  geom_bar(stat = "count") + ggtitle("No. of taxa from each kingdom")

ggsave("./output/non-fungal_taxa_assignment_plot.png")  

# Remove non-fungi ####
ps_FF = subset_taxa(ps_FF, Kingdom == "k__Fungi")


# Clean up taxonomy names
for(i in 1:7){
  tax_table(ps_FF)[,i] <- str_remove(tax_table(ps_FF)[,i],".__")  
}

# clean up empty samples and ESVs
summary(sample_sums(ps_FF))
summary(taxa_sums(ps_FF))
ps_FF <- subset_taxa(ps_FF,taxa_sums(ps_FF) > 1)
ps_FF <- subset_samples(ps_FF,sample_sums(ps_FF) > 1)



# Plot diversity ####
plot_richness(ps_FF,x="Set_ID",measures = "Shannon") + 
  stat_smooth() +
  facet_grid(~FireTreatment) + labs(y="Shannon diversity")
ggsave("./output/alpha_diversity_scatterplot.png")

ps_FF@sam_data$Shannon <- estimate_richness(ps_FF, measures="Shannon")$Shannon

sam = as.data.frame(ps_FF@sam_data)
sam = (as(sam, "data.frame"))

names(sample_data(ps_FF))
ps_FF@sam_data$Richness <- specnumber(otu_table(ps_FF))


# Add grouping factor for mean annual temp
temps <- ps_FF@sam_data$ANNUAL_MEAN_TEMP
tempgroups <- cut(temps,3)
tempgroups <- as.character(tempgroups) %>% str_replace("\\(-7.9,-6.4]","A: -7.9:-6.4") %>%
  str_replace("\\(-6.4,-4.9]","B: -6.4:-4.9") %>%
  str_replace("\\(-4.9,-3.4]","C: -4.9:-3.4") %>%
  factor(ordered = FALSE)
ps_FF@sam_data$TempGroup <- tempgroups



# export ps_FF
saveRDS(ps_FF,"./output/final_phyloseq_object.RDS")
#



# Basic plots of diversity for overall data
# stacked boxplots x-axis factor(years since burn), y-axis relative-abundance
ps_FF %>% merge_samples(group = "FireTreatment") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Phylum") + scale_fill_manual(values = pal.discrete) + labs(y="Relative abundance",x="Site type")
ggsave("./output/Barplot_Phylum_mean_by_burn_treatment.png",dpi=300,height = 6,width = 6)  


# Alph diversity vs fire treatment boxplot
ggplot(sam, aes(x=FireTreatment,y=Shannon,fill=FireTreatment)) +   geom_boxplot() + geom_jitter(width = .1) + labs(y="Shannon diversity") +
  scale_fill_manual(values = pal.discrete[c(2,6)])
ggsave("./output/Shannon_Div_by_treatment.png", dpi=300)  





mod1 = aov(Shannon ~ Set_ID * FireTreatment, data=sam)
summary(mod1)

sam = mutate(sam, YearsSinceBurn = 2019-BurnYear)
sam$Richness <- specnumber(otu_table(ps_FF))





ggplot(sam, aes(x=YearsSinceBurn,y=Shannon,color=FireTreatment)) +
  geom_jitter() + stat_smooth(method="lm", se = FALSE) + labs(x="Years Since Burn",y="Shannon Diversity") +
  scale_x_reverse()
ggsave("./output/Shannon_div_raw_over_time.png", dpi=300)

names(sam)
# richness vs time since fire (as a factor), colored by fire treatment
sam %>% 
  group_by(YearsSinceBurn, FireTreatment) %>%
  dplyr::summarize(N=n(),Mean_Shannon_Div=mean(Shannon),Upper_Shannon=Mean_Shannon_Div+sd(Shannon),Lower_Shannon=Mean_Shannon_Div-sd(Shannon)) %>%
  ggplot(aes(x=factor(YearsSinceBurn),y=Mean_Shannon_Div,ymin=Lower_Shannon,ymax=Upper_Shannon,color=FireTreatment)) +
  geom_errorbar() + geom_line()

summarize
sam %>%
ggplot(aes(x=factor(YearsSinceBurn),y=Shannon,color=FireTreatment)) + 
  geom_boxplot() + geom_line(aes(group=factor(YearsSinceBurn))) 

sam %>%
  ggplot(aes(x=factor(YearsSinceBurn),y=Richness,color=FireTreatment)) + 
  geom_boxplot() + geom_line(aes(group=factor(YearsSinceBurn))) 

  
  
  # Remove low-abundance taxa and empty samples ####


# quick barplot
ps_FF %>%
  merge_samples(group = "FireTreatment") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill = "Phylum") + theme_bw() + labs(y="Relative abundance") + scale_fill_manual(values=pal.discrete)
ggsave("./output/barplot_phylum_burntreatment.png", dpi = 300)

# merge samples by burn and group ####
newpastedvar = paste(sample_data(ps_FF)$FireTreatment, sample_data(ps_FF)$Set_ID, sep = "_")
sample_data(ps_FF)$NewPastedVar = newpastedvar
firelevels = levels(sample_data(ps_FF)$FireTreatment)
setlevels = levels(sample_data(ps_FF)$Set_ID)
lat = unique(sample_data(ps_FF)$Latitude)
lon = unique(sample_data(ps_FF)$Longitude)
latlon = as.data.frame(cbind(ps_FF@sam_data$NewPastedVar, as.numeric(ps_FF@sam_data$Latitude), as.numeric(ps_FF@sam_data$Longitude)))


psm = merge_samples(ps_FF, "NewPastedVar")

# repair values
psm@sam_data$FireTreatment <- c(rep(firelevels, each = 10),firelevels[2])
psm@sam_data$Set_ID <- setlevels[psm@sam_data$Set_ID]

# normalize (relabund) ####
psm_ra = transform_sample_counts(psm, function(OTU) OTU/sum(OTU) )

psm_ra = subset_taxa(psm_ra, Class != "c__NA|NA")
psm_ra = subset_taxa(psm_ra, Class != "c__NA")
plot_bar2(psm_ra, x="BurnYear",fill = "Class") + facet_wrap(~Class) + theme_bw()
ggsave("./output/Class_by_burnyear.png", dpi=300, width = 12,height = 10)


# Fix longitude
sample_data(psm_ra)$Longitude[which(sample_data(psm_ra)$Longitude > -100)] <- -111.33787



# Normalize full ps object ####
ps_ra = transform_sample_counts(ps_FF, function(OTU) OTU/sum(OTU) )


# Find community distance between burned and non-burned within each set ####

# find and remove any sets that don't have both burned and unburned samples remaining

sets = as.character(purrr::map(strsplit(unique(ps_ra@sam_data$NewPastedVar), split = "_"),3))
sets = as.data.frame(table(sets))
goodsets = as.character(sets[which(sets$Freq > 1),1])
goodsets = goodsets[order(as.numeric(goodsets))]
goodsets = paste0("Set_",goodsets)


# set vectors and counter
set = c()
meandist = c()
sddist = c()
x=1

# For-loop calculates mean and stdev between burned and unburned samples in each set ####
for(i in goodsets){
  ps = subset_samples(ps_ra, Set_ID == i)
  psb = subset_samples(ps, FireTreatment == "Burn")
  psn = subset_samples(ps, FireTreatment == "NonBurn")
  
  dist = as.data.frame(as.matrix(vegdist(otu_table(ps))))
  dist_fire = dist[sample_names(psb),sample_names(psn)]
  mean = mean(c(as.matrix(dist_fire)))
  sdev = sd(c(as.matrix(dist_fire)))
  
  set[x] <- i
  meandist[x] <- mean
  sddist[x] <- sdev
  x=x+1
}


setyear = as.data.frame(unique(cbind(as.character(ps_ra@sam_data$Set_ID),ps_ra@sam_data$BurnYear)))
setyear = setyear[setyear$V1 %in% set,]
setyear = setyear[c(3,4,5,6,7,8,9,10,1,2),]



# build data frame and plot ####
dist.df = data.frame(Set = set, MeanDist = meandist, StDevDist = sddist, BurnYear = as.numeric(as.character(setyear$V2)))
dist.df$upper = dist.df$MeanDist + dist.df$StDevDist
dist.df$lower = dist.df$MeanDist - dist.df$StDevDist

n=1
meantemps = c()
meanprecips = c()
meanslopes = c()
meancanopys = c()
for(i in levels(meta$Set_ID)[c(1,5,6,7,8,9,10,11,2,3)]) {
  df = meta[meta$Set_ID == i,]
  meantemps[n] = mean(df[,"ANNUAL_MEAN_TEMP"])
  meanprecips[n] = mean(df[,"ANNUAL_PRECIP"])
  meanslopes[n] = mean(df[,"SLOPE_AVG"])
  meancanopys[n] = mean(df[,"CANOPY_2001"])
n=n+1
}
dist.df$Temp <- meantemps 
dist.df$Precip <- meanprecips
dist.df$Slope <- meanslopes
dist.df$CanopyCover <- meancanopys


ggplot(dist.df, aes(x=(2019-BurnYear),y=MeanDist,ymin=lower,ymax=upper)) +
  geom_point() + 
  geom_errorbar() +
  stat_smooth(method = "lm") +
  theme_bw() +
  labs(y="Community Distance", x="Years Since Burn")
ggsave("./output/Community_Distance_vs_Burn_Year.png", dpi = 300)


# Add community distances to metadata
names(dist.df)[1] <- "Set_ID"
meta <- as(meta,"data.frame")
meta <- as.data.frame(left_join(meta,dist.df,by=c("Set_ID","BurnYear")))
glimpse(meta)
# names(meta)
# 
# from = (meta$Set_ID)[c(1,5,6,7,8,9,10,11,2,3,4)]
# to = c(dist.df$MeanDist,NA)
# distance = plyr::mapvalues(meta$Set_ID,from=from,to=to)
# meta$CommDist_Burn = as.numeric(as.character(distance))

# Models with GDC data
meta$YearsSinceBurn <-  2019 - meta$BurnYear
mod2 = glm(MeanDist ~ (ANNUAL_MEAN_TEMP * YearsSinceBurn), data = meta)
summary(mod2)
meta$BurnYear
mod3 = lmer(MeanDist ~ (ANNUAL_MEAN_TEMP) + (1|YearsSinceBurn), data = meta)
summary(mod3)

mod4 <- lmerTest::lmer(MeanDist ~ (ANNUAL_MEAN_TEMP) + (1|YearsSinceBurn), data = meta)
summary(mod4)

sink("./output/lmer_mod_of_Paired_Distance.txt")
summary(mod4)
sink(NULL)


ggplot(meta,aes(x=ANNUAL_MEAN_TEMP,y=MeanDist)) +
  geom_jitter(aes(color=factor(2019-BurnYear)),height = .01, width = 0,size=3,alpha=.75) + geom_smooth(color="Black",se=FALSE,method = "lm") + 
  # geom_boxplot(aes(group=(ANNUAL_MEAN_TEMP),fill=factor(2019-BurnYear)),alpha=.25)
  labs(x="Annual Mean Temperature",y="Fungal Community Distance",color = "Years Since Burn") + scale_color_viridis_d() +
  theme(axis.title = element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12,face="bold"),
        legend.title = element_text(size=14,face="bold"),
        legend.text = element_text(size=14,face="bold"),
        legend.position = "right") + scale_y_continuous(limits = c(0.5,1))

ggsave("./output/Community_Distance_vs_AnnualMeanTemp_partial-lims_lm.png",dpi=300,width = 10,height = 10)

cor(meta$ANNUAL_MEAN_TEMP,meta$BurnYear)

ggplot(meta,aes(x=ANNUAL_PRECIP,y=CommDist_Burn)) +
  geom_jitter(aes(color=2019-BurnYear),height = .01, width = 0,size=2,alpha=.5)+geom_smooth(method = "lm") +
  labs(x="Annual Average Precipitation (mm)",y="Fungal Community Distance",color = "Years Since Burn")
ggsave("./output/Community_Distance_vs_Precip_and_Burn_Year2.png",dpi=300,width = 10,height = 10)


ggplot(meta,aes(x=2019-BurnYear,y=CommDist_Burn)) +
  geom_jitter(aes(color=ANNUAL_MEAN_TEMP),height = .01, width = 0,size=2,alpha=.5)+geom_smooth(method="lm") +
  labs(x="Years Since Burn",y="Fungal Community Distance",color = "Annual Mean Temp (deg C)")

glimpse(meta)
# adonis(otu_table(ps_ra) ~ meta$BurnYear+meta$SLOPE_AVG+meta$FireTreatment+meta$ANNUAL_PRECIP)

# Ordinate ####
DCA = ordinate(ps_ra)
plot_ordination(ps_ra,DCA, color = "FireTreatment") + theme_bw()
ggsave("./output/DCA_ordination.png", dpi=300)

# PCoA 
pca = prcomp(as.matrix(otu_table(ps_ra)))

dim(sample_data(ps_ra))
dim(otu_table(ps_ra))

g <- ggbiplot(pca, groups = ps_ra@sam_data$Set_ID)
g

NMDS = ordinate(ps_ra, method = "NMDS")
plot_ordination(ps_ra,NMDS, color = "FireTreatment", type = "biplot") + theme_bw()




# PermANOVA ####
sink("./output/adonis_table.txt")
adonis(otu_table(ps_ra) ~ ps_ra@sam_data$BurnYear * ps_ra@sam_data$FireTreatment * ps_ra@sam_data$Set_ID)
adonis(otu_table(ps_ra) ~ (ps_ra@sam_data$BurnYear + ps_ra@sam_data$ANNUAL_MEAN_TEMP) * ps_ra@sam_data$Set_ID)


sink(NULL)

# Mantel Test and multiple regression on distance matrices ####
spatial.dist = dist(cbind(ps_ra@sam_data$Longitude, ps_ra@sam_data$Latitude))
comm.dist = vegdist(as.matrix(ps_ra@otu_table))

envir = cbind(ps_ra@sam_data$SLOPE_AVG,ps_ra@sam_data$CANOPY_2001,ps_ra@sam_data$ANNUAL_MEAN_TEMP,ps_ra@sam_data$ANNUAL_PRECIP)
envir.dist = vegdist(envir, method = "man")

mantel.test = mantel.rtest(spatial.dist, comm.dist, nrepet = 9999)

# MRM
dist_MRM <- MRM(comm.dist ~ spatial.dist + envir.dist, nperm = 9999)




sink("./output/mantel_test.txt")
print(mantel.test)
sink(NULL)

sink("./output/MRM_table.txt")
dist_MRM
sink(NULL)


# single groups ??? ####
sordaria = subset_taxa(psm_ra, Class == "c__Sordariomycetes")
plot_bar2(sordaria, x="FireTreatment", fill = "Family")



# Core members ####

# rename ASVs
ps_ra_rn <- ps_ra
taxa_count = length(colnames(ps_ra_rn@otu_table@.Data))

colnames(ps_ra_rn@otu_table@.Data) <- paste0("ASV_",1:taxa_count)



det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

plot_core(ps_ra_rn, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")

p <- plot_core(ps_ra_rn, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, horizontal = TRUE)
print(p)
ggsave("./output/Core_Taxa_Overall.png",dpi=300)

core_burned = subset_samples(ps_ra, FireTreatment == "Burn")
core_unburned = subset_samples(ps_ra, FireTreatment == "NonBurn")
core_burned = subset_taxa(core_burned,colSums(otu_table(core_burned)) > 0)
core_unburned = subset_taxa(core_unburned,colSums(otu_table(core_unburned)) > 0)

taxa_count = length(colnames(core_burned@otu_table@.Data))
colnames(core_burned@otu_table@.Data) <- paste0("ASV_",1:taxa_count)
taxa_count = length(colnames(core_unburned@otu_table@.Data))
colnames(core_unburned@otu_table@.Data) <- paste0("ASV_",1:taxa_count)

p <- plot_core(core_burned, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, horizontal = TRUE) + ggtitle("Burned Sites - Core Taxa")
print(p)
ggsave("./output/Core_Taxa_Burned.png", dpi=300)

p <- plot_core(core_unburned, plot.type = "heatmap", 
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, horizontal = TRUE) + ggtitle("Non-Burned Sites - Core Taxa")
print(p)
ggsave("./output/Core_Taxa_NonBurned.png", dpi=300)


core(otu_table(core_burned), detection = 0.01, prevalence = .2)
colSums(otu_table(core_burned))

# Mapping sites ####
shannon = diversity(t(otu_table(psm_ra)), "shannon")
simpson = diversity(t(otu_table(psm_ra)), "simpson")


world <- ne_countries(scale = "medium", returnclass = "sf")

sites = data.frame(latitude = sample_data(psm_ra)$Latitude, longitude = sample_data(psm_ra)$Longitude, 
                   diversity = diversity)

str(diversity)
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))
states$ID <- toTitleCase(states$ID)



ggplot(data = world) +
  geom_sf() +
  geom_sf(data = states, fill = "White") + 
  geom_text(data = states, aes(X, Y, label = ID), size = 5) +
  geom_point(data = sites, aes(x = longitude, y = latitude, color = diversity), 
             shape = 19, size=2) +
  coord_sf(xlim = c(-114.05, -108.95), ylim = c(42.05,36.95), expand = FALSE) +
  labs(color = "Shannon Diversity", x="Longtitude",y="Latitude") +
  scale_color_gradient(low="Blue",high= "Orange")
ggsave("./output/Site_Map_merged.png", dpi=300,height = 10,width = 12)


diversity = diversity(otu_table(ps_ra), "shannon")

sites = data.frame(latitude = sample_data(ps_ra)$Latitude, longitude = sample_data(ps_ra)$Longitude, 
                   diversity = diversity)


states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))
states$ID <- toTitleCase(states$ID)

theme_set(theme_bw())

ggplot(data = world) +
  geom_sf() +
  geom_sf(data = states, fill = "White") + 
  geom_text(data = states, aes(X, Y, label = ID), size = 5) +
  geom_point(data = sites, aes(x = longitude, y = latitude, color = diversity), 
             shape = 19, size=2.5) +
  coord_sf(xlim = c(-114.05, -108.95), ylim = c(42.05,36.95), expand = FALSE) +
  labs(color = "Shannon Diversity", x="Longtitude",y="Latitude") +
  scale_color_gradient(low="Blue",high= "Orange")
ggsave("./output/Site_Map.png", dpi=300,height = 10,width = 12)



# metadata figures ####

meta_figs = as.data.frame(sample_data(ps_FF))
meta_figs = meta_figs %>% select(Set_ID,FireTreatment,BurnYear,Shannon) 
meta_figs$ANNUAL_MEAN_TEMP
ggplot(meta_figs, aes(x=BurnYear,y=ANNUAL_MEAN_TEMP, color=FireTreatment)) +
  geom_jitter(size=2) + 
  scale_x_continuous(labels = c(unique(meta_figs$BurnYear)), breaks = c(unique(meta_figs$BurnYear))) +
  labs(x= "Year of burn",y="Annual Mean Temp (C)",color="Fire Treatment",caption = "Sample distribution through time and temperature") + 
  scale_color_manual(values=pal.discrete[c(2,6)]) +
  theme(axis.text.x = element_text(angle=60,hjust=1))
ggsave("./output/sample_distribution.png", dpi=300)


ggplot(meta_figs, mapping = aes(x=BurnYear, y=Shannon, color=FireTreatment)) +
  geom_point() + geom_smooth(method = "lm") + scale_color_manual(values=pal.discrete[c(2,6)])


# Alpha diversity between treatment pairs ####

sample_data(ps_FF)$Richness <- specnumber(otu_table(ps_FF))
alpha.df <- as(sample_data(ps_FF),"data.frame")


ggplot(alpha.df, aes(x=factor(BurnYear),y=Richness,fill=FireTreatment)) +
  geom_boxplot() + labs(x="Burn year", y="Fungal richness",fill= "Site status") +
  scale_fill_manual(values=pal.discrete[c(2,6)])
ggsave("./output/boxplot_richness_vs_burnyear.png")


# Boxplots between groups #
names(alpha.df)
comparisons <- list(c("in","out"))
alpha.df$InOrOut
ggboxplot(alpha.df, x = "BurnYear", y = "Richness",
                color = "InOrOut", palette =pal.discrete,
                add = "jitter") +
  stat_compare_means(comparisons = comparisons)


+
  theme(axis.title = element_text(face="bold",size=14),
        axis.text = element_text(face="bold")) +
  labs(x="")


# Site map ####

# ggmap::register_google(key = "??????????") # Key kept private
# 
# df2 = read.csv("Desktop/Fire_metadata.csv")
# 
# ggmap(get_googlemap(center = c(lon = -111.4, lat = 39.5),
#                     zoom = 7, scale = 2,
#                     maptype ='satellite')) + 
#   geom_point(aes(x = Longitude, y = Latitude, colour = BurnYear), data = df2, size = 4) +
#   theme(legend.position="right") +
#   scale_colour_gradient(low = 'orange', high = 'red') +
#   borders("state", colour = "dark blue", region = "utah", size = 2) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())

# Look for Rhizopogon and Wilcoxina differences between burned-unburned
# sum((ps_FF@tax_table[,6] == "g__Rhizopogon"),na.rm = TRUE) 
# sum((ps_FF@tax_table[,6] == "g__Wilcoxina"),na.rm = TRUE)
# tax_table(ps_FF)
# 
# plot_bar(ps_FF, fill="Order")
# plot_bar(psm_ra, fill="Genus")


# Find species, if any, that are more prevalent in burned sites than unburned, arrange by burn year

sample_data(ps_FF)[,"FireTreatment"]

ps_family <- tax_glom(ps_FF,"Family")

set.seed(123)
da_analysis <- differentialTest(formula = ~ FireTreatment, #abundance
                                phi.formula = ~ FireTreatment, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_family,
                                fdr_cutoff = 0.05)


da_analysis$significant_taxa

stax_1 <- paste(tax_table(ps_family)[da_analysis$significant_taxa,2:5][1],sep="_")
stax_1 <- paste(stax_1[1],stax_1[2],stax_1[3],stax_1[4],sep="_")
stax_2 <- paste(tax_table(ps_family)[da_analysis$significant_taxa,2:5][2],sep="_")
stax_2 <- paste(stax_2[1],stax_2[2],stax_2[3],stax_2[4],sep="_")
stax_3 <- paste(tax_table(ps_family)[da_analysis$significant_taxa,2:5][3],sep="_")
stax_3 <- paste(stax_3[1],stax_3[2],stax_3[3],stax_3[4],sep="_")
stax_4 <- paste(tax_table(ps_family)[da_analysis$significant_taxa,2:5][3],sep="_")
stax_4 <- paste(stax_3[1],stax_4[2],stax_4[3],stax_4[4],sep="_")

names(da_analysis$significant_models) <- c(stax_1,stax_2,stax_3,stax_4)

sink("./output/Differential_abundance_model_stats_tables.txt")
print("Family-level taxonomic comparisons...")
da_analysis$significant_models
sink(NULL)


daplot_1 <- plot(da_analysis) + labs(y="Differentially abundant taxa\n(relative abundance)") +
  theme(axis.text.y = element_text(face="bold.italic"),
        axis.title.y = element_text(face="bold",size=16),
        strip.text = element_text(face="bold",size=12))

ggsave(daplot_1, filename = "./output/Diff_Abund_Family_By_BurnTreatment.png",device = "png",width = 9,height = 5,dpi=300)

da_analysis$significant_taxa


set.seed(123)
corncob_da1 <- bbdml(formula = ACGAGTTACAAGTCGGTCGACCGTGCTGGCGGAAACGCACGTGCACGTCGGTCGCAAACCTCATCCACACACCTGTGAACGTATGGCCTTGGGTCTCTCGACCCGGGGCAAACCTTTTTTACCCACTCTGTTTGTAAAGGAATGTCATACGTGCGTAACGCATAAATGAA ~ FireTreatment,
                     phi.formula = ~ FireTreatment,
                     data = ps_family)

# pull out model results into df
corncob_da1_wald <- waldt(corncob_da1)
corncob_da1_wald <- corncob_da1_wald[grep("mu.",row.names(corncob_da1_wald)),]
corncob_da1_wald <- tidy(corncob_da1_wald)
corncob_da1_wald$OTU <- stax_1

p1 <- plot(corncob_da1,color="FireTreatment") + scale_color_manual(values = pal.discrete[c(2,6)] ) + ggtitle(stax_1)

###
set.seed(123)
corncob_da2 <- bbdml(formula = CTGAACTGTCAACACGAGTTGTTGCTGGTCCTCAAATGGGGGCATGTGCACGCTCTGTTTACATACCCACTCACACCCGTGCACCCTCTGTAGTTCTGTGGTGTGGGGGACTCTGTCCTCCCGCTGTGGTTCTATGTCTTTTACACACACACAGTCTCATAGAATGTATGTCGCGTTTAACGCAATACAATA ~ FireTreatment,
                     phi.formula = ~ FireTreatment,
                     data = ps_family)

# pull out model results into df
corncob_da2_wald <- waldt(corncob_da2)
corncob_da2_wald <- corncob_da2_wald[grep("mu.",row.names(corncob_da2_wald)),]
corncob_da2_wald <- tidy(corncob_da2_wald)
corncob_da2_wald$OTU <- stax_2

p2 <- plot(corncob_da2,color="FireTreatment") + scale_color_manual(values = pal.discrete[c(2,6)] ) + ggtitle(stax_2)

###
set.seed(123)
corncob_da3 <- bbdml(formula = AAGAGATAGGGTGCTCAGCGCCCGACCTCCAACCCTTTGTTGTTAAAACTACCTTGTTGCTTTGGCGGGACCGCTCGGTCTCGAGCCGCTGGGGATTCGTCCCAGGCGAGTGCCCGCCAGAGTTAAACCAAACTCTTGTTAATTAAACCGGTCGTCTGAGTTAAAATTTTGAATAAATCA ~ FireTreatment,
                     phi.formula = ~ FireTreatment,
                     data = ps_family)

# pull out model results into df
corncob_da3_wald <- waldt(corncob_da3)
corncob_da3_wald <- corncob_da3_wald[grep("mu.",row.names(corncob_da3_wald)),]
corncob_da3_wald <- tidy(corncob_da3_wald)
corncob_da3_wald$OTU <- stax_3

p3 <- plot(corncob_da3,color="FireTreatment") + scale_color_manual(values = pal.discrete[c(2,6)] ) + ggtitle(stax_3)

###
set.seed(123)
corncob_da4 <- bbdml(formula = CCGAAGTTACCTTCAAAACCCACTGTGAACCTTACCTCTTGCCGCGTTGTCTCGGCGGGAGGCGGTGGGCGTCGCGTGCCCTAGCGGGCCGTGCCGCTCCCGTCCCCGCCGGCGGCGCCAAACTCTAAATTTACAGCGGACTGTATGTTCTGATTTACAAAAAAAACAAGTTA ~ FireTreatment,
                     phi.formula = ~ FireTreatment,
                     data = ps_family)

# pull out model results into df
corncob_da4_wald <- waldt(corncob_da4)
corncob_da4_wald <- corncob_da4_wald[grep("mu.",row.names(corncob_da4_wald)),]
corncob_da4_wald <- tidy(corncob_da4_wald)
corncob_da4_wald$OTU <- stax_4

p4 <- plot(corncob_da4,color="FireTreatment") + scale_color_manual(values = pal.discrete[c(2,6)] ) + ggtitle(stax_4)



# combine plots
p1/p2/p3/p4
ggsave("./output/Diff_Abund_Family_By_BurnTreatment_Indiv_Taxa.png",height = 6,width = 12,dpi=300)

# join all 4 together
full_corncob <- rbind(corncob_da1_wald,corncob_da2_wald,corncob_da3_wald,corncob_da4_wald)

# plot
tidy_corncob <- full_corncob %>% select(FireTreatment = .rownames, Estimate, StdErr = Std..Error, t.value, P.val = Pr...t..,OTU) %>%
  filter(FireTreatment != "mu.(Intercept)") %>%
  # arrange(AgeGroup) %>%
  mutate(ymin = Estimate - StdErr, ymax=Estimate + StdErr)

tidy_corncob$FireTreatment <-  str_remove(tidy_corncob$FireTreatment,pattern = "mu.CoralAgeBinned") 


ggplot(tidy_corncob, aes(x=FireTreatment,y=Estimate)) +
  geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = .2) + theme_bw() +
  geom_hline(yintercept =0,linetype=2,alpha=.5) +
  coord_flip() +
  facet_grid(~OTU) + 
  labs(x="Fire Treatment", y= "Wald test estimate") + theme(strip.text = element_text(size=14,face="bold"),
                                                                     axis.title = element_text(size=14,face="bold"),
                                                                     axis.text = element_text(size=12,face = "bold"))






























taxa_names(ps_family)


sample_data(ps_FF)
sample_data(ps_FF)$DaysSinceBurn <- as.numeric(as.POSIXct("2019-01-01") - as.POSIXct(sample_data(ps_FF)$BurnDate,format='%Y-%m-%d'))

ggplot(sample_data(ps_FF),aes(x=DaysSinceBurn,y=Richness,color=FireTreatment,group=Set_ID)) +
  geom_point() + facet_wrap(~Set_ID)

tempgroups = kmeans(x = df2$ANNUAL_MEAN_TEMP,centers = 2)
tempgroups$cluster

mod3 = lmer(data = as(sample_data(ps_FF),"data.frame"), Richness ~ DaysSinceBurn + ANNUAL_MEAN_TEMP + (1|Set_ID)/(1|FireTreatment))
summary(mod3)

df2 = add_predictions(data = as(sample_data(ps_FF),"data.frame"),model = mod3,type = "response")
df2$tempgroup <- tempgroups$cluster
sink("./output/lmer_model_Richness.txt")
print("Richness as dependent var., DaysSinceBurn and ANNUAL_MEAN_TEMP as fixed vars., FireTreatment nested within SetID as random vars")
summary(mod3)
sink(NULL)

<<<<<<< HEAD
# # Site table ####
# df = as(sample_data(psm_ra),"data.frame")
# write.csv(df,"./output/site_info.csv",row.names = TRUE,quote = FALSE)
=======
ggplot(df2,aes(x=DaysSinceBurn,y=Shannon,color=factor(tempgroup))) +
  geom_point() + geom_smooth(method="lm",se=FALSE) + facet_wrap(~FireTreatment) +
  scale_x_reverse()

# difference between shannon diversity between burned/unburned
d = df2 %>% group_by(factor(NewPastedVar)) %>%
  dplyr::summarize(Shan = mean(Shannon))
d
d = data.frame(SetID = str_remove(d$`factor(NewPastedVar)`[1:10],pattern = "Burn_") ,DiffShannon = d[1:10,2] - d[c(11:13,15:21),2])
names(d)[2] <- "DiffShannon"
d
d$DiffShannon = as.character(d$DiffShannon * -1)
d

df2$ShannonDifference <- as.numeric(as.character(plyr::mapvalues(df2$Set_ID,from = d$SetID,to=d$DiffShannon)))


mod4 <- glm(data=df2,ShannonDifference ~ tempgroup * DaysSinceBurn)
summary(mod4)

ggplot(df2,aes(x=DaysSinceBurn,y=ShannonDifference,color=factor(tempgroup))) + geom_point() +
  geom_point(aes(y=predict(object = mod4,newdata = df2)),color="Red",size=3) +
  geom_smooth(method="lm",formula = y~x+color)
# this plot sucks. I hate it

ggplot(subset(df2,Set_ID != "Set_2"),aes(x=factor(BurnYear),y=Shannon,fill=FireTreatment)) +
  geom_boxplot(position="dodge")

ggplot(subset(df2,Set_ID != "Set_2"),aes(x=factor(round(ANNUAL_MEAN_TEMP)),y=Shannon,fill=FireTreatment)) +
  geom_boxplot(position="dodge")



df2$ANNUAL_MEAN_TEMP
# Site table ####
df = as(sample_data(psm_ra),"data.frame")
write.csv(df,"./output/site_info.csv",row.names = TRUE,quote = FALSE)
>>>>>>> 951dbbb763cdfb68db506575459e058ae34d8b10


# Import plant data -- not finished!!! ####
plants <- read.csv("./plant_metadata_by_site.csv")

plants$Unique.ID
plants$Family=as.character(plants$Family)
plants$Family[plants$Family == "Moss species*"] <- "Bryophyte"
plants$Family <- factor(plants$Family)

# make it match :(
<<<<<<< HEAD
plants <- plants[plants$Unique.ID %in% meta$ID,]

plantrichness <- table(plants$Unique.ID)
plantrichness <- as.data.frame(plantrichness)
plantrichness$Var1 <- as.character(plantrichness$Var1)
names(plantrichness) <- c("ID","Plant_Richness")

meta <- left_join(sam,plantrichness)


mod <- glm(data=meta,Richness ~ FireTreatment)
summary(mod)
# No noticable difference in alpha diversity based on burn year or fire treatment
# But... B-diversity is a different story!









=======
'%ni%' <- Negate('%in%')
# find any values that don't match
notinmeta <- which(meta$ID %ni% plants$Unique.ID)
nrow(meta[notinmeta,])


# find unique genus values in each site
x=1
y=list()
for(i in levels(plants$Unique.ID)){
  d <- plants %>% filter(Unique.ID == i)
  y[[x]] <- as.character(unique(d$Family))
  names(y)[x] <- as.character(i)
  x=x+1
}
y

# build into presence-absence table
# blank data frame
pa <- data.frame(row.names = names(y))
for(i in levels(d$Family)){
  pa[,i] <- NA
}

# fill, row by row
x=1
for(i in names(y)){
  pa[x,] <- colnames(pa) %in% y[[i]]
  x=x+1
}
pa

# convert to presence-absence for vegan
pa[pa==FALSE] <- 0
pa[pa==TRUE] <- 1

# NMDS
NMDS = vegan::metaMDS(pa,distance = "jaccard")

mds1=NMDS$points[,1]
mds2=NMDS$points[,2]

treatment=unlist(purrr::map(str_split(unique(as.character(plants$Unique.ID)),"_"),2))

plantord=data.frame(mds1=mds1,mds2=mds2,treatment=treatment,site=row.names(plantord))
plantord$treatment = as.character(plantord$treatment)
plantord$treatment[plantord$treatment == "in"] <- "Burned"
plantord$treatment[plantord$treatment == "out"] <- "Unburned"


plantord %>% filter(site != "46_out") %>%
  ggplot(plantord,mapping = aes(x=mds1,y=mds2,color=treatment,label=site)) +
  geom_point(size=3) + stat_ellipse() + labs(x="MDS1",y="MDS2",color="Fire Treatment",
                                             caption = "Plant community in sampling plot")
ggsave("./output/NMDS_plant_community_by_burn_treatment.png",dpi=300,device = "png")


# Build gen. linear models incorporating habitat info besides temp and precip (tree cover, etc)
# use lmer

# burnstatus nested within site
# sample ID nested within burnstatus for each site
>>>>>>> 951dbbb763cdfb68db506575459e058ae34d8b10

# site and burnstatus are random; temp, precip, etc are fixed
