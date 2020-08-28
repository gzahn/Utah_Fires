# This script is for importing and cleaning GeoDataCrawler variables for the selected Utah burn sites 


#### Libraries and aliases used: ####
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
"%ni%" <- Negate("%in%")

#### Loading ALL that GeoDataCrawler data ####

dat = read_excel("INIT20180405_100m_1.xlsx")
dat2 = read_excel("INIT20180405_100m_2.xlsx")
meta = read_excel("Geodata Crawler - Output Data Descriptions.xlsx")

#### Compare paired burned/un-burned sites across all variables for correlation (should be high) ####

  ## Subset burned and unburned sites
unique(dat$InOrOut)
burned = dat$InOrOut == "in"
unburned = dat$InOrOut == "out"

burned = dat[burned,]
unburned = dat[unburned,]

  ## Differing lengths ... find "pair" that is missing a partner
  
        # Find where mismatch begins
mismatches = unburned$ObjectID_1[order(unburned$ObjectID_1)] == burned$ObjectID_1[order(burned$ObjectID_1)]
first_mismatch = which(match(mismatches, FALSE) == 1)[1]
bad_object = dat[(first_mismatch-1),"ObjectID_1"]

        # re-subset to remove site with mismatch
burned = burned[burned$ObjectID_1 %ni% bad_object,]
unburned = unburned[unburned$ObjectID_1 %ni% bad_object,]
dat = dat[dat$ObjectID_1 %ni% bad_object,] # do this for main datasets too!
dat2 = dat2[dat2$ObjectID_1 %ni% bad_object,]         

        # reorder
burned = burned[order(burned$ObjectID_1),]
unburned = unburned[order(unburned$ObjectID_1),]
identical(burned$ObjectID_1,unburned$ObjectID_1) # check to see if subsets and orders are identical

  ## Plot paired Latitudes and elevations for sanity check .. should be very similar
plot(burned$Lat,unburned$Lat) # Looks good.
plot(burned$L100_ELEVATION_AVG,unburned$L100_ELEVATION_AVG) # Looks good.
plot(burned$L100_SOIL_DOM, unburned$L100_SOIL_DOM) # Looks good.


#### Clean up metadata and data and get them ordered logically ####
        
        # remove all those "L100_" tags from the column names
names(dat) = gsub("L100_","",names(dat), perl = TRUE) 

        # Convert M,D,Y to Dates
            # First, remove any rows that have a "0" in Y,M,orD
good_dates = !(dat$Year_ == 0 | dat$Month_ == 0 | dat$Day_ == 0)
dat = dat[good_dates,]
good_dates = !(dat2$Year_ == 0 | dat2$Month_ == 0 | dat2$Day_ == 0)
dat2 = dat2[good_dates,]

            # Now, extract dates from columns
year = dat$Year_
month = dat$Month_
day = dat$Day_

year2 = dat2$Year_
month2 = dat2$Month_
day2 = dat2$Day_

            # Add BurnDate column
burn_date = as.Date(paste(year,month,day, sep = "-"), format = "%Y-%m-%d")
dat$BurnDate = burn_date
burn_date2 = as.Date(paste(year2,month2,day2, sep = "-"), format = "%Y-%m-%d")
dat2$BurnDate = burn_date2

#### Merge dat and dat2 ####
dat = merge(dat,dat2)


# remove all those "L100_" tags from the column names (again!)
names(dat) = gsub("L100_","",names(dat), perl = TRUE) 


#### Trim down columns only to potentially important ones ####

# Let's get rid of all the OIL and GAS and WELLS, OTHER, FIRE, IMPERVIOUS stuff
good_cols = grep("^OIL_|^GAS_|^FIRE|^OTHER|^ROAD|^IMPERVIOUS|LIGHTS|^PESTICIDE", 
                 names(dat), invert = TRUE)

dat = (dat[,good_cols])

# Remove some other stuff
dat = dat %>% 
  select(-c(OBJECTID_12,Shape,OBJECTID,POINT_ID,POINT_X,POINT_Y,ORIG_X,ORIG_Y,DATUM,SKIP_ROW,MOVE_DIST,
            MOVE_ANGL,DUPLICATE,HOUSEHOLD_DEN_2000,HOUSEHOLD_2000,AQUIFER_DOM,AQUIFER_TYPE_DOM))


# Tidy up the column names and units ####
names(dat)

# Temperatures should be dived by 10 to get deg C
temps = names(select(dat,c(starts_with("TM"),contains("TMIN"),contains("TMAX"),contains("TEMP"))))
dat[,temps] = select(dat,starts_with("TM")) * 0.1

# Convert canopy cover to percentages
dat$CANOPY_2001 = dat$CANOPY_2001 * 0.01

# Convert soil code into soil names
soil_codes = read.csv("soil_codes.csv", stringsAsFactors = FALSE)
dictionary = (soil_codes[which(soil_codes$Code %in% dat$SOIL_DOM),])

codes = as.factor(dat$SOIL_DOM)
dictionary = dictionary[dictionary$Code %in% levels(codes),]
s_names = plyr::mapvalues(dat$SOIL_DOM, from=levels(as.factor(dictionary$Code)), 
                          to=levels(as.factor(dictionary$Soil.Names)))  
dat$Soil_Description = s_names
names(dat)


# Remove all the soil type columns now that plotting confirms good correlations ####
dat = dat %>% select(-starts_with("LC2011"))


# Remove duplicate ObjectIDs
sorted = dat[order(dat$ObjectID_1, dat$InOrOut,decreasing = TRUE),]
sorted$to_keep = c("yes","no")
cbind(sorted$to_keep,sorted$InOrOut)

sorted = sorted[!duplicated(sorted[,c("ObjectID_1","InOrOut","to_keep")]),]
sorted = sorted[!duplicated(sorted[,c("ObjectID_1","InOrOut")]),]

# order by date
dat = sorted[order(sorted$BurnDate),]

table(dat$Year_)

# Write intermediate file ####
write.csv(dat, "Cleaned_Site_Data.csv", row.names = FALSE, quote = FALSE)





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Beyond this point is extra crap ####


# Backup dat
dat2 = dat

# Subset original dat to just Burned areas
dat = subset(dat,InOrOut == "in")

# Plotting using Soil type % cover as NMDS ####

# Pull out dfs with Soil Info (percent cover 100m) only!

soil_info_L1 = dat %>%
  select(starts_with("LC2011_L1_")) 

soil_info_L2 = dat %>%
  select(starts_with("LC2011_L2_")) 


row.names(soil_info_L1) <- paste0(dat$Lat,"_",dat$Lon)
row.names(soil_info_L2) <- paste0(dat$Lat,"_",dat$Lon)

# Soil-based distance measures
L1_dist = as.matrix(dist(soil_info_L1))
L2_dist = as.matrix(dist(soil_info_L2))


# Heatmap
year = as.factor(dat$Year_)
year_colors = plyr::mapvalues(year, 
                              from = c(levels(year)),
                              to = gray.colors(length(levels(year))))

year_colors = as.character(year_colors)

h1 = heatmap(L1_dist, col = gray.colors(20),
        RowSideColors = year_colors,
        ColSideColors = year_colors)  ## Shows decent spread of years in most similar clusters

h2 = heatmap(L2_dist, col = gray.colors(20),
        RowSideColors = year_colors,
        ColSideColors = year_colors,
        cexRow = .25)  ## Shows decent spread of years in most similar clusters

?png
png("heatmap_in_only.png",res = 500, width = 5000, height = 5000)
heatmap(L2_dist, #col = gray.colors(20),
        RowSideColors = year_colors,
        ColSideColors = year_colors,
        cexRow = .25, cexCol = .25, distfun = vegdist)  ## Shows decent spread of years in most similar clusters

dev.off()

row_order = h2$rowInd
col_order = h2$colInd
row_similar = row_order[50:90]
col_similar = col_order[50:90]

plot(L2_dist)
plot(L2_dist[col_similar])


row.names(dat) <- paste0(dat$Lat,"_",dat$Lon)
soil_similar = row.names(L2_dist[row_similar,])
row.names(L2_dist)

table((dat[soil_similar,])[,"Soil_Description"])
table(dat$Soil_Description)
table(dat$SOIL_DOM)



names(dat)
# NMDS
soil_mds_1 = metaMDS(soil_info_L1)
ordiplot(soil_mds_1, type = "points")
stressplot(soil_mds_1)

soil_mds_2 = metaMDS(soil_info_L2)
ordiplot(soil_mds_2, type = "points")
stressplot(soil_mds_2)

# Make into dfs
mds1_L1 = soil_mds_1$points[,1]
mds2_L1 = soil_mds_1$points[,2]

mds1_L2 = soil_mds_2$points[,1]
mds2_L2 = soil_mds_2$points[,2]

# Colored plots of burn status influence on major soil types

ggplot(mapping = aes(x=mds1_L1,y=mds2_L1,color=dat$InOrOut)) +
  geom_point() + 
  stat_ellipse() +
  ggtitle("Soil type is independent of burn status")

ggplot(mapping = aes(x=mds1_L2,y=mds2_L2,color=dat$InOrOut)) +
  geom_point() + 
  stat_ellipse() +
  ggtitle("Soil type is independent of burn status")

# plots of year interaction with soil types

ggplot(mapping = aes(x=mds1_L2,y=mds2_L2,color=as.factor(dat$Year_))) +
  geom_point() + 
  stat_ellipse() +
  ggtitle("Pretty good overlap of burn years and soil types")

# Extract distance clusters
clustered_order = h2$rowInd
soil_similar = clustered_order[1:50]
plot(L2_dist)
plot(L2_dist[soil_similar])

soil_similar = row.names(L2_dist[soil_similar,soil_similar])

row.names(dat) <- paste0(dat$Lat,"_",dat$Lon)

# compare year spread in whole data set to hclust subset based on soil type similarity
plot(unique(dat$Year_)[order(unique(dat$Year_))])
plot(dat[soil_similar,"Year_"][order(dat[soil_similar,"Year_"])], col="Red",
     main = "Subset year distribution for hclust-similar soil types")

# Check in/out distribution
names(dat)
table(dat[soil_similar,"InOrOut"])
dat[soil_similar,which(dat$InOrOut == "in")]








# Plot yearly stuff to play around

plot(dat$BurnDate, dat$PRECIP_DRIEST_MONTH)
abline(lm(dat$PRECIP_DRIEST_MONTH ~dat$BurnDate))
plot(as.factor(dat$Year_), dat$PRECIP_SEASONALITY)



#### Remove Soil Type outliers


names(dat)




# L2_cols = grep("LC2011_L2",names(dat))
# dat[,L2_cols]

# 





# #### Select variables of interest
# c("SOIL_DOM","SLOPE_AVG", "ELEVATION_AVG","ASPECT_DOM","LC2011_L2_DOM")







dat$ObjectID_1 %in% burned$ObjectID_1
dat$ObjectID_1 %in% unburned$ObjectID_1

burned$ObjectID_1
unburned$ObjectID_1



  ## Compare a few obvious variables for sanity check
plot(burned$Lat, unburned$Lat)



