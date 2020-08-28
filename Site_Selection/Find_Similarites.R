# Libraries Used ####
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
"%ni%" <- Negate("%in%")

# Import cleaned Data ####
dat = read.csv("Cleaned_Site_Data.csv", stringsAsFactors = FALSE, row.names = "UniqID")
dat$ID = row.names(dat)

# Convert to Date class
dat$BurnDate = as.Date(dat$BurnDate, format = "%Y-%m-%d")

# Find all those numeric columns
classes = c()
i=1
for(col in names(dat)){
  classes[i] = class(dat[,col])
  i=i+1
}

# subset to just numeric environmental data
dat_matrix = dat[,classes == "numeric"]
dat_matrix = select(dat_matrix,-c(Lat,Lon))

# Compare burned sites only
dat_burned = dat_matrix[dat$InOrOut == "in",]

# NMDS on environmental data
NMDS = metaMDS(dat_matrix, distance = "euclidean")

MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
MDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Year = dat$Year_, ID = dat$ID, Soil = dat$SOIL_DOM)

# Plot NMDS
ggplot(MDS, aes(x=MDS1, y=MDS2, color=as.factor(Year))) +
  geom_point() +
  annotate("rect", xmin = 0, xmax = 30,ymin = -20,ymax = 20, fill = "Brown", alpha = 0.2) +
  ggtitle("Shaded rectangle shows selected sites")

  #xlim(c(0,30)) + ylim(c(-20,20))
ggsave("Environ_NMDS_Sample_Selection.png")


# Extract MDS points from selected range
mds1_range = MDS$MDS1 >= 0 & MDS$MDS1 <=30
mds2_range = MDS$MDS2 >= -20 & MDS$MDS2 <= 20

good_rows = which(mds1_range == TRUE & mds2_range == TRUE)
selection = as.character(MDS[good_rows,"ID"])

# Subset cleaned data set to selection
dat2 = dat[selection,]
IDs = unique(dat2$ObjectID_1)

dat = dat[dat$ObjectID_1 %in% IDs,]

#sanity check
ggplot(mapping = aes(x=dat$BurnDate,y=dat$PRECIP_DRIEST_MONTH)) +
  geom_point() + stat_smooth(method = "lm")


# Write sample subset
write.csv(dat,"Selected_Sites.csv", row.names = FALSE, quote = FALSE)

