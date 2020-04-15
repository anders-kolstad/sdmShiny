
# Get data

#See if there are matches without downloading...
dismo::gbif("Lynx", "pardinus", download = F)

# download as df
sp <- dismo::gbif("Lynx", "pardinus", download = T, geo=T, sp=F)

class(sp)  

# subset
w <- which(is.na(sp$lon))
sp <- sp[-w,]
which(is.na(sp$lat))

sp$species <- 1
sp <- sp[,c("lon", "lat","species")]
head(sp)

# define as spaTIAL POINTS
sp::coordinates(sp) <- ~lon + lat
head(sp)


# get env data
bio <- raster::getData('worldclim', var='bio',res=10)
bio


# check collinerity
library(usdm)
# two options
usdm::vifstep(bio)
v2 <- usdm::vifcor(bio, th=0.7)

biom <- usdm::exclude(bio, v2)

plot(biom[[1]]) # error
par(plt=c(0,1,0,1))
plot(biom[[1]])
points(sp)

#define projection and plot in mapview
proj4string(sp) <- proj4string(raster())
mapview(sp)
ras <- raster








# SPEED ####
setwd("/home/anders/Documents/SDMs/demoSDM")
prisca<-readRDS('data/prisca_obs.rds')
Priscamaxent <- readRDS("data/prisca_maxent.rds")
botlan<-readRDS('data/botlan_obs.rds')
Botlanmaxent <- readRDS("data/botlan_maxent.rds")

selectvarsshiny<-raster::stack('data/selectvars')

plot(selectvarsshiny)


p1<-dismo::predict(Priscamaxent,selectvarsshiny)
plot(p1) # current distribution of Primula scandinavica

#increase herbive biomass by 20%:
new_env <- selectvarsshiny
new_env[[1]] <- new_env[[1]]*((20+100)/100)
plot(new_env[[1]])
# and make new predictions:
p2<-dismo::predict(Priscamaxent,new_env)
plot(p2)


# find the areas that are different:
diff <- p2-p1
plot(diff)

#rasterVis::levelplot(diff)

# Use threshold to define presense (as opposed to probability of occurence)
pa1 <- raster(p1)
pa1[] <- ifelse(p1[] >= 0.2, 1, 0)
# selecting threshold = 0.2 (should be objective)

# now the same for the new env
pa2 <- raster(p2)
pa2[] <- ifelse(p2[] >= 0.2, 1, 0)

# and the differece
pa <- pa2-pa1
plot(pa)
c1 <- colorRampPalette(c("black", "grey", "green"))

plot(pa, col=c1(3))

# with three colours:
pa3 <- pa
pa3[pa3==0] <- ifelse(pa1[pa3==0] == 0, 2,0)
plot(pa3)

c2 <- colorRampPalette(c("black", "grey", "green", "white"))
plot(pa3, col=c2(4))

#outline norway, grey current distribution? 
ol <- p1 > -Inf
pp <- raster::rasterToPolygons(ol, dissolve=TRUE)
plot(pp, lwd=1, border='black', add=TRUE)


# compare three plots
par(mar=c(2,2,3,1), mfrow = c(2,3))
getwd()
tiff("c-e.tiff", width = 800, height = 800, units = "px")
plot(p1, main = "Est. prob. of occ.\nunder current cond.")
plot(p2, main = "Est. prob. of occ.\nunder new cond.")
plot(pa1, main = "Est. presence\nunder current cond.")
plot(pa2, main = "Est. presence\nunder new cond.")

plot(pa3, col=c2(4), main = "Est. collonisations (green)\n and extinctions (black)")
plot(pp, lwd=1, border='black', add=TRUE)
dev.off()
