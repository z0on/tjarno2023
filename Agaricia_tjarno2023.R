
library(tess3r)
library(vegan)
library(ggplot2)
library(sf)      
library(rnaturalearth)
library(rnaturalearthdata)
library(gradientForest)
library(dplyr)
library(raster)
library(fields)
library(maps)
coasts= ne_coastline(scale="large")
theme_set(theme_bw())

# ------------------------------
# function to run gradient forest on PCs of a vegan ordination object

makeGF=function(ordination,environment,ntrees=500) {
  Y=data.frame(scores(ordination,scaling=1,choices=c(1:length(ordination$CA$eig)))$sites)
  X=environment
  nSites = dim(Y)[1]
  nSpecs = dim(Y)[2]
  lev = floor(log2(nSites * 0.368/2))
  gf = gradientForest(cbind(X, Y),
                      predictor.vars = colnames(X), response.vars = colnames(Y),
                      ntree = ntrees, transform = NULL, trace=T,
                      maxLevel = lev, corr.threshold = 0.25)
  return(gf)
}


#--------------------------------

# edit this line to point to the path you downloaded the files to
setwd("~/Dropbox/methods_in_ecogeno_2023/Agaricia_RDA_gradientForest_practice/") 

# environmental parameters
load("Agaricia_env_nr2.RData")
load("Agaricia_latlon.RData")

# admixture proportions between 3 clusters
admix=read.csv("Agaricia_admix.csv")

# genetic distances
IBS=as.matrix(read.table("Agaricia.ibsMat"))
samples=paste("aa",sub(".bam","",row.names(env)),sep="")

# adding sample names to everything
dimnames(IBS)=list(samples,samples)
row.names(admix)=row.names(latlon)=row.names(env)=samples

# raster of variables, to predict adaptation
ll=load("Rasters/rasters_XY.RData")
ll
# "rasters" "XY"

# ------ hierarchical clustering tree to look for clones, wrong species collections

hc=hclust(as.dist(IBS),"ave")
plot(hc,cex=0.5) # clustering of samples by IBS (great to detect clones or closely related individuals)
# there are a few "low-hanging" groups of closely related individuals 
# Note: to make sure that two samples are clones, their distance must be compared to the distance
# between genotyping replicates (the same sample genotyped twice)
abline(h=0.05,col="red") # this seems like a reasonable  "low-hang" threshold for calling related groups
pheatmap::pheatmap(1-IBS)

#--------- retaining a single representative of each highly-related group

cuts=cutree(hc,h=0.05)
goods=c();i=1
for (i in unique(cuts)) {
  goods[i]=names(cuts[cuts==i])[1]
}
length(goods)  # how many samples are left?
# subsetting all data for only the retained samples
IBS=IBS[goods,goods] 
latlon=latlon[goods,]
admix=admix[goods,]
env=env[goods,]

#----- Isolation by distance: is there population structure?

xydist=dist(latlon)
gendist=as.dist(IBS)
plot(gendist~xydist,col=rgb(0,0,0,0.05),pch=19,cex=0.5)
# multiple clines! => multiple genetic lineages (not just gradual change dependiong on distance)

# -------- PCoA

ord.all=capscale(as.dist(IBS)~1)
plot(ord.all,scaling=1)

# eigenvectors (% of variation explained)
plot(100*ord.all$CA$eig/sum(ord.all$CA$eig)) 

# ----- plotting admixture bar plot (tess3Q functions)

# my.colors <- c("tomato", "lightblue", "wheat","olivedrab", "cyan3","hotpink","gold","orange")
my.colors <- c("red", "blue", "gold")
my.palette <- CreatePalette(my.colors, 7)

qm=as.qmatrix(admix[,c(1:3)])
bp=barplot(qm,border=NA,space=0,col.palette = my.palette,xaxt="null")
axis(1, at = 1:nrow(qm), labels = admix$Site[bp$order], las = 3, cex.axis = .4) 

# --- PCoA plot colored by admixture cluster assignment

cols=rep("grey80",nrow(admix))
cols[admix$Blue>0.75]="skyblue"
cols[admix$Yellow>0.75]="gold"
cols[admix$Indigo>0.75]="darkblue"

axes2plot=c(1,2)
scores=scores(ord.all,scaling=1,display="sites",choices=axes2plot)
plot(scores,pch=16,col=cols,asp=1)
abline(h=0,lty=3,col="grey60")
abline(v=0,lty=3,col="grey60")

# ======= admixture clusters on map (using some  tess3r functions) ======

# ------ creating mapping polygon

# plot the whole raster with map
plot(XY, pch=15,cex = 0.5, asp = 1, col = "grey95")
#points(sf1$long,sf1$lat,pch=15,cex=0.2)
map(coasts,add=T)

# click around the area of interest 
keys.polygon=locator(n=100, type="l")
keys.polygon=cbind("lon"=keys.polygon[["x"]],"lat"=keys.polygon[["y"]])
keys.polygon[nrow(keys.polygon),]=keys.polygon[1,]
# saving the polygon so we don't have to click on map every time
save(keys.polygon,file="keys.polygon.RData")
lines(keys.polygon)

# ---- interpolating admixture colors onto polygon

i2p=admix[,c("Sample","Site")]
names(i2p)=c("ind","pop")
meta=data.frame(cbind(i2p,env))
coordinates = latlon[,c(2,1)]
load("keys.polygon.RData")
names(keys.polygon)=c("lon","lat")
margin=.1

# sanity check: plotting map with sites and polygon to interpolate
plot(lat~lon,coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)",asp=1,xlim=c(min(keys.polygon$lon)-margin,max(keys.polygon$lon)+margin))

polygon(keys.polygon,col="grey90")
points(lat~lon,coordinates, pch = 19, cex = 0.6,col="red")
text(meta[,c("lon","lat")],labels=meta$pop,cex=0.7,col="red",pos=4)
map(coasts,add = T, interior = F)

# magic spells to format our simple polygon into maps object
Keys=Polygon(keys.polygon)
KL=Polygons(list(Keys),"Keys")
BLL=SpatialPolygons(list(KL))
keys=SpatialPolygonsDataFrame(BLL,data.frame(N=c("one"),row.names=c("Keys")))

# plotting interpolated admixture clusters (tess3Q)
plot(qm,coordinates,cex=0.4,map.polygon=keys,asp=1,col.palette = CreatePalette(my.colors, 20),xlab = "Longitude (°E)", ylab = "Latitude (°N)")
# adding population labels (may get messy quickly)
text(latlon[,2:1],labels=meta$pop,cex=0.4,col="cyan3",pos=4)

# ======== Gradient Forest analysis

ord.all=capscale(as.dist(IBS)~1)
gf=makeGF(ord.all,env,ntrees=500) # will take a bit

# plotting bar chart of importances:
ii=data.frame(importance(gf))
names(ii)="importance"
# reordering by decreasing importance
ii$var=factor(row.names(ii),levels=rev(row.names(ii)))
ggplot(ii[1:min(10,nrow(ii)),],aes(var,importance))+geom_bar(stat="identity")+theme_bw()+coord_flip()+xlab("")

# violin plot of lineages by depth to convince ourselves that
# the three lineages are really distributed unevenly across depths
ld=data.frame(cbind(lineage=cols,depth=env$Depth))
# removing individuals with ambiguous lineage assignment
ld=ld[cols!="grey80",]
ld$depth=as.numeric(ld$depth)
ggplot(ld,aes(lineage,depth))+geom_violin(scale="width")

# which of the MDSes are actually explained, and with which R2?
gf$result

# plotting turnover curves: overall importance (y axis is common scale to compare importances of predictors)
# how many top turnover curves to plot
nvars=6
#pdf(file=paste(pop,"_importance.pdf",sep=""),height=3.5,width=3.5)
plot(gf,plot.type="C",show.species=F,common.scale=T,
     imp.vars=names(importance(gf))[1:nvars],
     #     cex.axis=0.6,cex.lab=0.7,
     line.ylab = 0.9, par.args = list(
       mgp = c(1.5,0.5, 0), 
       mar = c(2.5, 1, 0.1, 0.5), 
       omi = c(0,0.3, 0, 0)
     )
)
#dev.off()

# projecting turnover curves to raster grid
#nvars=length(importance(gf)) # to keep all variables for prediction
nvars=6 # to choose top 6 variables
important=names(importance(gf))[1:nvars]
raster.vars=colnames(rasters)[colnames(rasters) %in% important]
# actual prediction line:
trans.grid = cbind(XY, predict(gf, rasters[,raster.vars]))

# rescaling projected turnover curves to 0:importance in the model 
# (because otherwise out-of-model-range values in the predicted raster will inflate the importance)
for (v in raster.vars){
  trans.grid[,v]=scales::rescale(trans.grid[,v])*importance(gf)[v]
}

# making unscaled PCA (so that variables play according to their importance)
# and setting up color scheme to visualize first three PCs
# (cannibalized from Ellis' Gradient Forest vignette)
pc <- prcomp(trans.grid[, raster.vars])
plot(pc$sdev)
pcs2show=c(1,2,3)
# color flippage flags - change between -1 and 1 to possibly improve color representation
flip.pc1=(1)
flip.pc2=(1)
flip.pc3=(1)

pc1 <- flip.pc1*pc$x[, pcs2show[1]]
pc2 <- flip.pc2*pc$x[, pcs2show[2]]
pc3 <- flip.pc3*pc$x[, pcs2show[3]]
b <- pc1 - pc2
g <- -pc1
r <- pc3 + pc2 - pc1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
nvs <- dim(pc$rotation)[pcs2show[1]]
lv <- length(raster.vars)
vind <- rownames(pc$rotation) %in% raster.vars
scal <- 40
xrng <- range(pc$x[, 1], pc$rotation[, pcs2show[1]]/scal) * 1.1
yrng <- range(pc$x[, 2], pc$rotation[, pcs2show[2]]/scal) * 1.1

# plotting PCA of the predicted adaptive communities
# (should look like a cloud-ish shape, not a long stick!)
plot((pc$x[, pcs2show[1:2]]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1)
points(pc$rotation[!vind, pcs2show[1:2]]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), pc$rotation[raster.vars, pcs2show[1]]/scal, pc$rotation[raster.vars, pcs2show[2]]/scal, length = 0.0625)
jit <- 0.0015
text(pc$rotation[raster.vars, 1]/scal + jit * sign(pc$rotation[raster.vars, pcs2show[1]]), pc$rotation[raster.vars, pcs2show[2]]/scal + jit * sign(pc$rotation[raster.vars, pcs2show[2]]), labels = raster.vars,cex=0.7)
man.colors=rgb(r, g, b, max = 255)

# map of adaptation gradients
#pdf(file=paste(pop,"_adaptation_map.pdf",sep=""),height=3.5,width=5)
plot(XY, pch=15,cex = 0.5, asp = 1, col = man.colors)
#points(sf1$long,sf1$lat,pch=15,cex=0.2)
map(coasts,add=T)
#dev.off()

# >>> ON YOUR OWN <<<< 
# - analyze each of the three lineages independently
# - compare adaptation maps: plot where there is difference in adaptation between lineages

# ======= some scrap code that will be useful:

# ------------ choosing subpop to work on
pop="Yellow"
choices=which(admix[,pop]>=0.75)
env.p=env[choices,]
IBS.p=IBS[choices,choices]
latlon.p=latlon[choices,]
admix.p=admix[choices,]

# ---- modification to ordination command when analyzing individual lineages,
#      to DISREGARD admixture with other lineages as a force shaping genetic distances
#      (conditional PCoA)

# removing the focal lineage from admixture table (covariates will be two others)
admix.cov=admix[choices,1:3]
admix.cov=admix.cov[,!(names(admix.cov) %in% pop)]

# computing conditional PCoA
ord.all=capscale(as.dist(IBS.p)~1+Condition(admix.cov[,1])+Condition(admix.cov[,2]))
plot(ord.all)

# ---- to compare how adaptation maps differ between lineages

yellow.raster=trans.grid.yellow
blue.raster=trans.grid.blue

# determining the scale for adaptation offset:
# value of 1 will be as much adaptation difference as the 90th percentile difference between points in the yellow raster
# (will take a bit)
dd=dist(trans.grid.yellow[,raster.vars])
qq=quantile(abs(dd),c(0.5,0.9),na.rm=T)
med=qq[1]
max=qq[2]

# calculating adaptation offset (sometimes called genomic offset)
#    function to compute distance between corresponding rows of two data frames
dist2dataframes=function(X,Y,method="euclidean"){
  if((nrow(X) == nrow(Y)) & (ncol(X) == ncol(Y))) {
    di=c()
    for (i in 1:nrow(X)) {
      #    message(i)
      dd=data.frame(rbind(X[i,],Y[i,]))
      di=c(di,dist(dd,method=method))
    }
    return(di)
  } else { 
    stop("dist2dataframes: dataframes are of different sizes!")
  }
}

# computing offset, rescaling to 90%th quantile of yellow raster distance
yb.offset=dist2dataframes(yellow.raster[,raster.vars],blue.raster[,raster.vars])
yb.offset.2max=yb.offset/qq[2]

ggplot(XY,aes(x,y,color=yb.offset.2max))+
  geom_point(cex=0.5)+
  theme_void()+
  coord_equal()




