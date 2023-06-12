
library(ggplot2)

# function to read and average (across times) one specified variable, while lowering resolution
rasterMeanDf=function(ncfile,varname,lowres.factor=4,align.with=0){
  require(raster)
#ncfile="~/Dropbox/cmems_mod_bal_bgc_anfc_P1M-m_1686434987026.nc";varname="chl";lowres.factor=4;align.with=ref
  raster_brick <- brick(ncfile,varname=varname)
  if (length(align.with)>10) { 
    raster_brick=resample(raster_brick,align.with,method="bilinear")
 #   lowres.factor=1
  }
  lores=aggregate(raster_brick, fact = lowres.factor, fun = mean)
#  plot(lores,asp=1)
  if(length(names(lores))>1){
    mean_raster <- calc(lores, fun = mean, na.rm = TRUE)
  } else {
    mean_raster=lores
  }
  pts = rasterToPoints(mean_raster, spatial = TRUE)
  pts.df  = data.frame(pts)[,1:3]
  names(pts.df)=c(varname,"lon","lat")
  pts.df[,2:3]=round(pts.df[,2:3],3)
  return(pts.df)
}

# function to read and average multi-variable nc files and convert them to dataframe,
# while aligning the raster to some reference (align.with argument)
ncread=function(filename,lowres.factor=4,align.with=10){
  require(ncdf4)
  nc=nc_open(filename)
  var_list <- names(nc$var)
  li=list()
  for(v in var_list){
    pts=rasterMeanDf(ncfile=filename,varname=v,lowres.factor,align.with)
    li[[v]]=pts[,1]
  }
  df=data.frame(do.call(cbind,li))
  df=cbind(pts[,2:3],df)
  return(df)
}

ref <- brick("~/Dropbox/t_sal_cmems_mod_bal_phy_anfc_P1M-m_1686562820182.nc",varname="so")
waterChem=ncread("~/Dropbox/waterChem.nc",lowres.factor=4,align.with = ref[[1]])
ts=ncread("~/Dropbox/t_sal_cmems_mod_bal_phy_anfc_P1M-m_1686562820182.nc",lowres.factor=4,align.with = ref[[1]])
names(ts)[3:4]=c("T","Sal")

dim(waterChem)
dim(ts)

table(waterChem$lon %in% ts$lon)
# shoudl be all TRUE because we requested alignment to ref[[1]]
# TRUE 
# 10426 

balt.raster=merge(waterChem,ts,by=c("lon","lat"))

# creating log-salinity to make it less skewed
balt.raster$logSal=log(balt.raster$Sal+2,2)

# plotting variables one by one, looking for ones showing reasonable variation
n=4
balt.raster[,4]=log(balt.raster[,4]+1,2)
plot(density(balt.raster[,n]))
XY=balt.raster[,1:2]
names(XY)=c("x","y")
ggplot(XY,aes(x,y,color=balt.raster[,n]))+geom_point(cex=1,pch=15)+theme_bw()+coord_equal()
goods=c(1,2,4,5,6,8,14,15,17)
balt.raster=balt.raster[,goods]
pairs(balt.raster[,-c(1,2)],pch=".",col=rgb(0,0,0,0.2))

names(balt.raster)
# "lon"    "lat"    "o2"     "chl"    "zooc"   "po4"    "dissic" "T"      "logSal"
# see baltic.chemVars.names to see what they mean
names(balt.raster)[4]="logChl"

rasters=balt.raster[,-c(1:2)]
save(rasters,XY,file="baltic.rasters.XY.RData")

#--------- finding values for the actual sampling locations

setwd('~/Dropbox/mega2019/idotea_2019') # change this to where your scp'd files are
meta=read.csv("Idotea environment 171005.csv")
head(meta)

new.meta=c();pop=3;i=1
for(pop in 1:nrow(meta)){
#  message(pop)
  closest=10;cl.point=0
  for (i in 1:nrow(balt.raster)){
    dpop=data.frame(rbind(meta[pop,2:3],balt.raster[i,1:2]))
    dd=dist(dpop)
    if(dd<closest){
      closest=dd
      cl.point=i
    }
  }
  message(paste(pop,cl.point,closest))
  new.meta=data.frame(rbind(new.meta,balt.raster[cl.point,]))
}

new.meta=data.frame(cbind(pop=meta$pop,new.meta))

plot(meta[,2:3])
points(new.meta[,2:3],col="red")

new.meta2=subset(new.meta,!(pop %in% c("SYL","LIL")))
dim(new.meta2)

meta=new.meta2
getwd()
save(meta,file="idotea_metadata_noSYLLIL.RData")


