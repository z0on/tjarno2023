
library(ggplot2)

?raster::calc

# function to read and average (across times) one specified variable, while lowering resolution
rasterMeanDf=function(ncfile,varname,lowres.factor=4,align.with=0,func="mean"){
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
#    stop("nothing to average")
    mean_raster=lores
  }
  pts = rasterToPoints(mean_raster, spatial = TRUE)
  pts.df  = data.frame(pts)[,1:3]
  names(pts.df)=c(varname,"lon","lat")
  pts.df[,2:3]=round(pts.df[,2:3],3)
  return(pts.df)
}

rasterSDDf=function(ncfile,varname,lowres.factor=4,align.with=0){
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
    mean_raster <- calc(lores, fun = sd, na.rm = TRUE)
  } else {
    stop("nothing to calculate SD of")
#    mean_raster=lores
  }
  pts = rasterToPoints(mean_raster, spatial = TRUE)
  pts.df  = data.frame(pts)[,1:3]
  names(pts.df)=c(varname,"lon","lat")
  pts.df[,2:3]=round(pts.df[,2:3],3)
  return(pts.df)
}

rasterLSD1Df=function(ncfile,varname,lowres.factor=4,align.with=0){
  require(raster)
  #ncfile="~/Dropbox/cmems_mod_bal_bgc_anfc_P1M-m_1686434987026.nc";varname="chl";lowres.factor=4;align.with=ref
  raster_brick <- brick(ncfile,varname=varname)
  if (length(align.with)>10) { 
    raster_brick=resample(raster_brick,align.with,method="bilinear")
    #   lowres.factor=1
  }
  lores=aggregate(raster_brick, fact = lowres.factor, fun = mean)
  #  plot(lores,asp=1)
  sdlog1=function(x,na.rm=TRUE){ return(sd(log(x+1,2))) }
  if(length(names(lores))>1){
    mean_raster <- calc(lores, fun = sdlog1, na.rm = TRUE)
  } else {
    stop("nothing to calculate log-SD of")
    #    mean_raster=lores
  }
  pts = rasterToPoints(mean_raster, spatial = TRUE)
  pts.df  = data.frame(pts)[,1:3]
  names(pts.df)=c(varname,"lon","lat")
  pts.df[,2:3]=round(pts.df[,2:3],3)
  return(pts.df)
}

# function to read, apply function 'func' ("mean","sd",or "lsd1" = log(x+1,2)) across timepoints, and convert them to dataframe
# for multi-variable nc files,
# while aligning the raster to some reference ('align.with' argument)
ncread=function(filename,lowres.factor=4,align.with=10,func="mean"){
  require(ncdf4)
  nc=nc_open(filename)
  var_list <- names(nc$var)
  li=list()
  if(func=="mean") { 
    for(v in var_list){
      pts=rasterMeanDf(ncfile=filename,varname=v,lowres.factor,align.with)
      li[[v]]=pts[,1]
    }
  } else {
    if(func=="sd"){
      for(v in var_list){
        pts=rasterSDDf(ncfile=filename,varname=v,lowres.factor,align.with)
        li[[v]]=pts[,1]
      }
    } else {
      if(func=="lsd1"){
        for(v in var_list){
          pts=rasterLSD1Df(ncfile=filename,varname=v,lowres.factor,align.with)
          li[[v]]=pts[,1]
        }
      } else { stop("unknown func, should be mean, sd, or lsd1")}
    }
  }
  df=data.frame(do.call(cbind,li))
  df=cbind(pts[,2:3],df)
  return(df)
}

#-------------------------

wcfile="~/Dropbox/waterChem_2000_cmems_mod_bal_bgc_my_P1M-m_1686581987299.nc"
tsfile="~/Dropbox/t_sal_2000_cmems_mod_bal_phy_my_P1M-m_1686581523152.nc"

ref <- brick(tsfile,varname="so")
waterChem=ncread(wcfile,lowres.factor=4,align.with = ref[[1]],func="mean")
waterChem.sd=ncread(wcfile,lowres.factor=4,align.with = ref[[1]],func="sd")
waterChem.lsd=ncread(wcfile,lowres.factor=4,align.with = ref[[1]],func="lsd1")
ts=ncread(tsfile,lowres.factor=4,align.with = ref[[1]],func="mean")
ts.sd=ncread(tsfile,lowres.factor=4,align.with = ref[[1]],func="sd")
ts.lsd=ncread(tsfile,lowres.factor=4,align.with = ref[[1]],func="lsd1")
names(ts)[3:4]=c("T","Sal")
names(ts.sd)[3:4]=c("T","Sal")
names(ts.lsd)[3:4]=c("T","Sal")

dim(waterChem)
dim(ts)
dim(waterChem.sd)
dim(ts.sd)

table(waterChem$lon %in% ts$lon)
# shoudl be all TRUE because we requested alignment to ref[[1]]
# TRUE 
# 10426 

balt.raster=merge(waterChem,ts,by=c("lon","lat"))
balt.raster.sd=merge(waterChem.sd,ts.sd,by=c("lon","lat"))
balt.raster.lsd=merge(waterChem.lsd,ts.lsd,by=c("lon","lat"))
dim(balt.raster)

logthese=c(5,7,8,10,11)
for (i in logthese){
  balt.raster[,i]=log(balt.raster[,i]+1,2)
  names(balt.raster)[i]=paste("log.",names(balt.raster)[i],sep="")
}
balt.raster.sd[,logthese]=balt.raster.lsd[,logthese]
sdnames=paste(names(balt.raster.sd[,-c(1,2)]),".sd",sep="")
sdnames[logthese]=paste(names(balt.raster.sd[,logthese]),".sdlog",sep="")
names(balt.raster.sd)[-c(1,2)]=sdnames

goods=c(1,2,4,5,6,7,8,9,10,11,12,13)
balt.raster=balt.raster[,goods]
balt.raster.sd=balt.raster.sd[,goods]

#pairs(balt.raster[,3:ncol(balt.raster)],pch=".")
#pairs(balt.raster.sd[,3:ncol(balt.raster.sd)],pch=".")


# # plotting variables one by one, looking for ones showing reasonable variation
# n=9
# names(balt.raster)[n]
# hist(balt.raster.sd[,n])
# hist(balt.raster.lsd[,n])
# #plot(density(balt.raster[,n]))
# XY=balt.raster[,1:2]
# names(XY)=c("x","y")
# ggplot(XY,aes(x,y,color=balt.raster.sd[,n]))+geom_point(cex=1,pch=15)+theme_bw()+coord_equal()+scale_color_viridis()
# ggplot(XY,aes(x,y,color=balt.raster.lsd[,n]))+geom_point(cex=1,pch=15)+theme_bw()+coord_equal()+scale_color_viridis()
# ggplot(XY,aes(x,y,color=log(balt.raster.sd[,n])))+geom_point(cex=1,pch=15)+theme_bw()+coord_equal()+scale_color_viridis()
# 

rasters=merge(balt.raster,balt.raster.sd,by=c("lat","lon"))

dim(rasters)
XY=rasters[,2:1]
head(XY)
#names(XY)=c("x","y")
rasters=rasters[,-c(1:2)]

names(rasters)
# [1] "o2"         "log.chl"    "zsd"        "log.no3"    "log.po4"    "ph"        
# [7] "log.nh4"    "log.nppv"   "T"          "Sal"        "o2.sd"      "chl.sd"    
# [13] "zsd.sd"     "chl.sdlog"  "po4.sd"     "no3.sdlog"  "po4.sdlog"  "nppv.sd"   
# [19] "nh4.sdlog"  "nppv.sdlog" 

save(rasters,XY,file="baltic.rasters.XY.2000.RData")

# goodvars=names(rasters)
# save(goodvars,logthese,file="baltic_environmentalVariableNames.RData")

getwd()
#--------- finding values for the actual sampling locations: idotea

setwd('~/Dropbox/mega2019/idotea_2019') # change this to where your scp'd files are
ll=load("idotea_37pops_rasters_2023_clean.RData")
head(meta)
load("baltic.rasters.XY.2014.RData")
head(XY)
new.meta=c();pop=3;i=1
for(pop in 1:nrow(meta)){
#  message(pop)
  closest=10;cl.point=0
  for (i in 1:nrow(XY)){
    dpop=data.frame(rbind(meta[pop,2:3],XY[i,1:2]))
    dd=dist(dpop)
    if(dd<closest){
      closest=dd
      cl.point=i
    }
  }
  message(paste(pop,cl.point,closest))
  new.meta=data.frame(rbind(new.meta,rasters[cl.point,]))
}

new.meta2=data.frame(cbind(meta[,1:3],new.meta))

plot(meta[,2:3])
points(new.meta2[,2:3],col="red")
head(new.meta2)
pop37=data.frame(pop=meta.xt$pop)
head(pop37)
meta.xt=merge(pop37,new.meta2,by="pop",all.X=T)
head(meta.xt)

meta=new.meta2
save(meta,bams,admix,ibs,meta.xt,file="idotea_37pops_env2014_clean.RData")
#load("idotea_37pops_env2014_clean.RData")

#--------- finding values for the actual sampling locations: Stefanie Ries eelgrass

setwd('~/Dropbox/mega2019/idotea_2019') # change this to where your scp'd files are
meta=read.table("Ries_eelgrass_sampling_sites_coordinates_Baltic_Sea.csv",sep=";",header=T)
names(meta)[5:6]=c("lat","lon")

new.meta=c();pop=3;i=1
for(pop in 1:nrow(meta)){
  #  message(pop)
  closest=10;cl.point=0
  for (i in 1:nrow(rasters)){
    dpop=data.frame(rbind(meta[pop,c(5,6)],XY[i,1:2]))
    dd=dist(dpop)
    if(dd<closest){
      closest=dd
      cl.point=i
    }
  }
  message(paste(pop,cl.point,closest))
  new.meta=data.frame(rbind(new.meta,rasters[cl.point,]))
}

new.meta2=data.frame(cbind(meta[,1:4],XY[row.names(new.meta),],new.meta))
head(new.meta2)
plot(meta[,c(6,5)])
points(new.meta2[,5:6],col="red")
map(coasts,add = T, interior = F)

eelgrass.meta=new.meta2
save(eelgrass.meta,file="eelgrass_metadata.RData")
