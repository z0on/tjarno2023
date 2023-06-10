
#------- making the table of environmental conditions *for each sample*

# sample names
samples=dimnames(ibs)[[1]]
# extracting first three characters of the sample name (code for sampling location)
pop=substr(samples,1,3)
popstab=data.frame(cbind(pop=pop))
env=read.csv("Idotea environment 171005.csv")
env.xt=merge(popstab,env,by="pop",all.X=T)
row.names(env.xt)=samples
# the resulting table is ordered like IBS and ADMIXTURE data and has sample names as row names

# ------------ setting up colors for samples and sampling locations (base R, RColorBrewer)
#              (uses popstab from the code bit above)

library(RColorBrewer)
popcols = brewer.pal(11, "RdYlBu") 
popcols = colorRampPalette(cols)(length(unique(pop)))
pp=data.frame(cbind(pop=unique(pop),cols=popcols))
pp2=merge(popstab,pp,by="pop",all.X=T)
cols=pp2$cols
# pop: sampling locations ("population") for each sample
# cols: colors for each sample (per-location)
# popcols: colors for each location

#-------------- PCoA with spiders and ellipses for sampling locations
#               (uses pop, cols and popcols from the code bit above)

pp0=capscale(ibs~1)
# which PCs we want to plot?
pcs2plot=c(1,2)
scores=scores(pp0,scaling=1,display="sites",choices=pcs2plot)
plot(scores,pch=16,col=cols,asp=1)
ordispider(scores,groups=pop,col=popcols,label=F)
ordiellipse(scores,groups=pop,label=T,draw="polygon",col=popcols,cex=0.4)


