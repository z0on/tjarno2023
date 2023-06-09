# Introduction to seascape/landscape genomics
#### - exploring and displaying population structure on a map
#### - finding genotype-environment associations
#### - mapping divergent adaptation
####  

## Installing stuff

These exercises are exclusively in R. If you are completely unfamiliar with R, consider working through chapters 1-7 here: http://swcarpentry.github.io/r-novice-gapminder/

First, [install R and Rstudio](https://rstudio-education.github.io/hopr/starting.html), unless you have it already.

The exercises will require quite a few R packages; so open the file **trarno2023_installations.R** in Rstudio and follow it. Some of the packages might require additional persuasion to install, so make sure to DO THIS BEFORE THE CLASS.

First, install **devtools**. 
```R
install.packages("devtools")
```
This may require additional installations outside R. Hopefully they will happen automatically, if not, see [here](https://www.r-project.org/nosvn/pandoc/devtools.html).
Then:
```R
# for multivariate analysis:
install.packages(c("vegan","dplyr","raster","fields","ggplot2","pheatmap"))

# for maps:
install.packages(c("remotes","cowplot", "googleway", "ggrepel", "ggspatial","libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
remotes::install_github("ropensci/rnaturalearthhires")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")

# tess3r: for making admixture maps
devtools::install_github("bcm-uga/TESS3_encho_sen")
```
This one is going to be the trickiest:
```R
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
```
If this one fails, chances are, you need to install **gfortran** first, FOR YOUR SYSTEM from here:
https://gcc.gnu.org/wiki/GFortranBinaries or, for Mac, https://github.com/fxcoudert/gfortran-for-macOS/releases
On a Mac you would also need to point your Rstudio compiler to the location **gfortran** is installed at, by creating/modifying the file *~/.R/Makevars*. The following spell in **Terminal** should work:
```sh
echo "FC = /usr/local/bin/gfortran
F77 = /usr/local/gfortran
FLIBS = -L/usr/local/gfortran/lib" >> ~/.R/Makevars
```
To check if everything was intalled correctly, do this in Rstudio and see if all packages are loaded without errors.

```R
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
```
If some packages are still missing, try *install.packages( "[package name]" )* command in Rstudio again. If it fails, on a Mac, try opening Termial, saying "R" in command-line, and repeating the *install.packages* commands after you see the R prompt. If that fails too, create a github issue here.

## Agaricia agaricites in the FL Keys
##
This dataset is from the seascape genomics project led by Kristina Black, who put together many of the analyses I am going to show you. It contains 250 samples of the stony coral *Agaricia agaricites* collected from 63 sites in Florida. The samples were sequenced with 2b-RAD, genotyped [de-novo](https://github.com/z0on/2bRAD_denovo) and genetic disstances between samples were computed in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) as identity-by-state (IBS) based on single-read resampling. The AD@MIXTURE proportions were computed by [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix).

>We recommend PCAngsd for ADMIXTURE analysis these days. 
>The method also computes SNP covariance matrix that 
>can be converted into correlations and used instead of the IBS matrix.
>We prefer IBS since it is the most assumption-free approach to computing genetic distances, 
>robust to variation in coverage across samples.

There are also environmental data for these same sites, as well as for the whole region. Some of these data come from satellite measurements (temperature, salinity, turbidity, chlorophyll), other data - from [SERC Water Quality Monitoring Network](http://serc.fiu.edu/wqmnetwork/).
###3
The script *Agaricia_tjarno2023.R* shows how to:
- explore the data and remove problematic samples (functions *hclust*, *pheatmap*, *cutree*)
- plot isolation-by-distance graph to see if there are multiple clines, signifying distinct genetic clusters
- explore genetic structure using principal coordinate analysis (PCoA, function *capscale*), plot ADMIXTURE bar chart using package *tess3R*, and color the PCoA ordination plots according to ADMIXTURE.
- project ADMIXTURE clusters on the map (package *trss3r*) to see where cluster representatives are predominantly found
- 
