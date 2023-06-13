# Introduction to seascape/landscape genomics
#### - exploring and displaying population structure on a map
#### - finding genotype-environment associations
#### - mapping divergent adaptation
####  

## Readings
- [short and sweet intro into decision trees and random forest](https://towardsdatascience.com/understanding-random-forest-58381e0602d2)
- `death00_classification_regression_trees.pdf` : great intro into regression trees, for ecological data
- `ellis_gradient_forest.pdf` : dense but very comprehensive outline of the gradient forest approach
- `fitzpatrick15_genomic_vulnerability.pdf` : introduces the idea of comparing current and future adaptation maps, to learn where our species will be the most in trouble
- `bay_songbird_science.pdf` : one of the best examples so far of putting all the above to work.

## Installing stuff

These exercises are exclusively in R. If you are completely unfamiliar with R, consider working through chapters 1-7 here: http://swcarpentry.github.io/r-novice-gapminder/

First, [install R and Rstudio](https://rstudio-education.github.io/hopr/starting.html), unless you have it already.

The exercises will require quite a few R packages; so open the file *`trarno2023_installations.R`* in Rstudio and follow it. Some of the packages might require additional persuasion to install, so make sure to DO THIS BEFORE THE CLASS.

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
cd
mkdir .R
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

## Seascape genomics of *Agaricia agaricites* in the FL Keys
###
This dataset is from the seascape genomics project led by Kristina Black, who put together many of the analyses I am going to show you. It contains 250 samples of the stony coral *Agaricia agaricites* collected from 63 sites in Florida. 
![Keys seascape](FL_seascape_agaricia.png)
The samples were sequenced with 2b-RAD, genotyped [de-novo](https://github.com/z0on/2bRAD_denovo) and genetic distances between samples were computed in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) as identity-by-state (IBS) based on single-read resampling. The ADMIXTURE clusters were delineated using [ngsAdmix](http://www.popgen.dk/software/index.php/NgsAdmix).

>We recommend *PCAngsd* for ADMIXTURE analysis these days. That method also computes SNP covariance matrix that can be converted into correlations and used instead of the IBS matrix. We prefer IBS since it is the most assumption-free approach to computing genetic distances, robust to variation in coverage across samples.

There are also environmental data for these same sites, as well as for the whole region. Some of these data come from satellite measurements (temperature, salinity, turbidity, chlorophyll), other data - from [SERC Water Quality Monitoring Network](http://serc.fiu.edu/wqmnetwork/).
###
The script *`Agaricia_tjarno2023.R`* shows how to:
- explore the data and remove problematic samples (functions *hclust*, *pheatmap*, *cutree*)
- plot isolation-by-distance graph to see if there are multiple clines, signifying distinct genetic clusters
- explore genetic structure using principal coordinate analysis (PCoA, function *capscale*), plot ADMIXTURE bar chart using package *tess3R*, and color the PCoA ordination plots according to ADMIXTURE cluster assignments.
- project ADMIXTURE clusters on the map (package *tess3r*) to see where cluster representatives are predominantly found
- perform a gradient forest analysis of genotype-environment associations (package *gradientForest*)
- predict and map areas of divergent adaptation across the whole region.

After going through *`Agaricia_tjarno2023.R`*, you will be asked to investigate adaptation *within* each of the three distinct *A.agaricites* lineages and compare them, all on your own. The script has some scrap code in the end that will be helpful.

## *Idotea baltica* dataset
###
*Idotea baltica* is a small, extremely cute isopod found on *Fucus* algae in the intertidal. The Baltic-wide 2b-RAD genotypoing dataset for this charismatic microfauna representative has been generated by our host, Pierre De Wit. 

![*Idotea baltica*](isopoda_idotea_balthica_01-10-15_1.jpg)

Unlike corals, *I. baltica* has a very narrow dispersal range (it is a marsupial brooding its young in a pouch, believe it or not). Because of this, it shows strong isolation-by-distance, and makes the most spectacular ADMIXTURE-on-map plot (let's make one!). But does it also adapt genetically to the variation in conditions across the Baltic sea? Gradient forest to the rescue! 

The environmental dataset used here (`meta` for sampling locations, `meta.xt` for individual samples, and `rasters` for the whole Baltic Sea) were donwloaded from [Copernicus](https://data.marine.copernicus.eu/products?q=baltic). I took temperature and salinity from [Baltic Sea Physics Analysis](https://data.marine.copernicus.eu/product/BALTICSEA_ANALYSISFORECAST_PHY_003_006/description) and the rest from [Biogeochemistry](https://data.marine.copernicus.eu/product/BALTICSEA_ANALYSISFORECAST_BGC_003_007/description). The donloaded `.nc` files (monthly averages and log2-SDs for two years, Jan 16, 2021 - Jan 16, 2023) were grid-aligned, averaged, and transformed into dataframes using script `baltic.variables.R`. I retained only the variables showing reasonable variation acoss Baltic, and log2-transformed mean chlorophyll and salinity measures (all SDs were log2-transformed). The doc with additional info on the variables is `baltic.chemVars.names.docx`.

> Hint: The GF analysis must be conditional on the strong genetic structure due to isolation-by-distance. This structure should be reasonably accounted for by the first two principal coordinates.

Data for *I. baltica*:
- `idotea_37pops_env2014_clean.RData`: a bundle of everything you need for PCoA, admixture, and initial gradient forest analysis: objects `bams`, `ibs`, `meta` (metadata for 37 sampling locations),  `meta.xt` (for 584 individual samples), `admix`. The files are already aligned and free of clonal duplicates, so go straight into PCoA.
- `baltic.rasters.XY.2014.RData`: data over 2013-2015 (`rasters`, `XY`) for projecting gradient forest model.
- `baltic.rasters.XY.2000.RData`: data over 1999-2001 (`rasters`, `XY`) for projecting gradient forest model.
- `Baltic_polygon.txt`: table containing longitude and latitude of points defining the Baltic polygon (for plotting ADMIXTURE clusters on map)
- `idotea_scrapCode.R` : some potentially helpful bits of code, for example for more advanced PCoA visualization

