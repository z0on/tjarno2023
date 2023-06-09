install.packages("devtools")
#this may require additional installations (see here: https://www.r-project.org/nosvn/pandoc/devtools.html)
#hopefully they will happen automatically

# for multivariate analysis:
install.packages(c("vegan","dplyr","raster","fields","ggplot2","pheatmap"))

# for maps:
install.packages(c("remotes","cowplot", "googleway", "ggrepel", "ggspatial","libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
remotes::install_github("ropensci/rnaturalearthhires")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")

devtools::install_github("bcm-uga/TESS3_encho_sen")

install.packages("gradientForest", repos="http://R-Forge.R-project.org")
#If this one fails, chances are, you need to install gfortran first, FOR YOUR SYSTEM from here:
#https://gcc.gnu.org/wiki/GFortranBinaries or, for Mac, https://github.com/fxcoudert/gfortran-for-macOS/releases
#On a Mac you would also need to point your Rstudio compiler to the location gfortran is installed at.
#The following spell in Terminal should work (just remove # symbols anad paste into Terminal)
#echo "FC = /usr/local/bin/gfortran
#F77 = /usr/local/gfortran
#FLIBS = -L/usr/local/gfortran/lib"  >>~/.R/Makevars

