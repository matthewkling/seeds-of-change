library(terra)
# library(gdalUtilities)
# library(raster)
library(tidyverse)


# SoilGrids v1 (2017) ##############

"https://files.isric.org/soilgrids/former/2017-03-10/data/AWCh1_M_sl1_250m_ll.tif"

vars <- tibble(var = c("AWCh1", "AWCh2", "AWCh3", "AWCtS", "BLDFIE", "CECSOL", "CLYPPT", "CRFVOL",
                       "OCDENS", "OCSTHA", "ORCDRC", "PHIHOX", "PHIKCL", "SLTPPT", "SNDPPT", "TEXMHT", "WWP")) %>%
      expand_grid(depth = c(1, 3, 5, 7))

dl <- function(var, depth){
      x <- paste0(var, "_M_sl", depth, "_250m_ll.tif")
      dest <- paste0("~/data/soilgrids/v1/250m/", basename(x))
      if(file.exists(dest)) return("skipped")
      base <- "https://files.isric.org/soilgrids/former/2017-03-10/data/"
      options(timeout = 100000)
      download.file(paste0(base, x), dest)
      download.file(paste0(base, x, ".xml"), paste(dest, ".xml"))
}

pmap(vars, dl)


# crop & project to CA

r <- rast("assets/climate/ensemble 1981-2010 NA AET.tif")[[1]] # template (chelsa)

f <- list.files("~/data/soilgrids/v1/250m", full.names = T)
f <- f[!grepl("xml", f)]

f %>% map(function(x){
      x %>% 
            rast() %>% 
            crop(r) %>% 
            project(y = r, method = "average") %>%
            mask(r) %>%
            writeRaster(paste0("~/data/soilgrids/v1/250m_California/", basename(x)))
})

