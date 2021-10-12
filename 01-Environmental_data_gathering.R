



library(tidyverse)
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing
library(heatwaveR) # For detecting MHWs
library(sf)
library(rnaturalearth)
library(fst)

# Function to detect marine heatwaves -------------------------------------

event_only <- function(df){
            require(heatwaveR)
            # First calculate the climatologies
            clim <- ts2clm(data = df, climatologyPeriod = c("1982-01-01", "2019-12-31"))
            # Then the events
            event <- detect_event(data = clim)
            # Return only the event metric dataframe of results
            return(event$event)
}

# function to convert longitude data from 180 to 360 degrees
conv = function(x) {
            ((360 + (x %% 360)) %% 360)
}

# Function to calculate slope and p-values

lin_fun <- function(ev) {
            mod1 <- glm(n ~ year, family = poisson(link = "log"), data = ev)
            # extract slope coefficient and its p-value
            tr <- data.frame(slope = summary(mod1)$coefficients[2,1],
                             p = summary(mod1)$coefficients[2,4])
            return(tr)
}



# Sea Surface Temperature -------------------------------------------------

# Loreto bounding box
xcoord_lor <- c(-111.5, -110.5)
xcoord2_lor <- conv(xcoord_lor)
ycoord_lor <- c(25.5, 26.5)

# La Paz bounding box
xcoord_lp <- c(-110.5, -109.5)
xcoord2_lp <- conv(xcoord_lp)
ycoord_lp <- c(24, 25)


# The information for the NOAA OISST data
rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")


# For the Mexican Pacific Marine Heatwaves

xcoord_mxPC <- c(-128, -98)
xcoord2_mxPC <- conv(xcoord_mxPC)
ycoord_mxPC <- c(15, 31)

# This function downloads and prepares data based on user provided start and end dates
OISST_sub_dl <- function(time_df){
            OISST_dat <- griddap(x = "ncdcOisst21Agg_LonPM180", 
                                 url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                                 time = c(time_df$start, time_df$end), 
                                 zlev = c(0, 0),
                                 latitude = ycoord_mxPC,
                                 longitude = xcoord_mxPC,
                                 fields = "sst")$data %>% 
                        mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
                        dplyr::rename(t = time, temp = sst) %>% 
                        select(lon, lat, t, temp) %>% 
                        na.omit()
}


# Date download range by start and end dates per year
dl_years <- data.frame(date_index = 1:5,
                       start = as.Date(c("1982-01-01", "1990-01-01", 
                                         "1998-01-01", "2006-01-01", "2014-01-01")),
                       end = as.Date(c("1989-12-31", "1997-12-31", 
                                       "2005-12-31", "2013-12-31", "2021-09-01")))

# Download all of the data with one nested request
# The time this takes will vary greatly based on connection speed
system.time(
            OISST_data <- dl_years %>% 
                        group_by(date_index) %>% 
                        group_modify(~OISST_sub_dl(.x)) %>% 
                        ungroup() %>% 
                        select(lon, lat, t, temp)
) # 636 seconds, ~127 seconds per batch

fst::write.fst(OISST_data, "data/env_data/full_oisst_pacific_data.fst")


OISST_data %>% 
            group_by(t) %>% 
            summarise(temp = mean(temp)) %>% 
            ggplot(aes(x = t, y = temp)) +
            geom_point()

# Detecting marine heatwaves in the Mexican Pacific ----------------------------------------------

# NB: One should never use ALL available cores, save at least 1 for other essential tasks

registerDoParallel(cores = 7)

system.time(
            MHW_result <- plyr::ddply(.data = OISST_data, .variables = c("lon", "lat"), .fun = event_only, .parallel = TRUE)
) # 46.50 seconds



# summarise the number of unique longitude, latitude and year combination for heatwaves events:
event_freq <- MHW_result %>% 
            mutate(year = lubridate::year(date_start)) %>% 
            group_by(lon, lat, year) %>% 
            summarise(n = n(), .groups = "drop")
head(event_freq)

# create complete grid for merging with:
sst_grid <- OISST_data %>% 
            select(lon, lat, t) %>% 
            mutate(t = lubridate::year(t)) %>% 
            dplyr::rename(year = t) %>% 
            distinct()

# and merge:
MHWs_n <- left_join(sst_grid, event_freq, by = c("lon", "lat", "year")) %>% 
            mutate(n = ifelse(is.na(n), 0, n))



# creating trend data 

MHWs_Trend <- plyr::ddply(MHWs_n, c("lon", "lat"), lin_fun, .parallel = T)
MHWs_Trend$pval <- cut(MHWs_Trend$p, breaks = c(0, 0.001, 0.01, 0.05, 1))
head(MHWs_Trend)

OISST_mean <- OISST_data %>% 
            mutate(year = lubridate::year(t)) %>% 
            group_by(lon, lat, year) %>% 
            summarise(temp = mean(temp))


saveRDS(MHWs_Trend, file = "data/env_data/MHWs_nTrend.Rds") # Mexican Pacific SST trends
saveRDS(MHWs_n, file = "data/env_data/MHWs_n.Rds") # Mexican Pacific SST heatwaves n


# For the Gulf of California study areas


# Central
OISST_dat_lor <- griddap(x = "ncdcOisst21Agg_LonPM180", 
                     url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                     time = c("1982-01-01", "2021-09-01"), 
                     zlev = c(0, 0),
                     latitude = ycoord_lor,
                     longitude = xcoord_lor,
                     fields = "sst")$data %>% 
            mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
            dplyr::rename(t = time, temp = sst) %>% 
            select(lon, lat, t, temp) %>% 
            na.omit()

# Transitional


OISST_dat_lp <- griddap(x = "ncdcOisst21Agg_LonPM180", 
                         url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                         time = c("1982-01-01", "2021-09-01"), 
                         zlev = c(0, 0),
                         latitude = ycoord_lp,
                         longitude = xcoord_lp,
                         fields = "sst")$data %>% 
            mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
            dplyr::rename(t = time, temp = sst) %>% 
            select(lon, lat, t, temp) %>% 
            na.omit()


fst::write_fst(OISST_dat_lor, "data/env_data/OISST_loreto.fst")
fst::write_fst(OISST_dat_lp, "data/env_data/OISST_lapaz.fst")








# Chlor_a -----------------------------------------------------------------

## Setting download parameters
ERDDAP_Node="https://oceanwatch.pifsc.noaa.gov/erddap/"


# seawifs 

dataInfo <- rerddap::info('sw_chla_monthly_2018_0', url = ERDDAP_Node) 
var = dataInfo$variable$variable_name

# This function downloads and prepares data based on user provided start and end dates

sw_dat <- griddap(x = "sw_chla_monthly_2018_0", 
                  url = ERDDAP_Node, 
                  time = c('1997-12-01', '2010-12-01'), 
                  latitude = ycoord_lor,
                  longitude = xcoord2_lor,
                  fields = var[1])$data %>% 
            mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
            dplyr::rename(t = time, chla = chlor_a) %>% 
            select(lon, lat, t, chla) %>% 
            mutate(lon = lon - 360)


# MODIS

dataInfo <- rerddap::info('aqua_chla_monthly_2018_0', url = ERDDAP_Node) 
var = dataInfo$variable$variable_name

mod_dat <- griddap(x = "aqua_chla_monthly_2018_0", 
                   url = ERDDAP_Node, 
                   time = c('2002-07-16', '2021-08-16'), 
                   latitude = ycoord_lor,
                   longitude = xcoord2_lor,
                   fields = var[1])$data %>% 
            mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
            dplyr::rename(t = time, chla = chlor_a) %>% 
            select(lon, lat, t, chla) %>% 
            mutate(lon = lon - 360)



# VIIRS 

dataInfo <- rerddap::info('noaa_snpp_chla_monthly', 
                          url = ERDDAP_Node) 
var=dataInfo$variable$variable_name

viirs_dat <- griddap(x = "noaa_snpp_chla_monthly", 
                     url = ERDDAP_Node, 
                     time = c('2012-01-02', '2021-08-16'), 
                     latitude = ycoord_lor,
                     longitude = xcoord2_lor,
                     fields = var[1])$data %>% 
            mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
            dplyr::rename(t = time, chla = chlor_a) %>% 
            select(lon, lat, t, chla) %>% 
            mutate(lon = lon - 360)



chla_data_lor <- rbind(sw_dat, mod_dat, viirs_dat) %>% 
            group_by(t) %>% 
            summarise(chla = mean(chla, na.rm = T))

## La Paz



# seawifs 

dataInfo <- rerddap::info('sw_chla_monthly_2018_0', url = ERDDAP_Node) 
var = dataInfo$variable$variable_name

# This function downloads and prepares data based on user provided start and end dates

sw_dat <- griddap(x = "sw_chla_monthly_2018_0", 
                  url = ERDDAP_Node, 
                  time = c('1997-12-01', '2010-12-01'), 
                  latitude = ycoord_lp,
                  longitude = xcoord2_lp,
                  fields = var[1])$data %>% 
            mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
            dplyr::rename(t = time, chla = chlor_a) %>% 
            select(lon, lat, t, chla) %>% 
            mutate(lon = lon - 360)


# MODIS

dataInfo <- rerddap::info('aqua_chla_monthly_2018_0', url = ERDDAP_Node) 
var = dataInfo$variable$variable_name

mod_dat <- griddap(x = "aqua_chla_monthly_2018_0", 
                   url = ERDDAP_Node, 
                   time = c('2002-07-16', '2021-08-16'), 
                   latitude = ycoord_lp,
                   longitude = xcoord2_lp,
                   fields = var[1])$data %>% 
            mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
            dplyr::rename(t = time, chla = chlor_a) %>% 
            select(lon, lat, t, chla) %>% 
            mutate(lon = lon - 360)



# VIIRS 

dataInfo <- rerddap::info('noaa_snpp_chla_monthly', 
                          url = ERDDAP_Node) 
var=dataInfo$variable$variable_name

viirs_dat <- griddap(x = "noaa_snpp_chla_monthly", 
                     url = ERDDAP_Node, 
                     time = c('2012-01-02', '2021-08-16'), 
                     latitude = ycoord_lp,
                     longitude = xcoord2_lp,
                     fields = var[1])$data %>% 
            mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
            dplyr::rename(t = time, chla = chlor_a) %>% 
            select(lon, lat, t, chla) %>% 
            mutate(lon = lon - 360)



chla_data_lp <- rbind(sw_dat, mod_dat, viirs_dat) %>% 
            group_by(t) %>% 
            summarise(chla = mean(chla, na.rm = T))



chla_data_lor$region <- "Loreto"
chla_data_lp$region <- "La Paz"


chla_data <- rbind(chla_data_lor, chla_data_lp)

chla_data %>% 
            ggplot(aes(x = t, y = chla)) +
            geom_point(col = "red") +
            geom_line(col = "red") +
            geom_smooth(method = "gam") +
            facet_grid(~region)



write_fst(chla_data, "data/env_data/chla_data.fst")


