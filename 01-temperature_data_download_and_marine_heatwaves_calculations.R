#'This scripts download OISST data and save them as outputs.
#'Once OISST data are downloaded, it calculates Marine Heatwaves and then store results in a file.  
#'There are portion with parallel processing so take care.

# Loading libraries -------------------------------------------------------

library(tidyverse)
library(lubridate) 
library(tidync) 
library(rerddap) 
library(doParallel) 
library(heatwaveR) 
library(fst)


# Loading custom functions ------------------------------------------------

event_only <- function(df){
            require(heatwaveR)
            # First calculate the climatologies
            clim <- ts2clm(data = df, climatologyPeriod = c("1982-01-01", "2019-12-31"))
            # Then the events
            event <- detect_event(data = clim)
            # Return only the event metric dataframe of results
            return(event$event)
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

# The information for the NOAA OISST data
rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")


# Coordinates for the Mexican Pacific (extended)
xcoord_mxPC <- c(-128, -98)
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

fst::write.fst(OISST_data, "data/full_oisst_pacific_data.fst")


OISST_data %>% 
            group_by(t) %>% 
            summarise(temp = mean(temp)) %>% 
            ggplot(aes(x = t, y = temp)) +
            geom_point()

# Detecting marine heatwaves in the Mexican Pacific ----------------------------------------------

# NB: One should never use ALL available cores, save at least 1 for other essential tasks


registerDoParallel(cores = detectCores() - 1)

system.time(
            MHW_result <- plyr::ddply(.data = OISST_data, .variables = c("lon", "lat"), .fun = event_only, .parallel = TRUE)
) # 46.50 seconds



# summaries the number of unique longitude, latitude and year combination for heatwaves events:
event_freq <- MHW_result %>% 
            mutate(year = lubridate::year(date_start)) %>% 
            group_by(lon, lat, year) %>% 
            summarise(n = n(), .groups = "drop")


# create complete grid for merging with:
sst_grid <- OISST_data %>% 
            select(lon, lat, t) %>% 
            mutate(t = lubridate::year(t)) %>% 
            dplyr::rename(year = t) %>% 
            distinct()

# and merge:
MHWs_n <- left_join(sst_grid, event_freq, by = c("lon", "lat", "year")) %>% 
            mutate(n = ifelse(is.na(n), 0, n))


## Saving count data

saveRDS(MHWs_n, file = "data/MHWs_n.Rds") # Mexican Pacific SST heatwaves n


# Creating trend data 

MHWs_Trend <- plyr::ddply(MHWs_n, c("lon", "lat"), lin_fun, .parallel = T)
MHWs_Trend$pval <- cut(MHWs_Trend$p, breaks = c(0, 0.001, 0.01, 0.05, 1))

## Saving trend data

saveRDS(MHWs_Trend, file = "data/MHWs_nTrend.Rds") # Mexican Pacific SST trends


# Extracting the sst for the study area of the rocky reefs ----------------


sst_data <-OISST_data %>% 
            filter(between(lon, -111.3, -109.4)) %>% 
            filter(between(lat, 22.8, 26.2)) 

sst_data_sf <- st_as_sf(sst_data, coords = c("lon", 'lat'), crs = 4326, remove = F)

st_area <- read_sf("shp/study_area.shp")

results <- st_intersection(sst_data_sf, st_area) %>% 
            as.data.frame() %>% 
            select(-geometry)

write_fst(results, "data/sst_extract_in_the_study_area.fst")

# End of script -----------------------------------------------------------




