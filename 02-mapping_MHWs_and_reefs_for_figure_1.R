
# Loading packages --------------------------------------------------------

library(tidyverse)
library(sf)
library(fst)
library(rnaturalearth)
library(raster)
library(patchwork)


# Loading data ------------------------------------------------------------

## Loading data 
sst_data <- read_fst("data/study_area_sst_extract.fst") 

# MHWs trends and p-values
OISST_n <- readRDS('data/MHWs_n.Rds')
OISST_nTrend <- readRDS('data/MHWs_nTrend.Rds')

# Ecoregion shapefile
ecoreg <- st_read('shp/GOC_ecoregions.shp')

# Names wrangling
ecoreg$region <- factor(c("Northern ecoregion", "Transitional ecoregion", "Subtropical ecoregion"), 
                        levels = c("Northern ecoregion", "Transitional ecoregion", "Subtropical ecoregion"))

# Mexico basemap 
spdf_mx <- st_transform(st_as_sf(ne_countries(country = 'mexico', scale = "large")), crs = 4326)

# Reef coordinates
reefs <- readRDS("data/coordinates_of_rocky_reefs_consistently_monitored.RDS")

# Converting data to sf and WGS84 projection
reefs_sf <- st_as_sf(reefs, coords = c("Longitude", "Latitude"), crs = 4326)



# Plotting ----------------------------------------------------------------
# The base map
map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>% 
            dplyr::rename(lon = long)


(
            map_slope <- ggplot() +
                        geom_rect(
                                    data = OISST_nTrend,
                                    size = 0.2,
                                    fill = NA,
                                    aes(
                                                x = lon,
                                                y = lat,
                                                xmin = lon - 0.1,
                                                xmax = lon + 0.1,
                                                ymin = lat - 0.1,
                                                ymax = lat + 0.1
                                    )
                        ) +
                        geom_raster(
                                    data = OISST_nTrend,
                                    aes(x = lon, y = lat, fill = slope),
                                    interpolate = FALSE,
                                    alpha = 0.9
                        ) +
                        scale_fill_gradient2( 
                                    name = "MHWs/year (slope)",
                                    high = "red",
                                    mid = "white",
                                    low = "darkblue",
                                    midpoint = 0,
                                    guide = guide_colourbar(direction = "horizontal",
                                                            title.position = "top",
                                                            ticks.colour = "black"
                                    )
                        ) +
                        geom_polygon(
                                    data = map_base,
                                    aes(x = lon, y = lat, group = group),
                                    colour = NA,
                                    fill = "grey80"
                        ) +
                        geom_sf(
                                    data = spdf_mx,
                                    colour = "black",
                                    fill = "gray99",
                                    alpha = .2
                        ) +
                        coord_sf(
                                    xlim = c(-116, -105),
                                    ylim = c(21, 31),
                                    clip = "on"
                        ) +
                        labs(x = "", y = "") +
                        theme_bw() +
                        theme(panel.background = element_rect(fill = "white"), 
                              axis.text.x = element_text(angle = 90),
                              legend.position = "right", 
                              legend.margin = margin(),
                              legend.text = element_text(angle = 90, vjust = .5))
)

(
            map_p <- ggplot() +
                        geom_raster(
                                    data = OISST_nTrend,
                                    aes(x = lon, y = lat, fill = pval),
                                    interpolate = FALSE
                        ) +
                        scale_fill_manual(
                                    breaks = c(
                                                "(0,0.001]",
                                                "(0.001,0.01]",
                                                "(0.01,0.05]",
                                                "(0.05,0.1]",
                                                "(0.1,0.5]",
                                                "(0.5,1]"
                                    ),
                                    values = c("black", "grey40", "grey60",
                                               "grey80", "grey95", "white"),
                                    name = "p-value"
                        ) +
                        geom_polygon(
                                    data = map_base,
                                    aes(x = lon, y = lat, group = group),
                                    colour = NA,
                                    fill = "grey80"
                        ) +
                        geom_sf(
                                    data = spdf_mx,
                                    colour = "black",
                                    fill = "gray99",
                                    alpha = .2
                        ) +
                        coord_sf(
                                    xlim = c(-116, -105),
                                    ylim = c(21, 31),
                                    clip = "on"
                        ) +
                        guides(fill = guide_legend(title.position = "top", title.hjust = .5)) +
                        labs(x = "", y = "") +
                        theme_bw() +
                        theme(panel.background = element_rect(fill = "white"), 
                              axis.text.x = element_text(angle = 90),
                              legend.position = "right", 
                              legend.margin = margin(),
                              legend.title = element_text(face = "italic"))
)


(p1 <- ggplot(reefs_sf) +
                        geom_sf(data = ecoreg, aes(fill = factor(region)), col = "black") +
                        geom_sf(data = spdf_mx, fill = "gray90", col = NA) +
                        scale_fill_viridis_d() +
                        coord_sf(xlim = c(-116, -105), ylim = c(21, 31), expand = TRUE) +
                        annotate("rect", 
                                 xmin = -112, 
                                 xmax = -109, 
                                 ymin = 22, 
                                 ymax = 26.5, 
                                 alpha = 0.2, 
                                 fill = "red",
                                 col = "black") +
                        theme_bw() +
                        theme(panel.background = element_rect(fill = "white"), 
                              axis.text.x = element_text(angle = 90),
                              legend.position = "right", 
                              legend.margin = margin(),
                              legend.title = element_blank()))

colors <- c("Study reefs" = "black")

(p2 <- ggplot(reefs_sf) +
                        geom_sf(data = filter(ecoreg, region != "Northern ecoregion"), 
                                aes(fill = factor(region)), col = "black") +
                        geom_sf(data = spdf_mx, fill = "gray90", col = NA) +
                        geom_sf(aes(color = "Study reefs")) +
                        coord_sf(xlim = c(-112, -109), ylim = c(22, 26.5), expand = TRUE) +
                        theme_bw() +
                        labs(color = "Study reefs") +
                        scale_color_manual(values = colors) +
                        scale_fill_manual(values = c("#21908C", "#FDE725")) +
                        guides(fill = FALSE) +
                        theme(panel.background = element_rect(fill = "white"), 
                              axis.text.x = element_text(angle = 90),
                              legend.position = "right", 
                              legend.margin = margin(),
                              legend.title = element_blank()))


overall_degree <- sst_data %>% 
            mutate(fDegree = factor(round(lat))) %>% 
            ggplot(aes(x = t, 
                       y = temp,
                       col = fDegree)) +
            geom_line() +
            scale_colour_manual(values = c("#b2182b", "#f4a582", "#053061", "#4393c3")) +
            facet_grid(~fDegree) +
            labs(x= "Time", y=expression(Temperature ~ (degree*C))) +
            scale_y_continuous(breaks = seq(18, 32, by = 1),
                               minor_breaks = NULL) +
            theme_bw() +
            theme(strip.background = element_blank(),
                  strip.text = element_text(face="bold"),
                  panel.grid = element_blank(), 
                  legend.position = "")


# Assembling plots --------------------------------------------------------

study_area <- p1 + p2 + 
            plot_layout(guides = "collect") & theme(legend.position = 'bottom')

heatwaves <- map_slope / map_p + 
            plot_layout(guides = "collect") & theme(legend.position = 'right')

maps <- study_area / heatwaves +
            plot_annotation(tag_levels = "A") 

((study_area/overall_degree)|heatwaves) +
            plot_annotation(tag_levels = "A") 

# Saving the final figure
ggsave('figs/Figure_1_map_of_study_area_sst_MHWs.png', device = "png", type = "cairo", dpi = 600, width = 15, height = 8)




