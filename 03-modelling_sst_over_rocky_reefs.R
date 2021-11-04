## Loading libraries

library(tidyverse)
library(lubridate)
library(fst)
library(mgcv)
library(gratia)
library(ggthemes)
library(tidymv)
library(patchwork)

## Loading data 
sst_data <- read_fst("data/study_area_sst_extract.fst") 


## Data wrangling 

sst_model <- sst_data %>% 
            mutate(Year = year(t), nMonth = month(t), fDegree = factor(round(lat, 0))) %>% 
            group_by(Year, nMonth, fDegree) %>% 
            summarise(Temperature = mean(temp))


sst_model %>% 
            ggplot(aes(x=as.Date(paste0(Year, "-", nMonth, "-01"), "%Y-%m-%d"), 
                       y = Temperature,
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

knots <- list(nMonth = c(0.5, seq(1,12, length = 10), 12.5))





sst_model_list <- split(sst_model, sst_model$fDegree)


m_23 <- gamm(Temperature ~ te(Year, nMonth, bs = c("cr", "cc"), k = c(10, 12)), 
             data = sst_model_list[[1]], 
             method = "REML", 
             correlation = corARMA(form = ~ 1 | Year, p = 1),
             knots = knots)


m_24 <- gamm(Temperature ~ te(Year, nMonth, bs = c("cr", "cc"), k = c(10, 12)), 
             data = sst_model_list[[2]], 
             method = "REML", 
             correlation = corARMA(form = ~ 1 | Year, p = 1),
             knots = knots)

m_25 <- gamm(Temperature ~ te(Year, nMonth, bs = c("cr", "cc"), k = c(10, 12)), 
             data = sst_model_list[[3]], 
             method = "REML", 
             correlation = corARMA(form = ~ 1 | Year, p = 1),
             knots = knots)


m_26 <- gamm(Temperature ~ te(Year, nMonth, bs = c("cr", "cc"), k = c(10, 12)), 
             data = sst_model_list[[4]], 
             method = "REML", 
             correlation = corARMA(form = ~ 1 | Year, p = 1),
             knots = knots)

pdat <- with(sst_model_list[[1]],
             data.frame(Year = rep(c(1985, 2021), each = 100),
                        nMonth = rep(seq(1, 12, length = 100), times = 2)))


plot_list <- list()
model_list <- list(m_23, m_24, m_25, m_26)

for(i in 1:length(sst_model_list)) {
            
            pred <- predict(model_list[[i]]$gam, newdata = pdat, se.fit = TRUE)
            crit <- qt(0.975, df = df.residual(model_list[[i]]$gam)) # ~95% interval critical t
            pdat <- transform(pdat, fitted = pred$fit, se = pred$se.fit, fYear = as.factor(Year))
            pdat <- transform(pdat,
                              upper = fitted + (crit * se),
                              lower = fitted - (crit * se))
            
            plot_list[[i]] <- ggplot(pdat, aes(x = nMonth, y = fitted, group = fYear)) +
                        geom_ribbon(mapping = aes(ymin = lower, ymax = upper,
                                                  fill = fYear), alpha = 0.2) + # confidence band
                        geom_line(aes(colour = fYear)) +    # predicted temperatures
                        theme_bw() +                        # minimal theme
                        theme(legend.position = "bottom") +    # push legend to the top
                        labs(y = expression(Temperature ~ (degree*C)), x = "Month", 
                             title = sst_model_list[[i]]$fDegree) +
                        # scale_color_viridis_d() +
                        # scale_fill_viridis_d() +
                        scale_fill_manual(name = "Year", values = c("#053061", "#b2182b")) + # correct legend name
                        scale_colour_manual(name = "Year", values = c("#053061", "#b2182b")) +
                        scale_x_continuous(breaks = 1:12,   # tweak where the x-axis ticks are
                                           labels = month.abb, # & with what labels
                                           minor_breaks = NULL) +
                        scale_y_continuous(breaks = seq(18, 32, by = 4), limits = c(18, 32),
                                           minor_breaks = NULL) 
}



Month_plot <- patchwork::wrap_plots(plot_list[[4]], plot_list[[3]], plot_list[[2]], plot_list[[1]], 
                            nrow = 4, ncol = 1, guides = "collect") & theme(legend.position = "bottom")

# 
# plot_list2 <- list()
# predicted_list <- list()
# for(i in 1:length(sst_model_list)) {
#             pdat2 <- with(model_list[[i]]$gam,
#                           data.frame(
#                                       Year = rep(1983:2021, each = 12),
#                                       nMonth = rep(1:12, times = 39)
#                           ))
#             pred2 <-
#                         predict(model_list[[i]]$gam,
#                                 newdata = pdat2,
#                                 se.fit = TRUE)
#             ## add predictions & SEs to the new data ready for plotting
#             pdat2 <- transform(
#                         pdat2,
#                         fitted = pred2$fit,
#                         # predicted values
#                         se = pred2$se.fit,
#                         # standard errors
#                         fYear = as.factor(Year),
#                         fMonth = factor(month.abb[nMonth], # month as a factor
#                                         levels = month.abb)
#             )
#             pdat2 <- transform(pdat2,
#                                upper = fitted + (crit * se),
#                                # upper and...
#                                lower = fitted - (crit * se)) # lower confidence bounds
#             
#             pdat2$fDegree <- names(sst_model_list)[i]
#             predicted_list[[i]] <- pdat2
#             
#             
#             plot_list2[[i]] <-
#                         ggplot(pdat2, aes(
#                                     x = fMonth,
#                                     y = fitted,
#                                     group = fYear
#                         )) +
#                         geom_line(aes(colour = fYear)) +   # draw trend lines
#                         theme_bw() +                        # minimal theme
#                         scale_color_viridis_d() +
#                         scale_fill_viridis_d() +
#                         theme(legend.position = "none") +   # no legend
#                         labs(
#                                     y = expression(SST ~ (degree * C)),
#                                     x = NULL,
#                                     title = sst_model_list[[i]]$fDegree
#                         ) +
#                         #facet_wrap(~ fMonth, ncol = 6) +    # facet on month
#                         scale_y_continuous(breaks = seq(18, 32, by = 2),
#                                            minor_breaks = NULL)
# }
# 
# 
# patchwork::wrap_plots(plot_list2[[4]], plot_list2[[3]], plot_list2[[2]], plot_list2[[1]], 
#                       nrow = 4, ncol = 1, guides = "collect") +
#             plot_annotation(tag_levels = "A")
# 
# 
# ggsave("figs/Monthly_predicted_change.png", dpi = 800, width = 5, height = 8)

 
mod_data_lat <- sst_data %>% 
            group_by(t, lat) %>% 
            summarise(Temperature = mean(temp)) %>% 
            mutate(Year = year(t), nMonth = month(t))

fullM_sst <- bam(Temperature ~ t2(Year, lat, bs = c("tp", "tp")), 
                 data = mod_data_lat, 
                 method = "REML")

appraise(fullM_sst)
gam.check(fullM_sst)
year_plot <- draw(fullM_sst) & labs(x="Year", y="Latitude", title = "SST variation in the study area") &
            theme_bw()

cowplot::plot_grid(year_plot, Month_plot, labels = "AUTO")

ggsave("figs/full_sst_image.png", height = 6, width = 8)

