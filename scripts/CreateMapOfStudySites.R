#Map of study locations

# Load packages -----------------------------------------------------------

library(sf)
library(ozmaps)
library(tidyverse)
library(cowplot)

# Load in site location data and map data ---------------------------------

site_locations <- read_sf(dsn = "rawdata/geographic", layer = "audio_sensors")

site_locations$site <- gsub("([a-z])\\..*", "\\1", site_locations$name)
#site_locations <- site_locations[,c(2,3,8)] %>% group_by(site) %>% summarize(across(.fns = mean))
#site_locations <- st_as_sf(site_locations, coords = c("lon", "lat"), crs = "EPSG:4326")

site_locations <- site_locations %>% group_by(site) %>% slice_head()

oz_states <- ozmaps::ozmap_states

oz_capitals <- st_as_sf(tibble::tribble( 
  ~city,           ~lat,     ~lon,
  "Sydney",    -33.8688, 151.2093,  
  "Melbourne", -37.8136, 144.9631, 
  "Brisbane",  -27.4698, 153.0251, 
  "Adelaide",  -34.9285, 138.6007, 
  "Perth",     -31.9505, 115.8605, 
  "Hobart",    -42.8821, 147.3272, 
  "Canberra",  -35.2809, 149.1300, 
  "Darwin",    -12.4634, 130.8456 
), coords = c("lon", "lat"), crs = "EPSG:4283")

# Plot map ----------------------------------------------------------------


Plot_Zoomed <- ggplot() +
  geom_sf(data = oz_states[c(1,3),]) +
  geom_point(data = site_locations, aes(lon, lat), colour = "red", size = 2) +
  geom_sf(data = oz_capitals[c(1,3),], colour = "black") +
  geom_sf_text(data = oz_capitals[c(1,3),], aes(label = city), nudge_y = 0.8, nudge_x = -0.8) +
  coord_sf() +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()

Plot_SiteNames <- ggplot() +
  geom_sf(data = oz_states[c(1,3),]) +
  geom_point(data = site_locations, aes(lon, lat), colour = "red", size = 2) +
  geom_text(data = site_locations, aes(x = lon, y = lat, label = str_to_title(site)), nudge_y = 0.5, nudge_x = -0.5, hjust = 1) +
  geom_sf(data = oz_capitals[c(1,3),], colour = "black") +
  #geom_sf_text(data = oz_capitals[c(1,3),], aes(label = city), nudge_y = 0.8, nudge_x = -0.8) +
  coord_sf() +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()

Zoomed_box <- list(xmin = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$x_range[1],
                   xmax = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$x_range[2],
                   ymin = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$y_range[1],
                   ymax = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$y_range[2])

Plot_Aus <- ggplot() +
  geom_sf(data = oz_states) +
  geom_rect(aes(xmin = Zoomed_box$xmin, xmax = Zoomed_box$xmax, 
                ymin = Zoomed_box$ymin, ymax = Zoomed_box$ymax),
            fill = "transparent", color = "red", size = 1.5) +
  coord_sf() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

Map_Inset <- ggdraw() +
  draw_plot(Plot_Zoomed, x = -0.2, y = 0, height = 1, width = 1) +
  draw_plot(Plot_Aus, x = 0.35, y = 0.65, width = 0.3, height = 0.3)

ggsave(filename = "outputs/figures/StudySiteLocations.png",
       Map_Inset,
       dpi = 800, width = 18, height = 18, units = "cm")

Map_Inset_Named <- ggdraw() +
  draw_plot(Plot_SiteNames, x = -0.2, y = 0, height = 1, width = 1) +
  draw_plot(Plot_Aus, x = 0.35, y = 0.65, width = 0.3, height = 0.3)

ggsave(filename = "outputs/figures/StudySiteLocations_Named.png",
       Map_Inset_Named,
       dpi = 800, width = 18, height = 18, units = "cm")


Site_schematic <- png::readPNG("Schematic_trapping_vert.png")
Site_schematic <- grid::rasterGrob(Site_schematic)
Site_schematic_plot <- qplot(geom="blank") +
  annotation_custom(Site_schematic, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point() + 
  theme_void()

CombinedMap_Schematic <- plot_grid(Map_Inset, Site_schematic_plot,
          rel_widths = c(0.7, 1),
          labels = c("a)", "b)"))

ggsave(filename = "outputs/figures/CombinedMap_Schematic.png",
       CombinedMap_Schematic,
       dpi = 800, width = 44, height = 18, units = "cm")

#New 1000dpi figure for publication (Ecological Indicators) ----

Plot_Zoomed <- ggplot() +
  geom_sf(data = oz_states[c(1,3),]) +
  geom_point(data = site_locations, aes(lon, lat), colour = "red", size = 1.5) +
  geom_sf(data = oz_capitals[c(1,3),], colour = "black", size = 1.5) +
  geom_sf_text(data = oz_capitals[c(1,3),], aes(label = city), nudge_y = 0.9, nudge_x = -0.9, size = 2) +
  coord_sf() +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()

Zoomed_box <- list(xmin = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$x_range[1],
                   xmax = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$x_range[2],
                   ymin = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$y_range[1],
                   ymax = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$y_range[2])

Plot_Aus <- ggplot() +
  geom_sf(data = oz_states, size = 0.2) +
  geom_rect(aes(xmin = Zoomed_box$xmin, xmax = Zoomed_box$xmax, 
                ymin = Zoomed_box$ymin, ymax = Zoomed_box$ymax),
            fill = "transparent", color = "red", size = 0.8) +
  coord_sf() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA))

Map_Inset <- ggdraw() +
  draw_plot(Plot_Zoomed, x = -0.2, y = 0, height = 1, width = 1) +
  draw_plot(Plot_Aus, x = 0.38, y = 0.65, width = 0.3, height = 0.3)

CombinedMap_Schematic <- plot_grid(Map_Inset, Site_schematic_plot,
                                   rel_widths = c(0.7, 1),
                                   labels = c("a)", "b)"), label_x = -0.015)

ggsave(filename = "outputs/figures/publication/Figure01.png",
       CombinedMap_Schematic,
       dpi = 1000, width = 190, height = 80, units = "mm")







# ----
Site_schematic <- ggdraw + draw_image("Schematic_trapping_vert.png")





Plot_Zoomed <- ggplot() +
  geom_sf(data = oz_states[c(1,3),]) +
  geom_point(data = site_locations, aes(lon, lat), colour = "red", size = 2) +
  geom_sf(data = oz_capitals[c(1,3),], colour = "black") +
  geom_sf_text(data = oz_capitals[c(1,3),], aes(label = city), nudge_y = 0.8, nudge_x = -0.8) +
  coord_sf() +
  theme_bw()

Zoomed_box <- list(xmin = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$x_range[1],
                   xmax = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$x_range[2],
                   ymin = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$y_range[1],
                   ymax = ggplot_build(Plot_Zoomed)$layout$panel_params[[1]]$y_range[2])

Plot_Aus <- ggplot() +
  geom_sf(data = oz_states) +
  geom_rect(aes(xmin = Zoomed_box$xmin, xmax = Zoomed_box$xmax, 
                ymin = Zoomed_box$ymin, ymax = Zoomed_box$ymax),
            fill = "transparent", color = "red", size = 1.5) +
  coord_sf() +
  theme_bw()