library(sf)
library(ggplot2)
library(terra)
library(patchwork)
library(rnaturalearth)
library(geodata)

setwd("PATH/05_aDNA/02_results/climatemaps")

# Load ice-sheet reconstructions
ice24 <- st_read("./DATED1_TimeSlices_shp/TS24_mc.shp", quiet = TRUE)
ice16 <- st_read("./DATED1_TimeSlices_shp/TS16_mc.shp", quiet = TRUE)

# http://www.paleoclim.org
# tp0: Pleistocene: late-Holocene, Meghalayan (4.2-0.3 ka), v1.0: http://sdmtoolbox.org/paleoclim.org/data/LH/LH_v1_2_5m.zip
# tp16: Pleistocene: Heinrich Stadial 1 (17.0-14.7 ka), v1.0: http://sdmtoolbox.org/paleoclim.org/data/HS1/HS1_v1_2_5m.zip
# tp24: Pleistocene: Last Glacial Maximum (ca. 21 ka), v1.2b**, NCAR CCSM4: http://sdmtoolbox.org/paleoclim.org/data/chelsa_LGM/chelsa_LGM_v1_2B_r2_5m.zip
# Bio1 = annual temperature x10
tp24_temp <- rast("./chelsa_LGM_v1_2B_r2_5m/bio_1.tif") /10
tp16_temp <- rast("./HS1_v1_2_5m/bio_1.tif") /10
tp0_temp <- rast("./LH_v1_2_5m/bio_1.tif") /10
tp24_temp <- project(tp24_temp, "EPSG:4326")
tp16_temp <- project(tp16_temp, "EPSG:4326")
tp0_temp  <- project(tp0_temp,  "EPSG:4326")

# Ensure all layers use the same projection
europe <- ne_countries(continent = "Europe",returnclass = "sf")
ice24 <- st_transform(ice24, 4326)
ice16 <- st_transform(ice16, 4326)
europe <- st_transform(europe, 4326)
#europe <- st_transform(europe, st_crs(ice24))
#ice16  <- st_transform(ice16, st_crs(ice24))

region_ext <- ext(-15, 90, 30, 72)
tp24_temp <- crop(tp24_temp, region_ext)
tp16_temp <- crop(tp16_temp, region_ext)
tp0_temp  <- crop(tp0_temp, region_ext)

df_tp24 <- as.data.frame(tp24_temp, xy = TRUE)
df_tp16 <- as.data.frame(tp16_temp, xy = TRUE)
df_tp0  <- as.data.frame(tp0_temp,  xy = TRUE)
names(df_tp24)[3] <- "bio1"
names(df_tp16)[3] <- "bio1"
names(df_tp0)[3]  <- "bio1"

dem <- geodata::elevation_global(res = 10,path = "dem/")
dem_eur <- crop(dem, region_ext)
dem_match24 <- resample(dem_eur,tp24_temp,method = "bilinear")
alps_mask24 <- dem_match24 > 500 & tp24_temp < -2
alps_poly24 <- as.polygons(alps_mask24, dissolve = TRUE)
alps_poly24 <- sf::st_as_sf(alps_poly24)
alps_poly24 <- alps_poly24[alps_poly24$wc2.1_10m_elev == 1, ]

dem_match24 <- resample(dem_eur,tp24_temp,method = "bilinear")
alps_mask24 <- dem_match24 > 500 & tp24_temp < -2
alps_poly24 <- as.polygons(alps_mask24, dissolve = TRUE)
alps_poly24 <- sf::st_as_sf(alps_poly24)
alps_poly24 <- alps_poly24[alps_poly24$wc2.1_10m_elev == 1, ]

dem_match16 <- resample(dem_eur,tp16_temp,method = "bilinear")
alps_mask16 <- dem_match16 > 500 & tp16_temp < -2
alps_poly16 <- as.polygons(alps_mask16, dissolve = TRUE)
alps_poly16 <- sf::st_as_sf(alps_poly16)
alps_poly16 <- alps_poly16[alps_poly16$wc2.1_10m_elev == 1, ]

dem_match2 <- resample(dem_eur,tp0_temp,method = "bilinear")
alps_mask2 <- dem_match2 > 500 & tp0_temp < -2
alps_poly2 <- as.polygons(alps_mask2, dissolve = TRUE)
alps_poly2 <- sf::st_as_sf(alps_poly2)
alps_poly2 <- alps_poly2[alps_poly2$wc2.1_10m_elev == 1, ]

all_vals <- c(values(tp24_temp),values(tp16_temp),values(tp0_temp))
lims <- range(all_vals, na.rm = TRUE)

bul_df <- data.frame(lon = 24.895556,lat = 43.236389)

p_lgm <- ggplot() + geom_tile(data = df_tp24, aes(x = x, y = y, fill = bio1)) +
  scale_fill_viridis_c(option = "C",limits = lims,name = "Mean annual\n temp (°C)") +
  geom_sf(data = europe,fill = NA,color = "grey40",linewidth = 0.3) +
  geom_sf(data = ice24,fill = "white",color = "steelblue",alpha = 0.6) +
  geom_sf(data = alps_poly24,fill = "white",color = "steelblue",alpha = 0.6) +
  geom_point(data = bul_df,aes(lon, lat),shape = 21,size = 3,fill = "red",color = "black") +
  coord_sf(crs = st_crs(4326),xlim = c(-15, 50), ylim = c(30, 70),expand = FALSE,clip = "on") +
  theme_minimal() + theme_bw(base_size=9) +
  theme(axis.text = element_text(),
    axis.title = element_text(), 
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
    legend.position = "none", plot.margin = margin(1, 1, 1, 1, "pt")) +
  labs(title = "", fill = "BIO1")

p_hs1 <- ggplot() + geom_tile(data = df_tp16, aes(x = x, y = y, fill = bio1)) +
  scale_fill_viridis_c(option = "C",limits = lims,name = "Mean annual\n temp (°C)") +
  geom_sf(data = europe,fill = NA,color = "grey40",linewidth = 0.3) +
  geom_sf(data = ice16,fill = "white",color = "steelblue",alpha = 0.6) +
  geom_point(data = bul_df,aes(lon, lat),shape = 21,size = 3,fill = "red",color = "black") +
  geom_sf(data = alps_poly,fill = "grey20",alpha = 0.15,color = NA) +
  geom_sf(data = alps_poly16,fill = "white",color = "steelblue",alpha = 0.6) +
  coord_sf(crs = st_crs(4326),xlim = c(-15, 50), ylim = c(30, 70),expand = FALSE,clip = "on") +
  theme_minimal() + theme_bw(base_size=9) +
  theme(axis.text = element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text()) +
  labs(title = "", fill = "BIO1")

p_holo <- ggplot() + geom_tile(data = df_tp0, aes(x = x, y = y, fill = bio1)) +
  scale_fill_viridis_c(option = "C",limits = lims,name = "Mean annual\n temp (°C)") +
  geom_sf(data = europe,fill = NA,color = "grey40",linewidth = 0.3) +
  geom_sf(data = alps_poly2,fill = "white",color = "steelblue",alpha = 0.6) +
  geom_point(data = bul_df,aes(lon, lat),shape = 21,size = 3,fill = "red",color = "black") +
  geom_sf(data = alps_poly,fill = "grey20",alpha = 0.15,color = NA) +
  coord_sf(crs = st_crs(4326),xlim = c(-15, 50), ylim = c(30, 70),expand = FALSE,clip = "on") +
  theme_minimal() + theme_bw(base_size=9) +
  theme(axis.text = element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(),
        legend.position = "top", plot.margin = margin(1, 1, 1, 1, "pt")) +
  labs(title = "", fill = "BIO1")

final_plot <- p_lgm + p_holo + plot_layout(guides = "collect") &
  theme(legend.position = "top")
final_plot #6x4 portrait 600x800 png/ tiff

ggsave("climateplot_updated2026JUN.tiff",
  final_plot,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)
