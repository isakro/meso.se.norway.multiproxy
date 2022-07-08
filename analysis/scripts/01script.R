library(tidyverse)
library(sf)
library(terra)
library(here)
library(spatstat)

surveyed <- st_read(here("analysis/data/raw_data/surveyed_sites.gpkg"))
excavated <- st_zm(read_sf(here("analysis/data/raw_data/excavated_sites.gpkg")))
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")

load(here("analysis/data/raw_data/displacement_data.RData"))

# Find elevation of surveyed sites
surveyed$elev <- terra::extract(dtm, vect(surveyed), fun = mean)[, -1]

# Surveyed sites with a quality score better than 4
survq <- filter(surveyed, quality < 4)

bmap <- st_read(here("analysis/data/raw_data/naturalearth_countries.gpkg"))

bboxsites <- st_bbox(surveyed)
bboxsites[1] <- bboxsites[1] - 15000
bboxsites[3] <- bboxsites[3] + 15000
bboxsites[2] <- bboxsites[2] - 5000
bboxsites[4] <- bboxsites[4] + 5000
bboxsitespoly <- st_as_sf(st_as_sfc(bboxsites))

bound_reproj <- st_transform(bboxsitespoly, st_crs(bmap))
bmap2 <- bmap %>%
  dplyr::filter(st_intersects(., bound_reproj, sparse = FALSE))
bmap_reproj <- st_transform(bmap2, st_crs(surveyed))

anc <- as.numeric(c(bboxsites$ymin, bboxsites$xmax))
ggplot() +
  geom_sf(data = bmap_reproj, fill = "grey", colour = NA) +
  geom_sf(data = st_centroid(surveyed), fill = "#ff8c00",
          size = 2.25, shape = 21,
          colour = "black") +
  geom_sf(data = st_centroid(survq), fill = "#00BA38",
          size = 2.25, shape = 21,
          colour = "black") +
  geom_sf(data = st_centroid(excavated), fill = "white",
          size = 2.25, shape = 21,
          colour = "black") +
  ggsn::scalebar(data = surveyed, dist = 20, dist_unit = "km",
                 transform = FALSE, st.size = 4, height = 0.02,
                 border.size = 0.1, st.dist = 0.03,
                 anchor = c(x = anc[2] - 15500, y = anc[1]) + 8000) +
  coord_sf(xlim = c(bboxsites[1], bboxsites[3]),
           ylim = c(bboxsites[2], bboxsites[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "none")

ggplot() + geom_histogram(data = surveyed, aes(x = elev), bins = 50)
ggplot() + geom_histogram(data = survq, aes(x = elev), bins = 50)

#Spatstat

ssites  <- as.ppp(st_centroid(survq))
marks(ssites) <- NULL
ssites <- rescale(ssites, 1000)
Window(ssites) <- as.owin(ssites)

plot(ssites)

K1 <- density(ssites, sigma = 7) # Using the default bandwidth
plot(K1, main=NULL, las=1)
contour(K1, add=TRUE)
points(as.ppp(sample(K1)), add = TRUE)

class(K1)

shorelinedates <- list()
for(i in 1:nrow(survq)){
  print(paste(i, survq$askeladden_id[i]))
  shorelinedates[[i]] <- shoreline_date(survq[i,],
                                      dtm,
                                      displacement_curves,
                                      isobases,
                                      expratio = 0.168,
                                      reso = 0.001)
}

sdates <- bind_rows(shorelinedates) %>% group_by(site_name)

# save(sdates,
#       file = here::here("analysis/data/derived_data/sdates.RData"))

ggplot(sdates) + geom_line(aes(x = years, y = probability, col = site_name)) +
  theme(legend.position="none")

sdates %>% group_by(years) %>%
  filter(!is.na(probability)) %>% # Excludes Askeladden ID 108682-1 at 209masl
  summarise(prob_sum = sum(probability)) %>%
  ggplot() + geom_line(aes(x = years, y = prob_sum)) +
  scale_x_continuous(breaks = seq(-10000, 2500, 500))

