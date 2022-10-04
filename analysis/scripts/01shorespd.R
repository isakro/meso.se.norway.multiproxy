library(tidyverse)
library(sf)
library(terra)
library(here)
library(shoredate)
library(patchwork)

surveyed <- st_read(here("analysis/data/raw_data/surveyed_sites.gpkg"))
excavated <- st_zm(read_sf(here("analysis/data/raw_data/excavated_sites.gpkg")))
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")
site_dat <- read.csv(here("analysis/data/raw_data/sites.csv"))
munc <- read_sf(here("analysis/data/raw_data/municipalities.gpkg"))

# Exclude sites that haven't been excavated and which haven't been shoreline
# dated
shore_sites <- site_dat %>% filter(investigated == "t") %>%
  filter(dating_method == "shore_typo" & !is.na(reported_earliest_bce)
         & !is.na(reported_min_elev))

# Mørland D11, Solum 2, Kvastad A4 west, Ragnhildrød, Marum mellom 4
# Join shoreline dated excavated sites with spatial data
excshore <- st_as_sf(left_join(shore_sites, excavated,
                                  by = c("name" = "site_name", "ask_id"))) %>%
  filter(!(is.na(st_dimension(geom))))

step_reso <- 0.1

# Find mean reported elevation for the sites
excshore$elev <- rowMeans(subset(st_drop_geometry(excshore),
                                 select = c(reported_min_elev,
                                            reported_max_elev)))

# # Print sites that are excluded due to lacking site limits for inspection
# excshore %>% filter(st_is_empty(excshore)) %>%  pull(name)

load(here("analysis/data/raw_data/displacement_data.RData"))
load(here("analysis/data/derived_data/incpolys.RData"))

# Find elevation of surveyed sites
surveyed$elev <- terra::extract(dtm, vect(surveyed), fun = mean)[, -1]

# Surveyed sites with a quality score better than 4
survq <- filter(surveyed, quality < 4)


surveyed <- surveyed %>%  mutate(lab = ifelse(quality < 4, "t", "f"))
excavated$lab <- "y"

# Select and rename columns necessary for shoreline dating
rcolexcshore <- excshore %>% select(ask_id, elev)
rcolsurvq <- survq  %>% select(askeladden_id, elev) %>%
  rename(ask_id = askeladden_id)

# Combine the data
ssites <-rbind(rcolexcshore, rcolsurvq)

bmap <- st_read(here("analysis/data/raw_data/naturalearth_countries.gpkg"))

# Create a bounding box around site points to limit the overview map.
sitbbox <- st_bbox(surveyed)
sitbbox[1:2] <- sitbbox[1:2] - 1000000
sitbbox[3:4] <- sitbbox[3:4] + 1000000
boundingpoly <- st_as_sf(st_as_sfc(sitbbox))

# Reproject the bounding box to match world map, and crop the world map
# with the bounding box polygon.
bound_reproj <- st_transform(boundingpoly, st_crs(bmap))
mapcountries <- bmap %>%
  filter(st_intersects(., bound_reproj, sparse = FALSE))
count_reproj <- st_transform(mapcountries, st_crs(surveyed))

bboxsites <- st_bbox(surveyed)
bboxsites[1] <- bboxsites[1] - 15000
bboxsites[3] <- bboxsites[3] + 15000
bboxsites[2] <- bboxsites[2] - 5000
bboxsites[4] <- bboxsites[4] + 5000
bboxsitespoly <- st_as_sf(st_as_sfc(bboxsites))

overview <-
  ggplot() +
  geom_sf(data = count_reproj, fill = "grey", colour = NA) +
  geom_sf(data = bboxsitespoly,
          fill = NA, colour = "black", size = 0.5) +
  coord_sf(xlim = c(sitbbox[1], sitbbox[3]), ylim = c(sitbbox[2], sitbbox[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank())

bound_reproj <- st_transform(bboxsitespoly, st_crs(bmap))
bmap2 <- bmap %>%
  filter(st_intersects(., bound_reproj, sparse = FALSE))
bmap_reproj <- st_transform(bmap2, st_crs(surveyed))


anc <- as.numeric(c(bboxsites$ymin, bboxsites$xmax))
sa <- ggplot() +
  geom_sf(data = bmap_reproj, fill = "grey", colour = NA) +
  geom_sf(data = st_transform(munc, st_crs(surveyed)), fill = NA,
          colour = "black", lwd = 0.25) +
  geom_sf(data = st_centroid(surveyed), aes(fill = lab),
          size = 2, shape = 21,
          colour = "black", show.legend = "point") +
  geom_sf(data = st_centroid(filter(surveyed, quality < 4)),
          aes(fill = lab),
          size = 2, shape = 21,
          colour = "black", show.legend = "point") +
  geom_sf(data = st_centroid(excavated), aes(fill = lab),
          size = 2, shape = 21,
          colour = "black", show.legend = "point") +
  scale_fill_manual(labels = c(paste0("Surveyed sites,\nincluded (n = ",
                                      nrow(filter(surveyed, quality < 4)) ,")"),
                               paste0("Surveyed sites,\nexcluded (n = ",
                                      nrow(filter(surveyed, quality > 3)) ,")"),
                               paste0("Excavated sites (n = ", nrow(excavated) ,
                                      ")")),
                    values = c("t" = "#00BA38", "f" = "black", "y" = "white"),
                    name = "") +
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
                     legend.position = "bottom",
                     legend.text=element_text(size = 11))

cowplot::ggdraw() +
  cowplot::draw_plot(sa) +
  cowplot::draw_plot(overview, x = 0.126,
                     y = 0.65, width = 0.35, height = 0.35)


# overview + sa +   plot_layout(widths = c(1, 2.5))
ggsave(here::here("analysis/figures/map.png"),
       units = "px", width = 2194*1.7, height = 2380)


ggplot() + geom_histogram(data = surveyed, aes(x = elev), bins = 50)
ggplot() + geom_histogram(data = survq, aes(x = elev), bins = 50)

step_reso <- 0.1

shorelinedates <- list()
for(i in 1:nrow(ssites)){
  print(paste(i, ssites$ask_id[i]))
  shorelinedates[[i]] <- shoreline_date(ssites[i,],
                                      dtm,
                                      displacement_curves,
                                      isobases,
                                      expratio = 0.168,
                                      reso = step_reso)


}

sdates <- bind_rows(shorelinedates)

# save(sdates, sdates_001,
#       file = here("analysis/data/derived_data/sdates_0.RData"))

# save(sdates,
#       file = here("analysis/data/derived_data/sdates.RData"))

load("analysis/data/derived_data/sdates.RData")

spd <- sdates %>% filter(years < -2500)
spdn <- spd %>% filter(probability != 0)
spdn <- length(unique(spdn$site_name)) # Couldn't figure out n_distinct() here

spd <- spd %>% group_by(years) %>%
  filter(!is.na(probability)) %>% # Excludes Askeladden ID 108682-1 at 209masl
  summarise(prob_sum = sum(probability))

exnoprob <- spd %>% filter(prob_sum > 0)
lower_lim <- min(exnoprob$years)
upper_lim <- max(exnoprob$years)

exp <- TRUE

if(exp){
  # Code repurposed from modelTest() from rcarbon
  # Exponential
  print("Exponential")
  fit <- nls(y ~ (exp(a + b * x)),
             data = data.frame(x = exnoprob$years, y = exnoprob$prob_sum),
             start = list(a = 0, b = 0))
  est <- predict(fit, list(x = exnoprob$years))
  predgrid <- data.frame(bce = exnoprob$years, prob_dens = est)
  predgrid <- filter(predgrid, bce >= lower_lim & bce <= upper_lim)
} else{
  print("Uniform")
  # # Uniform
  predgrid <- data.frame(bce = exnoprob$years,
                         prob_dens = mean(exnoprob$prob_sum))
  predgrid <- filter(predgrid, bce > lower_lim & bce < upper_lim)
}

ggplot(spd) + geom_line(aes(x = years, y = prob_sum)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens)) +
  scale_x_continuous(breaks = seq(-10000, 2500, 500))

incpolys$dens <- lengths(st_intersects(incpolys, ssites)) /
  sum(lengths(st_intersects(incpolys, ssites)))

bboxpolys <- st_bbox(surveyed)
bboxpolys[1] <- bboxpolys[1] - 900
bboxpolys[3] <- bboxpolys[3] + 900
bboxpolys[2] <- bboxpolys[2] - 900
bboxpolys[4] <- bboxpolys[4] + 900
bboxpolyspoly <- st_as_sf(st_as_sfc(bboxpolys))

ggplot() +
  geom_sf(data = incpolys, aes(fill = dens), colour = NA, alpha = 0.7) +
  geom_sf(data = bmap_reproj, fill = NA, colour = "black") +
  geom_sf(data = st_centroid(ssites), size = 0.5, colour = "black") +
  # ggsn::scalebar(data = surveyed, dist = 20, dist_unit = "km",
  #                transform = FALSE, st.size = 4, height = 0.02,
  #                border.size = 0.1, st.dist = 0.03,
  #                anchor = c(x = anc[2] - 15500, y = anc[1]) + 8000) +
  scale_fill_gradient(low = "white", high = "black", name = "Site density") +
  coord_sf(xlim = c(bboxpolys[1], bboxpolys[3]),
           ylim = c(bboxpolys[2], bboxpolys[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "right")

ggsave(here::here("analysis/figures/incpolys.png"),
       units = "px", width = 2650, height = 1600)

ssize = spdn
nsim <- 100
random_dates <- data.frame(matrix(nrow = ssize*nsim, ncol = 3))
names(random_dates) <- c("sample", "simn", "displacement_curve")


random_dates$sample = sample(predgrid$bce, replace = TRUE,
                            size = ssize*nsim, prob = predgrid$prob_dens)
random_dates$simn <- rep(1:nsim, each = nrow(random_dates)/nsim)
random_dates$displacement_curve <- incpolys[sample(incpolys$id, replace = TRUE, size = ssize*nsim,
                              prob = incpolys$dens),]$disp[[1]]

# do.call(reverse_shoredate(sample, incpolys[random_dates$incpoly,]$disp[[1]]), random_dates)

random_dates$reverse_elevation <- mapply(reverse_shoredate,
                                         random_dates$sample,
                                         random_dates$displacement_curve)

shorelinedates <- list()
# start_time <- Sys.time()
# shorelinedates <- mapply(shoredate::shoreline_date,
#                          site = random_dates$simn,
#                          reso = step_reso,
#                          elevation = random_dates$reverse_elevation,
#                          interpolated_curve = random_dates$displacement_curve)
# end_time <- Sys.time()
# end_time - start_time
#
# simdates <- lapply(shorelinedates, '[[', 'date') %>% bind_rows() %>%
#   group_by(bce, site_name) %>%
#   filter(!is.na(probability)) %>%
#   summarise(prob_sum = sum(probability))

# 299

# shorelinedates <- list()
# start_time <- Sys.time()
# for (i in 1:nrow(random_dates)){
#   print(paste(i, "/", nrow(random_dates)))
#
#   shorelinedates[[i]] <- shoreline_date(
#                  site = random_dates$simn[i],
#                  elev = NA,
#                  displacement_curves = displacement_curves,
#                  isobases = isobases,
#                  reso = step_reso,
#                  expratio = 0.168,
#                  siteelev = "mean",
#                  specified_elev = TRUE,
#                  elevation = random_dates$reverse_elevation[i],
#                  interpolate_curve = FALSE,
#                  interpolated_curve = random_dates$displacement_curve[i])
#
  # # Save data externally and overwrite shorelinedate for memory purposes
  # if(i %% ssize == 0){
  #   tmp <- bind_rows(shorelinedates) %>% group_by(bce, site_name) %>%
  #     filter(!is.na(probability)) %>%
  #     summarise(prob_sum = sum(probability))
  #
  #   save(tmp,
  #        file = paste0(here::here("../external_data/shorespd/"),
  #                      i, "shore.Rdata"))
  #
  #   shorelinedates <- list()
  # }

  # if(i %% divinc == 0){
  #   tmp <- shorelinedates[(i+1-divinc):i]
  #   save(tmp,
  #        file = paste0(here::here("../external_data/shorespd/"),
  #                      i, "shore.Rdata"))
  #   shorelinedates <- list()
  # }
  # if(i == nrow(random_dates)){
  #   tmp <- shorelinedates[(i+1-(i %% divinc)):i]
  #   save(tmp,
  #        file = paste0(here::here("../external_data/shorespd/"),
  #                      i, "shore.Rdata"))
  # }
# }
# end_time <- Sys.time()
# end_time - start_time
#
# tst <- list()
# for(i in 1:1000){
#   tst[[i]] <- i
#   if (i %% 100 == 0){
#     print(tst[(i+1-100):i])
#     tst <- list()
#   }
# }


shorelinedates <- list()
start_time <- Sys.time()
for (i in 1:nrow(random_dates)){
  print(paste(i, "/", nrow(random_dates)))

  shorelinedates[[i]] <- shoredate::shoreline_date(
    site = random_dates$simn[i],
    reso = step_reso,
    elevation = random_dates$reverse_elevation[i],
    interpolated_curve = random_dates$displacement_curve[i],
    sparse = TRUE)

  if(i %% ssize == 0){
  # Save data externally and overwrite shorelinedate for memory purposes
    if(exp){
      tmp <- bind_rows(shorelinedates) %>% group_by(bce, site_name) %>%
        filter(!is.na(probability)) %>%
        summarise(prob_sum = sum(probability))

      save(tmp,
           file = paste0(here::here("../external_data/shorespd/exp/"),
                         i, "shore.rds"))

      shorelinedates <- list()
    }  else{
      tmp <- bind_rows(shorelinedates) %>% group_by(bce, site_name) %>%
        filter(!is.na(probability)) %>%
        summarise(prob_sum = sum(probability))

      save(tmp,
           file = paste0(here::here("../external_data/shorespd/uni/"),
                         i, "shore.rds"))

      shorelinedates <- list()
    }
  }
  # # Save data externally and overwrite shorelinedate for memory purposes
  # if(i %% divinc == 0){
  #   tmp <- shorelinedates[(i+1-divinc):i]
  #   save(tmp,
  #        file = paste0(here::here("../external_data/shorespd/"),
  #                      i, "shore.rds"))
  #   shorelinedates <- list()
  # }
  # if(i == nrow(random_dates)){
  #   tmp <- shorelinedates[(i+1-(i %% divinc)):i]
  #   save(tmp,
  #        file = paste0(here::here("../external_data/shorespd/"),
  #                      i, "shore.rds"))
  # }

}
end_time <- Sys.time()
end_time - start_time

if(exp){
  resfiles <- list.files(here::here("../external_data/shorespd/exp/"),
                         full.names = TRUE)
} else{
    resfiles <- list.files(here::here("../external_data/shorespd/uni/"),
                          full.names = TRUE)
}


results <- list()
# seq_along(resfiles)
for(i in 1:length(resfiles)){
  load(resfiles[i])
  results[[i]] <- bind_rows(tmp)
  rm(tmp)
}

simdates <- do.call(rbind.data.frame, results)

simdates <- simdates %>% group_by(bce) %>%
  mutate(low = quantile(prob_sum, prob = 0.025),
         high = quantile(prob_sum, prob = 0.975),
         mean = mean(prob_sum)) %>%
  distinct(bce, .keep_all = TRUE)

sspdu2 <- ggplot() +
  ggplot2::geom_ribbon(data = simdates, aes(x = bce, ymin = low, ymax = high),
              fill = "grey", alpha = 0.9) +
  # geom_line(data = simdates, aes(x = bce, y = mean), col = "red", lwd = 0.5) +
  geom_line(data = spd, aes(x = years, y = prob_sum)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500)) +
  scale_y_continuous(limits = c(0, 0.7)) +
  labs(x = "BCE", y = "Summed probability") +
  theme_bw()

(sspdu1 + sspdu2)  /
(rspdu1 + rspdu2) + plot_annotation(tag_levels = 'A')
ggsave("analysis/figures/spds_uniform2.png",
       units = "px", width = 2194*1.5, height = 2380/1.25)


vctrs::vec_as_names(names(tst), repair = "universal")

simdates <- do.call(rbind.data.frame,
        do.call(rbind, unlist(results, recursive = FALSE))) %>%
  group_by(bce, site_name) %>%
  filter(!is.na(probability)) %>%
  summarise(prob_sum = sum(probability))


simdates <- bind_rows(shorelinedates) %>% group_by(bce, site_name) %>%
  filter(!is.na(probability)) %>%
  summarise(prob_sum = sum(probability))

ggplot() +
  # geom_line(data = spd, aes(x = years, y = prob_sum)) +
  geom_line(data = simdates, aes(x = bce, y = prob_sum, group = site_name),
            col = "red") +
  # geom_line(data = predgrid, aes(x = bce, y = prob_dens)) +
  scale_x_continuous(breaks = seq(-10000, 2500, 500),
                     limits = c(-10000, -2500))
scale_y_continuous(limits = c(0,0.5))


simdates <- bind_rows(shorelinedates, .id = "site_name")

simdates %>% mutate(simn = strsplit(site_name, split ="_"))

strsplit(simdates[1,]$site_name, split = "_")

group_by(years, ) %>%
  filter(!is.na(probability)) %>%
  summarise(prob_sum = sum(probability))

# Shoredate cKDE
sdats <- sdates  %>% group_by(site_name) %>% filter(!is.na(probability))

sampledats <- sampleshoredates(sdats, nsim = 10)
bootdats <- sampleshoredates(sdats, nsim = 100, boot = TRUE)

ckdedates <- ckdeshore(sampledats, bw = 50)
ckdeboot <- ckdeshore(bootdats, bw = 50)
ggplot(ckdedates) + geom_line(aes(x = years, y = prob, group = n), col = "grey") + xlim(c(-9500, -2500)) + theme_classic()
ggplot(ckdeboot) + geom_line(aes(x = years, y = prob, group = n), col = "grey") + xlim(c(-9500, -2500)) + theme_classic()


# Plot example date

target_pt <- sf::st_sfc(sf::st_point(c(567929, 6559112)), crs = 32632)
target_date <- shoredate::shoreline_date(site = target_pt, elevation = 56)
plt1 <- shoredate::shoredate_plot(target_date, hdr_label = FALSE, greyscale = TRUE)

target_curve <- shoredate::interpolate_curve(target_pt)
plt2 <- shoredate::displacement_plot(target_curve, greyscale = TRUE)

isobases <- st_read(here("analysis/data/raw_data/isobases.gpkg"))


plt3 <- ggplot() +
  geom_sf(data = bmap_reproj, fill = "grey", colour = NA) +
  geom_sf(data = target_pt, size = 2.5, shape = 21,
          fill = "black") +
  geom_sf(data = isobases, aes(linetype = name)) +
  geom_sf(data = target_isobase) +
  ggsn::scalebar(data = surveyed, dist = 20, dist_unit = "km",
                 transform = FALSE, st.size = 4, height = 0.02,
                 border.size = 0.1, st.dist = 0.03,
                 anchor = c(x = anc[2] - 18000, y = anc[1]) + 8000) +
  scale_linetype_manual(values = c("Horten" = "twodash",
                                  "Larvik" = "dashed",
                                  "Tvedestrand" = "dotted",
                                  "Arendal" = "longdash"), name = "")+
  coord_sf(xlim = c(bboxsites[1], bboxsites[3]),
           ylim = c(bboxsites[2], bboxsites[4]),
           expand = FALSE) +
  theme_bw() + theme(axis.title=element_blank(),
                     axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     rect = element_rect(),
                     axis.ticks = element_blank(),
                     panel.grid.major = element_blank(),
                     legend.position = "bottom")

# Retrieve the legend from the site plot
legend <- retrieve_legend(plt2)



plt2.2 <- plt2 + theme(legend.position = "none")
  # theme(legend.text=element_text(size = 12))
plt3.2 <- plt3 + theme(legend.position = "none")

((plt3.2 + plt2.2) / legend/ (plot_spacer() + plt1 + plot_spacer())) +
  plot_layout(guides = "collect", heights = c(3,0.5,3)) +
  plot_annotation(tag_levels = 'A')
ggsave(here::here("analysis/figures/shoredate.png"),
       units = "px", width = 2205*1.4, height = 2330)

plt2.2 <- plt2 + guides(linetype = guide_legend(nrow=2, byrow=TRUE),
              colour = guide_legend(nrow=2, byrow=TRUE)) +
  theme(legend.text=element_text(size = 12))
plt3.2 <- plt3 + theme(legend.position = "none")


plt3.2 + plt2.2  + plt1 + plot_layout(ncol = 3) +
  plot_annotation(tag_levels = 'A')
ggsave(here::here("analysis/figures/shoredate.png"),
       units = "px", width = 2205*2.2, height = 2380*1.2)


deg <- 327

# Create empty sf object to hold isobases for the current direction
isobases_temp <- sf::st_sf(id = 1:nrow(centrepoints),
                           name = NA, direction = NA,
                           crs = 32632,
                           geometry = sf::st_sfc(lapply(1:nrow(centrepoints),
                                                        function(x) sf::st_linestring())))

# Loop over crentre points and create isobases
isobase_length = 9000000
  # Find x and y coords
  x <- sf::st_coordinates(target_pt)[1]
  y <- sf::st_coordinates(target_pt)[2]

  # Find coords at the specified distance from the point at deg degree angle
  # adding isobase_length
  xx <- x + isobase_length * (cos(deg * pi / 180))
  yy <- y + isobase_length * (sin(deg * pi / 180))

  # Find coords at the specified distance from the point at deg degree angle
  # subtracting isobase_length
  xx2 <- x - isobase_length * (cos(deg * pi / 180))
  yy2 <- y - isobase_length * (sin(deg * pi / 180))

  # Create points from the identified coordinates
  pts <- sf::st_sfc(
    sf::st_multipoint(
      rbind(
        sf::st_coordinates(target_pt)[1,],
        c(xx, yy), c(xx2, yy2))))
  pts <-  sf::st_sf(sf::st_set_crs(pts, 32632)) # WGS 84 / UTM 32N

  # Combine these into a line and update the temporary data frame
  target_isobase <- sf::st_cast(pts, to = 'LINESTRING')
