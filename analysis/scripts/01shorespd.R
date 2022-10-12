library(tidyverse)
library(sf)
library(terra)
library(here)
library(shoredate)
library(patchwork)
library(ADMUR)
library(DEoptimR)

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
  # shorelinedates[[i]] <- shoreline_date(,
  #                                     dtm,
  #                                     displacement_curves,
  #                                     isobases,
  #                                     expratio = 0.168,
  #                                     reso = step_reso)

    shorelinedates[[i]] <- shoredate::shoreline_date(
    site = ssites[i,],
    elev_raster = dtm,
    elev_reso = step_reso,
    cal_reso = 5,
    sparse = TRUE)


}

sdates <- bind_rows(shorelinedates)

# save(sdates,
#       file = here("analysis/data/derived_data/sdates_5.rds"))

# save(sdates,
#       file = here("analysis/data/derived_data/sdates.RData"))

load("analysis/data/derived_data/sdates.RData")

spd <- sdates %>% rename("bce" = "years") %>%  filter(bce < -2500)
spdn <- spd %>% filter(probability != 0)
spdn <- length(unique(spdn$site_name)) # Couldn't figure out n_distinct() here

spd <- spd %>% group_by(bce) %>%
  filter(!is.na(probability)) %>% # Excludes Askeladden ID 108682-1 at 209masl
  group_by(site_name) %>%
  mutate(prob_sum = sum(probability))

exnoprob <- spd %>% filter(prob_sum > 0) %>%
  group_by(bce) %>%
  filter(probability != 0) %>%
  summarise(prob_sum = sum(probability))
# Excludes "159969"   "94301-1"  "144111-1" "144113-1" "144114-1",
# "138457-1" "230584-0" "265781-0" "265784-0"
lower_lim <- min(exnoprob$bce)
upper_lim <- max(exnoprob$bce)

exp <- TRUE

if(exp){
  # Code repurposed from modelTest() from rcarbon
  # Exponential
  print("Exponential")
  fit <- nls(y ~ (exp(a + b * x)),
             data = data.frame(x = exnoprob$bce, y = exnoprob$prob_sum),
             start = list(a = 0, b = 0))
  est <- predict(fit, list(x = exnoprob$bce))
  predgrid <- data.frame(bce = exnoprob$bce, prob_dens = est)
  predgrid <- filter(predgrid, bce >= lower_lim & bce <= upper_lim)
} else{
  print("Uniform")
  # # Uniform
  predgrid <- data.frame(bce = exnoprob$bce,
                         prob_dens = mean(exnoprob$prob_sum))
  predgrid <- filter(predgrid, bce > lower_lim & bce < upper_lim)
}

ggplot(exnoprob) + geom_line(aes(x = bce, y = prob_sum)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens)) +
  scale_x_continuous(breaks = seq(-10000, 2500, 500))

ggplot(spd) + geom_line(aes(x = bce, y = prob_sum)) +
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
random_dates$displacement_curve <- incpolys[sample(incpolys$id, replace = TRUE,
                                                   size = ssize*nsim,
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
    elev_reso = step_reso,
    cal_reso = 1,
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
  geom_line(data = exnoprob, aes(x = bce, y = prob_sum)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens),
            linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500)) +
  labs(x = "BCE", y = "Summed probability") +
  theme_bw()

ggplot(exnoprob) + geom_line(aes(x = bce, y = prob_sum)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens)) +

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


library(cowplot)

toprow <- plot_grid(plt3.2, plt2.2, labels = c('A', 'B'))
bottomrow <- plot_grid(NULL, plt1, NULL, ncol = 3,
                       rel_widths = c(1.2, 2, 1.2),
                       labels = c('', 'C', ""))

cowplot::plot_grid(toprow, legend, bottomrow, nrow = 3,
                   rel_heights = c(5/11, 1/11, 5/11))

ggsave(here::here("analysis/figures/shoredate.png"), width = 30, height = 30,
       units = 'cm')


gridExtra::grid.arrange(gridExtra::arrangeGrob(plt3.2 , plt2.2, nrow = 1),
                        legend, plt1, nrow = 3, widths = c(2, 2, 1),
                        heights = c(10,1, 10))

island_hist <- grid.arrange(arrangeGrob(hull_isl +
                             theme(legend.position = 'none'),
                             buff_isl + theme(legend.position = 'none'),
                             sites_isl + theme(legend.position = 'none'), nrow = 1),
                            legend, nrow = 2, heights = c(10,1))




ggsave('../figures/island_hist.png', island_hist, width = 15, height = 8,
       units = 'cm', dpi = 600)


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


# Prepare shoreline spd for admur

tst <- exnoprob[1:500,]

spd <- as.data.frame(tst)
spd$prob_sum <- spd$prob_sum/sum(spd$prob_sum)

spd$years <- (spd$years-1950)*-1
maxbp <- max(spd$years)
minbp <- min(spd$years)


row.names(spd) <- spd$years
spd[1] <- NULL
names(spd) <- NA


spd$years <- (spd$years-1950)*-1
maxbp <- max(spd$years)
minbp <- min(spd$years)

load("analysis/data/derived_data/sdates.rds")

spd <- sdates %>% filter(bce < -2500) %>%
  group_by(bce) %>%
  filter(!is.na(probability)) %>% # Excludes Askeladden ID 108682-1 at 209masl
  group_by(site_name) %>%
  mutate(prob_sum_site = sum(probability)) %>%
  filter(prob_sum_site > 0) %>%
  group_by(bce) %>%
  mutate(prob_sum = sum(probability)) %>%
  filter(prob_sum > 0)

# Prep data for JDEoptim
pd <- sdates %>% filter(bce < -2500) %>%
  group_by(bce) %>%
  filter(!is.na(probability)) %>%
  mutate(prob_sum = sum(probability)) %>%
  # ungroup() %>%
  filter(prob_sum > 0) %>%
  group_by(site_name) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = site_name, values_from = probability)

pd <- as.data.frame(pd[!(is.na(pd[,4])),])
row.names(pd) <- pd[,"bce"]
# yrnames <- pd[,"bce"]
pd <- select(pd, -bce, -prob_sum, -row)
pd <- pd[, which(colSums(pd) != 0)]

# scaledpd <- as.data.frame(scale(pd, center = FALSE, scale = colSums(pd)))
# row.names(scaledpd) <- yrnames


spd <- spd %>%
  group_by(bce) %>%
  summarise(prob_sum = sum(probability)) %>%
  column_to_rownames(var = "bce")

maxbce <- max(as.numeric(row.names(spd)))
minbce <- min(as.numeric(row.names(spd)))

cpl1 <- JDEoptim(lower = 0, upper = 1, fn = objectiveFunction,
                 PDarray = pd, type = 'CPL', NP = 20, trace = TRUE)

start_time <- Sys.time()
cpl2 <- JDEoptim(lower = rep(0,3), upper = rep(1,3), fn = objectiveFunction,
                 PDarray = pd, type = 'CPL', NP = 60, trace = TRUE)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
cpl3 <- JDEoptim(lower = rep(0,5), upper = rep(1,5), fn = objectiveFunction,
                 PDarray = pd, type = 'CPL', NP = 100, trace = TRUE)
end_time <- Sys.time()
(cpl3t <- end_time - start_time)

start_time <- Sys.time()
cpl4 <- JDEoptim(lower = rep(0,7), upper = rep(1,7), fn = objectiveFunction,
                  PDarray = pd, type = 'CPL', NP = 140, trace =TRUE)
end_time <- Sys.time()
end_time - start_time
(cpl4t <- end_time - start_time)

exp <- JDEoptim(lower = -0.01, upper = 0.01, fn=objectiveFunction,
                PDarray = pd, type = 'exp', NP = 20)
uniform <- objectiveFunction(NULL, pd, type='uniform')

save(cpl1, cpl2, cpl3, cpl4, file = here::here("analysis/data/derived_data/cpl2.RData"))
# save(cpl1, cpl2, file = here::here("analysis/data/derived_data/cpl.RData"))

h1 <- CPLparsToHinges(pars=cpl1$par, years=minbce:maxbce)
h2 <- CPLparsToHinges(pars=cpl2$par, years=minbce:maxbce)
h3 <- CPLparsToHinges(pars=cpl3$par, years=minbp:maxbp)
h4 <- CPLparsToHinges(pars=cpl4$par, years=minbce:maxbce)

ggplot() + geom_point(data = h4, aes(x = year, y = pdf))

cols <- c('firebrick','orchid2','coral2','steelblue','goldenrod3')

plotPD(spd/(sum(spd) * 5))
lines(h1$year, h1$pdf, lwd=2, col=cols[1])
lines(h2$year, h2$pdf, lwd=2, col=cols[2])
lines(h3$year, h3$pdf, lwd=2, col=cols[3])
lines(h4$year, h4$pdf, lwd=2, col=cols[4])


modelpd1 <- convertPars(pars = cpl1$par, years = minbce:maxbce, type='CPL')
modelpd2 <- convertPars(pars = cpl2$par, years = minbce:maxbce, type='CPL')
modelpd3 <- convertPars(pars = cpl3$par, years = minbce:maxbce, type='CPL')
modelpd4 <- convertPars(pars = cpl4$par, years = minbce:maxbce, type='CPL')
modelexp <- convertPars(pars = exp$par, years = minbce:maxbce, type='CPL')

plotPD(spd/(sum(spd) * 5)) # Remember to multiply by calendar resolution
lines(modelpd1$year, modelpd1$pdf, col = cols[1], lwd=2)
lines(modelpd2$year, modelpd2$pdf, col = cols[2], lwd=2)
lines(modelpd3$year, modelpd3$pdf, col = cols[3], lwd=2)
lines(modelpd4$year, modelpd4$pdf, col = cols[4], lwd=2)

ggplot() +
  geom_line(data = spd/(sum(spd)*5), aes(x = as.numeric(row.names(spd)), y = prob_sum)) +
  geom_line(data = modelpd1, aes(x = year, y = pdf), col = cols[1]) +
  geom_line(data = modelpd2, aes(x = year, y = pdf), col = cols[2])+
  geom_line(data = modelpd3, aes(x = year, y = pdf), col = cols[3]) +
  # geom_line(data = modelpd4, aes(x = year, y = pdf), col = cols[4]) +
  geom_line(data = modelexp, aes(x = year, y = pdf))

data.frame(L1= -cpl1$value,
           L2= -cpl2$value,
           L3= -cpl3$value,
           L4= -cpl4$value,
           Lexp= -exp$value,
           Lunif= -uniform)

BIC.1 <- 1*log(303) - 2*(-cpl1$value)
BIC.2 <- 3*log(303) - 2*(-cpl2$value)
BIC.3 <- 5*log(303) - 2*(-cpl3$value)
BIC.4 <- 7*log(303) - 2*(-cpl4$value)
BIC.exp <- 1*log(303) - 2*(-exp$value)
BIC.uniform <- 0 - 2*(-uniform)
data.frame(BIC.1,BIC.2,BIC.3,BIC.4,BIC.exp,BIC.uniform)
