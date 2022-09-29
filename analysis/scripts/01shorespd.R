library(tidyverse)
library(sf)
library(terra)
library(here)
library(shoredate)

surveyed <- st_read(here("analysis/data/raw_data/surveyed_sites.gpkg"))
excavated <- st_zm(read_sf(here("analysis/data/raw_data/excavated_sites.gpkg")))
dtm <- rast("/home/isak/phd/eaa_presentation/dtm10/dtm10.tif")
site_dat <- read.csv(here("analysis/data/raw_data/sites.csv"))

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

# Select and rename columns necessary for shoreline dating
rcolexcshore <- excshore %>% select(ask_id, elev)
rcolsurvq <- survq  %>% select(askeladden_id, elev) %>%
  rename(ask_id = askeladden_id)

# Combine the data
ssites <-rbind(rcolexcshore, rcolsurvq)

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

# Code repurposed from modelTest() from rcarbon
# Exponential
# fit <- nls(y ~ exp(a + b * x),
#            data = data.frame(x = spd$years, y = spd$prob_sum),
#            start = list(a = 0, b = 0))
# est <- predict(fit, list(x = spd$years))
# predgrid <- data.frame(bce = spd$years, prob_dens = est)
# predgrid <- filter(predgrid, bce >= lower_lim & bce <= upper_lim)

# # Uniform
predgrid <- data.frame(bce = spd$years, prob_dens = mean(spd$prob_sum))
predgrid <- filter(predgrid, bce > lower_lim & bce < upper_lim)

ggplot(spd) + geom_line(aes(x = years, y = prob_sum)) +
  geom_line(data = predgrid, aes(x = bce, y = prob_dens)) +
  scale_x_continuous(breaks = seq(-10000, 2500, 500))

incpolys$dens <- lengths(st_intersects(incpolys, survq)) /
  sum(lengths(st_intersects(incpolys, survq)))

ssize = spdn
nsim <- 10
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
start_time <- Sys.time()
shorelinedates <- mapply(shoredate::shoreline_date,
                         site = random_dates$simn,
                         reso = step_reso,
                         elevation = random_dates$reverse_elevation,
                         interpolated_curve = random_dates$displacement_curve)
end_time <- Sys.time()
end_time - start_time

simdates <- lapply(shorelinedates, '[[', 'date') %>% bind_rows() %>%
  group_by(bce, site_name) %>%
  filter(!is.na(probability)) %>%
  summarise(prob_sum = sum(probability))

# 299

shorelinedates <- list()
start_time <- Sys.time()
for (i in 1:nrow(random_dates)){
  print(paste(i, "/", nrow(random_dates)))

  shorelinedates[[i]] <- shoreline_date(
                 site = random_dates$simn[i],
                 elev = NA,
                 displacement_curves = displacement_curves,
                 isobases = isobases,
                 reso = step_reso,
                 expratio = 0.168,
                 siteelev = "mean",
                 specified_elev = TRUE,
                 elevation = random_dates$reverse_elevation[i],
                 interpolate_curve = FALSE,
                 interpolated_curve = random_dates$displacement_curve[i])

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
}
end_time <- Sys.time()
end_time - start_time

tst <- list()
for(i in 1:1000){
  tst[[i]] <- i
  if (i %% 100 == 0){
    print(tst[(i+1-100):i])
    tst <- list()
  }
}


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

  # Save data externally and overwrite shorelinedate for memory purposes
  if(i %% ssize == 0){
    tmp <- bind_rows(shorelinedates) %>% group_by(bce, site_name) %>%
      filter(!is.na(probability)) %>%
      summarise(prob_sum = sum(probability))

    save(tmp,
         file = paste0(here::here("../external_data/shorespd/"),
                       i, "shore.rds"))

    shorelinedates <- list()
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

results <- list()
resfiles <- list.files(here::here("../external_data/shorespd/"),
                       full.names = TRUE)
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
         mean = mean(prob_sum))

ggplot() +
  geom_ribbon(data = simdates, aes(x = bce, ymin = low, ymax = high),
              fill = "grey", alpha = 0.9) +
  geom_line(data = simdates, aes(x = bce, y = mean), col = "red", lwd = 0.5) +
  geom_line(data = spd, aes(x = years, y = prob_sum)) +
  # geom_line(data = predgrid, aes(x = bce, y = prob_dens)) +
  scale_x_continuous(breaks = seq(-10000, 2000, 1000),
                     limits = c(-10000, -2500)) +
  labs(x = "BCE", y = "Summed probability") +
  theme_bw()





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

