library(tidyverse)


sites <- read.csv(here::here("analysis/data/raw_data/sites.csv"))

curation_sites <- sites %>%
  # Exclude multiphased sites and sites that have not been manually excavated
  filter(multiphase == "f" | !is.na(excavated_m3)) %>%
  filter(dating_method == "c14" | dating_method == "shore_typo"
         & !is.na(reported_min_elev)) %>%
  rowwise() %>%
  mutate(
    # Find total artefact count and flint count for each site
    mean_date = mean(c(reported_earliest_bce, reported_latest_bce),
                     na.rm = TRUE),
    artefact_count = sum(c_across(core:chip_retouch_nonflint)),
                     na.rm = TRUE,
    flint_count = sum(c_across(core:flake_retouch), na.rm = TRUE),
    nonflint_count = sum(c_across(axe_nonflint:chip_retouch_nonflint),
                         na.rm = TRUE)) %>%
  ungroup() %>%
  # Exclude sites with a artefact count < 200
  filter(artefact_count > 200) %>%
  # Prepare for curation index
  select(-beach_flint_flint_nodule, -hammerstone) %>%
  mutate(
    nonflint_prop = nonflint_count / artefact_count,
    secondary_prop =
      sum(axe_debitage, axe, chisel,
          microlith, tanged_point, transverse_point, single_edged_point,
          høgnipen_point, pressure_flaked_point, nøklegaard_point,
          projectile_undefined, sickle, dagger_spear, microburin, scraper,
          knife, burin, drill, blade_retouch, microblade_retouch,
          chip_retouch, fragment_retouch, flake_retouch, axe_nonflint,
          axe_debitage_nonflint, chisel_nonflint, club_hatchet_nonflint,
          grinding_slab_stone, projectile_nonflint, knife_nonflint,
          scraper_nonflint, drill_nonflint, burin_nonflint,
          microblade_retouch_nonflint, blade_retouch_nonflint,
          fragment_retouch_nonflint,flake_retouch_nonflint,
          chip_retouch_nonflint) / artefact_count,
    secondary_prop_flint =
      sum(axe_debitage, axe, chisel,
          microlith, tanged_point, transverse_point, single_edged_point,
          høgnipen_point, pressure_flaked_point, nøklegaard_point,
          projectile_undefined, sickle, dagger_spear, microburin, scraper,
          knife, burin, drill, blade_retouch, microblade_retouch,
          chip_retouch, fragment_retouch, flake_retouch) / flint_count,
    lithic_dens = artefact_count / excavated_m3,
    flint_dens = flint_count / excavated_m3) %>%
    mutate(minmaxdens =  (lithic_dens - min(lithic_dens, na.rm = TRUE)) /
             (max(lithic_dens, na.rm = TRUE) - min(lithic_dens, na.rm = TRUE)),
       minmaxsecondary = (secondary_prop - min(secondary_prop, na.rm = TRUE)) /
    (max(secondary_prop, na.rm = TRUE) - min(secondary_prop, na.rm = TRUE)),
    minmaxdens_flint =  (flint_dens - min(flint_dens, na.rm = TRUE)) /
      (max(flint_dens, na.rm = TRUE) - min(flint_dens, na.rm = TRUE)),
    minmaxsecondary_flint = (secondary_prop_flint - min(secondary_prop_flint,
                                                        na.rm = TRUE)) /
      (max(secondary_prop_flint, na.rm = TRUE) - min(secondary_prop_flint,
                                                     na.rm = TRUE))) %>%
  rowwise() %>%
  mutate(curation_index = mean(c(minmaxdens * -1, minmaxsecondary))) %>%
  mutate(curation_index_flint = mean(c(minmaxdens_flint * -1,
                                   minmaxsecondary_flint)))

curation_sites$normalised_index <-
  (curation_sites$curation_index - min(curation_sites$curation_index)) /
  (max(curation_sites$curation_index) - min(curation_sites$curation_index))




curation_sites <- curation_sites %>% mutate(weighted_index =
  ifelse(nonflint_prop > 0,
         normalised_index / (nonflint_prop^0.0005), normalised_index)
) %>% filter(nonflint_prop < 0.9)

ggplot(curation_sites, aes(mean_date, weighted_index)) +
  geom_point() +
  scale_x_continuous(breaks = c(seq(-9000,-4000, 1000)),
                     limits = c(-9200, -3500)) +
  geom_smooth(method = "loess")

ggplot(curation_sites, aes(mean_date, normalised_index)) +
  geom_point() +
  scale_x_continuous(breaks = c(seq(-9000,-4000, 1000)),
                     limits = c(-9200, -3500))+
  geom_smooth(method = "loess")

ggplot(curation_sites) +
  geom_point(aes(lithic_dens, secondary_prop, col = mean_date)) +
  scale_y_log10() +
  scale_x_log10()

curation_sites %>% select(name, mean_date) %>% View()
curation_sites %>% select(name, lithic_dens) %>% View()
curation_sites %>% select(name, curation_index) %>% View()
curation_sites %>% select(name, curation_index_flint) %>% View()

ggplot(curation_sites) +
  geom_point(aes(mean_date, curation_index)) +
  scale_x_continuous(breaks = c(seq(-9000,-4000, 1000)))

ggplot(curation_sites) +
  geom_point(aes(mean_date, normalised_index)) +
  scale_x_continuous(breaks = c(seq(-9000,-4000, 1000)))

ggplot(curation_sites) +
  geom_point(aes(mean_date, weighted_index)) +
  scale_x_continuous(breaks = c(seq(-9000,-4000, 1000)))


curation_sites %>%
  filter(nonflint_prop < 0.9) %>%
ggplot(aes(mean_date, curation_index_flint)) +
  geom_point() +
  scale_x_continuous(breaks = c(seq(-9000,-4000, 1000)),
                     limits = c(-9200, -3500)) +
  geom_smooth(method = "loess")

curation_sites %>% select(name, mean_date, curation_index_flint) %>%  View()

