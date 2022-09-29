library(rcarbon)
library(ADMUR)

c14 <- read.csv(here::here("analysis/data/raw_data/radiocarbon.csv"))

caldates <- calibrate(x = c14$c14_bp, normalised = FALSE,
                      errors = c14$error, calCurves = "intcal20")
c14spd <- spd(caldates, timeRange = c(10000, 2000))


plot(c14spd, calendar = "BCAD")



c14bins = binPrep(sites = c14$site_name, ages = c14$c14_bp, h = 200)
c14spd2 <- spd(caldates, bins = c14bins, timeRange = c(11000, 2000))

par(mfrow = c(2,1))
plot(c14spd, calendar = "BCAD")
plot(c14spd2, ylim = c(0, 0.11), xlim = c(-9300, -2000), calendar = "BCAD")

caldates <- calibrate(x = c14$c14_bp, normalised = FALSE,
                      errors = c14$error, calCurves = "intcal20")
c14randates <- sampleDates(caldates, bins = c14bins, nsim = 1000,
                           boot = TRUE, verbose = FALSE)
c14ckde <- ckde(c14randates, timeRange = c(10000, 4500), bw = 100,
                normalised = FALSE)

par(mfrow = c(1,1))
plot(c14ckde, calendar = "BCAD", xlim=c(-8000, -2500))

nsim = 1000
expnull <- modelTest(caldates, errors = c14$error, bins = c14bins, nsim = nsim,
                     timeRange = c(12000, 4500), model = "exponential",
                     runm = 100)

plot(expnull, calendar = "BCAD")

unfnull <- modelTest(caldates, errors = c14$error, bins = c14bins, nsim = nsim,
                     timeRange = c(12000, 4000), model = "uniform",
                     runm = 100)
plot(unfnull, calendar = "BCAD", xlim = c(-9300, -2000))

ggplot(data = unfnull$result, aes(x = (calBP-1950)*-1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              colour = "lightgrey", alpha = 0.4) +
  labs(x = "BCE", y = "Summed probability") +
  xlim(c(-9300, -2500)) +
  geom_line(aes(y = PrDens)) +
  theme_bw()

# save(unfnull,
#     file = here("analysis/data/derived_data/unfnull.rds"))

ggplot(data = expnull$result, aes(x = (calBP-1950)*-1)) +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                colour = "lightgrey", alpha = 0.5) +
  geom_line(aes(y = PrDens)) +
    theme_bw()


# ADMUR
c14admur <- c14 %>% rename("site" = "site_name",
                           "age" = "c14_bp", "sd" = "error")
calar <- makeCalArray(calcurve = intcal20, calrange = c(11000, 4000))
cal <- summedCalibrator(c14admur, calar)
plotPD(cal)

# Phase the data
x <- phaseCalibrator(data = c14admur, CalArray = calar)
plotPD(x)

# Sum by site
spdadmur <- as.data.frame(rowSums(x))
# Normalise
spdadmur <- spdadmur/(sum(spdadmur) * calar$inc)

# Wrapper function should do the same
spdwrap <- summedPhaseCalibrator(data = c14admur, calcurve = intcal20,
                                 calrange =  c(11000, 4000))

# These should be identical
par(mfrow = c(2,1))
plotPD(spdadmur)
plotPD(spdwrap)

