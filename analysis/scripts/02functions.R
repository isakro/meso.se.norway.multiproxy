# Function to interpolate displacement curve using IDW
interpolate_curve <- function(years, target, dispdat, isodat){

  dists <- as.data.frame(st_distance(target, isobases))
  names(dists) <- isobases$name

  values <- data.frame(matrix(ncol = 3, nrow = length(years)))
  names(values) <- c("years", "lowerelev", "upperelev")

  # In the case that a site is on the isobase of a
  # displacement curve (e.g. Alveberget 8), simply return that displacement
  # curve
  if(any(as.numeric(dists) == 0)){
    values <-
      dispdat[dispdat$name == names(dists)[which(as.numeric(dists) == 0)],]

  } else { for(i in 1:length(years)){
    for(j in 1:ncol(dists)){
      le <- dispdat[which(dispdat$name == names(dists)[j] &
                            dispdat$years == years[i]),
                    "lowerelev"]

      ue <- dispdat[which(dispdat$name == names(dists)[j] &
                            dispdat$years == years[i]),
                    "upperelev"]

      dists[2, j] <- le
      dists[3, j] <- ue
    }
    distdat <- as.data.frame(t(dists))
    names(distdat) <- c("distance", "lower", "upper")

    # No sites are older than the lowest limit of any displacement curve
    # so in case of NA values, simply assign NA
    if(any(is.na(distdat))){
      lowerval <- upperval <- NA
    } else {
      # Inverse distance weighting
      lowerval <- sum(apply(distdat, 1,
                            function(x) x["lower"] * x["distance"]^-2)) /
        sum(apply(distdat, 1, function(x) x["distance"] ^-2))
      upperval <- sum(apply(distdat, 1,
                            function(x) x["upper"] * x["distance"]^-2)) /
        sum(apply(distdat, 1, function(x) x["distance"] ^-2))

    }
    values[i,] <- c(years[i], lowerval, upperval)
  }
  }
  return(values)
}

shoreline_date <- function(site,
                           elev,
                           displacement_curves,
                           isobases,
                           reso = 0.1,
                           expratio,
                           siteelev = "mean",
                           specified_elev = NA){

  sitecurve <- interpolate_curve(years = xvals,
                                 target = site,
                                 dispdat = displacement_curves,
                                 isodat = isobases)

  if(!(is.na(specified_elev))){
    siteelev <- specified_elev
  } else{
    siteelev <- terra::extract(elev, vect(site), fun = siteelev)[,-1]
  }

  inc <- seq(0, siteelev, reso)

  expdat <- data.frame(
    offset = inc,
    px = pexp(inc, rate = expratio)) %>%
    mutate(probs = px - lag(px, default =  dplyr::first(px))) %>%
    tail(-1) %>%
    filter(px < 0.99999) # Probability cut-off

  dategrid <- data.frame(
    years = seq(-10000, 2000, 1),
    probability = 0)

  for(i in 1:nrow(expdat)){
    # Subtract offset
    adjusted_elev <- as.numeric(siteelev - expdat$offset[i])

    # Find lower date
    lowerd <- round(approx(sitecurve[,"lowerelev"],
                           xvals, xout = adjusted_elev)[['y']])

    # Find upper date
    upperd <- round(approx(sitecurve[,"upperelev"],
                           xvals, xout = adjusted_elev)[['y']])

    # Find youngest and oldest date
    earliest <- min(c(lowerd, upperd))
    latest <- max(c(lowerd, upperd))

    # Add probability to each year in range
    if(!is.na(latest) && !is.na(earliest)){

      year_range <- seq(earliest, latest, 1)
      prob <- 1/length(year_range)*expdat$probs[i]

      dategrid[dategrid$years %in% year_range, "probability"] <-
        dategrid[dategrid$years %in% year_range, "probability"] + prob
    }
  }
  dategrid$site_name <- as.character(st_drop_geometry(site[1]))

  # Normalise to sum to unity
  dategrid$probability <- dategrid$probability / sum(dategrid$probability)

  return(dategrid)
}
