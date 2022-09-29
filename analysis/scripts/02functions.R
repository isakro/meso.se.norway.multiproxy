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
                           specified_elev = FALSE,
                           elevation = NA,
                           interpolate_curve = TRUE,
                           interpolated_curve = NA){

  if(interpolate_curve == TRUE){
    sitecurve <- interpolate_curve(bce = xvals,
                                   target = site,
                                   dispdat = displacement_curves,
                                   isodat = isobases)
  } else{
    sitecurve <- interpolated_curve
  }

  if(class(sitecurve) == "list"){
    sitecurve <- rbind.data.frame(sitecurve)
  }

  if(specified_elev == TRUE){
    if(is.na(elevation)){
      return(NA)
    } else{
    siteelev <- elevation
    }
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
    bce = seq(-10000, 2000, 1),
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

      dategrid[dategrid$bce %in% year_range, "probability"] <-
        dategrid[dategrid$bce %in% year_range, "probability"] + prob
    }
  }

  if(class(site)[1] == "sf") {
    dategrid$site_name <- as.character(st_drop_geometry(site[1]))
  } else{
    dategrid$site_name <- site
  }

  # Normalise to sum to unity
  dategrid$probability <- dategrid$probability / sum(dategrid$probability)

  return(dategrid)
}

# Reverse a single calendar date to a corresponding elevation
# Repurposed from Bchron::unCalibrate()
reverse_shoredate <- function(shoreline_date,
                              displacement_curve){

  # if(class(displacement_curve) == "list"){
  #   displacement_curve <- rbind.data.frame(displacement_curve)
  # }

  lwelev <- round(stats::approx(displacement_curve$bce,
                displacement_curve$lowerelev,
                xout = shoreline_date,
                rule = 1
  )$y, 1)

  upelev <- round(stats::approx(displacement_curve$bce,
                displacement_curve$upperelev,
                xout = shoreline_date,
                rule = 1
  )$y, 1)

  if(!(is.na(lwelev) | is.na(upelev))){
    # In case the displacement curves intersect
    minelev <- min(lwelev, upelev)
    maxelev <- max(lwelev, upelev)

    # Make sure elevation is never 0 (i.e. present day sea-level)
    if(minelev <= 0){
      minelev <- 0.05
    }

    return(sample(seq(minelev, maxelev, 0.05), 1))
  } else{
    return(NA)
  }
}

# Repurposed from ckde() from rcarbon
ckdeshore <- function(sampledates, bw){

  # timerange <- c(max(sampledates$years), min(sampledates$years))
  # nsim <- length(unique(sampledates$site_name))
  # rawkde <- apply(as.data.frame(sampledates$years), 1, density, bw = bw, na.rm = TRUE,
  #                from = timerange[1], to = timerange[2])
  #
  # # resmatrix = data.frame(matrix(NA, nrow = length(timerange[1]:timerange[2]), ncol = nsim))
  # resdf <- data.frame(matrix(NA, nrow = nsim*length(timerange[1]:timerange[2]), ncol = 3))
  #
  # res <- list(nsim)
  # for (i in 1:nsim){
  #
  #
  #   tmpres <- data.frame(matrix(NA, nrow = length(timerange[1]:timerange[2]), ncol = 3))
  #   names(tmpres) <- c("years", "prob", "n")
  #   tmpres[,1] <- timerange[1]:timerange[2]
  #   tmpres[,2] <- approx(x = rawkde[[i]]$x, xout = timerange[1]:timerange[2],
  #                           y = rawkde[[i]]$y)$y
  #   tmpres[,3] <- i
  #   res[[i]] <- tmpres
  # }

  lsim <- length(unique(sampledates$nsim))

  res <- list(lsim)
  timerange <- c(max(sampledates$years), min(sampledates$years))
  for(i in 1:lsim){
    datedat <- sampledates[sampledates$nsim == i,]
    rawkde <- density(datedat$years, bw = bw, na.rm = TRUE, from = timerange[1], to = timerange[2],
                       weights= datedat$probability)

    tmpres <- data.frame(matrix(NA, nrow = length(timerange[1]:timerange[2]), ncol = 3))
    tmpres[,1] <- timerange[1]:timerange[2]
    tmpres[,2] <- approx(x = rawkde$x, xout = timerange[1]:timerange[2],
                         y = rawkde$y)$y
    tmpres[,3] <- i
    names(tmpres) <- c("years", "prob", "n")
    res[[i]] <- tmpres
  }

  return(bind_rows(res))
}

# Randomly sample dates from shoreline dates
sampleshoredates <- function(dates, nsim, boot = FALSE){
  nsites <- length(unique(dates$site_name))

  result <- list(nsites)
  if(!boot){
    for(i in 1:nsites){

      cdate <- dates[dates$site_name == unique(dates$site_name)[i],]

      index <- sample(1:nrow(cdate), size = nsim,
                    prob = cdate$probability, replace = TRUE)
      res <- cdate[index,]
      res$nsim <- 1:nsim
      result[[i]] <- res
    }
  }

  if(boot){
    bsindex = sample(nsites, replace=TRUE)
    for(i in 1:nsites){
      cdate <- dates[dates$site_name == unique(dates$site_name)[bsindex[i]],]
      # res <- sample(cdate$years, size = nsim,
      #               prob = cdate$probability, replace = TRUE)
      # result <- c(result, res)
      index <- sample(1:nrow(cdate), size = nsim,
                      prob = cdate$probability, replace = TRUE)
      res <- cdate[index,]
      res$nsim <- 1:nsim
      result[[i]] <- res
    }
  }
  return(bind_rows(result))
}


