# Utility function to retrieve a legend for having one legend with
# multiple ggplots (used for the island histograms)
retrieve_legend <- function(plot){
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
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


