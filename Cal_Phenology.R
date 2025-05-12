suppressMessages({
  library(data.table)
  library(magrittr)
  library(lubridate)
  library(purrr)
  library(plyr)
  library(ggplot2)
  library(phenofit)
  library(tidyverse)
  library(terra)
  library(zoo)
})

# Main function
main <- function() {
  # Configuration parameters
  config <- list(
    input_dir = "./data/",          # Input data directory
    output_file = "./GOME2-SIF.tif",  # Output file
    pattern = ".tif$",              # File pattern
    n_cores = 12                    
  )
  
  run_phenology_analysis(config)
}

# Main analysis pipeline
run_phenology_analysis <- function(config) {
  # 1. Load input data
  raster_data <- load_input_data(config$input_dir, config$pattern)
  
  # 2. Calculate phenology metrics
  cat("Starting phenology metrics calculation...\n")
  t1 <- proc.time()
  pheno_result <- terra::app(
    raster_data, 
    calculate_phenology, 
    cores = config$n_cores
  )
  t2 <- proc.time()
  cat(sprintf('Calculation completed, time elapsed: %.2f seconds\n', (t2 - t1)[3]))
  
  # 3. Save results
  cat("Saving results...\n")
  terra::writeRaster(pheno_result, filename = config$output_file)
  cat("Analysis completed, results saved to:", config$output_file, "\n")
}

# Load input data
load_input_data <- function(input_dir, pattern) {
  if (!dir.exists(input_dir)) {
    stop(paste("Input directory does not exist:", input_dir))
  }
  
  tif_files <- list.files(path = input_dir, pattern = pattern, full.names = TRUE)
  raster_data <- terra::rast(tif_files)
  
  return(raster_data)
}

# Generate time series data (2007-01-01 to 2018-12-27 with 8-day intervals)
generate_time_series <- function() {
  start_date <- as.Date("2007-01-01")
  end_date <- as.Date("2018-12-31")

  time_series <- seq(from = start_date, to = end_date, by = "8 days")
  
  return(time_series)
}

# Phenology calculation function
calculate_phenology <- function(x) {
  if (length(na.omit(x[1:552])) < 276) return(rep(NA, 48))
  
  tryCatch({
    SOS_TRS3 <- rep(NA, 12)
    EOS_TRS3 <- rep(NA, 12)
    POS_TRS3 <- rep(NA, 12)
    POP_TRS3 <- rep(NA, 12)
    
    params <- list(
      wFUN = phenofit::wBisquare,
      ymax_min = 0.1,
      rymin_less = 0.8,
      nptperyear = 46,
      w_critical = 0.2,
      south = FALSE
    )
    
    time_series <- generate_time_series()
    
    sif_data <- process_sif_data(x[1:552], time_series, params)
    
    pheno_metrics <- compute_pheno_metrics(sif_data, params)
    
    if (!is.null(pheno_metrics)) {
      SOS_TRS3 <- pheno_metrics$SOS
      EOS_TRS3 <- pheno_metrics$EOS
      POS_TRS3 <- pheno_metrics$POS
      POP_TRS3 <- pheno_metrics$POP
    }
    
    return(c(SOS_TRS3, EOS_TRS3, POS_TRS3, POP_TRS3))
    
  }, error = function(e) {
    return(rep(NA, 48))
  })
}

process_sif_data <- function(sif_site, t, params) {
  n <- length(sif_site)
  w <- rep(1, n)
  
  w[sif_site < 0] <- 0.5
  sif_site[sif_site < 0] <- 0
  
  w[is.na(sif_site)] <- 0.01
  sif_site[is.na(sif_site)] <- 0
  
  for (i in seq_len(params$nptperyear)) {
    pts <- seq(i, length(sif_site), by = params$nptperyear)
    mean_val <- mean(sif_site[pts])
    sd_val <- sd(sif_site[pts])
    lower <- mean_val - 3 * sd_val
    upper <- mean_val + 3 * sd_val
    
    for (j in seq_along(pts)) {
      if (sif_site[pts[j]] < lower | sif_site[pts[j]] > upper) {
        w[pts[j]] <- params$w_critical
        sif_site[pts[j]] <- mean_val
      }
    }
  }
  
  return(list(
    t = t,
    y = sif_site,
    w = w,
    QC_flag = "good"
  ))
}

compute_pheno_metrics <- function(data, params) {
  INPUT_SIF <- phenofit::check_input(
    data$t, data$y, data$w, data$QC_flag,
    params$nptperyear, params$south,
    maxgap = params$nptperyear / 4, alpha = 0.02,
    ymin = 0.05
  )
  
  brks2_sif <- phenofit::season_mov(
    INPUT_SIF,
    list(
      rFUN = "smooth_wHANTS",
      wFUN = params$wFUN,
      maxExtendMonth = 12,
      minpeakdistance = 1 / 3 * params$nptperyear,
      r_min = 0.2, r_max = 0.6,
      rtrough_max = 0.6,
      calendarYear = TRUE
    )
  )
  
  if (is.null(brks2_sif)) return(NULL)
  
  fit_sif <- phenofit::curvefits(
    INPUT_SIF, brks2_sif,
    options = list(
      methods = c("Elmore"),
      iters = 3,
      wFUN = params$wFUN,
      nextend = 2,
      maxExtendMonth = 3,
      minExtendMonth = 1,
      minPercValid = 0.2
    )
  )
  
  pheno_sif <- phenofit::get_pheno(fit_sif, "Elmore", TRS = 0.3, IsPlot = FALSE)
  day_sif <- pheno_sif$doy$Elmore
  
  sif_fit <- purrr::map2_dfr(
    fit_sif[1:12], 1:12,
    ~ {
      par <- .x$model$Elmore$par
      tout <- .x$model$Elmore$tout
      y_sif <- phenofit::doubleLog.Elmore(par, tout)
      y0_sif <- rep(NA, 365)
      sifmax <- max(y_sif)
      pos_sif <- which.max(y_sif)
      
      if (!is.na(day_sif$DER.pos[.y])) {
        interal_sif <- day_sif$DER.pos[.y] - pos_sif
        if (length(y_sif) > 365 - interal_sif) {
          y0_sif[interal_sif:365] <- y_sif[1:(365 - interal_sif)]
        } else {
          y0_sif[interal_sif:(length(y_sif) + interal_sif)] <- y_sif
        }
      }
      data.frame(doys = 1:365, sif_y0 = y0_sif)
    }
  )
  
  sif_yr_doy <- stats::aggregate(sif_fit$sif_y0, by = list(type = sif_fit$doys), mean, na.rm = TRUE)
  sif_yr_doy <- zoo::na.locf(sif_yr_doy$x)
  sif_yr_doy <- zoo::na.locf(sif_yr_doy, fromLast = TRUE)
  
  fit_doy <- phenofit::curvefit(sif_yr_doy, seq(1, 365, 1), seq(1, 365, 1), methods = c("Elmore"))
  l_pheno <- phenofit::get_pheno(fit_doy, "Elmore", TRS = 0.3, IsPlot = FALSE)
  threshold <- sif_yr_doy[l_pheno$TRS3[1]]
  
  # Calculate SOS and EOS
  for (num_yr in 1:12) {
    sif_test <- sif_fit$sif_y0[((num_yr - 1) * 365 + 1):(num_yr * 365)]
    if (length(na.omit(sif_test)) > 300) {
      SOS_TRS3[num_yr] <- min(which(sif_test > threshold))
      EOS_TRS3[num_yr] <- max(which(sif_test > threshold))
    }
  }
  
  # Get POS and POP
  POS_TRS3 <- pheno_sif$doy$Elmore$DER.pos
  
  gpp <- purrr::map2_dfr(
    fit_sif[1:12], 1:12,
    ~ {
      par <- .x$model$Elmore$par
      tout <- .x$model$Elmore$tout
      y <- phenofit::doubleLog.Elmore(par, tout)
      data.frame(gppmax = max(y))
    }
  )
  POP_TRS3 <- gpp$gppmax
  
  return(list(
    SOS = SOS_TRS3,
    EOS = EOS_TRS3,
    POS = POS_TRS3,
    POP = POP_TRS3
  ))
}

# Execute main function
if (!interactive()) {
  main()
}