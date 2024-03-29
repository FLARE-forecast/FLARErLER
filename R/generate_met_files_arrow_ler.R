##' @title Convert historical meteorology and NOAA forecasts to GLM format via arrow connection
##' @details Function combines historical meteorology and NOAA forecasts to create meteorology input files in the GLM format.  A file is generated for each ensemble member.
##' @param obs_met_file string; full path to netcdf that is observed historical meteorology
##' @param out_dir string; full path to directory where the converted files will be saved
##' @param forecast_dir string; full path to arrow file system with the NOAA forecast netcdf files
##' @return list; vector of full path for the converted files and boolean flag if issues with historical meteorology files
##' @export
##' @import dplyr
##' @import ncdf4
##' @importFrom stringr str_sub str_split str_detect
##' @importFrom tibble tibble
##' @importFrom lubridate as_datetime days hours ymd_hm
##' @author Quinn Thomas
##' @examples
##' \dontrun{
##' met_out <- FLAREr::generate_glm_met_files(obs_met_file = observed_met_file, out_dir = config$file_path$execute_directory, forecast_dir = config$file_path$noaa_directory, config)
##' }
generate_met_files_arrow_ler <- function(obs_met_file = NULL,
                                         out_dir,
                                         forecast_dir = NULL,
                                         config,
                                         use_s3 = FALSE,
                                         bucket = NULL,
                                         endpoint = NULL,
                                         local_directory = NULL){


  start_datetime <- lubridate::as_datetime(config$run_config$start_datetime)
  if(is.na(config$run_config$forecast_start_datetime)){
    end_datetime <- lubridate::as_datetime(config$run_config$end_datetime) #- lubridate::hours(1)
    forecast_start_datetime <- end_datetime
  }else{
    forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime)
    end_datetime <- forecast_start_datetime + lubridate::days(config$run_config$forecast_horizon) #- lubridate::hours(1)
  }

  if(!is.na(config$run_config$forecast_start_datetime)){

    forecast_date <- lubridate::as_date(config$run_config$forecast_start_datetime)
    forecast_hour <- lubridate::hour(config$run_config$forecast_start_datetime)

    if(forecast_hour != 0){
      stop("Only forecasts that start at 00:00:00 UTC are currently supported")
    }

    if(use_s3){
      if(is.null(bucket) | is.null(endpoint)){
        stop("inflow forecast function needs bucket and endpoint if use_s3=TRUE")
      }
      vars <- FLAREr:::arrow_env_vars()
      forecast_dir <- arrow::s3_bucket(bucket = file.path(bucket, "stage2/parquet", forecast_hour,forecast_date),
                                       endpoint_override =  endpoint, anonymous = TRUE)
      past_dir <- arrow::s3_bucket(bucket = file.path(bucket, "stage3/parquet", config$location$site_id),
                                   endpoint_override =  endpoint, anonymous = TRUE)
      FLAREr:::unset_arrow_vars(vars)
    }else{
      if(is.null(local_directory)){
        stop("inflow forecast function needs local_directory if use_s3=FALSE")
      }
      forecast_dir <- arrow::SubTreeFileSystem$create(file.path(local_directory, "stage2/parquet", forecast_hour,forecast_date))
      past_dir <-  arrow::SubTreeFileSystem$create(file.path(local_directory, "stage3/parquet", config$location$site_id))
    }
  }

  if(is.null(obs_met_file) & is.null(forecast_dir)){
    stop("missing files to convert")
  }



  full_time <- seq(start_datetime, end_datetime, by = "1 hour")
  if(config$met$use_forecasted_met){
    if(forecast_start_datetime > start_datetime){
      full_time_hist <- seq(start_datetime, forecast_start_datetime, by = "1 hour")
    }else{
      full_time_hist <- NULL
    }
  }else{
    full_time_hist <- seq(start_datetime, end_datetime, by = "1 hour")
  }

  if(!is.null(obs_met_file)){

    target <- readr::read_csv(obs_met_file, show_col_types = FALSE) |>
      tidyr::pivot_wider(names_from = variable, values_from = observation) |>
      dplyr::arrange(datetime) |>
      dplyr::mutate(WindSpeed = wind_speed) |>
      dplyr::rename(AirTemp = air_temperature,
                    ShortWave = surface_downwelling_shortwave_flux_in_air,
                    LongWave = surface_downwelling_longwave_flux_in_air,
                    RelHum = relative_humidity,
                    Rain = precipitation_flux,
                    time = datetime) |>
      dplyr::mutate(AirTemp = AirTemp - 273.15,
                    RelHum = RelHum * 100,
                    Rain = Rain * (60 * 60), #This is different from GLM
                    Snow = 0.0) |>
      dplyr::mutate_at(dplyr::vars(all_of(c("AirTemp", "ShortWave","LongWave","RelHum","WindSpeed"))), list(~round(., 2))) |>
      dplyr::mutate(Rain = round(Rain, 5),
                    time = strftime(time, format="%Y-%m-%d %H:%M", tz = "UTC")) |>
      dplyr::select(ensemble, time, AirTemp,ShortWave, LongWave, RelHum, WindSpeed,Rain) |>
      dplyr::rename(Air_Temperature_celsius = AirTemp,
                    Shortwave_Radiation_Downwelling_wattPerMeterSquared = ShortWave,
                    Longwave_Radiation_Downwelling_wattPerMeterSquared = LongWave,
                    Relative_Humidity_percent = RelHum,
                    Ten_Meter_Elevation_Wind_Speed_meterPerSecond = WindSpeed,
                    Precipitation_millimeterPerHour = Rain,
                    Snowfall_millimeterPerHour = Snow)
  }else{
    target <- arrow::open_dataset(past_dir) |>
      dplyr::filter(site_id == config$location$site_id) |>
      dplyr::select(datetime, parameter,variable,prediction) |>
      dplyr::collect() |>
      dplyr::filter(datetime %in% full_time_hist) |>
      tidyr::pivot_wider(names_from = variable, values_from = prediction) |>
      dplyr::arrange(parameter, datetime) |>
      dplyr::mutate(WindSpeed = sqrt(eastward_wind^2 + northward_wind^2)) |>
      dplyr::rename(AirTemp = air_temperature,
                    ShortWave = surface_downwelling_shortwave_flux_in_air,
                    LongWave = surface_downwelling_longwave_flux_in_air,
                    RelHum = relative_humidity,
                    Rain = precipitation_flux,
                    ensemble = parameter,
                    time = datetime) |>
      dplyr::mutate(AirTemp = AirTemp - 273.15,
                    RelHum = RelHum * 100,
                    Rain = Rain * (60 * 60), #This is different from GLM
                    Snow = 0.0) |>
      dplyr::mutate_at(dplyr::vars(all_of(c("AirTemp", "ShortWave","LongWave","RelHum","WindSpeed"))), list(~round(., 2))) |>
      dplyr::mutate(Rain = round(Rain, 5),
                    time = strftime(time, format="%Y-%m-%d %H:%M", tz = "UTC")) |>
      dplyr::select(ensemble, time, AirTemp,ShortWave, LongWave, RelHum, WindSpeed,Rain) |>
      dplyr::rename(Air_Temperature_celsius = AirTemp,
                    Shortwave_Radiation_Downwelling_wattPerMeterSquared = ShortWave,
                    Longwave_Radiation_Downwelling_wattPerMeterSquared = LongWave,
                    Relative_Humidity_percent = RelHum,
                    Ten_Meter_Elevation_Wind_Speed_meterPerSecond = WindSpeed,
                    Precipitation_millimeterPerHour = Rain,
                    Snowfall_millimeterPerHour = Snow) |>
      dplyr::group_by(ensemble) |>
      dplyr::slice(-dplyr::n()) |>
      dplyr::ungroup()
  }


  if(is.null(forecast_dir)){

    if(!is.null(obs_met_file)){

      current_filename <- "met.csv"
      current_filename <- file.path(out_dir, current_filename)


      write.csv(target, file = current_filename, quote = FALSE, row.names = FALSE)
    }else{

      ensemble_members <- unique(target$ensemble)

      current_filename <- purrr::map_chr(ensemble_members, function(ens, out_dir, target){
        df <- target |>
          dplyr::filter(ensemble == ens) |>
          dplyr::select(-ensemble) |>
          dplyr::arrange(time)

        if( any(df$Relative_Humidity_percent <= 0.0)) {
          idx <- which(df$Relative_Humidity_percent <= 0.0)
          df$Relative_Humidity_percent[idx] <- NA
          df$Relative_Humidity_percent <- zoo::na.approx(df$Relative_Humidity_percent, rule = 2)
        }

        df$Relative_Humidity_percent[df$Relative_Humidity_percent > 100] <- 100

        fn <- paste0("met_",stringr::str_pad(ens, width = 2, side = "left", pad = "0"),".csv")
        fn <- file.path(out_dir, fn)
        write.csv(df, file = fn, quote = FALSE, row.names = FALSE)
        return(fn)
      },
      out_dir = out_dir,
      forecast,
      target)
    }
  }else{

    forecast <- arrow::open_dataset(forecast_dir) |>
      dplyr::filter(site_id == config$location$site_id) |>
      dplyr::select(datetime, parameter,variable,prediction) |>
      dplyr::collect() |>
      tidyr::pivot_wider(names_from = variable, values_from = prediction) |>
      dplyr::arrange(parameter, datetime) |>
      dplyr::mutate(WindSpeed = sqrt(eastward_wind^2 + northward_wind^2)) |>
      dplyr::rename(AirTemp = air_temperature,
                    ShortWave = surface_downwelling_shortwave_flux_in_air,
                    LongWave = surface_downwelling_longwave_flux_in_air,
                    RelHum = relative_humidity,
                    Rain = precipitation_flux,
                    ensemble = parameter,
                    time = datetime) |>
      dplyr::mutate(AirTemp = AirTemp - 273.15,
                    RelHum = RelHum * 100,
                    Rain = Rain * (60 * 60), #This is different from GLM
                    Snow = 0.0) |>
      dplyr::mutate_at(dplyr::vars(all_of(c("AirTemp", "ShortWave","LongWave","RelHum","WindSpeed"))), list(~round(., 2))) |>
      dplyr::mutate(Rain = round(Rain, 5),
                    time = strftime(time, format="%Y-%m-%d %H:%M", tz = "UTC")) |>
      dplyr::select(ensemble, time, AirTemp,ShortWave, LongWave, RelHum, WindSpeed,Rain) |>
      dplyr::rename(Air_Temperature_celsius = AirTemp,
                    Shortwave_Radiation_Downwelling_wattPerMeterSquared = ShortWave,
                    Longwave_Radiation_Downwelling_wattPerMeterSquared = LongWave,
                    Relative_Humidity_percent = RelHum,
                    Ten_Meter_Elevation_Wind_Speed_meterPerSecond = WindSpeed,
                    Precipitation_millimeterPerHour = Rain,
                    Snowfall_millimeterPerHour = Snow) |>
      dplyr::group_by(ensemble) |>
      dplyr::slice(-dplyr::n()) |>
      dplyr::ungroup()

    ensemble_members <- unique(forecast$ensemble)

    current_filename <- purrr::map_chr(ensemble_members, function(ens, out_dir, forecast, target){

      if("ensemble" %in% names(target)){
        target <- target |>
          dplyr::filter(ensemble == ens) |>
          dplyr::select(-ensemble)
      }


      df <- forecast |>
        dplyr::filter(ensemble == ens) |>
        dplyr::select(-ensemble) |>
        dplyr::bind_rows(target) |>
        dplyr::arrange(time)

      if( any(df$Relative_Humidity_percent <= 0.0)) {
        idx <- which(df$Relative_Humidity_percent <= 0.0)
        df$Relative_Humidity_percent[idx] <- NA
        df$Relative_Humidity_percent <- zoo::na.approx(df$Relative_Humidity_percent, rule = 2)
      }

      df$Relative_Humidity_percent[df$Relative_Humidity_percent > 100] <- 100

      fn <- paste0("met_",stringr::str_pad(ens, width = 2, side = "left", pad = "0"),".csv")
      fn <- file.path(out_dir, fn)
      write.csv(df, file = fn, quote = FALSE, row.names = FALSE)
      return(fn)
    },
    out_dir = out_dir,
    forecast,
    target)
  }
  return(list(filenames = current_filename, historical_met_error = FALSE))
}
