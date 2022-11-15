##' @title Generate csv output file from arrow
##' @details Function generates a netcdf file from the object that is returned by run_da_forecast()
##' @param da_forecast_output list; object that is returned by run_da_forecast()
##' @param forecast_output_directory string; full path of directory where the csv file will be written
##' @return None
##' @export
##' @importFrom lubridate with_tz
##' @author Quinn Thomas
##' @examples
##' \dontrun{
##' write_forecast_csv(da_forecast_output = da_forecast_output, forecast_output_directory = config$file_path$forecast_output_directory, use_short_filename = TRUE)
##' }
##'

write_forecast_arrow_ler <- function(da_forecast_output,
                                     use_s3 = FALSE,
                                     bucket = NULL,
                                     endpoint = NULL,
                                     local_directory = NULL){


  if(use_s3){
    if(is.null(bucket) | is.null(endpoint)){
      stop("scoring function needs bucket and endpoint if use_s3=TRUE")
    }
    vars <- FLAREr:::arrow_env_vars()
    output_directory <- arrow::s3_bucket(bucket = bucket,
                                         endpoint_override =  endpoint)
    FLAREr:::unset_arrow_vars(vars)
  }else{
    if(is.null(local_directory)){
      stop("scoring function needs local_directory if use_s3=FALSE")
    }
    output_directory <- arrow::SubTreeFileSystem$create(local_directory)
  }

  x <- da_forecast_output$x
  pars <- da_forecast_output$pars
  lake_depth <- da_forecast_output$lake_depth
  snow_ice_thickness <- da_forecast_output$snow_ice_thickness
  data_assimilation_flag <- da_forecast_output$data_assimilation_flag
  forecast_flag <- da_forecast_output$forecast_flag
  full_time <- da_forecast_output$full_time
  forecast_start_datetime <- da_forecast_output$forecast_start_datetime
  config <- da_forecast_output$config
  states_config <- da_forecast_output$states_config
  obs_config <- da_forecast_output$obs_config
  pars_config <- da_forecast_output$pars_config
  diagnostics <- da_forecast_output$diagnostics

  forecast_flag[which(is.na(forecast_flag))] <- 0

  output_list <- NULL
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      for(k in 1:dim(x)[3]){
        tmp <- tibble(predicted = x[i, j, k, ],
                      time  = full_time[i],
                      depth = config$model_settings$modeled_depths[k],
                      variable = states_config$state_names[j],
                      forecast = forecast_flag[i],
                      ensemble = 1:dim(x)[4],
                      variable_type = "state")
        output_list <- rbind(output_list, tmp)
      }
    }
  }

  if(length(config$output_settings$diagnostics_names) > 0){
    for(i in 1:dim(diagnostics)[1]){
      for(j in 1:dim(diagnostics)[2]){
        for(k in 1:dim(diagnostics)[3]){
          tmp <- tibble(predicted = diagnostics[i, j, k, ],
                        time  = full_time[j],
                        variable = config$output_settings$diagnostics_names[i],
                        depth = config$model_settings$modeled_depths[k],
                        forecast = forecast_flag[j],
                        ensemble = 1:dim(diagnostics)[4],
                        variable_type = "diagnostic")
          output_list <- rbind(output_list, tmp)
        }
      }
    }
  }

  for(i in 1:dim(pars)[1]){
    for(j in 1:dim(pars)[2]){
      tmp <- tibble(predicted = pars[i, j, ],
                    time  = full_time[i],
                    variable = pars_config$par_names_save[j],
                    depth = NA,
                    forecast = forecast_flag[i],
                    ensemble = 1:dim(pars)[3],
                    variable_type = "parameter")
      output_list <- rbind(output_list, tmp)
    }
  }

  if(!is.null(da_forecast_output$restart_list)){
    lake_depth <- da_forecast_output$restart_list$lake_depth
  }
  for(i in 1:dim(lake_depth)[1]){
    tmp <- tibble(predicted = lake_depth[i, ],
                  time  = full_time[i],
                  variable = "depth",
                  depth = NA,
                  forecast = forecast_flag[i],
                  ensemble = 1:dim(lake_depth)[2],
                  variable_type = "state")
    output_list <- rbind(output_list, tmp)
  }


  if(length(config$output_settings$diagnostics_names) > 0){
    for(i in 1:dim(diagnostics)[2]){
      tmp <- tibble(predicted = 1.7 / diagnostics[1, i, which.min(abs(config$model_settings$modeled_depths-1.0)), ],
                    time = full_time[i],
                    variable = "secchi",
                    depth = NA,
                    forecast = forecast_flag[i],
                    ensemble = 1:dim(diagnostics)[4],
                    variable_type = "state")
      output_list <- rbind(output_list, tmp)
    }

  }

  # for(i in 1:dim(snow_ice_thickness)[1]){
  #   tmp <- tibble(predicted = apply(snow_ice_thickness[2:3, i, ], 2, sum),
  #                 time = full_time[i],
  #                 variable = "ice_thickness",
  #                 depth = NA,
  #                 forecast = forecast_flag[i],
  #                 ensemble = 1:dim(snow_ice_thickness)[3],
  #                 variable_type = "state")
  #   output_list <- rbind(output_list, tmp)
  # }

  time_of_forecast <- lubridate::with_tz(da_forecast_output$time_of_forecast, tzone = "UTC")

  output_list <- output_list %>%
    dplyr::mutate(pub_time = time_of_forecast,
                  reference_datetime = forecast_start_datetime,
                  site_id = config$location$site_id,
                  model_id = config$run_config$sim_name,
                  family = "ensemble") %>%
    rename(datetime = time,
           parameter = ensemble,
           prediction = predicted) %>%
    dplyr::select(reference_datetime, datetime, pub_time, model_id, site_id, depth, family, parameter, variable, prediction, forecast, variable_type)

  for(i in 1:length(states_config$state_names)){
    if(length(which(obs_config$state_names_obs == states_config$state_names[i])) >0){
      obs_name <- obs_config$target_variable[which(obs_config$state_names_obs == states_config$state_names[i])]
      output_list <- output_list %>%
        mutate(variable = ifelse(variable == states_config$state_names[i], obs_name, variable))
    }
  }

  reference_datetime_format <- "%Y-%m-%d %H:%M:%S"

  output_list <- output_list |> mutate(reference_datetime = strftime(lubridate::as_datetime(reference_datetime),format=reference_datetime_format,tz = "UTC"))
  arrow::write_dataset(dataset = output_list,
                       path = output_directory,
                       partitioning = c("site_id","model_id","reference_datetime"))

  return(output_list)

}
