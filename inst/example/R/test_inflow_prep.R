
template_folder <- system.file("example", package = "FLARErLER")
temp_dir <- tempdir()
# dir.create("example")
file.copy(from = template_folder, to = temp_dir, recursive = TRUE)

test_directory <- file.path(temp_dir, "example")

lake_directory <- test_directory

configure_run_file <- "configure_run.yml"
config_set_name <- "ler"

config <- FLAREr::set_configuration(configure_run_file,
                                    lake_directory,
                                    config_set_name = config_set_name)
config$run_config$forecast_start_datetime <- "2018-10-05 00:00:00"
config$run_config$forecast_horizon <- 1

config$da_setup$ensemble_size <- 21
config$model_settings$ncore <- 2
config$model_settings$use_ler <- TRUE

noaa_forecast_path <- FLAREr::get_driver_forecast_path(config,
                                                       forecast_model = config$met$forecast_met_model)

inflow_forecast_path <- FLAREr::get_driver_forecast_path(config,
                                                         forecast_model = config$inflow$forecast_inflow_model)

if(!is.null(noaa_forecast_path)){
  # FLAREr::get_driver_forecast(lake_directory, forecast_path = noaa_forecast_path)
  forecast_dir <- file.path(config$file_path$noaa_directory, noaa_forecast_path)
}else{
  forecast_dir <- NULL
}

if(!is.null(inflow_forecast_path)){
  # FLAREr::get_driver_forecast(lake_directory, forecast_path = inflow_forecast_path)
  inflow_file_dir <- file.path(config$file_path$noaa_directory,inflow_forecast_path)
}else{
  inflow_file_dir <- NULL
}

cleaned_inflow_file <- file.path(config$file_path$qaqc_data_directory, "/inflow_postQAQC.csv")


