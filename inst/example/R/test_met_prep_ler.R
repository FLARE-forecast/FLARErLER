
template_folder <- system.file("example", package = "FLARErLER")
temp_dir <- tempdir()
# dir.create("example")
file.copy(from = template_folder, to = temp_dir, recursive = TRUE)

test_directory <- file.path(temp_dir, "example")

lake_directory <- test_directory

configure_run_file <- "configure_run.yml"
config_set_name <- "ler"

config <- FLARErLER::set_configuration(configure_run_file,
                                    lake_directory,
                                    config_set_name = config_set_name)
config$run_config$forecast_start_datetime <- "2018-10-05 00:00:00"
config$run_config$forecast_horizon <- 1

config$da_setup$ensemble_size <- 21
config$model_settings$ncore <- 2
config$model_settings$use_ler <- TRUE

noaa_forecast_path <- FLARErLER::get_driver_forecast_path(config,
                                                       forecast_model = config$met$forecast_met_model)

inflow_forecast_path <- FLARErLER::get_driver_forecast_path(config,
                                                         forecast_model = config$inflow$forecast_inflow_model)

if(!is.null(noaa_forecast_path)){
  # FLARErLER::get_driver_forecast(lake_directory, forecast_path = noaa_forecast_path)
  forecast_dir <- file.path(config$file_path$noaa_directory, noaa_forecast_path)
}else{
  forecast_dir <- NULL
}

if(!is.null(inflow_forecast_path)){
  # FLARErLER::get_driver_forecast(lake_directory, forecast_path = inflow_forecast_path)
  inflow_file_dir <- file.path(config$file_path$noaa_directory,inflow_forecast_path)
}else{
  inflow_file_dir <- NULL
}


# file.copy(file.path(configuration_directory, "glm3.nml"), execute_directory)

pars_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$par_config_file), col_types = readr::cols())
obs_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$obs_config_file), col_types = readr::cols())
states_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$states_config_file), col_types = readr::cols())

#Download and process observations (already done)

cleaned_observations_file_long <- file.path(config$file_path$qaqc_data_directory,"observations_postQAQC_long.csv")
cleaned_inflow_file <- file.path(config$file_path$qaqc_data_directory, "/inflow_postQAQC.csv")
observed_met_file <- file.path(config$file_path$qaqc_data_directory,
                               paste0("observed-met_",config$location$site_id,".nc"))
