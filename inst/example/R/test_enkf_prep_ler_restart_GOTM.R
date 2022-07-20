
template_folder <- system.file("example", package = "FLARErLER")
temp_dir <- tempdir()
# dir.create("example")
file.copy(from = template_folder, to = temp_dir, recursive = TRUE)

test_directory <- file.path(temp_dir, "example")

lake_directory <- test_directory

configure_run_file <- "configure_run.yml"
config_set_name <- "ler"
model <- "GOTM"

config <- FLAREr::set_configuration(configure_run_file,
                                    lake_directory,
                                    config_set_name = config_set_name)
config$run_config$forecast_start_datetime <- "2018-10-05 00:00:00"
config$run_config$forecast_horizon <- 1

config$da_setup$ensemble_size <- 21
config$model_settings$ncore <- 2
config$model_settings$use_ler <- TRUE
config$model_settings$model <- model

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


# file.copy(file.path(configuration_directory, "glm3.nml"), execute_directory)

pars_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$par_config_file), col_types = readr::cols())
obs_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$obs_config_file), col_types = readr::cols())
states_config <- readr::read_csv(file.path(config$file_path$configuration_directory, config$model_settings$states_config_file), col_types = readr::cols())

#Download and process observations (already done)

cleaned_observations_file_long <- file.path(config$file_path$qaqc_data_directory,"observations_postQAQC_long.csv")
cleaned_inflow_file <- file.path(config$file_path$qaqc_data_directory, "/inflow_postQAQC.csv")
observed_met_file <- file.path(config$file_path$qaqc_data_directory,
                               paste0("observed-met_",config$location$site_id,".nc"))

#Step up Drivers
met_out <- FLARErLER::generate_met_files_ler(obs_met_file = observed_met_file,
                                             out_dir = config$file_path$execute_directory,
                                             forecast_dir = forecast_dir,
                                             config = config)
met_file_names <- met_out$filenames

historical_met_error <- met_out$historical_met_error

inflow_forecast_path <- config$file_path$inflow_directory

inflow_outflow_files <- FLARErLER::create_inflow_outflow_files_ler(inflow_file_dir = inflow_file_dir,
                                                                   inflow_obs = cleaned_inflow_file,
                                                                   working_directory = config$file_path$execute_directory,
                                                                   config,
                                                                   state_names = NULL)

inflow_file_names <- inflow_outflow_files$inflow_file_name
outflow_file_names <- inflow_outflow_files$outflow_file_name

obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long,
                                 obs_config,
                                 config)

states_config <- FLAREr::generate_states_to_obs_mapping(states_config, obs_config)

model_sd <- FLAREr::initiate_model_error(config, states_config)
init <- FLARErLER::generate_initial_conditions_ler(states_config,
                                                   obs_config,
                                                   pars_config,
                                                   obs,
                                                   config)
