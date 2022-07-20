test_that("LER met files are generated", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R/test_met_prep_ler.R"))
  config$model_settings$use_ler <- TRUE
  config$model_settings$model <- "GLM"

  met_out <- FLARErLER::generate_met_files_ler(obs_met_file = observed_met_file,
                                       out_dir = config$file_path$execute_directory,
                                       forecast_dir = forecast_dir,
                                       config)

  met_file_names <- met_out$filenames
  testthat::expect_equal(file.exists(met_file_names), expected = rep(TRUE, 21))
})


test_that("LER inflow & outflow files are generated", {

  template_folder <- system.file("example", package = "FLARErLER")
  source(file.path(template_folder, "R/test_inflow_prep.R"))

  config$model_settings$use_ler <- TRUE
  config$model_settings$model <- "GLM"

  inflow_outflow_files <- FLARErLER::create_inflow_outflow_files_ler(inflow_file_dir = inflow_file_dir,
                                                              inflow_obs = cleaned_inflow_file,
                                                              working_directory = config$file_path$execute_directory,
                                                              config,
                                                              state_names = NULL)

  inflow_file_names <- inflow_outflow_files$inflow_file_name
  outflow_file_names <- inflow_outflow_files$outflow_file_name

  testthat::expect_equal(file.exists(inflow_outflow_files[[1]]), expected = rep(TRUE, 21))
  testthat::expect_equal(file.exists(inflow_outflow_files[[2]]), expected = rep(TRUE, 21))
})

#Set initial conditions
test_that("LER-GLM initial conditions are generated", {

  template_folder <- system.file("example", package = "FLARErLER")
  source(file.path(template_folder, "R/test_met_prep_ler.R"))

  config$model_settings$model <- "GLM"

  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  init <- FLARErLER::generate_initial_conditions_ler(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config)

  testthat::expect_true(is.list(init))
  chk <- lapply(init, is.array)
  testthat::expect_true(any(unlist(chk)))
})

#Set initial conditions
test_that("LER-GOTM initial conditions are generated", {

  template_folder <- system.file("example", package = "FLARErLER")
  source(file.path(template_folder, "R/test_met_prep_ler.R"))

  config$model_settings$model <- "GOTM"

  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  init <- FLARErLER::generate_initial_conditions_ler(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config)

  testthat::expect_true(is.list(init))
  chk <- lapply(init, is.array)
  testthat::expect_true(any(unlist(chk)))
})

# LER-GLM-EnKF Tests ----
test_that("LER-GLM-EnKF can be run", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler.R"))
  config$model_settings$model <- "GLM"

  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  init <- FLARErLER::generate_initial_conditions_ler(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config)

  # setwd(config$file_path$execute_directory)


  # states_init = init$states
  # pars_init = init$pars
  # aux_states_init = init$aux_states_init
  # obs = obs
  # obs_sd = obs_config$obs_sd
  # model_sd = model_sd
  # working_directory = config$file_path$execute_directory
  # met_file_names = (met_file_names)
  # inflow_file_names = (inflow_file_names)
  # outflow_file_names = (outflow_file_names)
  # config = config
  # pars_config = pars_config
  # states_config = states_config
  # obs_config = obs_config
  # management = NULL
  # da_method = "enkf"
  # par_fit_method = "inflate"
  # debug = TRUE

  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                         pars_init = init$pars,
                                         aux_states_init = init$aux_states_init,
                                         obs = obs,
                                         obs_sd = obs_config$obs_sd,
                                         model_sd = model_sd,
                                         working_directory = config$file_path$execute_directory,
                                         met_file_names = met_file_names,
                                         inflow_file_names = inflow_file_names,
                                         outflow_file_names = outflow_file_names,
                                         config = config,
                                         pars_config = pars_config,
                                         states_config = states_config,
                                         obs_config = obs_config,
                                         management = NULL,
                                         da_method = "enkf",
                                         par_fit_method = "inflate",
                                         debug = TRUE
  )

  # saveRDS(object = enkf_output, file = "inst/example/benchmark_data/sampenkf_output_GLM.RDS")
  samp_enkf_output <- readRDS(file.path(template_folder, "benchmark_data", "sampenkf_output_GLM.RDS"))

  testthat::expect_true(is.list(enkf_output))
  chk <- lapply(1:length(enkf_output), function(x) {
    class(enkf_output[[x]]) == class(samp_enkf_output[[x]])
  })

  testthat::expect_true(all(unlist(chk)))

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                              forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)
})

# LER-GLM-EnKF Tests ----
test_that("LER-GLM-EnKF can be run with no inflows", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler.R"))
  config$model_settings$model <- "GLM"

  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long,
                                      obs_config,
                                      config)

  init <- FLARErLER::generate_initial_conditions_ler(states_config,
                                                 obs_config,
                                                 pars_config,
                                                 obs,
                                                 config)

  # setwd(config$file_path$execute_directory)


  # states_init = init$states
  # pars_init = init$pars
  # aux_states_init = init$aux_states_init
  # obs = obs
  # obs_sd = obs_config$obs_sd
  # model_sd = model_sd
  # working_directory = config$file_path$execute_directory
  # met_file_names = (met_file_names)
  # inflow_file_names = NULL
  # outflow_file_names = NULL
  # config = config
  # pars_config = pars_config
  # states_config = states_config
  # obs_config = obs_config
  # management = NULL
  # da_method = "enkf"
  # par_fit_method = "inflate"
  # debug = TRUE

  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                                pars_init = init$pars,
                                                aux_states_init = init$aux_states_init,
                                                obs = obs,
                                                obs_sd = obs_config$obs_sd,
                                                model_sd = model_sd,
                                                working_directory = config$file_path$execute_directory,
                                                met_file_names = met_file_names,
                                                inflow_file_names = NULL,
                                                outflow_file_names = NULL,
                                                config = config,
                                                pars_config = pars_config,
                                                states_config = states_config,
                                                obs_config = obs_config,
                                                management = NULL,
                                                da_method = "enkf",
                                                par_fit_method = "inflate",
                                                debug = TRUE
  )

  samp_enkf_output <- readRDS(file.path(template_folder, "benchmark_data", "sampenkf_output_GLM.RDS"))

  testthat::expect_true(is.list(enkf_output))
  chk <- lapply(1:length(enkf_output), function(x) {
    class(enkf_output[[x]]) == class(samp_enkf_output[[x]])
  })

  testthat::expect_true(all(unlist(chk)))

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                                 forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                   da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                              qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)
})

# LER-GOTM-EnKF Tests ----
test_that("LER-GOTM-EnKF can be run", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler.R"))
  config$model_settings$model <- "GOTM"
  config$output_settings$diagnostics_names <- NULL

  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  init <- FLARErLER::generate_initial_conditions_ler(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config)

  # setwd(config$file_path$execute_directory)


  # states_init = init$states
  # pars_init = init$pars
  # aux_states_init = init$aux_states_init
  # obs = obs
  # obs_sd = obs_config$obs_sd
  # model_sd = model_sd
  # working_directory = config$file_path$execute_directory
  # met_file_names = (met_file_names)
  # inflow_file_names = (inflow_file_names)
  # outflow_file_names = (outflow_file_names)
  # config = config
  # pars_config = pars_config
  # states_config = states_config
  # obs_config = obs_config
  # management = NULL
  # da_method = "enkf"
  # par_fit_method = "inflate"
  # debug = FALSE

  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                             pars_init = init$pars,
                                             aux_states_init = init$aux_states_init,
                                             obs = obs,
                                             obs_sd = obs_config$obs_sd,
                                             model_sd = model_sd,
                                             working_directory = config$file_path$execute_directory,
                                             met_file_names = met_file_names,
                                             inflow_file_names = inflow_file_names,
                                             outflow_file_names = outflow_file_names,
                                             config = config,
                                             pars_config = pars_config,
                                             states_config = states_config,
                                             obs_config = obs_config,
                                             management = NULL,
                                             da_method = "enkf",
                                             par_fit_method = "inflate",
                                             debug = TRUE
  )

  # saveRDS(object = enkf_output, file = "inst/example/benchmark_data/sampenkf_output_GOTM.RDS")
  samp_enkf_output <- readRDS(file.path(template_folder, "benchmark_data", "sampenkf_output_GOTM.RDS"))

  testthat::expect_true(is.list(enkf_output))
  chk <- lapply(1:length(enkf_output), function(x) {
    class(enkf_output[[x]]) == class(samp_enkf_output[[x]])
  })

  testthat::expect_true(all(unlist(chk)))

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                                     forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)
})

# LER-Simstrat-EnKF Tests ----
test_that("LER-Simstrat-EnKF can be run", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler.R"))
  config$model_settings$model <- "Simstrat"
  config$output_settings$diagnostics_names <- NULL

  obs <- FLAREr::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  init <- FLARErLER::generate_initial_conditions_ler(states_config,
                                                     obs_config,
                                                     pars_config,
                                                     obs,
                                                     config)

  # setwd(config$file_path$execute_directory)


  # states_init = init$states
  # pars_init = init$pars
  # aux_states_init = init$aux_states_init
  # obs = obs
  # obs_sd = obs_config$obs_sd
  # model_sd = model_sd
  # working_directory = config$file_path$execute_directory
  # met_file_names = (met_file_names)
  # inflow_file_names = (inflow_file_names)
  # outflow_file_names = (outflow_file_names)
  # config = config
  # pars_config = pars_config
  # states_config = states_config
  # obs_config = obs_config
  # management = NULL
  # da_method = "enkf"
  # par_fit_method = "inflate"
  # debug = TRUE

  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                             pars_init = init$pars,
                                             aux_states_init = init$aux_states_init,
                                             obs = obs,
                                             obs_sd = obs_config$obs_sd,
                                             model_sd = model_sd,
                                             working_directory = config$file_path$execute_directory,
                                             met_file_names = met_file_names,
                                             inflow_file_names = inflow_file_names,
                                             outflow_file_names = outflow_file_names,
                                             config = config,
                                             pars_config = pars_config,
                                             states_config = states_config,
                                             obs_config = obs_config,
                                             management = NULL,
                                             da_method = "enkf",
                                             par_fit_method = "inflate",
                                             debug = TRUE
  )

  # saveRDS(object = enkf_output, file = "inst/example/benchmark_data/sampenkf_output_Simstrat.RDS")
  samp_enkf_output <- readRDS(file.path(template_folder, "benchmark_data", "sampenkf_output_Simstrat.RDS"))

  testthat::expect_true(is.list(enkf_output))
  chk <- lapply(1:length(enkf_output), function(x) {
    class(enkf_output[[x]]) == class(samp_enkf_output[[x]])
  })

  testthat::expect_true(all(unlist(chk)))

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                                     forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)
})

# LER-GLM-PF Tests ----
test_that("LER-GLM-PF can be run", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler.R"))
  config$model_settings$model <- "GLM"

  obs <- FLARErLER::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  init <- generate_initial_conditions(states_config,
                                      obs_config,
                                      pars_config,
                                      obs,
                                      config)


  # states_init = init$states
  # pars_init = init$pars
  # aux_states_init = init$aux_states_init
  # obs = obs
  # obs_sd = obs_config$obs_sd
  # model_sd = model_sd
  # working_directory = config$file_path$execute_directory
  # met_file_names = (met_file_names)
  # inflow_file_names = (inflow_file_names)
  # outflow_file_names = (outflow_file_names)
  # config = config
  # pars_config = pars_config
  # states_config = states_config
  # obs_config = obs_config
  # management = NULL
  # da_method = "pf"
  # par_fit_method = "perturb"
  # debug = FALSE

  #Run EnKF
  enkf_output <- run_da_forecast_ler(states_init = init$states,
                                     pars_init = init$pars,
                                     aux_states_init = init$aux_states_init,
                                     obs = obs,
                                     obs_sd = obs_config$obs_sd,
                                     model_sd = model_sd,
                                     working_directory = config$file_path$execute_directory,
                                     met_file_names = met_file_names,
                                     inflow_file_names = inflow_file_names,
                                     outflow_file_names = outflow_file_names,
                                     config = config,
                                     pars_config = pars_config,
                                     states_config = states_config,
                                     obs_config = obs_config,
                                     management = NULL,
                                     da_method = "pf",
                                     par_fit_method = "perturb",
                                     debug = FALSE
  )

  samp_enkf_output <- readRDS(file.path(template_folder, "benchmark_data", "sampenkf_output_GLM.RDS"))

  testthat::expect_true(is.list(enkf_output))
  chk <- lapply(1:length(enkf_output), function(x) {
    class(enkf_output[[x]]) == class(samp_enkf_output[[x]])
  })

  testthat::expect_true(all(unlist(chk)))

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                                     forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)
})

# LER-GLM can be restarted ----
test_that("LER-GLM can be restarted", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler_restart_GLM.R"))
  config$output_settings$diagnostics_names <- NULL


  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                             pars_init = init$pars,
                                             aux_states_init = init$aux_states_init,
                                             obs = obs,
                                             obs_sd = obs_config$obs_sd,
                                             model_sd = model_sd,
                                             working_directory = config$file_path$execute_directory,
                                             met_file_names = met_file_names,
                                             inflow_file_names = inflow_file_names,
                                             outflow_file_names = outflow_file_names,
                                             config = config,
                                             pars_config = pars_config,
                                             states_config = states_config,
                                             obs_config = obs_config,
                                             management = NULL,
                                             da_method = "enkf",
                                             par_fit_method = "inflate",
                                             debug = FALSE
  )

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf(da_forecast_output = enkf_output,
                                              forecast_output_directory = config$file_path$forecast_output_directory)
  FLARErLER::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  FLARErLER::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)


  # Begin restart ----

  config <- FLARErLER::update_run_config(config = config, lake_directory =  lake_directory, configure_run_file =  configure_run_file, saved_file = saved_file, new_horizon = 2)

  config$run_config
  # Rename restart dir
  config$run_config$sim_name <- "test-restart"
  config$file_path$execute_directory <- file.path(lake_directory, "flare_tempdir",
                                                  config$location$site_id,
                                                  config$run_config$sim_name)
  dir.create(config$file_path$execute_directory, recursive = TRUE)

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

  #Step up Drivers
  met_out <- FLARErLER::generate_met_files(obs_met_file = observed_met_file,
                                        out_dir = config$file_path$execute_directory,
                                        forecast_dir = forecast_dir,
                                        config = config)
  met_file_names <- met_out$filenames

  historical_met_error <- met_out$historical_met_error

  inflow_forecast_path <- config$file_path$inflow_directory

  inflow_outflow_files <- FLARErLER::create_inflow_outflow_files(inflow_file_dir = inflow_file_dir,
                                                              inflow_obs = cleaned_inflow_file,
                                                              working_directory = config$file_path$execute_directory,
                                                              config,
                                                              state_names = NULL)

  inflow_file_names <- inflow_outflow_files$inflow_file_name
  outflow_file_names <- inflow_outflow_files$outflow_file_name

  obs <- FLARErLER::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  states_config <- FLARErLER::generate_states_to_obs_mapping(states_config, obs_config)

  model_sd <- FLARErLER::initiate_model_error(config, states_config)


  init <- FLARErLER::generate_initial_conditions(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config)


  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                             pars_init = init$pars,
                                             aux_states_init = init$aux_states_init,
                                             obs = obs,
                                             obs_sd = obs_config$obs_sd,
                                             model_sd = model_sd,
                                             working_directory = config$file_path$execute_directory,
                                             met_file_names = met_file_names,
                                             inflow_file_names = inflow_file_names,
                                             outflow_file_names = outflow_file_names,
                                             config = config,
                                             pars_config = pars_config,
                                             states_config = states_config,
                                             obs_config = obs_config,
                                             management = NULL,
                                             da_method = "enkf",
                                             par_fit_method = "inflate",
                                             debug = FALSE
  )

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                                     forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)
})

# LER-GOTM can be restarted ----
test_that("LER-GOTM can be restarted", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler_restart_GOTM.R"))
  config$output_settings$diagnostics_names <- NULL


  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                             pars_init = init$pars,
                                             aux_states_init = init$aux_states_init,
                                             obs = obs,
                                             obs_sd = obs_config$obs_sd,
                                             model_sd = model_sd,
                                             working_directory = config$file_path$execute_directory,
                                             met_file_names = met_file_names,
                                             inflow_file_names = inflow_file_names,
                                             outflow_file_names = outflow_file_names,
                                             config = config,
                                             pars_config = pars_config,
                                             states_config = states_config,
                                             obs_config = obs_config,
                                             management = NULL,
                                             da_method = "enkf",
                                             par_fit_method = "inflate",
                                             debug = FALSE
  )

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf(da_forecast_output = enkf_output,
                                              forecast_output_directory = config$file_path$forecast_output_directory)
  FLARErLER::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  FLARErLER::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)


  # Begin restart ----

  config <- FLARErLER::update_run_config(config = config, lake_directory =  lake_directory, configure_run_file =  configure_run_file, saved_file = saved_file, new_horizon = 2)

  config$run_config
  # Rename restart dir
  config$run_config$sim_name <- "test-restart"
  config$file_path$execute_directory <- file.path(lake_directory, "flare_tempdir",
                                                  config$location$site_id,
                                                  config$run_config$sim_name)
  dir.create(config$file_path$execute_directory, recursive = TRUE)

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

  #Step up Drivers
  met_out <- FLARErLER::generate_met_files(obs_met_file = observed_met_file,
                                        out_dir = config$file_path$execute_directory,
                                        forecast_dir = forecast_dir,
                                        config = config)
  met_file_names <- met_out$filenames

  historical_met_error <- met_out$historical_met_error

  inflow_forecast_path <- config$file_path$inflow_directory

  inflow_outflow_files <- FLARErLER::create_inflow_outflow_files(inflow_file_dir = inflow_file_dir,
                                                              inflow_obs = cleaned_inflow_file,
                                                              working_directory = config$file_path$execute_directory,
                                                              config,
                                                              state_names = NULL)

  inflow_file_names <- inflow_outflow_files$inflow_file_name
  outflow_file_names <- inflow_outflow_files$outflow_file_name

  obs <- FLARErLER::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  states_config <- FLARErLER::generate_states_to_obs_mapping(states_config, obs_config)

  model_sd <- FLARErLER::initiate_model_error(config, states_config)


  init <- FLARErLER::generate_initial_conditions(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config)


  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                             pars_init = init$pars,
                                             aux_states_init = init$aux_states_init,
                                             obs = obs,
                                             obs_sd = obs_config$obs_sd,
                                             model_sd = model_sd,
                                             working_directory = config$file_path$execute_directory,
                                             met_file_names = met_file_names,
                                             inflow_file_names = inflow_file_names,
                                             outflow_file_names = outflow_file_names,
                                             config = config,
                                             pars_config = pars_config,
                                             states_config = states_config,
                                             obs_config = obs_config,
                                             management = NULL,
                                             da_method = "enkf",
                                             par_fit_method = "inflate",
                                             debug = FALSE
  )

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                                     forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)
})

# LER-Simstrat can be restarted ----
test_that("LER-Simstrat can be restarted", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler_restart_Simstrat.R"))
  config$output_settings$diagnostics_names <- NULL


  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                     pars_init = init$pars,
                                     aux_states_init = init$aux_states_init,
                                     obs = obs,
                                     obs_sd = obs_config$obs_sd,
                                     model_sd = model_sd,
                                     working_directory = config$file_path$execute_directory,
                                     met_file_names = met_file_names,
                                     inflow_file_names = inflow_file_names,
                                     outflow_file_names = outflow_file_names,
                                     config = config,
                                     pars_config = pars_config,
                                     states_config = states_config,
                                     obs_config = obs_config,
                                     management = NULL,
                                     da_method = "enkf",
                                     par_fit_method = "inflate",
                                     debug = FALSE
  )

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf(da_forecast_output = enkf_output,
                                              forecast_output_directory = config$file_path$forecast_output_directory)
  FLARErLER::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  FLARErLER::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)


  # Begin restart ----
  config <- FLARErLER::update_run_config(config = config, lake_directory =  lake_directory, configure_run_file =  configure_run_file, saved_file = saved_file, new_horizon = 2)

  config$run_config
  # Rename restart dir
  config$run_config$sim_name <- "test-restart"
  config$file_path$execute_directory <- file.path(lake_directory, "flare_tempdir",
                                                  config$location$site_id,
                                                  config$run_config$sim_name)
  dir.create(config$file_path$execute_directory, recursive = TRUE)

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

  #Step up Drivers
  met_out <- FLARErLER::generate_met_files(obs_met_file = observed_met_file,
                                        out_dir = config$file_path$execute_directory,
                                        forecast_dir = forecast_dir,
                                        config = config)
  met_file_names <- met_out$filenames

  historical_met_error <- met_out$historical_met_error

  inflow_forecast_path <- config$file_path$inflow_directory

  inflow_outflow_files <- FLARErLER::create_inflow_outflow_files(inflow_file_dir = inflow_file_dir,
                                                              inflow_obs = cleaned_inflow_file,
                                                              working_directory = config$file_path$execute_directory,
                                                              config,
                                                              state_names = NULL)

  inflow_file_names <- inflow_outflow_files$inflow_file_name
  outflow_file_names <- inflow_outflow_files$outflow_file_name

  obs <- FLARErLER::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  states_config <- FLARErLER::generate_states_to_obs_mapping(states_config, obs_config)

  model_sd <- FLARErLER::initiate_model_error(config, states_config)


  init <- FLARErLER::generate_initial_conditions(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config)


  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                             pars_init = init$pars,
                                             aux_states_init = init$aux_states_init,
                                             obs = obs,
                                             obs_sd = obs_config$obs_sd,
                                             model_sd = model_sd,
                                             working_directory = config$file_path$execute_directory,
                                             met_file_names = met_file_names,
                                             inflow_file_names = inflow_file_names,
                                             outflow_file_names = outflow_file_names,
                                             config = config,
                                             pars_config = pars_config,
                                             states_config = states_config,
                                             obs_config = obs_config,
                                             management = NULL,
                                             da_method = "enkf",
                                             par_fit_method = "inflate",
                                             debug = FALSE
  )

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                                     forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)
})

# LER-GOTM-EnKF with fail ----
test_that("LER-GOTM-EnKF can fail and write output", {

  template_folder <- system.file("example", package = "FLARErLER")

  source(file.path(template_folder, "R", "test_enkf_prep_ler.R"))
  config$model_settings$model <- "GOTM"
  config$output_settings$diagnostics_names <- NULL
  config$run_config$restart_file <- NA

  obs <- FLARErLER::create_obs_matrix(cleaned_observations_file_long,
                                   obs_config,
                                   config)

  init <- FLARErLER::generate_initial_conditions(states_config,
                                              obs_config,
                                              pars_config,
                                              obs,
                                              config)

  inf <- read.csv(inflow_outflow_files$inflow_file_names[1])
  inf[3:5, 2] <- NA
  write.csv(inf, inflow_outflow_files$inflow_file_names[1], row.names = FALSE, quote = FALSE)

  # states_init = init$states
  # pars_init = init$pars
  # aux_states_init = init$aux_states_init
  # obs = obs
  # obs_sd = obs_config$obs_sd
  # model_sd = model_sd
  # working_directory = config$file_path$execute_directory
  # met_file_names = (met_file_names)
  # inflow_file_names = (inflow_file_names)
  # outflow_file_names = (outflow_file_names)
  # config = config
  # pars_config = pars_config
  # states_config = states_config
  # obs_config = obs_config
  # management = NULL
  # da_method = "enkf"
  # par_fit_method = "inflate"
  # debug = FALSE

  #Run EnKF
  enkf_output <- FLARErLER::run_da_forecast_ler(states_init = init$states,
                                             pars_init = init$pars,
                                             aux_states_init = init$aux_states_init,
                                             obs = obs,
                                             obs_sd = obs_config$obs_sd,
                                             model_sd = model_sd,
                                             working_directory = config$file_path$execute_directory,
                                             met_file_names = met_file_names,
                                             inflow_file_names = inflow_file_names,
                                             outflow_file_names = outflow_file_names,
                                             config = config,
                                             pars_config = pars_config,
                                             states_config = states_config,
                                             obs_config = obs_config,
                                             management = NULL,
                                             da_method = "enkf",
                                             par_fit_method = "inflate",
                                             debug = FALSE
  )

  # saveRDS(object = enkf_output, file = "inst/example/benchmark_data/sampenkf_output_Simstrat.RDS")
  samp_enkf_output <- readRDS(file.path(template_folder, "benchmark_data", "sampenkf_output_GOTM.RDS"))

  testthat::expect_true(is.list(enkf_output))
  chk <- lapply(1:length(enkf_output), function(x) {
    class(enkf_output[[x]]) == class(samp_enkf_output[[x]])
  })

  testthat::expect_true(all(unlist(chk)))

  # Save forecast
  saved_file <- FLARErLER::write_forecast_netcdf_ler(da_forecast_output = enkf_output,
                                                     forecast_output_directory = config$file_path$forecast_output_directory)
  testthat::expect_true(file.exists(saved_file))

  #Create EML Metadata
  FLAREr::create_flare_metadata(file_name = saved_file,
                                da_forecast_output = enkf_output)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".xml")
  testthat::expect_true(length(file_chk) > 0)

  FLAREr::plotting_general(file_name = saved_file,
                           qaqc_data_directory = config$file_path$qaqc_data_directory)
  file_chk <- list.files(config$file_path$forecast_output_directory, pattern = ".pdf")
  testthat::expect_true(length(file_chk) > 0)

})

# end
