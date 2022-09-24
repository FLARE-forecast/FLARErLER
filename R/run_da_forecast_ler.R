#' @title Run ensemble data assimilation and/or produce forecasts
#'
#' @details Uses the ensemble data assimilation to predict water quality for a lake
#' or reservoir.  The function requires the initial conditions (`states_init`) for each
#' state and ensemble member using an array with the following dimension order:
#' states, depth, ensembles member.  If you are fitting parameters, it also requires
#' initial conditions for each parameter and ensemble member using an array (`par_init`) with the
#' following dimension order: parameters, ensemble member.  The arrays for states_init
#' and pars_init can be created using the `generate_initial_conditions()` function, if
#' starting from initial conditions in the  `states_config` data frame or from observations
#' in first time column of the `obs` array.
#'
#' @param states_init array of the initial states.  Required dimensions are `[states, depths, ensemble]`
#' @param pars_init array of the initial states.  Required dimensions are `[pars, depths, ensemble]`.  (Default = NULL)
#' @param aux_states_init list of initial conditions for auxillary states.  These are states in the GLM that
#' are require for restarting the model but are not included in data assimilation.  These are states that are not associated
#' with a value in `model_sd`.
#' @param obs array; array of the observations. Required dimensions are `[nobs, time, depth]`
#' @param obs_sd vector; vector of standard deviation for observation
#' @param model_sd vector vector of standard deviations describing the model error for each state
#' @param working_directory string; full path to directory where model executes
#' @param met_file_names vector; vector of full path meteorology file names
#' @param inflow_file_names vector or matrix;; vector of inflow file names
#' @param outflow_file_names vector or matrix; vector of outflow file names
#' @param config list; list of configurations
#' @param pars_config list; list of parameter configurations  (Default = NULL)
#' @param states_config list; list of state configurations
#' @param obs_config list; list of observation configurations
#' @param management list; list of management inputs and configuration  (Default = NULL)
#' @param da_method string; data assimilation method (enkf or pf; Default = enkf)
#' @param par_fit_method string; method for adding noise to parameters during calibration
#' @param debug boolean; add extra diagnostics for debugging (Default = FALSE)
#' @return a list is passed to `write_forecast_netcdf()` to write the
#' netcdf output and `create_flare_eml()` to generate the EML metadata
#' @export
#' @importFrom parallel clusterExport detectCores clusterEvalQ parLapply stopCluster
#' @importFrom GLM3r glm_version
#' @examples
##' \dontrun{
#' da_forecast_output <- FLAREr::run_da_forecast(states_init = init$states, pars_init = init$pars, aux_states_init = init$aux_states_init, obs = obs, obs_sd = obs_config$obs_sd, model_sd = model_sd, working_directory = config$file_path$execute_directory, met_file_names = met_file_names, inflow_file_names = inflow_file_names, outflow_file_names = outflow_file_names, config = config, pars_config = pars_config, states_config = states_config, obs_config = obs_config)
#' }

run_da_forecast_ler <- function(states_init,
                                pars_init = NULL,
                                aux_states_init,
                                obs,
                                obs_sd,
                                model_sd,
                                working_directory,
                                met_file_names,
                                inflow_file_names = NULL,
                                outflow_file_names = NULL,
                                config,
                                pars_config = NULL,
                                states_config,
                                obs_config,
                                management = NULL,
                                da_method = "enkf",
                                par_fit_method = "inflate",
                                debug = FALSE,
                                log_wq = FALSE,
                                obs_secchi = NULL,
                                obs_depth = NULL){

  if(length(states_config$state_names) > 2){
    config$include_wq <- TRUE
  }else{
    config$include_wq <- FALSE
  }

  nstates <- dim(states_init)[1]
  ndepths_modeled <- dim(states_init)[2]
  nmembers <- dim(states_init)[3]
  n_met_members <- length(met_file_names)
  model <- config$model_settings$model
  if(!is.null(pars_config)){
    if("model" %in% names(pars_config)){
      pars_config <- pars_config[pars_config$model == model, ]
    }
    npars <- nrow(pars_config)
    par_names <- pars_config$par_names
    par_file <- pars_config$par_file
  }else{
    npars <- 0
    par_names <- NA
    par_file <- NA
  }

  start_datetime <- lubridate::as_datetime(config$run_config$start_datetime)
  if(is.na(config$run_config$forecast_start_datetime)){
    end_datetime <- lubridate::as_datetime(config$run_config$end_datetime)
    forecast_start_datetime <- end_datetime
  }else{
    forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime)
    end_datetime <- forecast_start_datetime + lubridate::days(config$run_config$forecast_horizon)
  }

  hist_days <- as.numeric(forecast_start_datetime - start_datetime)
  start_forecast_step <- 1 + hist_days
  full_time <- seq(start_datetime, end_datetime, by = "1 day")
  forecast_days <- as.numeric(end_datetime - forecast_start_datetime)
  nsteps <- length(full_time)

  data_assimilation_flag <- rep(NA, nsteps)
  forecast_flag <- rep(NA, nsteps)
  da_qc_flag <- rep(NA, nsteps)

  x <- array(NA, dim=c(nsteps, nstates, ndepths_modeled, nmembers))
  x[1, , , ]  <- states_init

  if(npars > 0){
    pars <- array(NA, dim=c(nsteps, npars, nmembers))
    pars[1, , ] <- pars_init
  }else{
    pars <- NULL
  }

  q_v <- rep(NA, ndepths_modeled)
  w <- rep(NA, ndepths_modeled)
  w_new <- rep(NA, ndepths_modeled)

  alpha_v <- 1 - exp(-states_config$vert_decorr_length)

  output_vars <- states_config$state_names

  if(config$include_wq){
    num_wq_vars <- dim(x)[2] - 2
  }else{
    num_wq_vars <- 0
  }

  if(length(config$output_settings$diagnostics_names) > 0){
    diagnostics <- array(NA, dim=c(length(config$output_settings$diagnostics_names), nsteps, ndepths_modeled, nmembers))
  }else{
    diagnostics <- NA
  }

  full_time_char <- strftime(full_time,
                             format="%Y-%m-%d %H:%M",
                             tz = "UTC")

  file.copy(from = file.path(config$file_path$configuration_directory, config$model_settings$ler_bathymetry_file),
            to = file.path(working_directory, config$model_settings$ler_bathymetry_file), overwrite = TRUE)

  if(!is.null(inflow_file_names)){
    inflow_file_names <- as.matrix(inflow_file_names)
    outflow_file_names <- as.matrix(outflow_file_names)
  }else{
    inflow_file_names <- NULL
    outflow_file_names <- NULL
  }

  config$model_settings$ncore <- min(c(config$model_settings$ncore, parallel::detectCores()))
  if(config$model_settings$ncore == 1) {
    if(!dir.exists(file.path(working_directory, "1"))) {
      dir.create(file.path(working_directory, "1"), showWarnings = FALSE)
    } else {
      unlink(file.path(working_directory, "1"), recursive = TRUE)
      dir.create(file.path(working_directory, "1"), showWarnings = FALSE)
    }
    FLARErLER:::set_up_model_ler(model,
                                 config,
                                 working_directory = working_directory,
                                 state_names = states_config$state_names,
                                 met_file_names = met_file_names,
                                 inflow_file_names = inflow_file_names,
                                 outflow_file_names = outflow_file_names,
                                 member = 1,
                                 start_datetime,
                                 end_datetime)
  } else {
    purrr::walk(1:nmembers, function(m){
      if(!dir.exists(file.path(working_directory, m))) {
        dir.create(file.path(working_directory, m), showWarnings = FALSE)
      } else {
        unlink(file.path(working_directory, m), recursive = TRUE)
        dir.create(file.path(working_directory, m), showWarnings = FALSE)
      }
      FLARErLER:::set_up_model_ler(model,
                                   config,
                                   working_directory = working_directory,
                                   state_names = states_config$state_names,
                                   met_file_names = met_file_names,
                                   inflow_file_names = inflow_file_names,
                                   outflow_file_names = outflow_file_names,
                                   member = m,
                                   start_datetime,
                                   end_datetime)
    })
  }

  lake_depth <- array(NA, dim = c(nsteps, nmembers))

  restart_list <- NULL
  if(model == "GLM") {

    lake_depth <- array(NA, dim = c(nsteps, nmembers))
    the_depths <- array(NA, dim = c(nsteps, ndepths_modeled, nmembers))
    model_internal_depths <- array(NA, dim = c(nsteps, 500, nmembers))
    the_sals <- array(NA, dim = c(nsteps, ndepths_modeled, nmembers))
    snow_thickness <- array(NA, dim = c(nsteps, nmembers))
    white_ice_thickness <- array(NA, dim = c(nsteps, nmembers))
    blue_ice_thickness <- array(NA, dim = c(nsteps, nmembers))
    avg_surf_temp <- array(NA, dim = c(nsteps, nmembers))
    restart_variables <- array(NA, dim = c(17, nsteps, nmembers))

    the_depths[1, ,] <- aux_states_init$the_depths
    lake_depth[1, ] <- aux_states_init$lake_depth
    model_internal_depths[1, ,] <- aux_states_init$model_internal_depths
    snow_thickness[1, ] <- aux_states_init$snow_thickness
    white_ice_thickness[1, ] <- aux_states_init$white_ice_thickness
    blue_ice_thickness[1, ] <- aux_states_init$blue_ice_thickness
    the_sals[1, , ] <- aux_states_init$the_sals
    avg_surf_temp[1, ] <- aux_states_init$avg_surf_temp
    restart_variables[, 1, ] <- 0

    restart_list <- list(lake_depth = lake_depth,
                         model_internal_depths = model_internal_depths,
                         the_depths = the_depths,
                         the_sals = the_sals,
                         snow_thickness = snow_thickness,
                         white_ice_thickness = white_ice_thickness,
                         blue_ice_thickness = blue_ice_thickness,
                         avg_surf_temp = avg_surf_temp,
                         restart_variables = restart_variables)

  } else if(model == "GOTM") {
    yaml <- gotmtools::read_yaml(file.path(working_directory, "1/GOTM", "gotm.yaml"))
    nz <- yaml$grid$nlev
    nzi <- nz + 1

    # z vars
    z <- array(NA, dim = c(nsteps, nz, nmembers))
    temp <- array(NA, dim = c(nsteps, nz, nmembers))
    salt <- array(NA, dim = c(nsteps, nz, nmembers))
    u <- array(NA, dim = c(nsteps, nz, nmembers))
    uo <- array(NA, dim = c(nsteps, nz, nmembers))
    v <- array(NA, dim = c(nsteps, nz, nmembers))
    vo <- array(NA, dim = c(nsteps, nz, nmembers))
    xP <- array(NA, dim = c(nsteps, nz, nmembers))
    h <- array(NA, dim = c(nsteps, nz, nmembers))
    ho <- array(NA, dim = c(nsteps, nz, nmembers))

    #zi vars
    tke <- array(NA, dim = c(nsteps, nzi, nmembers))
    zi <- array(NA, dim = c(nsteps, nzi, nmembers))
    tkeo <- array(NA, dim = c(nsteps, nzi, nmembers))
    eps <- array(NA, dim = c(nsteps, nzi, nmembers))
    num <- array(NA, dim = c(nsteps, nzi, nmembers))
    nuh <- array(NA, dim = c(nsteps, nzi, nmembers))
    nus <- array(NA, dim = c(nsteps, nzi, nmembers))

    if(length(aux_states_init) > 0) {
      # z vars
      z[1, , ] <- aux_states_init$z_vars$z
      lake_depth[1, ] <- max(aux_states_init$z_vars$z)
      # temp[1, , ] <- aux_states_init$z_vars$temp
      #salt[1, , ] <- aux_states_init$z_vars$salt
      u[1, , ] <- aux_states_init$z_vars$u
      uo[1, , ] <- aux_states_init$z_vars$uo
      v[1, , ] <- aux_states_init$z_vars$v
      vo[1, , ] <- aux_states_init$z_vars$uo
      xP[1, , ] <- aux_states_init$z_vars$xP
      h[1, , ] <- aux_states_init$z_vars$h
      ho[1, , ] <- aux_states_init$z_vars$ho

      #zi vars
      tke[1, , ] <- aux_states_init$zi_vars$tke
      zi[1, , ] <- aux_states_init$zi_vars$zi
      tkeo[1, , ] <- aux_states_init$zi_vars$tkeo
      eps[1, , ] <- aux_states_init$zi_vars$eps
      num[1, , ] <- aux_states_init$zi_vars$num
      nuh[1, , ] <- aux_states_init$zi_vars$nuh
      nus[1, , ] <- aux_states_init$zi_vars$nus
    }

    restart_list <- list(lake_depth = lake_depth,
                         z_vars = list(z = z,
                                       temp = temp,
                                       salt = salt,
                                       u = u,
                                       uo = uo,
                                       v = v,
                                       vo = vo,
                                       xP = xP,
                                       h = h,
                                       ho = ho),
                         zi_vars = list(tke = tke,
                                        zi = zi,
                                        tkeo = tkeo,
                                        eps = eps,
                                        num = num,
                                        nuh = nuh,
                                        nus = nus))
  } else if(model == "Simstrat") {

    ngrid <- LakeEnsemblR::get_json_value(file.path(working_directory, "1", model,
                                                    "simstrat.par"), label = "Input", key = "Grid")
    if(is.character(ngrid)) {
      nlev <- read.delim(file.path(working_directory, "1", model, ngrid))
      ngrid <- nlev[[1]]
    }
    nzi <- (ngrid * 2) + 1

    #zi vars
    zi <- array(NA, dim = c(nsteps, nzi, nmembers))
    u <- array(NA, dim = c(nsteps, nzi, nmembers))
    v <- array(NA, dim = c(nsteps, nzi, nmembers))
    temp <- array(NA, dim = c(nsteps, nzi, nmembers))
    S <- array(NA, dim = c(nsteps, nzi, nmembers))
    k <- array(NA, dim = c(nsteps, nzi, nmembers))
    eps <- array(NA, dim = c(nsteps, nzi, nmembers))
    num <- array(NA, dim = c(nsteps, nzi, nmembers))
    nuh <- array(NA, dim = c(nsteps, nzi, nmembers))

    seicheE <- array(NA, dim = c(nsteps, nmembers))
    b_ice <- array(NA, dim = c(nsteps, nmembers))
    w_ice <- array(NA, dim = c(nsteps, nmembers))
    snow <- array(NA, dim = c(nsteps, nmembers))

    if(length(aux_states_init) > 0) {
      lake_depth[1, ] <- max(aux_states_init$zi)
      zi[1, , ] <- aux_states_init$zi
      u[1, , ] <- aux_states_init$u
      v[1, , ] <- aux_states_init$v
      # temp[1, , ] <- aux_states_init$temp
      #S[1, , ] <- aux_states_init$S
      k[1, , ] <- aux_states_init$k
      eps[1, , ] <- aux_states_init$eps
      num[1, , ] <- aux_states_init$num
      nuh[1, , ] <- aux_states_init$nuh

      seicheE[1, ] <- aux_states_init$seicheE
      b_ice[1, ] <- aux_states_init$b_ice
      w_ice[1, ] <- aux_states_init$w_ice
      snow[1, ] <- aux_states_init$snow
    }

    restart_list <- list(lake_depth = lake_depth,
                         zi = zi,
                         u = u,
                         v = v,
                         temp = temp,
                         S = S,
                         k = k,
                         eps = eps,
                         num = num,
                         nuh = nuh,
                         seicheE = seicheE,
                         b_ice = b_ice,
                         w_ice = w_ice,
                         snow = snow)
  }

  if(config$da_setup$assimilate_first_step){
    start_step <- 1
  }else{
    start_step <- 2
  }


  if(model == "GLM") {
    # Print GLM version
    glm_v <- GLM3r::glm_version()
    glm_v <- substr(glm_v[3], 35, 58)
    message("Using GLM ", glm_v)
    config$metadata$model_description$version <- substr(glm_v, 9, 16)
  }
  ###START EnKF

  for(i in start_step:nsteps){

    #setwd("/Users/quinn/workfiles/Research/SSC_forecasting/ler/FCRE_LER-forecast-code/flare_tempdir/fcre/ms1_ler_flare_Simstrat")
    #load("../../../../workspace_inENKF.RData")

    if(i > 1){
      curr_start <- strftime(full_time[i - 1],
                             format="%Y-%m-%d %H:%M:%S",
                             tz = "UTC")
    }else{
      curr_start <- "Restart"
    }
    curr_stop <- strftime(full_time[i],
                          format="%Y-%m-%d %H:%M:%S",
                          tz = "UTC")

    message(paste0("Running time step ", i-1, "/", (nsteps - 1), " : ",
                   curr_start, " - ",
                   curr_stop, " [", Sys.time(), "]"))

    setwd(working_directory)

    met_index <- rep(1:length(met_file_names), times = nmembers)
    if(!is.null(ncol(inflow_file_names))) {
      inflow_outflow_index <- rep(1:nrow(inflow_file_names), times = nmembers)
    } else {
      inflow_outflow_index <- NULL
    }

    #Create array to hold GLM predictions for each ensemble
    x_star <- array(NA, dim = c(nstates, ndepths_modeled, nmembers))
    x_corr <- array(NA, dim = c(nstates, ndepths_modeled, nmembers))
    curr_pars <- array(NA, dim = c(npars, nmembers))

    if(i == start_step) {
      restart = FALSE
    } else {
      restart = TRUE
    }

    # if i = start_step set up cluster for parallelization
    # Switch for
    switch(Sys.info() [["sysname"]],
           Linux = { machine <- "unix" },
           Darwin = { machine <- "mac" },
           Windows = { machine <- "windows"})

    #If i == 1 then assimilate the first time step without running the process
    #model (i.e., use yesterday's forecast of today as initial conditions and
    #assimilate new observations)
    if(i > 1){

      if(config$model_settings$ncore == 1){
        future::plan("future::sequential", workers = config$model_settings$ncore)
      }else{
        future::plan("future::multisession", workers = config$model_settings$ncore)
      }

      out <- furrr::future_map(1:nmembers, function(m) {

        ens_dir_index <- m
        if(config$model_settings$ncore == 1) ens_dir_index <- 1

        setwd(file.path(working_directory, ens_dir_index))

        curr_met_file <- met_file_names[met_index[m]]

        if(npars > 0){
          if(par_fit_method == "inflate" & da_method == "enkf"){
            curr_pars_ens <-  pars[i-1, , m]
          }else if(par_fit_method %in% c("perturb","perturb_const") & da_method != "none"){
            if(par_fit_method == "perturb_const"){
              if(npars > 1){
                par_mean <- apply(pars[i-1, , ], 1, mean)
                par_sd <- apply(pars[i-1, , ], 1, sd)
              }else{
                par_mean <- mean(pars[i-1, , ])
                par_sd <- sd(pars[i-1, , ])
              }


              par_z <- (pars[i-1, ,m] - par_mean)/par_sd

              curr_pars_ens <- par_z * pars_config$perturb_par + par_mean

            }else{
              if(i < (hist_days + 1)){
                curr_pars_ens <- pars[i-1, , m] + rnorm(npars, mean = rep(0, npars), sd = pars_config$perturb_par)
              }else{
                curr_pars_ens <- pars[i-1, , m]
              }

            }
          }else if(da_method == "none" | par_fit_method == "perturb_init"){
            curr_pars_ens <- pars[i-1, , m]
          }else{
            message("parameter fitting method not supported.  inflate or perturb are supported. only inflate is supported for enkf")

          }
        }else{
          curr_pars_ens <- NULL
        }

        if(!is.null(ncol(inflow_file_names))){
          inflow_file_name <- inflow_file_names[inflow_outflow_index[m], ]
          outflow_file_name <- outflow_file_names[inflow_outflow_index[m], ]
        }else{
          inflow_file_name <- NULL
          outflow_file_name <- NULL
        }

        if(config$model_settings$ncore == 1) {
          wdir <- file.path(working_directory, 1)
        } else {
          wdir <- file.path(working_directory, m)
        }

        out <- suppressMessages({FLARErLER:::run_model_ler(model,
                             i,
                             m,
                             curr_start,
                             curr_stop,
                             par_names,
                             curr_pars = curr_pars_ens,
                             working_directory = wdir,
                             par_file,
                             x_start = x[i-1, , ,m ],
                             full_time,
                             management = management,
                             hist_days,
                             modeled_depths = config$model_settings$modeled_depths,
                             ndepths_modeled,
                             curr_met_file,
                             inflow_file_name = inflow_file_name,
                             outflow_file_name = outflow_file_name,
                             output_vars = output_vars,
                             diagnostics_names = config$output_settings$diagnostics_names,
                             npars,
                             num_wq_vars,
                             nstates,
                             state_names = states_config$state_names,
                             include_wq = config$include_wq,
                             restart = restart,
                             restart_list = restart_list)
          })

      }, .options = furrr::furrr_options(seed = TRUE))

      # Loop through output and assign to matrix
      # Loop through output and assign to matrix
      for(m in 1:nmembers) {

        x_star[, , m] <- out[[m]]$x_star_end
        curr_pars[, m] <- out[[m]]$curr_pars
        # lake_depth[i ,m ] <- out[[m]]$lake_depth_end
        # snow_ice_thickness[,i ,m] <- out[[m]]$snow_ice_thickness_end
        if(length(config$output_settings$diagnostics_names) > 0) {
          diagnostics[, i, , m] <- out[[m]]$diagnostics_end
        }
        # model_internal_depths[i, ,m] <- out[[m]]$model_internal_depths
        # salt[i, , m]  <- out[[m]]$salt_end
        if(model == "GLM") {
          restart_list$lake_depth[i, m] <- out[[m]]$lake_depth_end
          restart_list$model_internal_depths[i, , m] <- out[[m]]$model_internal_depths
          # restart_list$the_depths[i, , m] <- out[[m]]$the_depths
          #restart_list$the_sals[i, , m] <- approx(out[[m]]$restart_vars$the_depths, out[[m]]$restart_vars$the_sals, config$model_settings$modeled_depths, rule = 2)$y
          restart_list$snow_thickness[i, m] <- out[[m]]$restart_vars$snow_thickness
          restart_list$white_ice_thickness[i, m] <- out[[m]]$restart_vars$white_ice_thickness
          restart_list$blue_ice_thickness[i, m] <- out[[m]]$restart_vars$blue_ice_thickness
          restart_list$avg_surf_temp[i , m] <- out[[m]]$restart_vars$avg_surf_temp
          restart_list$restart_variables[, i, m] <- out[[m]]$restart_vars$restart_variables

        } else if(model == "GOTM") {
          # z vars
          restart_list$lake_depth[i, m] <- max(-out[[m]]$restart_vars$zi_vars$zi)
          restart_list$z_vars$z[i, , m] <- out[[m]]$restart_vars$z_vars$z
          restart_list$z_vars$temp[i, , m] <- out[[m]]$restart_vars$z_vars$temp
          restart_list$z_vars$salt[i, , m] <- out[[m]]$restart_vars$z_vars$salt
          restart_list$z_vars$u[i, , m] <- out[[m]]$restart_vars$z_vars$u
          restart_list$z_vars$uo[i, , m] <- out[[m]]$restart_vars$z_vars$uo
          restart_list$z_vars$v[i, , m] <- out[[m]]$restart_vars$z_vars$v
          restart_list$z_vars$vo[i, , m] <- out[[m]]$restart_vars$z_vars$vo
          restart_list$z_vars$xP[i, , m] <- out[[m]]$restart_vars$z_vars$xP
          restart_list$z_vars$h[i, , m] <- out[[m]]$restart_vars$z_vars$h
          restart_list$z_vars$ho[i, , m] <- out[[m]]$restart_vars$z_vars$ho

          # zi vars
          restart_list$zi_vars$tke[i, , m] <- out[[m]]$restart_vars$zi_vars$tke
          restart_list$zi_vars$zi[i, , m] <- out[[m]]$restart_vars$zi_vars$zi
          restart_list$zi_vars$tkeo[i, , m] <- out[[m]]$restart_vars$zi_vars$tkeo
          restart_list$zi_vars$eps[i, , m] <- out[[m]]$restart_vars$zi_vars$eps
          restart_list$zi_vars$num[i, , m] <- out[[m]]$restart_vars$zi_vars$num
          restart_list$zi_vars$nuh[i, , m] <- out[[m]]$restart_vars$zi_vars$nuh
          restart_list$zi_vars$nus[i, , m] <- out[[m]]$restart_vars$zi_vars$nus
        } else if(model == "Simstrat") {
          restart_list$lake_depth[i, m] <- out[[m]]$lake_depth_end
          restart_list$zi[i, , m] <- out[[m]]$restart_vars$zi
          restart_list$u[i, , m] <- out[[m]]$restart_vars$u
          restart_list$v[i, , m] <- out[[m]]$restart_vars$v
          restart_list$temp[i, , m] <- out[[m]]$restart_vars$temp
          restart_list$S[i, , m] <- out[[m]]$restart_vars$S
          restart_list$k[i, , m] <- out[[m]]$restart_vars$k
          restart_list$eps[i, , m] <- out[[m]]$restart_vars$eps
          restart_list$num[i, , m] <- out[[m]]$restart_vars$num
          restart_list$nuh[i, , m] <- out[[m]]$restart_vars$nuh
          restart_list$seicheE[i , m] <- out[[m]]$restart_vars$seicheE
          restart_list$b_ice[i , m] <- out[[m]]$restart_vars$b_ice
          restart_list$w_ice[i , m] <- out[[m]]$restart_vars$w_ice
          restart_list$snow[i , m] <- out[[m]]$restart_vars$snow
        }

        #Add process noise
        q_v[] <- NA
        w[] <- NA
        w_new[] <- NA
        for(jj in 1:nrow(model_sd)){
          w[] <- rnorm(ndepths_modeled, 0, 1)
          for(kk in 1:ndepths_modeled){
            #q_v[kk] <- alpha_v * q_v[kk-1] + sqrt(1 - alpha_v^2) * model_sd[jj, kk] * w[kk]
            if(kk == 1){
              w_new[kk] <- w[kk]
            }else{
              w_new[kk] <- (alpha_v[jj] * w_new[kk-1] + sqrt(1 - alpha_v[jj]^2) * w[kk])
            }
            q_v[kk] <- w_new[kk] * model_sd[jj, kk] #sqrt(log(1 + model_sd[jj, kk] ^ 2))
            x_corr[jj, kk, m] <-
              x_star[jj, kk, m] + q_v[kk] #- 0.5*log(1 + model_sd[jj, kk]^2)
          }
        }
        #if(log_wq){
        #  index <- which(x_corr[m, ] <=  0.0000001)
        #  x_corr[m, index[which(index > ndepths_modeled)]] <- 0.0000001
        #}

      } # END ENSEMBLE LOOP

      #Correct any negative water quality states
      if(length(states_config$state_names) > 1 & !log_wq){
        for(s in 2:nstates){
          for(k in 1:ndepths_modeled){
            index <- which(x_corr[s, k, ] < 0.0)
            x_corr[s, k, index] <- 0.0
          }
        }
      }

      if(npars > 0){
        pars_corr <- curr_pars
        if(npars == 1){
          pars_corr <- matrix(pars_corr,nrow = 1 ,ncol = length(pars_corr))
        }
        pars_star <- pars_corr
      }

    }else{
      x_star <- x[i, , , ]
      x_corr <- x_star
      if(npars > 0){
        pars_corr <- curr_pars
        if(npars == 1){
          pars_corr <- matrix(pars_corr,nrow = length(pars_corr),ncol = 1)
        }
        pars_star <- pars_corr
      }
    }

    if(dim(obs)[1] > 1){
      z_index <- which(!is.na(c(aperm(obs[,i , ], perm = c(2,1)))))
      max_index <- length(c(aperm(obs[,i , ], perm = c(2,1)))) + 1
    }else{
      z_index <- which(!is.na(c(obs[1,i , ])))
      max_index <- length(c(obs[1,i , ])) + 1
    }

    if(!is.null(obs_secchi)){
      if(!is.na(obs_secchi[i])){
        z_index <- c(z_index,max_index)
      }
    }

    #if(!is.null(obs_depth)){
    #  if(!is.na(obs_depth[i])){
    #    if(!is.na(obs_secchi[i])){
    #    z_index <- c(z_index,max_index +1)
    #  }else{
    #    z_index <- c(z_index,max_index)
    #  }
    #  }
    #}

    #if no observations at a time step then just propogate model uncertainity

    if(length(z_index) == 0 | config$da_setup$da_method == "none" | !config$da_setup$use_obs_constraint){

      if(i > (hist_days + 1)){
        data_assimilation_flag[i] <- 0
        forecast_flag[i] <- 1
        da_qc_flag[i] <- 0
      }else if(i <= (hist_days + 1) & config$da_setup$use_obs_constraint){
        data_assimilation_flag[i] <- 1
        forecast_flag[i] <- 0
        da_qc_flag[i] <- 1
      }else{
        data_assimilation_flag[i] <- 0
        forecast_flag[i] <- 0
        da_qc_flag[i] <- 0
      }

      if(npars > 0){

        x[i, , , ] <- x_corr
        if(npars > 0) pars[i, , ] <- pars_star

        if(config$uncertainty$process == FALSE & i > (hist_days + 1)){
          #don't add process noise if process uncertainty is false (x_star doesn't have noise)
          #don't add the noise to parameters in future forecast mode ()
          x[i, , , ] <- x_star
          if(npars > 0) pars[i, , ] <- pars_star
        }

        if(i == (hist_days + 1) & config$uncertainty$initial_condition == FALSE){
          pars[i, , ] <- pars_star
          for(s in 1:nstates){
            for(k in 1:ndepths_modeled){
              x[i, s, k , ] <- mean(x_star[s, k, ])
            }
          }
        }
      }

      for(s in 1:nstates){
        for(m in 1:nmembers){
          depth_index <- which(config$model_settings$modeled_depths > restart_list$lake_depth[i, m])
          x[i, s, depth_index, m ] <- NA
        }
      }
      if(length(config$output_settings$diagnostics_names) > 0){
        for(d in 1:dim(diagnostics)[1]){
          for(m in 1:nmembers){
            depth_index <- which(config$model_settings$modeled_depths > restart_list$lake_depth[i, m])
            diagnostics[d,i, depth_index, m] <- NA
          }
        }
      }
      #if(log_wq){
      #  x[i, , (ndepths_modeled+1):nstates] <- exp(x[i, , (ndepths_modeled+1):nstates])
      #}
    }else{


      x_matrix <- apply(aperm(x_corr[,1:ndepths_modeled,], perm = c(2,1,3)), 3, rbind)
      if(length(config$model_settings$diagnostics_names) > 0){
        modeled_secchi <- 1.7 / diagnostics[1, i, which.min(abs(config$model_settings$modeled_depths-1.0)), ]
        if(!is.null(obs_secchi)){
          if(!is.na(obs_secchi[i])){
            x_matrix <- rbind(x_matrix, modeled_secchi)
          }
        }
      }

      max_obs_index <- 0
      for(s in 1:dim(obs)[1]){
        if(length(which(!is.na(obs[s,i ,]))) > 0){
          max_obs_index <- max(max_obs_index, max(max_obs_index, max(which(!is.na(obs[s,i ,])))))
        }
      }

      for(m in 1:nmembers){
        if(config$model_settings$modeled_depths[max_obs_index] > restart_list$lake_depth[i, m]){
          restart_list$lake_depth[i, m] = config$model_settings$modeled_depths[max_obs_index]
        }
      }

      data_assimilation_flag[i] <- 1
      forecast_flag[i] <- 0
      da_qc_flag[i] <- 0

      curr_obs <- obs[,i,]

      #if(log_wq){
      #  for(kk in 2:dim(curr_obs)[1]){
      #    curr_obs[kk, which(curr_obs[kk, ] <= 0)] <- 0.0000001
      #    curr_obs[kk, ] <- log(curr_obs[kk, ])
      #  }
      #}

      #if observation then calucate Kalman adjustment
      if(dim(obs)[1] > 1){
        zt <- c(aperm(curr_obs, perm = c(2,1)))
        #zt[(ndepths_modeled+1):nstates] <- zt[(ndepths_modeled+1):nstates]  - 0.5*log(1 + psi[(ndepths_modeled+1):nstates]^2)
      }else{
        zt <- curr_obs
      }

      zt <- zt[which(!is.na(zt))]

      secchi_index <- 0
      if(!is.null(obs_secchi)){
        secchi_index <- 1
        if(!is.na(obs_secchi[i])){
          zt <- c(zt, obs_secchi[i])

        }
      }

      depth_index <- 0
      #if(!is.null(obs_depth)){
      #  depth_index <- 1
      #  if(!is.na(obs_depth[i])){
      #    zt <- c(zt, obs_depth[i])
      #  }
      #}


      #Assign which states have obs in the time step
      h <- matrix(0, nrow = length(obs_sd) * ndepths_modeled + secchi_index + depth_index, ncol = nstates * ndepths_modeled + secchi_index + depth_index)

      index <- 0
      for(k in 1:nstates){
        for(j in 1:ndepths_modeled){
          index <- index + 1
          if(!is.na(dplyr::first(states_config$states_to_obs[[k]]))){
            for(jj in 1:length(states_config$states_to_obs[[k]])){
              if(!is.na((obs[states_config$states_to_obs[[k]][jj], i, j]))){
                states_to_obs_index <- states_config$states_to_obs[[k]][jj]
                index2 <- (states_to_obs_index - 1) * ndepths_modeled + j
                h[index2,index] <- states_config$states_to_obs_mapping[[k]][jj]
              }
            }
          }
        }
      }

      if(!is.null(obs_secchi)){
        if(!is.na(obs_secchi[i])){
          h[dim(h)[1],dim(h)[2]] <- 1
        }
      }

      z_index <- c()
      for(j in 1:nrow(h)){
        if(sum(h[j, ]) > 0){
          z_index <- c(z_index, j)
        }
      }

      h <- h[z_index, ]

      if(!is.matrix(h)){
        h <- t(as.matrix(h))
      }

      if(da_method == "enkf"){

        #Extract the data uncertainty for the data
        #types present during the time-step

        psi <- rep(NA, length(obs_sd) * ndepths_modeled + secchi_index)
        index <- 0
        for(k in 1:length(obs_sd)){
          for(j in 1:ndepths_modeled){
            index <- index + 1
            if(k == 1){
              psi[index] <- obs_sd[k]
            }else{
              psi[index] <- obs_sd[k] #sqrt(log(1 + obs_sd[i] ^ 2))
            }
          }
        }

        if(!is.null(obs_secchi)){
          psi[length(psi)] <- 0.2
        }

        curr_psi <- psi[z_index] ^ 2

        if(length(z_index) > 1){
          psi_t <- diag(curr_psi)
        }else{
          #Special case where there is only one data
          #type during the time-step
          psi_t <- curr_psi
        }

        d_mat <- t(mvtnorm::rmvnorm(n = nmembers, mean = zt, sigma=as.matrix(psi_t)))

        #Set any negative observations of water quality variables to zero
        if(!log_wq){
          d_mat[which(z_index > ndepths_modeled & d_mat < 0.0)] <- 0.0
        }

        #Ensemble mean
        ens_mean <- apply(x_matrix, 1, mean)

        if(npars > 0){
          par_mean <- apply(pars_corr, 1, mean)
          if(par_fit_method == "inflate"){
            for(m in 1:nmembers){
              pars_corr[, m] <- pars_config$inflat_pars * (pars_corr[, m] - par_mean) + par_mean
            }
            par_mean <- apply(pars_corr, 1, mean)
          }
        }

        dit <- matrix(NA, nrow = nmembers, ncol = dim(x_matrix)[1])

        if(npars > 0) dit_pars<- array(NA, dim = c(nmembers, npars))

        #Loop through ensemble members
        for(m in 1:nmembers){
          #  #Ensemble specific deviation
          dit[m, ] <- x_matrix[, m] - ens_mean
          if(npars > 0){
            dit_pars[m, ] <- pars_corr[, m] - par_mean
          }
          if(m == 1){
            p_it <- dit[m, ] %*% t(dit[m, ])
            if(npars > 0){
              p_it_pars <- dit_pars[m, ] %*% t(dit[m, ])
            }
          }else{
            p_it <- dit[m, ] %*% t(dit[m, ]) +  p_it
            if(npars > 0){
              p_it_pars <- dit_pars[m, ] %*% t(dit[m, ]) + p_it_pars
            }
          }
        }

        #estimate covariance
        p_t <- p_it / (nmembers - 1)
        if(npars > 0){
          p_t_pars <- p_it_pars / (nmembers - 1)
        }

        #if(!is.na(config$da_setup$localization_distance)){
        #  p_t <- localization(p_t,
        #                      nstates,
        #                      modeled_depths = config$model_settings$modeled_depths,
        #                      localization_distance = config$da_setup$localization_distance)
        #}
        #Kalman gain
        k_t <- p_t %*% t(h) %*% solve(h %*% p_t %*% t(h) + psi_t, tol = 1e-17)
        if(npars > 0){
          k_t_pars <- p_t_pars %*% t(h) %*% solve(h %*% p_t %*% t(h) + psi_t, tol = 1e-17)
        }

        #Update states array (transposes are necessary to convert
        #between the dims here and the dims in the EnKF formulations)
        update <-  x_matrix + k_t %*% (d_mat - h %*% x_matrix)
        update <- update[1:(ndepths_modeled*nstates), ]
        update <- aperm(array(c(update), dim = c(ndepths_modeled, nstates, nmembers)), perm = c(2,1,3))

        if(!is.null(obs_depth)){
          if(!is.na(obs_depth[i])){
            restart_list$lake_depth[i, ] <- rnorm(nmembers, obs_depth[i], sd = 0.05)
            for(m in 1:nmembers){
              depth_index <- which(model_internal_depths[i, , m] > restart_list$lake_depth[i, m])
              model_internal_depths[i,depth_index , m] <- NA
            }
          }
        }

        for(s in 1:nstates){
          for(m in 1:nmembers){
            depth_index <- which(config$model_settings$modeled_depths <= restart_list$lake_depth[i, m])
            x[i, s, depth_index, m ] <- update[s,depth_index , m]
          }
        }

        if(length(config$model_settings$diagnostics_names) > 0){
          for(d in 1:dim(diagnostics)[1]){
            for(m in 1:nmembers){
              depth_index <- which(config$model_settings$modeled_depths > restart_list$lake_depth[i, m])
              diagnostics[d,i, depth_index, m] <- NA
            }
          }
        }

        #  if(max_depth_index < ndepths_modeled){
        #    x[i, s,(max_depth_index+1):ndepths_modeled, ] <- NA
        #  }
        #  if(max_depth_index < max_obs_index){
        #    for(m in 1:nmembers){
        #      if(!is.na(states_config$states_to_obs[[s]]) &
        #         length(which(!is.na(obs[states_config$states_to_obs[[s]],i,(max_depth_index+1):ndepths_modeled]))) > 0){
        #        new_depth <- config$model_settings$modeled_depths[max_depth_index:ndepths_modeled]
        #        modeled <- c(x[i, s,max_depth_index, m ], obs[states_config$states_to_obs[[s]],i,(max_depth_index+1):ndepths_modeled])
        #        x[i, s,(max_depth_index+1):ndepths_modeled, m] <- approx(new_depth,modeled , xout = new_depth[-1], rule = 2)$y
        #      }else{
        #        x[i, s,(max_depth_index+1):ndepths_modeled, m] <- x[i, s,max_depth_index, m ]
        #      }
        #      lake_depth[i, m] <- max(lake_depth[i, m], onfig$model_settings$modeled_depths[max_obs_index])
        #    }
        #  }
        #}

        if(npars > 0){
          if(par_fit_method != "perturb_init"){
            pars[i, , ] <- pars_corr + k_t_pars %*% (d_mat - h %*% x_matrix)
          }else{
            pars[i, , ]  <- pars[i-1, , ]
          }
        }

        #if(log_wq){
        #  x[i, , (ndepths_modeled+1):nstates] <- exp(x[i, , (ndepths_modeled+1):nstates])
        #}

      }else if(da_method == "pf"){

        obs_states <- t(h %*% t(x_matrix))

        LL <- rep(NA, length(nmembers))
        for(m in 1:nmembers){
          LL[m] <- sum(dnorm(zt, mean = obs_states[m, ], sd = psi[z_index], log = TRUE))
        }

        sample <- sample.int(nmembers, replace = TRUE, prob = exp(LL))

        update <- x_matrix[, sample]

        if(npars > 0){
          pars[i] <- pars_star[, sample]
        }

        if(model == "GLM") {
          restart_list$lake_depth[i, ] <- restart_list$lake_depth[i, sample]
          restart_list$model_internal_depths[i, , ] <- restart_list$model_internal_depths[i, , sample]
          # restart_list$the_depths[i, , ] <- restart_list$the_depths[i, , sample]
          restart_list$the_sals[i, , ] <- restart_list$the_sals[i, , sample]
          restart_list$snow_thickness[i, ] <- restart_list$snow_thickness[i, sample]
          restart_list$white_ice_thickness[i, ] <- restart_list$white_ice_thickness[i, sample]
          restart_list$blue_ice_thickness[i, ] <- restart_list$blue_ice_thickness[i, sample]
          restart_list$avg_surf_temp[i , ] <- restart_list$avg_surf_temp[i , sample]
          restart_list$restart_variables[, i, ] <- restart_list$restart_variables[, i, sample]
        } else if(model == "GOTM") {
          # z vars
          restart_list$z_vars$z[i, , ] <- restart_list$z_vars$z[i, , sample]
          restart_list$z_vars$temp[i, , ] <- restart_list$z_vars$temp[i, , sample]
          restart_list$z_vars$salt[i, , ] <- restart_list$z_vars$salt[i, , sample]
          restart_list$z_vars$u[i, , ] <- restart_list$z_vars$u[i, , sample]
          restart_list$z_vars$uo[i, , ] <- restart_list$z_vars$uo[i, , sample]
          restart_list$z_vars$v[i, , ] <- restart_list$z_vars$v[i, , sample]
          restart_list$z_vars$vo[i, , ] <- restart_list$z_vars$vo[i, , sample]
          restart_list$z_vars$xP[i, , ] <- restart_list$z_vars$xP[i, , sample]
          restart_list$z_vars$h[i, , ] <- restart_list$z_vars$h[i, , sample]
          restart_list$z_vars$ho[i, , ] <- restart_list$z_vars$ho[i, , sample]

          # zi vars
          restart_list$zi_vars$tke[i, , ] <- restart_list$zi_vars$tke[i, , sample]
          restart_list$zi_vars$zi[i, , ] <- restart_list$zi_vars$zi[i, , sample]
          restart_list$zi_vars$tkeo[i, , ] <- restart_list$zi_vars$tkeo[i, , sample]
          restart_list$zi_vars$eps[i, , ] <- restart_list$zi_vars$eps[i, , sample]
          restart_list$zi_vars$num[i, , ] <- restart_list$zi_vars$num[i, , sample]
          restart_list$zi_vars$nuh[i, , ] <- restart_list$zi_vars$nuh[i, , sample]
          restart_list$zi_vars$nus[i, , ] <- restart_list$zi_vars$nus[i, , sample]
        }
        if(model == "Simstrat") {
          restart_list$zi[i, , ] <- restart_list$zi[i, , sample]
          restart_list$u[i, , ] <- restart_list$u[i, , sample]
          restart_list$v[i, , ] <- restart_list$v[i, , sample]
          restart_list$temp[i, , ] <- restart_list$temp[i, , sample]
          restart_list$S[i, , ] <- restart_list$S[i, , sample]
          restart_list$k[i, , ] <- restart_list$k[i, , sample]
          restart_list$eps[i, , ] <- restart_list$eps[i, , sample]
          restart_list$num[i, , ] <- restart_list$num[i, , sample]
          restart_list$nuh[i, , ] <- restart_list$nuh[i, , sample]
          restart_list$seicheE[i , ] <- restart_list$seicheE[i , sample]
          restart_list$b_ice[i , ] <- restart_list$b_ice[i , sample]
          restart_list$w_ice[i , ] <- restart_list$w_ice[i , sample]
          restart_list$snow[i , ] <- restart_list$snow[i , sample]
        }

        for(s in 1:nstates){
          for(m in 1:nmembers){
            depth_index <- which(config$model_settings$modeled_depths > lake_depth[i, m])
            if(length(depth_index) > 0){
              x[i, s, depth_index, m ] <- NA
            }
          }
        }
        if(length(config$model_settings$diagnostics_names) > 0){
          for(d in 1:dim(diagnostics)[1]){
            for(m in 1:nmembers){
              depth_index <- which(config$model_settings$modeled_depths > lake_depth[i, m])
              if(length(depth_index > 0)){
                diagnostics[d,i, depth_index, m] <- NA
              }
            }
          }
        }

      }else{
        message("da_method not supported; select enkf or pf or none")
      }
    }

    #IF NO INITIAL CONDITION UNCERTAINITY THEN SET EACH ENSEMBLE MEMBER TO THE MEAN
    #AT THE INITIATION OF ThE FUTURE FORECAST
    # if(i == (hist_days + 1)){
    #   if(config$uncertainty$initial_condition == FALSE){
    #     state_means <- colMeans(x[i, ,1:nstates])
    #     for(m in 1:nmembers){
    #       x[i, m, 1:nstates]  <- state_means
    #     }
    #   }
    #   if(npars > 0){
    #     if(config$uncertainty$parameter == FALSE){
    #       par_means <- colMeans(x[i, ,(nstates + 1):(nstates + npars)])
    #       for(m in 1:nmembers){
    #         x[i, m, (nstates + 1):(nstates + npars)] <- par_means
    #       }
    #     }
    #   }
    # }

    ###################
    ## Quality Control Step
    ##################

    #Correct any negative water quality states
    if(length(states_config$state_names) > 1 & !log_wq){
      for(s in 2:nstates){
        for(k in 1:ndepths_modeled){
          index <- which(x[i, s, k, ] < 0.0)
          x[i, s, k, index] <- 0.0
        }
      }
    }

    #Correct any parameter values outside bounds
    if(npars > 0){
      for(par in 1:npars){
        low_index <- which(pars[i,par ,] < pars_config$par_lowerbound[par])
        high_index <- which(pars[i,par ,] > pars_config$par_upperbound[par])
        pars[i,par, low_index] <- pars_config$par_lowerbound[par]
        pars[i,par, high_index]  <- pars_config$par_upperbound[par]
      }
    }

    ###############

    #Print parameters to screen
    if(npars > 0){
      for(par in 1:npars){
        message(paste0(pars_config$par_names_save[par],": mean ",
                       round(mean(pars_corr[par,]),4)," sd ",
                       round(sd(pars_corr[par,]),4)))
      }
    }
  }

  if(lubridate::day(full_time[1]) < 10){
    file_name_H_day <- paste0("0",lubridate::day(full_time[1]))
  }else{
    file_name_H_day <- lubridate::day(full_time[1])
  }
  if(lubridate::day(full_time[hist_days+1]) < 10){
    file_name_H_end_day <- paste0("0",lubridate::day(full_time[hist_days+1]))
  }else{
    file_name_H_end_day <- lubridate::day(full_time[hist_days+1])
  }
  if(lubridate::month(full_time[1]) < 10){
    file_name_H_month <- paste0("0",lubridate::month(full_time[1]))
  }else{
    file_name_H_month <- lubridate::month(full_time[1])
  }
  if(lubridate::month(full_time[hist_days+1]) < 10){
    file_name_H_end_month <- paste0("0",lubridate::month(full_time[hist_days+1]))
  }else{
    file_name_H_end_month <- lubridate::month(full_time[hist_days+1])
  }

  time_of_forecast <- Sys.time()
  curr_day <- lubridate::day(time_of_forecast)
  curr_month <- lubridate::month(time_of_forecast)
  curr_year <- lubridate::year(time_of_forecast)
  curr_hour <- lubridate::hour(time_of_forecast)
  curr_minute <- lubridate::minute(time_of_forecast)
  curr_second <- round(lubridate::second(time_of_forecast),0)
  if(curr_day < 10){curr_day <- paste0("0",curr_day)}
  if(curr_month < 10){curr_month <- paste0("0",curr_month)}
  if(curr_hour < 10){curr_hour <- paste0("0",curr_hour)}
  if(curr_minute < 10){curr_minute <- paste0("0",curr_minute)}
  if(curr_second < 10){curr_second <- paste0("0",curr_second)}

  forecast_iteration_id <- paste0(curr_year,
                                  curr_month,
                                  curr_day,
                                  "T",
                                  curr_hour,
                                  curr_minute,
                                  curr_second)

  save_file_name <- paste0(config$run_config$sim_name, "_H_",
                           (lubridate::year(full_time[1])),"_",
                           file_name_H_month,"_",
                           file_name_H_day,"_",
                           (lubridate::year(full_time[hist_days+1])),"_",
                           file_name_H_end_month,"_",
                           file_name_H_end_day,"_F_",
                           forecast_days,"_",
                           forecast_iteration_id)



  if(lubridate::day(full_time[hist_days+1]) < 10){
    file_name_F_day <- paste0("0",lubridate::day(full_time[hist_days+1]))
  }else{
    file_name_F_day <- lubridate::day(full_time[hist_days+1])
  }
  if(lubridate::month(full_time[hist_days+1]) < 10){
    file_name_F_month <- paste0("0",lubridate::month(full_time[hist_days+1]))
  }else{
    file_name_F_month <- lubridate::month(full_time[hist_days+1])
  }

  if(length(full_time) >= hist_days+1){
    save_file_name_short <- paste0(config$location$site_id, "-",
                                   (lubridate::year(full_time[hist_days+1])),"-",
                                   file_name_F_month,"-",
                                   file_name_F_day,"-",
                                   config$run_config$sim_name)
  }else{
    save_file_name_short <- paste0(config$location$site_id, "-",
                                   (lubridate::year(full_time[hist_days+1])),"-",
                                   file_name_F_month,"-",
                                   file_name_F_day,"-",
                                   paste0(config$run_config$sim_name,"_spinup"))
  }

  #for(m in 1:nmembers){
  #  unlink(file.path(working_directory, m), recursive = TRUE)
  #}


  return(list(full_time = full_time,
              forecast_start_datetime = forecast_start_datetime,
              x = x,
              pars = pars,
              obs = obs,
              save_file_name_short = save_file_name_short,
              forecast_iteration_id = forecast_iteration_id,
              forecast_project_id = config$run_config$sim_name,
              time_of_forecast = time_of_forecast,
              restart_list =  restart_list,
              diagnostics = diagnostics,
              data_assimilation_flag = data_assimilation_flag,
              forecast_flag = forecast_flag,
              da_qc_flag = da_qc_flag,
              config = config,
              states_config = states_config,
              pars_config = pars_config,
              obs_config = obs_config,
              met_file_names = met_file_names))
}
