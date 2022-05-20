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
#' da_forecast_output <- FLARErLER::run_da_forecast(states_init = init$states, pars_init = init$pars, aux_states_init = init$aux_states_init, obs = obs, obs_sd = obs_config$obs_sd, model_sd = model_sd, working_directory = config$file_path$execute_directory, met_file_names = met_file_names, inflow_file_names = inflow_file_names, outflow_file_names = outflow_file_names, config = config, pars_config = pars_config, states_config = states_config, obs_config = obs_config)
#' }
#' @noRd
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
                                debug = FALSE){

  if(length(states_config$state_names) > 1) {
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
    pars_config <- pars_config[pars_config$model == model, ]
    npars <- nrow(pars_config)
    par_names <- pars_config$par_names
    par_file <- pars_config$par_file
  }else{
    npars <- 0
    par_names <- NA
    par_file <- NA
  }

  x_init <- array(NA, dim = c(nmembers, nstates * ndepths_modeled + npars))
  for(m in 1:nmembers){
    if(nstates > 1){
      x_init[m,1:(nstates * ndepths_modeled)] <- c(aperm(states_init[, ,m ], perm = c(2,1)))
    }else{
      x_init[m,1:(nstates * ndepths_modeled)] <- states_init[1, ,m]
    }
    if(!is.null(pars_init) | npars == 0){
      x_init[m,(nstates * ndepths_modeled + 1):(nstates * ndepths_modeled + npars)] <- pars_init[, m]
    }
  }

  psi <- rep(NA, length(obs_sd) * ndepths_modeled)
  index <- 0
  for(i in 1:length(obs_sd)){
    for(j in 1:ndepths_modeled){
      index <- index + 1
      psi[index] <- obs_sd[i]
    }
  }

  wq_start <- rep(NA, nstates)
  wq_end <- rep(NA, nstates)
  for(wq in 1:nstates){
    if(wq == 1){
      wq_start[wq] <- 1
      wq_end[wq] <- ndepths_modeled
    }else{
      wq_start[wq] <- wq_end[wq-1]+1
      wq_end[wq] <- wq_end[wq-1] + (ndepths_modeled)
    }
  }

  states_config$wq_start <- wq_start
  states_config$wq_end <- wq_end

  FLARErLER:::check_enkf_inputs(states_init,
                            pars_init,
                            obs,
                            psi,
                            model_sd,
                            config,
                            pars_config,
                            states_config,
                            obs_config)

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
  save_filenames <- FLARErLER:::get_savefile_name(full_time, hist_days, forecast_days, config)

  nstates <- dim(x_init)[2] -  npars
  nsteps <- length(full_time)
  nmembers <- dim(x_init)[1]
  n_met_members <- length(met_file_names)
  ndepths_modeled <- length(config$model_settings$modeled_depths)

  data_assimilation_flag <- rep(NA, nsteps)
  forecast_flag <- rep(NA, nsteps)
  da_qc_flag <- rep(NA, nsteps)

  x <- array(NA, dim=c(nsteps, nmembers, nstates + npars))

  x[1, ,]  <- x_init

  q_v <- rep(NA, ndepths_modeled)
  w <- rep(NA, ndepths_modeled)
  w_new <- rep(NA, ndepths_modeled)

  alpha_v <- 1 - exp(-states_config$vert_decorr_length)


  output_vars <- states_config$state_names


  if(config$include_wq){
    num_wq_vars <- length(states_config$wq_start) - 1
  }else{
    num_wq_vars <- 0
  }

  if(length(config$diagnostics_names) > 0){
    diagnostics <- array(NA, dim=c(length(config$diagnostics_names), nsteps, ndepths_modeled, nmembers))
  }else{
    diagnostics <- array(NA, dim=c(1, nsteps, ndepths_modeled, nmembers))
  }

  num_phytos <- length(which(stringr::str_detect(states_config$state_names,"PHY_") & !stringr::str_detect(states_config$state_names,"_IP") & !stringr::str_detect(states_config$state_names,"_IN")))

  full_time_char <- strftime(full_time,
                             format="%Y-%m-%d %H:%M",
                             tz = "UTC")

  x_prior <- array(NA, dim = c(nsteps, nmembers, nstates + npars))

  # if(!is.null(ncol(inflow_file_names))) {
  #   inflow_file_names <- as.matrix(inflow_file_names)
  #   outflow_file_names <- as.matrix(outflow_file_names)
  # }else{
  #   inflow_file_names <- NULL
  #   outflow_file_names <- NULL
  # }

  file.copy(from = file.path(config$file_path$configuration_directory, config$model_settings$ler_bathymetry_file),
            to = file.path(working_directory, config$model_settings$ler_bathymetry_file), overwrite = TRUE)

  config$model_settings$ncore <- min(c(config$model_settings$ncore, parallel::detectCores()))
  # if(config$model_settings$ncore == 1) {
  #   if(!dir.exists(file.path(working_directory, "1"))) {
  #     dir.create(file.path(working_directory, "1"), showWarnings = FALSE)
  #   } else {
  #     unlink(file.path(working_directory, "1"), recursive = TRUE)
  #     dir.create(file.path(working_directory, "1"), showWarnings = FALSE)
  #   }
  #   FLARErLER:::set_up_model_ler(model,
  #                            config,
  #                            working_directory = working_directory,
  #                            state_names = states_config$state_names,
  #                            inflow_file_names = inflow_file_names,
  #                            outflow_file_names = outflow_file_names,
  #                            member = "1")
  # } else {
    lapply(1:nmembers, function(m) {
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
  # }

  # model_internal_depths
  # lake_depth <- array(NA, dim = c(nsteps, nmembers))
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
    yaml <- gotmtools::read_yaml(file.path(working_directory, "1", "LakeEnsemblR.yaml"))
    nz <- round(yaml$location$depth / yaml$output$depths)
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
      # temp[1, , ] <- aux_states_init$z_vars$temp
      salt[1, , ] <- aux_states_init$z_vars$salt
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

    restart_list <- list(z_vars = list(z = z,
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

      zi[1, , ] <- aux_states_init$zi
      u[1, , ] <- aux_states_init$u
      v[1, , ] <- aux_states_init$v
      # temp[1, , ] <- aux_states_init$temp
      S[1, , ] <- aux_states_init$S
      k[1, , ] <- aux_states_init$k
      eps[1, , ] <- aux_states_init$eps
      num[1, , ] <- aux_states_init$num
      nuh[1, , ] <- aux_states_init$nuh

      seicheE[1, ] <- aux_states_init$seicheE
      b_ice[1, ] <- aux_states_init$b_ice
      w_ice[1, ] <- aux_states_init$w_ice
      snow[1, ] <- aux_states_init$snow
    }

    restart_list <- list(zi = zi,
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

  ###START EnKF
  for(i in start_step:nsteps) {

    curr_start <- strftime(full_time[i - 1],
                           format="%Y-%m-%d %H:%M:%S",
                           tz = config$local_tzone)
    curr_stop <- strftime(full_time[i],
                          format="%Y-%m-%d %H:%M:%S",
                          tz = config$local_tzone)

    message(paste0("Running time step ", i-1, "/", (nsteps - 1), " : ",
                   curr_start, " - ",
                   curr_stop, " [", Sys.time(), "]"))

    setwd(working_directory)

    met_index <- rep(1:length(met_file_names), times = nmembers)
    inflow_outflow_index <- rep(1:length(inflow_file_names), times = nmembers)

    #Create array to hold GLM predictions for each ensemble
    x_star <- array(NA, dim = c(nmembers, nstates))
    x_corr <- array(NA, dim = c(nmembers, nstates))

    #Matrix to store calculated ensemble specific deviations and innovations
    dit <- array(NA, dim = c(nmembers, nstates))

    if(npars > 0){
      pars_corr <-  array(NA, dim = c(nmembers, npars))
      dit_pars<- array(NA, dim = c(nmembers, npars))
    }

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
    if(i == start_step) {
      if(machine == "windows") {
        cl <- parallel::makeCluster(config$model_settings$ncore, setup_strategy = "sequential")
        parallel::clusterEvalQ(cl, library(FLARErLER))
      } else {
        cl <- parallel::makeCluster(config$model_settings$ncore, setup_strategy = "sequential")
      }

      parallel::clusterExport(cl, varlist = list("working_directory", "met_file_names", "met_index",
                                                 "par_fit_method", "da_method", "nstates", "npars",
                                                 "pars_config", "inflow_file_names", "inflow_outflow_index",
                                                 "outflow_file_names", "i", "curr_start",
                                                 "curr_stop", "par_names", "par_file",
                                                 "num_phytos", "full_time", "management",
                                                 "hist_days", "config", "states_config",
                                                 "ndepths_modeled", "output_vars", "num_wq_vars", "model", "restart"),
                              envir = environment())
    }

	  # Variables that need to be exported at each timestep
    parallel::clusterExport(cl, varlist = list("x", "restart_list"),
                            envir = environment())

    #If i == 1 then assimilate the first time step without running the process
    #model (i.e., use yesterday's forecast of today as initial conditions and
    #assimilate new observations)
    if(i > 1) {

      # Start loop through ensemble members

      out <- tryCatch({

        parallel::parLapply(cl, 1:nmembers, function(m) {
          # lapply(1:nmembers, function(m) { # Commented out for debugging
          # print(m)

          curr_met_file <- met_file_names[met_index[m]]

          if(npars > 0){
            if(par_fit_method == "inflate" & da_method == "enkf"){
              curr_pars <- x[i - 1, m , (nstates+1):(nstates+ npars)]
            }else if(par_fit_method == "perturb"){
              if(i > (hist_days + 1)){
                curr_pars <- x[i - 1, m , (nstates+1):(nstates+ npars)] + rnorm(npars, mean = rep(0, npars), sd = pars_config$perturb_par)
              }else{
                curr_pars <- x[i - 1, m , (nstates+1):(nstates+ npars)]
              }
            }else{
              message("parameter fitting method not supported.  inflate or perturb are supported. only inflate is supported for enkf")
            }
          }

          if(length(inflow_file_names) > 0){
            inflow_file_name <- inflow_file_names[inflow_outflow_index[m]]
            outflow_file_name <- outflow_file_names[inflow_outflow_index[m]]
          }else{
            inflow_file_name <- NULL
            outflow_file_name <- NULL
          }

          if(config$model_settings$ncore == 1) {
            wdir <- file.path(working_directory, 1)
          } else {
            wdir <- file.path(working_directory, m)
          }

          # model
          # i
          # m
          # curr_start
          # curr_stop
          # par_names
          # curr_pars
          # working_directory = wdir
          # par_file
          # num_phytos
          # x_start = x[i-1, m, ]
          # full_time
          # wq_start = states_config$wq_start
          # wq_end = states_config$wq_end
          # management = management
          # hist_days
          # modeled_depths = config$model_settings$modeled_depths
          # ndepths_modeled
          # curr_met_file
          # inflow_file_name = inflow_file_name
          # outflow_file_name = outflow_file_name
          # output_vars = output_vars
          # diagnostics_names = config$output_settings$diagnostics_names
          # npars
          # num_wq_vars
          # nstates
          # state_names = states_config$state_names
          # include_wq = config$include_wq
          # restart = restart
          # restart_list = restart_list

          out <- FLARErLER:::run_model_ler(model,
                                        i,
                                        m,
                                        curr_start,
                                        curr_stop,
                                        par_names,
                                        curr_pars,
                                        working_directory = wdir,
                                        par_file,
                                        num_phytos,
                                        x_start = x[i-1, m, ],
                                        full_time,
                                        wq_start = states_config$wq_start,
                                        wq_end = states_config$wq_end,
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
          }
        )

      },
      error = function(e) {

        message("Forecast failed on time step ", i-1, "/", (nsteps - 1), " : ",
                curr_start, " - ",
                curr_stop)
        if(i == start_step) {
          stop("Failed on the first timestep.")
        }
        message("Returning output from ", curr_start)


        # Subsetting of times
        # full_time <- full_time[1:(i-1)]
        #
        # if(model == "GLM") {
        #   restart_list <- list(lake_depth = restart_list$lake_depth,
        #                        model_internal_depths = restart_list$model_internal_depths[1:(i-1), , ],
        #                        the_depths = restart_list$the_depths[1:(i-1), , ],
        #                        the_sals = restart_list$the_sals[1:(i-1), , ],
        #                        snow_thickness = restart_list$snow_thickness[1:(i-1), ],
        #                        white_ice_thickness = restart_list$white_ice_thickness[1:(i-1), ],
        #                        blue_ice_thickness = restart_list$blue_ice_thickness[1:(i-1), ],
        #                        avg_surf_temp = restart_list$avg_surf_temp[1:(i-1), ],
        #                        restart_variables = restart_list$restart_variables[, 1:(i-1), ])
        # } else if(model == "GOTM") {
        #   restart_list <- list(z_vars = list(z = restart_list$z_vars$z[1:(i-1), , ],
        #                                      temp = restart_list$z_vars$temp[1:(i-1), , ],
        #                                      salt = restart_list$z_vars$salt[1:(i-1), , ],
        #                                      u = restart_list$z_vars$u[1:(i-1), , ],
        #                                      uo = restart_list$z_vars$uo[1:(i-1), , ],
        #                                      v = restart_list$z_vars$v[1:(i-1), , ],
        #                                      vo = restart_list$z_vars$vo[1:(i-1), , ],
        #                                      xP = restart_list$z_vars$xP[1:(i-1), , ],
        #                                      h = restart_list$z_vars$h[1:(i-1), , ],
        #                                      ho = restart_list$z_vars$ho[1:(i-1), , ]),
        #                        zi_vars = list(tke = restart_list$zi_vars$tke[1:(i-1), , ],
        #                                       zi = restart_list$zi_vars$zi[1:(i-1), , ],
        #                                       tkeo = restart_list$zi_vars$tkeo[1:(i-1), , ],
        #                                       eps = restart_list$zi_vars$eps[1:(i-1), , ],
        #                                       num = restart_list$zi_vars$num[1:(i-1), , ],
        #                                       nuh = restart_list$zi_vars$nuh[1:(i-1), , ],
        #                                       nus = restart_list$zi_vars$nus[1:(i-1), , ]))
        # } else if(model == "Simstrat") {
        #   restart_list <- list(zi = restart_list$zi[1:(i-1), , ],
        #                        u = restart_list$u[1:(i-1), , ],
        #                        v = restart_list$v[1:(i-1), , ],
        #                        temp = restart_list$temp[1:(i-1), , ],
        #                        S = restart_list$S[1:(i-1), , ],
        #                        k = restart_list$k[1:(i-1), , ],
        #                        eps = restart_list$eps[1:(i-1), , ],
        #                        num = restart_list$num[1:(i-1), , ],
        #                        nuh = restart_list$nuh[1:(i-1), , ],
        #                        seicheE = restart_list$seicheE[1:(i-1), ])
        # }


        list(full_time = full_time,
            forecast_start_datetime = forecast_start_datetime,
            x = x,
            obs = obs,
            save_file_name = save_filenames$save_file_name,
            save_file_name_short = save_filenames$save_file_name_short,
            forecast_iteration_id = save_filenames$forecast_iteration_id,
            forecast_project_id = config$run_config$sim_name,
            time_of_forecast = save_filenames$time_of_forecast,
            restart_list =  restart_list,
            diagnostics = diagnostics,
            data_assimilation_flag = data_assimilation_flag,
            forecast_flag = forecast_flag,
            da_qc_flag = da_qc_flag,
            config = config,
            states_config = states_config,
            pars_config = pars_config,
            obs_config = obs_config,
            met_file_names = met_file_names)
      })

      if("data_assimilation_flag" %in% names(out)) {
        return(out)
      }



	  # Loop through output and assign to matrix
	  for(m in 1:nmembers) {

	    x_star[m, ] <- out[[m]]$x_star_end
	    # lake_depth[i ,m ] <- out[[m]]$lake_depth_end
	    # snow_ice_thickness[,i ,m] <- out[[m]]$snow_ice_thickness_end
	    if(length(config$diagnostics_names) > 0) {
	      diagnostics[, i, , m] <- out[[m]]$diagnostics_end
	    }
	    # model_internal_depths[i, ,m] <- out[[m]]$model_internal_depths
	    # salt[i, , m]  <- out[[m]]$salt_end
	    if(model == "GLM") {
	      restart_list$lake_depth[i, m] <- out[[m]]$lake_depth_end
	      restart_list$model_internal_depths[i, , m] <- out[[m]]$model_internal_depths
	      # restart_list$the_depths[i, , m] <- out[[m]]$the_depths
	      restart_list$the_sals[i, , m] <- approx(out[[m]]$restart_vars$the_depths, out[[m]]$restart_vars$the_sals, config$model_settings$modeled_depths, rule = 2)$y
	      restart_list$snow_thickness[i, m] <- out[[m]]$restart_vars$snow_thickness
	      restart_list$white_ice_thickness[i, m] <- out[[m]]$restart_vars$white_ice_thickness
	      restart_list$blue_ice_thickness[i, m] <- out[[m]]$restart_vars$blue_ice_thickness
	      restart_list$avg_surf_temp[i , m] <- out[[m]]$restart_vars$avg_surf_temp
	      restart_list$restart_variables[, i, m] <- out[[m]]$restart_vars$restart_variables
	    } else if(model == "GOTM") {
	      # z vars
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
	      w_new[1] <- w[1]
	      q_v[1] <- model_sd[jj, 1] * w_new[1]
	      for(kk in 2:ndepths_modeled){
	        #q_v[kk] <- alpha_v * q_v[kk-1] + sqrt(1 - alpha_v^2) * model_sd[jj, kk] * w[kk]

	        w_new[kk] <- (alpha_v[jj] * w_new[kk-1] + sqrt(1 - alpha_v[jj]^2) * w[kk])
	        q_v[kk] <- w_new[kk] * model_sd[jj, kk]
	      }

	      x_corr[m, (((jj-1)*ndepths_modeled)+1):(jj*ndepths_modeled)] <-
	        x_star[m, (((jj-1)*ndepths_modeled)+1):(jj*ndepths_modeled)] + q_v
	    }

	  } # END ENSEMBLE LOOP

	  # Close clusters at the last time step
	  if(i == nsteps) {
	    parallel::stopCluster(cl)
	  }

	  #Correct any negative water quality states
	  if(config$include_wq & config$da_setup$no_negative_states){
	    for(m in 1:nmembers){
	      index <- which(x_corr[m,] < 0.0)
	      x_corr[m, index[which(index <= states_config$wq_end[num_wq_vars + 1] & index >= states_config$wq_start[2])]] <- 0.0
	    }
	  }

    if(npars > 0) {
      pars_corr <- x[i - 1, , (nstates + 1):(nstates+ npars)]
      if(npars == 1) {
        pars_corr <- matrix(pars_corr,nrow = length(pars_corr),ncol = 1)
      }
      pars_star <- pars_corr
      }
	  } else {
      x_star <- x[i, ,1:nstates]
      x_corr <- x_star
      if(npars > 0){
        pars_corr <- x[i, , (nstates + 1):(nstates+ npars)]
        if(npars == 1){
          pars_corr <- matrix(pars_corr,nrow = length(pars_corr),ncol = 1)
        }
        pars_star <- pars_corr
      }
    }

    if(npars > 0){
      x_prior[i, , ] <- cbind(x_corr, pars_corr)
    }else{
      x_prior[i, , ] <- x_corr
    }

    if(dim(obs)[1] > 1){
      z_index <- which(!is.na(c(aperm(obs[,i , ], perm = c(2,1)))))
    }else{
      z_index <- which(!is.na(c(obs[1,i , ])))
    }

    #if no observations at a time step then just propogate model uncertainity

    if(length(z_index) == 0){

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

        x[i, , ] <- cbind(x_corr, pars_star)

        if(config$uncertainty$process == FALSE & i > (hist_days + 1)){
          #don't add process noise if process uncertainty is false (x_star doesn't have noise)
          #don't add the noise to parameters in future forecast mode ()
          x[i, , ] <- cbind(x_star, pars_star)
        }

        if(i == (hist_days + 1) & config$uncertainty$initial_condition == FALSE){
          for(m in 1:nmembers){
            x[i, m, ] <- c(colMeans(x_star), pars_star[m, ])
          }
        }

      }else{
        x[i, , ] <- cbind(x_corr)

        if(config$uncertainty$process == FALSE & i > (hist_days + 1)){
          x[i, , ] <- x_star
        }

      }
    } else {

      data_assimilation_flag[i] <- 1
      forecast_flag[i] <- 0
      da_qc_flag[i] <- 0

      #if observation then calucate Kalman adjustment
      if(dim(obs)[1] > 1){
        zt <- c(aperm(obs[,i , ], perm = c(2,1)))
      }else{
        zt <- c(obs[1,i , ])
      }
      zt <- zt[which(!is.na(zt))]

      if(da_method == "enkf") {

        #Assign which states have obs in the time step
        h <- matrix(0, nrow = length(obs_sd) * ndepths_modeled, ncol = nstates)

        index <- 0
        for(k in 1:((nstates/ndepths_modeled))){
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

        #Extract the data uncertainity for the data



        #types present during the time-step



        curr_psi <- psi[z_index]  ^ 2

        if(length(z_index) > 1){
          psi_t <- diag(curr_psi)
        }else{
          #Special case where there is only one data
          #type during the time-step
          psi_t <- curr_psi
        }

        d_mat <- t(mvtnorm::rmvnorm(n = nmembers, mean = zt, sigma=as.matrix(psi_t)))


        #Set any negative observations of water quality variables to zero
        d_mat[which(z_index > length(config$modeled_depths) & d_mat < 0.0)] <- 0.0






        #Ensemble mean
        ens_mean <- apply(x_corr[,], 2, mean)





        if(npars > 0){
          par_mean <- apply(pars_corr, 2, mean)
          if(par_fit_method == "inflate"){
            for(m in 1:nmembers){
              pars_corr[m, ] <- pars_config$inflat_pars * (pars_corr[m,] - par_mean) + par_mean
            }
            par_mean <- apply(pars_corr, 2, mean)
          }
        }


        #Loop through ensemble members
        for(m in 1:nmembers){
          #  #Ensemble specific deviation
          dit[m, ] <- x_corr[m, ] - ens_mean
          if(npars > 0){
            dit_pars[m, ] <- pars_corr[m, ] - par_mean
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

        if(!is.na(config$da_setup$localization_distance)){
          p_t <- localization(p_t,
                              nstates,
                              modeled_depths = config$model_settings$modeled_depths,
                              localization_distance = config$da_setup$localization_distance)
        }
        #Kalman gain
        k_t <- p_t %*% t(h) %*% solve(h %*% p_t %*% t(h) + psi_t, tol = 1e-17)
        if(npars > 0){
          k_t_pars <- p_t_pars %*% t(h) %*% solve(h %*% p_t %*% t(h) + psi_t, tol = 1e-17)
        }

        #Update states array (transposes are necessary to convert
        #between the dims here and the dims in the EnKF formulations)
        update_increment <-  k_t %*% (d_mat - h %*% t(x_corr))
        x[i, , 1:nstates] <- t(t(x_corr) + update_increment)
        if(npars > 0){
          x[i, , (nstates+1):(nstates+npars)] <- t(t(pars_corr) +
                                                     k_t_pars %*% (d_mat - h %*% t(x_corr)))
        }

      } else if(da_method == "pf"){

        LL <- rep(NA, length(nmembers))
        for(m in 1:nmembers){
          LL[m] <- sum(dnorm(zt, mean = x_corr[m,z_index], sd = psi[z_index], log = TRUE))
        }

        sample <- sample.int(nmembers, replace = TRUE, prob = exp(LL))

        x[i, , ] <- cbind(x_corr, pars_star)[sample, ]

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
        if(length(config$diagnostics_names) > 0){
          diagnostics[ ,i, , ] <- diagnostics[ ,i, ,sample]
        }

      } else {
        message("da_method not supported; select enkf or pf")
      }
    }

    #IF NO INITIAL CONDITION UNCERTAINITY THEN SET EACH ENSEMBLE MEMBER TO THE MEAN
    #AT THE INITIATION OF ThE FUTURE FORECAST
    if(i == (hist_days + 1)) {

      if(config$uncertainty$initial_condition == FALSE){
        state_means <- colMeans(x[i, ,1:nstates])
        for(m in 1:nmembers){
          x[i, m, 1:nstates]  <- state_means
        }
      }
      if(npars > 0){
        if(config$uncertainty$parameter == FALSE){
          par_means <- colMeans(x[i, ,(nstates + 1):(nstates + npars)])
          for(m in 1:nmembers){
            x[i, m, (nstates + 1):(nstates + npars)] <- par_means
          }
        }
      }
    }

    ###################
    ## Quality Control Step
    ##################

    #Correct any negative water quality states
    if(config$include_wq & config$da_setup$no_negative_states){
      for(m in 1:nmembers){
        index <- which(x[i,m,] < 0.0)
        x[i, m, index[which(index <= states_config$wq_end[num_wq_vars + 1] & index >= states_config$wq_start[2])]] <- 0.0
      }
    }

    #Correct any parameter values outside bounds
    if(npars > 0){
      for(par in 1:npars){
        low_index <- which(x[i, ,nstates + par] < pars_config$par_lowerbound[par])
        high_index <- which(x[i, ,nstates + par] > pars_config$par_upperbound[par])
        x[i,low_index ,nstates + par] <- pars_config$par_lowerbound[par]
        x[i,high_index ,nstates + par] <- pars_config$par_upperbound[par]
      }
    }

    ###############

    if(debug) {

      da_forecast_output <- list(full_time = full_time,
                                 forecast_start_datetime = forecast_start_datetime,
                                 x = x,
                                 obs = obs,
                                 save_file_name = save_filenames$save_file_name,
                                 save_file_name_short = save_filenames$save_file_name_short,
                                 forecast_iteration_id = save_filenames$forecast_iteration_id,
                                 forecast_project_id = config$run_config$sim_name,
                                 time_of_forecast = save_filenames$time_of_forecast,
                                 restart_list =  restart_list,
                                 diagnostics = diagnostics,
                                 data_assimilation_flag = data_assimilation_flag,
                                 forecast_flag = forecast_flag,
                                 da_qc_flag = da_qc_flag,
                                 config = config,
                                 states_config = states_config,
                                 pars_config = pars_config,
                                 obs_config = obs_config,
                                 met_file_names = met_file_names)
      if(i > 2) {
        add <- TRUE
      } else {
        add <- FALSE
      }

       FLARErLER::write_forecast_netcdf(da_forecast_output = da_forecast_output,
                                     forecast_output_directory = config$file_path$forecast_output_directory,
                                     add = add)

       # Fix weird restart_list fill values -  ncdf4 seems to re-write them to missvals
       if(config$model_settings$model ==  "GLM") {
         restart_list$lake_depth[restart_list$lake_depth == -99] <- NA
         restart_list$snow_thickness[restart_list$snow_thickness == -99] <- NA
         restart_list$white_ice_thickness[restart_list$white_ice_thickness == -99] <- NA
         restart_list$blue_ice_thickness[restart_list$blue_ice_thickness == -99] <- NA
         restart_list$avg_surf_temp[restart_list$avg_surf_temp == -99] <- NA
         restart_list$model_internal_depths[restart_list$model_internal_depths > 1e31] <- NA
         restart_list$the_sals[restart_list$the_sals > 1e31] <- NA
         restart_list$restart_variables[restart_list$restart_variables > 1e31] <- NA
       }


    }

    #Print parameters to screen
    if(npars > 0){
      for(par in 1:npars){
        message(paste0(pars_config$par_names_save[par],": mean ",
                       signif(mean(pars_corr[,par]),4)," sd ",
                       signif(sd(pars_corr[,par]),4)))
      }
    }

    message(paste0("surface_temp: mean ", signif(mean(x[i, , 1])), " sd ", signif(sd(x[i, , 1]))))
    message(paste0("bottom_temp: mean ", signif(mean(x[i, , ndepths_modeled]), 4), " sd ", signif(sd(x[i, , ndepths_modeled], 4))))

  }

  #for(m in 1:nmembers){
  #  unlink(file.path(working_directory, m), recursive = TRUE)
  #}


  return(list(full_time = full_time,
              forecast_start_datetime = forecast_start_datetime,
              x = x,
              obs = obs,
              save_file_name = save_filenames$save_file_name,
              save_file_name_short = save_filenames$save_file_name_short,
              forecast_iteration_id = save_filenames$forecast_iteration_id,
              forecast_project_id = config$run_config$sim_name,
              time_of_forecast = save_filenames$time_of_forecast,
              restart_list =  restart_list,
              # snow_ice_thickness = snow_ice_thickness,
              # lake_depth = lake_depth,
              # salt = salt,
              # model_internal_depths = model_internal_depths,
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
