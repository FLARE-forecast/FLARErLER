#' @title Generate initial conditions for FLARE
#' @details Function to generate initial conditions from either default values in the states_config, observations (if available), or a previous run using the output as a restart file.
#' @param states_config list; list of state configurations
#' @param obs_config list; list of observation configurations
#' @param pars_config list; list of parameter configurations  (Default = NULL)
#' @param obs array; array of the observations. Required dimensions are `[nobs, time, depth]`
#' @param config list; list of configurations
#' @param historical_met_error boolean; producted by generate_glm_met_files()
#' @import ncdf4
#' @return list; list contains the initial conditions objects required by run_da_forecast()
#' @export
#' @author Quinn Thomas
#' @examples
#' \dontrun{
#'   init <- generate_initial_conditions(states_config, obs_config, pars_config, obs, config, restart_file = config$run_config$restart_file, historical_met_error = met_out$historical_met_error)
#'   }
generate_initial_conditions_ler <- function(states_config,
                                        obs_config,
                                        pars_config = NULL,
                                        obs,
                                        config,
                                        historical_met_error = FALSE){

  pars_config <- pars_config[pars_config$model == config$model_settings$model, ]

  if(is.na(config$run_config$restart_file)){

    init <- list()
    if(!is.null(pars_config) & any(pars_config$model == config$model_settings$model)){
      pars_config <- pars_config[pars_config$model == config$model_settings$model, ]
      npars <- nrow(pars_config)
    }else{
      npars <- 0
    }

    ndepths_modeled <- length(config$model_settings$modeled_depths)
    nmembers <- config$da_setup$ensemble_size
    nstates <- length(states_config$state_names)

    init$states <- array(NA, dim=c(nstates, ndepths_modeled, nmembers))
    init$pars <- array(NA, dim=c(npars, nmembers))
    if(config$model_settings$model == "GLM") {
      init$lake_depth = array(NA, dim=c(nmembers))
      init$model_internal_depths <- array(NA, dim = c(500, nmembers)) # Original FLARE
      init$the_depths <- array(NA, dim = c(ndepths_modeled, nmembers))
      init$the_sals <- array(NA, dim = c(ndepths_modeled, nmembers))
      init$snow_thickness <- array(NA, dim = c(nmembers))
      init$white_ice_thickness <- array(NA, dim = c(nmembers))
      init$blue_ice_thickness <- array(NA, dim = c(nmembers))
      init$avg_surf_temp <- array(NA, dim = c(nmembers))
      init$restart_variables <- array(NA, dim=c(17, nmembers))
    } else if(config$model_settings$model == "Simstrat") {
      init$U <- array(0, dim = c(ndepths_modeled, nmembers))
      init$V <- array(0, dim = c(ndepths_modeled, nmembers))
      init$k <- array(3e-6, dim = c(ndepths_modeled, nmembers))
      init$eps <- array(5e-10, dim = c(ndepths_modeled, nmembers))
    }

    alpha_v <- 1 - exp(-states_config$vert_decorr_length)

    q_v <- rep(NA ,ndepths_modeled)
    w <- rep(NA, ndepths_modeled)
    w_new <- rep(NA, ndepths_modeled)

    init_depth <- array(NA, dim = c(nrow(states_config),ndepths_modeled))
    for(i in 1:nrow(states_config)){
      if(!is.na(states_config$init_obs_name[i])){
        obs_index <- which(obs_config$state_names_obs == states_config$init_obs_name[i])
        init_obs <- obs[obs_index, 1, ] * (1/states_config$states_to_obs_mapping_1[i]) * states_config$init_obs_mapping[i]
        if(length(which(!is.na(init_obs))) == 0){
          init_depth[i, ] <- rep(states_config$initial_conditions[i], ndepths_modeled)
          if(states_config$init_obs_name[i] == "temp"){
            init_depth[i, ] <- approx(x = config$default_init$temp_depths, y = config$default_init$temp, xout = config$model_settings$modeled_depths, rule=2)$y
          }
        }else if(length(which(!is.na(init_obs))) == 1){
          init_depth[i, ]  <- rep(init_obs[!is.na(init_obs)], ndepths_modeled)
        }else{
          init_depth[i, ]  <- approx(x = config$model_settings$modeled_depths[!is.na(init_obs)], y = init_obs[!is.na(init_obs)], xout = config$model_settings$modeled_depths, rule=2)$y
        }
      }else{
        init_depth[i, ]  <- rep(states_config$initial_conditions[i], ndepths_modeled)
      }
    }

    for(m in 1:nmembers){
      q_v[] <- NA
      w[] <- NA
      w_new[] <- NA
      for(jj in 1:nstates){
        w[] <- rnorm(ndepths_modeled, 0, 1)
        w_new[1] <- w[1]
        q_v[1] <- states_config$initial_model_sd[jj] * w[1]
        for(kk in 2:ndepths_modeled){
          w_new[kk] <- (alpha_v[jj] * w_new[kk-1] + sqrt(1 - alpha_v[jj]^2) * w[kk])
          q_v[kk] <- w_new[kk] * states_config$initial_model_sd[jj]
          #q_v[kk] <- alpha_v * q_v[kk-1] + sqrt(1 - alpha_v^2) * states_config$initial_model_sd[jj] * w[kk]
        }

        if(config$uncertainty$initial_condition == FALSE){
          init$states[jj, , m] <- init_depth[jj, ]
        }else{
          init$states[jj, , m] <- init_depth[jj, ] + q_v
        }
        if(jj > 1){
          init$states[jj,which(init$states[jj, , m] < 0.0) , m] <- 0.0
        }
      }
    }

    if(npars > 0){
      for(par in 1:npars){
        init$pars[par, ] <- runif(n=nmembers,pars_config$par_init_lowerbound[par], pars_config$par_init_upperbound[par])
      }
    }

    if(config$model_settings$model == "GLM") {
      init$lake_depth[] <- round(config$default_init$lake_depth, 3)
      init$snow_thickness[] <- config$default_init$snow_thickness
      init$white_ice_thickness[] <- config$default_init$white_ice_thickness
      init$blue_ice_thickness[] <- config$default_init$blue_ice_thickness
      init$avg_surf_temp[] <- init$states[1 , 1, ]
      init$restart_variables[, ] <- 0.0
      init$the_sals[, ] <- config$default_init$salinity
      for(m in 1:nmembers){
        init$the_depths[1:ndepths_modeled, m] <- config$model_settings$modeled_depths
        init$model_internal_depths[1:ndepths_modeled, m] <- config$model_settings$modeled_depths
      }
    }

    aux_states_init <- list()
    if(config$model_settings$model == "GLM") {
      aux_states_init$snow_thickness <- init$snow_thickness
      aux_states_init$white_ice_thickness <- init$white_ice_thickness
      aux_states_init$blue_ice_thickness <- init$blue_ice_thickness
      aux_states_init$the_sals <- init$the_sals
      aux_states_init$the_depths <- init$the_depths
      aux_states_init$model_internal_depths <- init$model_internal_depths
      aux_states_init$lake_depth <- init$lake_depth
      aux_states_init$avg_surf_temp <- init$avg_surf_temp
      aux_states_init$mixing_vars <- init$mixing_vars
    }

    init <- list(states = init$states,
                 pars = init$pars,
                 aux_states_init = aux_states_init)

  } else {
    nc <- ncdf4::nc_open(config$run_config$restart_file)
    #forecast <- ncdf4::ncvar_get(nc, "forecast")
    t <- ncdf4::ncvar_get(nc,'time')
    local_tzone <- ncdf4::ncatt_get(nc, 0)$local_time_zone_of_simulation
    datetime <- as.POSIXct(t,
                           origin = '1970-01-01 00:00.00 UTC',
                           tz = "UTC")
    ncdf4::nc_close(nc)

    restart_index <- which(datetime == lubridate::as_datetime(config$run_config$start_datetime))

    if(length(restart_index) != 1){
      warning("start_datetime for this simulation is missing from restart file")
    }

    out <- FLARErLER:::generate_restart_initial_conditions_ler(
      restart_file = config$run_config$restart_file,
      state_names = states_config$state_names,
      par_names = pars_config$par_names_save,
      restart_index = restart_index,
      config = config)

    aux_states_init <- list()
    if(config$model_settings$model == "GLM") {
      aux_states_init$snow_thickness <- out$restart_list$snow_thickness
      aux_states_init$white_ice_thickness <- out$restart_list$white_ice_thickness
      aux_states_init$blue_ice_thickness <- out$restart_list$blue_ice_thickness
      aux_states_init$the_sals <- out$restart_list$the_sals
      aux_states_init$the_depths <- out$restart_list$the_depths
      aux_states_init$model_internal_depths <- out$restart_list$model_internal_depths
      aux_states_init$lake_depth <- out$restart_list$lake_depth
      aux_states_init$avg_surf_temp <- out$restart_list$avg_surf_temp
      aux_states_init$restart_variables <- out$restart_list$restart_variables
    } else if(config$model_settings$model == "GOTM") {

      # z vars
      aux_states_init$z_vars$z <- out$restart_list$z_vars$z
      aux_states_init$z_vars$temp <- out$restart_list$z_vars$temp
      aux_states_init$z_vars$salt <- out$restart_list$z_vars$salt
      aux_states_init$z_vars$u <- out$restart_list$z_vars$u
      aux_states_init$z_vars$uo <- out$restart_list$z_vars$uo
      aux_states_init$z_vars$v <- out$restart_list$z_vars$v
      aux_states_init$z_vars$vo <- out$restart_list$z_vars$vo
      aux_states_init$z_vars$xP <- out$restart_list$z_vars$xP
      aux_states_init$z_vars$h <- out$restart_list$z_vars$h
      aux_states_init$z_vars$ho <- out$restart_list$z_vars$ho

      # zi vars
      aux_states_init$zi_vars$tke <- out$restart_list$zi_vars$tke
      aux_states_init$zi_vars$zi <- out$restart_list$zi_vars$zi
      aux_states_init$zi_vars$tkeo <- out$restart_list$zi_vars$tkeo
      aux_states_init$zi_vars$eps <- out$restart_list$zi_vars$eps
      aux_states_init$zi_vars$num <- out$restart_list$zi_vars$num
      aux_states_init$zi_vars$nuh <- out$restart_list$zi_vars$nuh
      aux_states_init$zi_vars$nus <- out$restart_list$zi_vars$nus

    } else if(config$model_settings$model == "Simstrat") {
      aux_states_init$zi <- out$restart_list$zi
      aux_states_init$u <- out$restart_list$u
      aux_states_init$v <- out$restart_list$v
      aux_states_init$temp <- out$restart_list$temp
      aux_states_init$S <- out$restart_list$S
      aux_states_init$k <- out$restart_list$k
      aux_states_init$eps <- out$restart_list$eps
      aux_states_init$num <- out$restart_list$num
      aux_states_init$nuh <- out$restart_list$nuh
      aux_states_init$seicheE <- out$restart_list$seicheE
      aux_states_init$b_ice <- out$restart_list$b_ice
      aux_states_init$w_ice <- out$restart_list$w_ice
      aux_states_init$snow <- out$restart_list$snow
    }



    init <- list(states = out$states,
                 pars = out$pars,
                 aux_states_init = aux_states_init)

  }
  return(init)
}
