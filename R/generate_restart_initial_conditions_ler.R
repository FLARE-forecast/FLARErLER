
#' Generate initial conditions from existing output file
#'
#' @param restart_file file name of FLARE output csv
#' @param state_names names of states that are initialized
#' @param par_names (optional) names of parameters that are initialized
#' @param restart_index (optional) time index in restart file used for initialization
#' @param config list; list of configurations
#' @noRd
#' @return list of initial conditions
#'
generate_restart_initial_conditions_ler <- function(restart_file,
                                                state_names,
                                                par_names = NULL,
                                                restart_index = NULL,
                                                config){



  nc <- ncdf4::nc_open(restart_file)
  on.exit({
    ncdf4::nc_close(nc)
  })
  restart_nmembers <- length(ncdf4::ncvar_get(nc, "ensemble"))
  forecast <- ncdf4::ncvar_get(nc, "forecast")
  if(is.null(restart_index)){
    restart_index <- max(which(forecast == 0))
    if(is.na(restart_index)){
      restart_index <- length(data_assimilation)
    }
  }

  message(paste0("Using restart file with restart index of ", restart_index))

  modeled_depths <- ncdf4::ncvar_get(nc, "depth")

  if(config$model_settings$model ==  "GLM") {
    the_depths <- ncdf4::ncvar_get(nc, "depth")
    lake_depth_restart <- ncdf4::ncvar_get(nc, "lake_depth")[restart_index, ]
    snow_thickness_restart <- ncdf4::ncvar_get(nc, "snow_thickness")[restart_index, ]
    white_ice_thickness_restart <- ncdf4::ncvar_get(nc, "white_ice_thickness")[restart_index, ]
    blue_ice_thickness_restart <- ncdf4::ncvar_get(nc, "blue_ice_thickness")[restart_index, ]
    avg_surf_temp_restart <- ncdf4::ncvar_get(nc, "avg_surf_temp")[restart_index, ]
    the_sals_restart <- ncdf4::ncvar_get(nc, "salt")[restart_index, , ]

    restart_variables_restart <- ncdf4::ncvar_get(nc, "restart_variables")[ ,restart_index, ]
    model_internal_depths  <- ncdf4::ncvar_get(nc, "model_internal_depths")[restart_index, , ]

    restart_list <- list(lake_depth = lake_depth_restart,
                         the_depths = the_depths,
                         the_sals = the_sals_restart,
                         snow_thickness = snow_thickness_restart,
                         white_ice_thickness = white_ice_thickness_restart,
                         blue_ice_thickness = blue_ice_thickness_restart,
                         avg_surf_temp = avg_surf_temp_restart,
                         model_internal_depths = model_internal_depths,
                         restart_variables = restart_variables_restart)
  } else if(config$model_settings$model == "GOTM") {
    # z variables
    z <- ncdf4::ncvar_get(nc, "z")
    temp <- ncdf4::ncvar_get(nc, "temp")[restart_index, , ]
    salt <- ncdf4::ncvar_get(nc, "salt")[restart_index, , ]
    u <- ncdf4::ncvar_get(nc, "u")[restart_index, , ]
    uo <- ncdf4::ncvar_get(nc, "uo")[restart_index, , ]
    v <- ncdf4::ncvar_get(nc, "v")[restart_index, , ]
    vo <- ncdf4::ncvar_get(nc, "vo")[restart_index, , ]
    xP <- ncdf4::ncvar_get(nc, "xP")[restart_index, , ]
    h <- ncdf4::ncvar_get(nc, "h")[restart_index, , ]
    ho <- ncdf4::ncvar_get(nc, "ho")[restart_index, , ]

    # zi variables
    zi <- ncdf4::ncvar_get(nc, "zi")
    tke <- ncdf4::ncvar_get(nc, "tke")[restart_index, , ]
    tkeo <- ncdf4::ncvar_get(nc, "tkeo")[restart_index, , ]
    eps <- ncdf4::ncvar_get(nc, "eps")[restart_index, , ]
    num <- ncdf4::ncvar_get(nc, "num")[restart_index, , ]
    nuh <- ncdf4::ncvar_get(nc, "nuh")[restart_index, , ]
    nus <- ncdf4::ncvar_get(nc, "nus")[restart_index, , ]

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
  } else if(config$model_settings$model == "Simstrat") {

    zi <- ncdf4::ncvar_get(nc, "zi")
    u <- ncdf4::ncvar_get(nc, "u")[restart_index, , ]
    v <- ncdf4::ncvar_get(nc, "v")[restart_index, , ]
    temp <- ncdf4::ncvar_get(nc, "temp")[restart_index, , ]
    S <- ncdf4::ncvar_get(nc, "salt")[restart_index, , ]
    k <- ncdf4::ncvar_get(nc, "k")[restart_index, , ]
    eps <- ncdf4::ncvar_get(nc, "eps")[restart_index, , ]
    num <- ncdf4::ncvar_get(nc, "num")[restart_index, , ]
    nuh <- ncdf4::ncvar_get(nc, "nuh")[restart_index, , ]
    seicheE <- ncdf4::ncvar_get(nc, "seicheE")[restart_index, ]
    b_ice <- ncdf4::ncvar_get(nc, "b_ice")[restart_index, ]
    w_ice <- ncdf4::ncvar_get(nc, "w_ice")[restart_index, ]
    snow <- ncdf4::ncvar_get(nc, "snow")[restart_index, ]

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



  states_restart <- array(NA, dim = c(length(state_names), length(modeled_depths), restart_nmembers))
  for(i in 1:length(state_names)){
    states_restart[i, , ] <- ncdf4::ncvar_get(nc,state_names[i])[restart_index, , ]
  }

  if(!is.null(par_names)){
    pars_restart <- array(NA, dim = c(length(par_names), restart_nmembers))
    for(i in 1:length(par_names)){
      pars_restart[i, ] <- ncdf4::ncvar_get(nc, par_names[i])[restart_index, ]
    }
  }else{
    pars_restart = NULL
  }

  return(list(states = states_restart,
              pars = pars_restart,
              restart_list = restart_list)
  )
}
