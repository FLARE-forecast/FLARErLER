#' @title Download and Downscale NOAA GEFS for a single site
#' @return None
#'
#' @param site_index, index of site_list, lat_list, lon_list to be downloaded
#' @param lat_list, vector of latitudes that correspond to site codes
#' @param lon_list, vector of longitudes that correspond to site codes
#' @param site_list, vector of site codes, used in directory and file name generation
#' @param downscale, logical specifying whether to downscale from 6-hr to 1-hr
#' @param overwrite, logical stating to overwrite any existing output_file
#' @param model_name, directory name for the 6-hr forecast, this will be used in directory and file name generation
#' @param model_name_ds, directory name for the 1-hr forecast, this will be used in directory and file name generation
#' @param output_directory, directory where the model output will be save
#' @noRd
#' @importFrom LakeEnsemblR export_config
#' @importFrom gotmtools write_yaml
#' @importFrom yaml read_yaml
#'
#' @author Quinn Thomas
#'


run_model_ler <- function(model,
                      i,
                      m,
                      curr_start,
                      curr_stop,
                      par_names,
                      curr_pars,
                      working_directory,
                      par_file,
                      x_start,
                      full_time,
                      management = NULL,
                      hist_days,
                      modeled_depths,
                      ndepths_modeled,
                      curr_met_file,
                      inflow_file_name,
                      outflow_file_name,
                      output_vars,
                      diagnostics_names,
                      npars,
                      num_wq_vars,
                      nstates,
                      state_names,
                      include_wq,
                      restart,
                      restart_list){

  if(is.null(management)) {
    simulate_sss <- FALSE
  }else{
    simulate_sss <- management$simulate_sss
  }

  yml <- yaml::read_yaml(file.path(working_directory, config$model_settings$base_ler_yaml))

  # Remove all parameters - updated ones will be added
  #yml$model_parameters[[model]] <- NULL


  # Restart with GLM depths
  if(model == "GLM") {

    model_depths_start <- restart_list$model_internal_depths[i-1, , m]
    #
    model_depths_tmp <- model_depths_start[!is.na(model_depths_start)]

    # Initial temperature
    the_temps_enkf_tmp <- x_start[1, ]

    the_temps <- approx(modeled_depths, the_temps_enkf_tmp, model_depths_tmp, rule = 2)$y
    the_sals <- approx(modeled_depths, round(x_start[2, ]), model_depths_tmp, rule = 2)$y

    init_prof <- data.frame(Depth_meter = round(model_depths_tmp, 4),
                            Water_Temperature_celsius = round(the_temps, 4))
  } else if(model == "GOTM"){


    model_depths_start <- restart_list$z_vars$z[i-1, , m]

    if(length(!is.na(model_depths_start))){
      model_depths_tmp <- modeled_depths
    }else{
    model_depths_tmp <- model_depths_start[!is.na(model_depths_start)]
    }

    # Initial temperature
    the_temps_enkf_tmp <- x_start[1, ]

    the_temps <- approx(modeled_depths, the_temps_enkf_tmp, model_depths_tmp, rule = 2)$y
    the_sals <- approx(modeled_depths, round(x_start[2, ]), model_depths_tmp, rule = 2)$y

    init_prof <- data.frame(Depth_meter = round(model_depths_tmp, 4),
                            Water_Temperature_celsius = round(the_temps, 4))


  }else if(model == "Simstrat"){

    model_depths_start <- restart_list$zi[i-1, , m]
    #
    if(length(!is.na(model_depths_start))){
      model_depths_tmp <- modeled_depths
    }else{
      model_depths_tmp <- model_depths_start[!is.na(model_depths_start)]
    }


    # Initial temperature
    the_temps_enkf_tmp <- x_start[1, ]

    the_temps <- approx(modeled_depths, the_temps_enkf_tmp, model_depths_tmp, rule = 2)$y
    the_sals <- approx(modeled_depths, round(x_start[2, ]), model_depths_tmp, rule = 2)$y

    init_prof <- data.frame(Depth_meter = round(model_depths_tmp, 4),
                            Water_Temperature_celsius = round(the_temps, 4))


  }

  write.csv(init_prof, file.path(working_directory, "initial_profile.csv"),
            row.names = FALSE, quote = FALSE)


  # GLM ----
  if( model == "GLM") {

    model_depths_start = restart_list$the_depths[i-1, , m]
    lake_depth_start = restart_list$lake_depth[i-1, m]
    yml$location$depth <- lake_depth_start
    yml$location$init_depth <- lake_depth_start

    update_aed_nml_list <- list()
    update_phyto_nml_list <- list()
    update_aed_nml_names <- c()
    update_phyto_nml_names <- c()
    list_index_aed <- 1
    list_index_phyto <- 1


    diagnostics <- array(NA, dim = c(length(diagnostics_names),ndepths_modeled))

    x_star_end <- array(NA, dim =c(nstates, ndepths_modeled))

    if(npars > 0){

      unique_pars <- unique(par_names)

      for(par in 1:length(unique_pars)){

        curr_par_set <- which(par_names == unique_pars[par])
        curr_nml <- par_file[curr_par_set[1]]
        if(curr_nml == "glm3.nml"){
          yml$model_parameters$GLM[[unique_pars[par]]] <- round(curr_pars[curr_par_set], 4)
        }else if(curr_nml == "aed2.nml"){
          update_aed_nml_list[[list_index_aed]] <- round(curr_pars[curr_par_set], 4)
          update_aed_nml_names[list_index_aed] <- unique_pars[par]
          list_index_aed <- list_index_aed + 1
        }else if(curr_nml == "aed2_phyto_pars.nml"){
          update_phyto_nml_list[[list_index_phyto]] <- rep(round(curr_pars[curr_par_set],4), num_phytos)
          update_phyto_nml_names[list_index_phyto] <- unique_pars[par]
          list_index_phyto <- list_index_phyto + 1
        }
      }
    }

    if(include_wq){

      wq_init_vals <- c()
      start_index <- 2

      for(wq in 1:num_wq_vars){
        wq_enkf_tmp <- x_start[start_index + wq, ]
        wq_enkf_tmp[wq_enkf_tmp < 0] <- 0
        wq_init_vals <- c(wq_init_vals,
                          approx(modeled_depths,wq_enkf_tmp, model_depths_mid, rule = 2)$y)
      }
      update_glm_nml_list[[list_index]] <- round(wq_init_vals, 4)
      update_glm_nml_names[list_index] <- "wq_init_vals"
      list_index <- list_index + 1

      if(simulate_sss){
        if(is.na(management$specified_sss_inflow_file)){
          FLARErLER:::create_sss_input_output(x = x_start,
                                          i,
                                          m,
                                          full_time,
                                          working_directory,
                                          wq_start,
                                          management$management_input,
                                          hist_days,
                                          management$forecast_sss_on,
                                          management$sss_depth,
                                          management$use_specified_sss,
                                          state_names,
                                          modeled_depths = modeled_depths,
                                          forecast_sss_flow = management$forecast_sss_flow,
                                          forecast_sss_oxy = management$forecast_sss_oxy,
                                          salt = salt_start)
        }else{
          file.copy(file.path(working_directory, management$specified_sss_inflow_file), paste0(working_directory,"/sss_inflow.csv"))
          if(!is.na(management$specified_sss_outflow_file)){
            file.copy(file.path(working_directory, management$specified_sss_outflow_file), paste0(working_directory,"/sss_outflow.csv"))
          }
        }
      }
    }

    glm_depths_end <- rep(NA, length(init_prof$Depth_meter))

    inp_list <- list(lake_depth = round(restart_list$lake_depth[i-1, m], 4),
                     the_depths = init_prof$Depth_meter,
                     the_temps = init_prof$Water_Temperature_celsius,
                     the_sals = the_sals,
                     snow_thickness = round(restart_list$snow_thickness[i-1, m]),
                     white_ice_thickness = round(restart_list$white_ice_thickness[i-1, m], 4),
                     blue_ice_thickness = round(restart_list$blue_ice_thickness[i-1, m], 4),
                     avg_surf_temp = round(restart_list$avg_surf_temp[i-1, m], 4),
                     restart_variables = restart_list$restart_variables[, i-1, m]
                     )


    yml[["model_parameters"]][[model]][["init_profiles/wq_names"]] <- "''"
    yml[["model_parameters"]][[model]][["init_profiles/wq_init_vals"]] <- 0

  }

  # GOTM ----
  if( model == "GOTM") {

    diagnostics <- array(NA, dim = c(length(diagnostics_names), ndepths_modeled))

    x_star_end <- array(NA, dim =c(nstates, ndepths_modeled))

    if(npars > 0){

      unique_pars <- unique(par_names)

      for(par in 1:length(unique_pars)){

        curr_par_set <- which(par_names == unique_pars[par])
        curr_nml <- par_file[curr_par_set[1]]
        yml[["model_parameters"]][[model]][[unique_pars[par]]] <- signif(curr_pars[curr_par_set], 4)

      }
    }

    # if(restart) {
    #   yml$model_parameters$GOTM$`restart/load` <- TRUE
    # } else {
    #   yml$model_parameters$GOTM$`restart/load` <- FALSE
    # }

    got_deps <- abs(restart_list$z_vars$z[i-1, , m])
    # got_deps <- got_deps[order(got_deps)]

    got_temps <- approx(modeled_depths, x_start[1, ], got_deps, rule = 2)$y
    got_salts <- approx(modeled_depths, x_start[2, ], got_deps, rule = 2)$y


    inp_list <- list(z_vars = list(z = restart_list$z_vars$z[i-1, , m],
                                   temp = got_temps,
                                   salt = got_salts,
                                   u = restart_list$z_vars$u[i-1, , m],
                                   uo = restart_list$z_vars$uo[i-1, , m],
                                   v = restart_list$z_vars$v[i-1, , m],
                                   vo = restart_list$z_vars$vo[i-1, , m],
                                   xP = restart_list$z_vars$xP[i-1, , m],
                                   h = restart_list$z_vars$h[i-1, , m],
                                   ho = restart_list$z_vars$ho[i-1, , m]),
                     zi_vars = list(tke = restart_list$zi_vars$tke[i-1, , m],
                                    zi = restart_list$zi_vars$zi[i-1, , m],
                                    tkeo = restart_list$zi_vars$tkeo[i-1, , m],
                                    eps = restart_list$zi_vars$eps[i-1, , m],
                                    num = restart_list$zi_vars$num[i-1, , m],
                                    nuh = restart_list$zi_vars$nuh[i-1, , m],
                                    nus = restart_list$zi_vars$nus[i-1, , m]))
  }

  # Simstrat ----
  if(model == "Simstrat") {

    diagnostics <- array(NA, dim = c(length(diagnostics_names), ndepths_modeled))

    x_star_end <- array(NA, dim =c(nstates, ndepths_modeled))

    if(npars > 0){

      unique_pars <- unique(par_names)

      for(par in 1:length(unique_pars)){

        curr_par_set <- which(par_names == unique_pars[par])
        curr_nml <- par_file[curr_par_set[1]]
        yml$model_parameters$Simstrat[[par_names[par]]] <- signif(curr_pars[curr_par_set], 4)

      }
    }

    if(i > 1){
      the_temps_enkf_tmp <- x_start[1, ]
      sim_deps <- abs(restart_list$zi[i-1, , m])
      sim_temp <- approx(modeled_depths, x_start[1, ], sim_deps, rule = 2)$y
      sim_salts <- approx(modeled_depths, x_start[2, ], sim_deps, rule = 2)$y
    }else{
      sim_temp <- x_start[1, ]
      sim_salts <- x_start[2, ]
    }


    inp_list <- list(zi = restart_list$zi[i-1, , m],
                     u = restart_list$u[i-1, , m],
                     v = restart_list$v[i-1, , m],
                     temp = sim_temp,
                     S = sim_salts,
                     k = restart_list$k[i-1, , m],
                     eps = restart_list$eps[i-1, , m],
                     num = restart_list$num[i-1, , m],
                     nuh = restart_list$nuh[i-1, , m],
                     seicheE = restart_list$seicheE[i-1, m],
                     b_ice = restart_list$b_ice[i-1, m],
                     w_ice = restart_list$w_ice[i-1, m],
                     snow = restart_list$snow[i-1, m])

    #yml[["model_parameters"]][[model]][["Simulation/Timestep s"]] <- 300

  }

  #ALLOWS THE LOOPING THROUGH NOAA ENSEMBLES
  yml$input$init_temp_profile$file <- "initial_profile.csv"
  # yml$input$meteo$file <- paste0("../", basename(curr_met_file))
  # yml$inflows$file <- paste0("../", basename(unlist(inflow_file_name)))

  yml$time$start <- curr_start
  yml$time$stop <- curr_stop

  gotmtools::write_yaml(yml, file.path(working_directory, config$model_settings$base_ler_yaml))

  if(model == "Simstrat") {
    l_ext <- TRUE
  } else {
    l_ext <- FALSE
  }

  suppressMessages({
    LakeEnsemblR::export_config(config_file = config$model_settings$base_ler_yaml, model = model, dirs = FALSE,
                                time = TRUE, location = TRUE, output_settings = TRUE,
                                meteo = TRUE, init_cond = TRUE, extinction = l_ext,
                                inflow = FALSE, model_parameters = TRUE,
                                folder = working_directory, print = FALSE)
  })

  # Don't write restart on the first time step for GOTM & Simstrat
  if((model == "GOTM" & !all(is.na(inp_list$z_vars$u))) |
     (model == "Simstrat" & !all(is.na(inp_list$temp)))) {
    LakeEnsemblR::write_restart(folder = working_directory, model = model,
                                restart_list = inp_list)
  } else if(model == "GLM") {
    LakeEnsemblR::write_restart(folder = working_directory, model = model,
                                restart_list = inp_list)
  }



  if(model == "GLM" & include_wq){
    FLARErLER:::update_nml(update_aed_nml_list,
                       update_aed_nml_names,
                       working_directory,
                       "aed2.nml")
    FLARErLER:::update_nml(update_phyto_nml_list,
                       update_phyto_nml_names,
                       working_directory,
                       "aed2_phyto_pars.nml")
  }


  #Use GLM NML files to run GLM for a day
  # Only allow simulations without NaN values in the output to proceed.
  #Necessary due to random Nan in AED output
  pass <- FALSE
  num_reruns <- 0
  setwd(working_directory)

  if(i == 2 & m == 1){
    file.copy(from = paste0(working_directory, "/", "LakeEnsemblR.yaml"), #GLM SPECIFIC
              to = paste0(working_directory, "/", "LakeEnsemblR_initial.yaml"),
              overwrite = TRUE) #GLM SPECIFIC
  }

  # From LakeEnsemblR::run_ensemble()
  time_step <- yml$output$time_step
  time_unit <- yml$output$time_unit
  if(time_unit == "second"){
    # Needed to create out_time vector
    time_unit <- "sec"
  }
  # Create output time vector
  out_time <- data.frame(datetime = seq.POSIXt(as.POSIXct(curr_start),
                                               as.POSIXct(curr_stop),
                                               by = paste(time_step, time_unit)))
  out_vars <- yml$output$variables

  # Create output folder
  dir.create(file.path(working_directory, model, "output"), showWarnings = FALSE)

  while(!pass) {

    #Delete previous output
    old_output <- list.files(file.path(working_directory, model, "output"))
    unlink(file.path(working_directory, model, "output", old_output), recursive = TRUE)

    model_states <- FLARErLER:::run_models_ler(model = model,
                          folder = working_directory,
                          verbose = FALSE,
                          restart = restart,
                          member = m,
                          the_temps = the_temps,
                          model_depths = modeled_depths)

    fils <- list.files(file.path(working_directory, model, "output"))
    run_success <- FLARErLER:::check_model_output(folder = working_directory, model = model)


    if(length(fils) != 0 & run_success) { #&


      if(run_success){
          # nc_close(nc)

        output_vars_multi_depth <- state_names
        output_vars_no_depth <- NA

        # LakeEnsemblR Output
        ler_temp_out <- FLARErLER:::get_ler_var_all(model = model,
                                                  working_dir = working_directory,
                                                  z_out = modeled_depths,
                                                  vars_depth = output_vars_multi_depth,
                                                  vars_no_depth = output_vars_no_depth,
                                                  diagnostic_vars = diagnostics_names,
                                                  ler_yaml = yml,
                                                  run_success = run_success)

        num_model_depths <- length(ler_temp_out$depths_enkf)

        if(model == "GLM") {
          num_glm_depths <- length(ler_temp_out$depths_enkf)
          glm_temps <- rev(ler_temp_out$output[ ,1])
          glm_depths_end[1:num_glm_depths] <- ler_temp_out$depths_enkf
          glm_depths_tmp <- c(ler_temp_out$depths_enkf,ler_temp_out$lake_depth)
          glm_depths_mid <- glm_depths_tmp[1:(length(glm_depths_tmp)-1)] + diff(glm_depths_tmp)/2

          x_star_end[1, ] <- approx(glm_depths_mid,glm_temps, modeled_depths, rule = 2)$y

          glm_salt <- rev(ler_temp_out$output[ ,2])
          x_star_end[2, ] <- approx(glm_depths_mid, glm_salt, modeled_depths, rule = 2)$y

        } else if(model == "GOTM") {
          #x_star_end[1, ] <- ler_temp_out$output[ ,1]
          #x_star_end[2, ] <- ler_temp_out$output[ ,2]

          num_glm_depths <- length(ler_temp_out$depths_enkf)
          gotm_temps <-ler_temp_out$output[ ,1]
          #gotm_depths_end  <- ler_temp_out$depths_enkf
          gotm_depths_tmp <- ler_temp_out$depths_enkf

          x_star_end[1, ] <- approx(gotm_depths_tmp,gotm_temps, modeled_depths, rule = 2)$y

          gotm_salt <- ler_temp_out$output[ ,2]
          x_star_end[2, ] <- approx(gotm_depths_tmp, gotm_salt, modeled_depths, rule = 2)$y


        } else if (model == "Simstrat") {

          num_glm_depths <- length(ler_temp_out$depths_enkf)
          simstrat_temps <-ler_temp_out$output[ ,1]
          simstrat_depths_tmp <- ler_temp_out$depths_enkf
          simstrat_depths_mid <-  simstrat_depths_tmp[1:(length( simstrat_depths_tmp)-1)] + diff( simstrat_depths_tmp)/2

          x_star_end[1, ] <- approx(simstrat_depths_tmp,simstrat_temps, modeled_depths, rule = 2)$y

          simstrat_salt <- ler_temp_out$output[ ,2]
          x_star_end[2, ] <- approx(simstrat_depths_tmp, simstrat_salt, modeled_depths, rule = 2)$y

        }

        model_depths_end <- rep(NA, 500)
        model_depths_end[1:num_model_depths] <- ler_temp_out$depths_enkf

        #
        # model_depths_tmp <- c(ler_temp_out$depths_enkf, ler_temp_out$lake_depth)
        #
        # model_depths_mid <- model_depths_tmp[1:(length(model_depths_tmp)-1)] + diff(model_depths_tmp)/2



#
#         x_star_end[1:ndepths_modeled] <- approx(model_depths_mid, temps,
#                                                 modeled_depths, rule = 2)$y

        # salt_end <- approx(model_depths_mid, ler_temp_out$salt,
        #                    modeled_depths, rule = 2)$y

        if(include_wq){
          for(wq in 1:num_wq_vars){
            glm_wq <-  rev(ler_temp_out$output[ ,1+wq])
            x_star_end[wq_start[1 + wq]:wq_end[1 + wq]] <- approx(model_depths_mid,glm_wq, modeled_depths, rule = 2)$y
          }
        }

        if(length(diagnostics_names) > 0){
          for(wq in 1:length(diagnostics_names)){
            glm_wq <-  rev(ler_temp_out$diagnostics_output[ , wq])
            diagnostics[wq , ] <- wq # approx(model_depths_mid, glm_wq, modeled_depths, rule = 2)$y
          }
        } else {
          diagnostics <- rep(NA, length(modeled_depths))
        }

        if(length(which(is.na(x_star_end))) == 0){
          pass = TRUE
        }else{
          message("Re-running model setup in ", working_directory, " due to errors in model output.")
          num_reruns <- num_reruns + 1
        }
      }else{
        message("Re-running model setup in ", working_directory, " due to errors in model output.")
        num_reruns <- num_reruns + 1
      }
    }else{
      message("Re-running model setup in ", working_directory, " due to errors in model output.")
      num_reruns <- num_reruns + 1
    }
    if(num_reruns > 100){
      stop("Too many re-runs (> 100) due to NaN values in output")
    }
  }

  return(list(x_star_end  = x_star_end,
                lake_depth_end  = ler_temp_out$lake_depth,
                restart_vars = ler_temp_out$restart_vars,
                diagnostics_end  = diagnostics,
                model_internal_depths = model_depths_end,
                curr_pars = curr_pars))
}
