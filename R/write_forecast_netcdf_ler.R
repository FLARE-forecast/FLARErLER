##' @title Generate netcdf output file
##' @details Function generates a netcdf file from the object that is returned by run_da_forecast()
##' @param da_forecast_output list; object that is returned by run_da_forecast()
##' @param forecast_output_directory string; full path of directory where the netcdf file will be written
##' @param use_short_filename use shortened file name; this results in less informatoin in the file name and potentially overwriting existing files
##' @param add logical; Add output to netCDF file
##' @return None
##' @export
##' @import ncdf4
##' @import ggplot2
##' @importFrom lubridate with_tz
##' @author Quinn Thomas
##' @examples
##' \dontrun{
##' write_forecast_netcdf(da_forecast_output = da_forecast_output, forecast_output_directory = config$file_path$forecast_output_directory, use_short_filename = TRUE)
##' }
##'

write_forecast_netcdf_ler <- function(da_forecast_output,
                                  forecast_output_directory,
                                  use_short_filename = TRUE,
                                  add = FALSE){

  dir.create(forecast_output_directory, recursive = TRUE, showWarnings = FALSE)

  x <- da_forecast_output$x
  pars <- da_forecast_output$pars
  data_assimilation_flag <- da_forecast_output$data_assimilation_flag
  forecast_flag <- da_forecast_output$forecast_flag
  da_qc_flag <- da_forecast_output$da_qc_flag
  full_time <- da_forecast_output$full_time
  forecast_start_datetime <- da_forecast_output$forecast_start_datetime
  config <- da_forecast_output$config
  if(config$model_settings$model ==  "GLM") {
    lake_depth <- da_forecast_output$restart_list$lake_depth
    snow_thickness <- da_forecast_output$restart_list$snow_thickness
    white_ice_thickness <- da_forecast_output$restart_list$white_ice_thickness
    blue_ice_thickness <- da_forecast_output$restart_list$blue_ice_thickness

    snow_ice_thickness <- array(c(snow_thickness, white_ice_thickness, blue_ice_thickness),
                                dim = c(3, dim(snow_thickness)[1], dim(snow_thickness)[2]))

    avg_surf_temp <- da_forecast_output$restart_list$avg_surf_temp
    restart_variables <- da_forecast_output$restart_list$restart_variables
    model_internal_depths <- da_forecast_output$restart_list$model_internal_depths

  } else if(config$model_settings$model == "GOTM") {
    z <- da_forecast_output$restart_list$z_vars$z
    u <- da_forecast_output$restart_list$z_vars$u
    uo <- da_forecast_output$restart_list$z_vars$uo
    v <- da_forecast_output$restart_list$z_vars$v
    vo <- da_forecast_output$restart_list$z_vars$vo
    xP <- da_forecast_output$restart_list$z_vars$xP
    h <- da_forecast_output$restart_list$z_vars$h
    ho <- da_forecast_output$restart_list$z_vars$ho

    tke <- da_forecast_output$restart_list$zi_vars$tke
    zi <- da_forecast_output$restart_list$zi_vars$zi
    tkeo <- da_forecast_output$restart_list$zi_vars$tkeo
    eps <- da_forecast_output$restart_list$zi_vars$eps
    num <- da_forecast_output$restart_list$zi_vars$num
    nuh <- da_forecast_output$restart_list$zi_vars$nuh
    nus <- da_forecast_output$restart_list$zi_vars$nus

  } else if(config$model_settings$model == "Simstrat") {
    zi <- da_forecast_output$restart_list$zi
    u <- da_forecast_output$restart_list$u
    v <- da_forecast_output$restart_list$v
    temp <- da_forecast_output$restart_list$temp
    k <- da_forecast_output$restart_list$k
    eps <- da_forecast_output$restart_list$eps
    num <- da_forecast_output$restart_list$num
    nuh <- da_forecast_output$restart_list$nuh
    seicheE <- da_forecast_output$restart_list$seicheE
    b_ice <- da_forecast_output$restart_list$b_ice
    w_ice <- da_forecast_output$restart_list$w_ice
    snow <- da_forecast_output$restart_list$snow
  }

  states_config <- da_forecast_output$states_config
  obs_config <- da_forecast_output$obs_config
  pars_config <- da_forecast_output$pars_config
  obs <- da_forecast_output$obs

  diagnostics <- da_forecast_output$diagnostics

  if(forecast_start_datetime %in% full_time) {
    hist_days <- as.numeric(forecast_start_datetime - full_time[1])
    start_forecast_step <- 1 + hist_days
  } else {
    hist_days <- length(full_time)
    start_forecast_step <- NA
  }


  if(!is.null(pars_config)){
    npars <- nrow(pars_config)
  }else{
    npars <- 0
  }

  x_efi <- aperm(x, c(1,3,4,2))
  diagnostics_efi <- diagnostics

  #Set dimensionsda_forecast_output
  ens <- seq(1,dim(x)[4],1)
  depths <- config$model_settings$modeled_depths
  t <- as.numeric(as.POSIXct(lubridate::with_tz(full_time),origin = '1970-01-01 00:00.00 UTC'))
  #obs_states <- seq(1,dim(obs)[3],1)

  #Set variable that states whether value is forecasted
  forecasted <- rep(1, length(t))
  forecasted[1:(hist_days + 1)] <- 0

  if(!use_short_filename | is.na(da_forecast_output$save_file_name_short) | length(which(forecasted == 1)) == 0){
    ncfname <- file.path(forecast_output_directory, paste0(da_forecast_output$save_file_name,".nc"))
  }else{
    ncfname <- file.path(forecast_output_directory, paste0(da_forecast_output$save_file_name_short,".nc"))
  }

  #Define dims
  ensdim <- ncdf4::ncdim_def("ensemble",units = "-",vals = ens, longname = 'ensemble member')
  depthdim <- ncdf4::ncdim_def("depth",units = "meters",vals = as.double(depths), longname = 'Depth from surface')
  timedim <- ncdf4::ncdim_def("time",units = "seconds since 1970-01-01 00:00.00 UTC", longname = "",vals = t)
  if(config$model_settings$model ==  "GLM") {
    # snow_ice_dim <- ncdf4::ncdim_def("snow_ice_dim",units = "meters", vals = c(1, 2, 3), longname = 'snow ice dims')
    # h_dim <- ncdf4::ncdim_def("h",units = "meters", vals = as.double(1), longname = 'height')
    restart_variables_dim <- ncdf4::ncdim_def("restart_variables_dim",units = '', vals = seq(1, dim(restart_variables)[1], 1), longname = 'number of mixing restart variables')
    internal_model_depths_dim <- ncdf4::ncdim_def("internal_model_depths_dim",units = '', vals = seq(1, dim(model_internal_depths)[2]), longname = 'number of possible depths that are simulated in GLM')
  } else if(config$model_settings$model ==  "GOTM") {
    zdim <- ncdf4::ncdim_def("z",units = "meters", vals = as.double(z[2, , 1]), longname = 'Depth (GOTM internal)')
    zidim <- ncdf4::ncdim_def("zi",units = "meters", vals = as.double(zi[2, , 1]), longname = 'interfaces (GOTM internal)')
  } else if(config$model_settings$model ==  "Simstrat") {
    zidim <- ncdf4::ncdim_def("zi",units = "meters", vals = zi[2, , 1], longname = 'depth interfaces (Simstrat internal)')
  }


  #Define variables
  fillvalue <- 1e32

  def_list <- list()
  index <- 1
  def_list[[index]] <- ncdf4::ncvar_def("temp","degC",list(timedim,depthdim, ensdim),fillvalue,'state:temperature',prec="single")
  index <- index + 1
  def_list[[index]] <- ncdf4::ncvar_def("data_assimilation","dimensionless",list(timedim),missval = -99,longname = '0 = no data assimilation; 1 = data assimilated',prec="integer")
  index <- index + 1
  def_list[[index]] <- ncdf4::ncvar_def("forecast","dimensionless",list(timedim),missval = -99,longname = '0 = historical; 1 = forecasted',prec="integer")
  index <- index + 1
  def_list[[index]] <- ncdf4::ncvar_def("da_qc","dimensionless",list(timedim),missval = -99,longname = '0 = successful DA; 1 = no data assimilated',prec="integer")
  index <- index + 1
  if(config$model_settings$model ==  "GLM") {
    # def_list[[index]] <- ncdf4::ncvar_def("snow_ice_thickness","meter", list(snow_ice_dim, timedim, ensdim),missval = -99,longname = 'Ice Thickness',prec="single")
    # index <- index + 1

    def_list[[index]] <- ncdf4::ncvar_def("snow_thickness","meter", list(timedim, ensdim),missval = -99,longname = 'Snow Thickness',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("white_ice_thickness","meter", list(timedim, ensdim),missval = -99,longname = 'White Ice Thickness',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("blue_ice_thickness","meter", list(timedim, ensdim),missval = -99,longname = 'Blue Ice Thickness',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("lake_depth","meter",list(timedim,ensdim),missval = -99,longname = 'Depth of lake',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("avg_surf_temp","degC",list(timedim, ensdim),missval = -99,longname ='Running Average of Surface Temperature',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("restart_variables","dimensionless",list(restart_variables_dim, timedim, ensdim),fillvalue,longname = "variables required to restart mixing",prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("model_internal_depths","meter",list(timedim, internal_model_depths_dim, ensdim), fillvalue, longname = "depths simulated by glm that are required to restart ",prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("salt","g_kg",list(timedim, depthdim, ensdim),fillvalue,longname = "state:salt",prec="single")
  } else if(config$model_settings$model ==  "GOTM") {
    # z vars
    def_list[[index]] <- ncdf4::ncvar_def("salt","g_kg", list(timedim, depthdim, ensdim), missval = -99,longname = 'state:salt',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("u","m/s", list(timedim, zdim, ensdim), missval = -99, longname = 'x-velocity',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("uo","m/s", list(timedim, zdim, ensdim), missval = -99, longname = 'x-velocity - old time step', prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("v","m/s", list(timedim, zdim, ensdim), missval = -99,longname = 'y-velocity',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("vo","m/s", list(timedim, zdim, ensdim), missval = -99, longname = 'y-velocity - old time step', prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("xP","m2/s3", list(timedim, zdim, ensdim), missval = -99,longname = 'extra turbulence production',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("h","m", list(timedim, zdim, ensdim), missval = -99,longname = 'layer thickness',prec="single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("ho","m", list(timedim, zdim, ensdim), missval = -99,longname = 'layer thickness - old time step', prec = "single")
    index <- index + 1


    # zi vars
    def_list[[index]] <- ncdf4::ncvar_def("tke","m2/s2", list(timedim, zidim, ensdim), missval = -99, longname = "turbulent kinetic energy", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("tkeo","m2/s2", list(timedim, zidim, ensdim), missval = -99, longname = "turbulent kinetic energy - old time step", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("eps","m2/s3", list(timedim, zidim, ensdim), missval = -99, longname = "energy dissipation rate", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("num","m2/s", list(timedim, zidim, ensdim), missval = -99, longname = "turbulent diffusivity of momentum", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("nuh","m2/s", list(timedim, zidim, ensdim), missval = -99, longname = "turbulent diffusivity of heat", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("nus","m2/s", list(timedim, zidim, ensdim), missval = -99, longname = "turbulent diffusivity of salt", prec = "single")

  } else if(config$model_settings$model ==  "Simstrat") {
    def_list[[index]] <- ncdf4::ncvar_def("u", "m/s", list(timedim, zidim, ensdim), missval = -99, longname = "x-velocity", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("v", "m/s", list(timedim, zidim, ensdim), missval = -99, longname = "y-velocity", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("simstrat_temp", "degC", list(timedim, zidim, ensdim), missval = -99, longname = "temperature (Simstrat)", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("salt", "g/kg", list(timedim, depthdim, ensdim), missval = -99, longname = "state:salt", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("k", "J/kg", list(timedim, zidim, ensdim), missval = -99, longname = "turbulent kinetic energy", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("eps", "W/kg", list(timedim, zidim, ensdim), missval = -99, longname = "energy dissipation rate", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("num", "m2/s", list(timedim, zidim, ensdim), missval = -99, longname = "turbulent diffusivity of momentum", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("nuh", "m2/s", list(timedim, zidim, ensdim), missval = -99, longname = "turbulent diffusivity of heat", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("seicheE","J", list(timedim, ensdim), missval = -99, longname = "seiche energy", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("b_ice","m", list(timedim, ensdim), missval = -99, longname = "black ice", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("w_ice","m", list(timedim, ensdim), missval = -99, longname = "white ice", prec = "single")
    index <- index + 1
    def_list[[index]] <- ncdf4::ncvar_def("snow","m", list(timedim, ensdim), missval = -99, longname = "snow", prec = "single")
  }


  if(npars > 0){
    for(par in 1:npars){
      def_list[[index+par]] <-ncdf4::ncvar_def(pars_config$par_names_save[par],pars_config$par_units[par],list(timedim,ensdim),fillvalue,paste0("parameter:",pars_config$par_names_save[par]),prec="single")
    }
  }

  if(config$include_wq){
    for(s in 2:length(states_config$state_names)){
      if(states_config$state_names[s] %in% obs_config$state_names_obs){
        tmp_index <- which(obs_config$state_names_obs == states_config$state_names[s])
        long_name <- paste0("state:",obs_config$target_variable[tmp_index])
      }else{
        long_name <- "state"
      }
      def_list[[index+npars+s-1]]<- ncdf4::ncvar_def(states_config$state_names[s],"mmol m-3",list(timedim,depthdim, ensdim),fillvalue,long_name,prec="single")
    }
  }

  if(length(config$output_settings$diagnostics_names) > 0){
    for(s in 1:length(config$output_settings$diagnostics_names)){
      def_list[[index+npars+length(states_config$state_names)-2 + s]]<- ncdf4::ncvar_def(config$output_settings$diagnostics_names[s],"-",list(timedim,depthdim, ensdim),fillvalue,paste0("diagnostic:",config$output_settings$diagnostics_names[s]),prec="single")
    }
  }

  tmp_index <- index+npars+length(states_config$state_names)+length(config$output_settings$diagnostics_names)-1
  for(s in 1:length(obs_config$state_names_obs)){
    if(!obs_config$state_names_obs[s] %in% states_config$state_names){
      tmp_index <- tmp_index + 1
      longname <- paste0("state:",obs_config$target_variable[s])
      def_list[[tmp_index]] <- ncdf4::ncvar_def(obs_config$state_names_obs[s],obs_config$obs_units[s],list(timedim,depthdim, ensdim),fillvalue,longname,prec="single")
    }
  }

  if(add) {
    ncout <- ncdf4::nc_open(ncfname, write = TRUE)
  } else {
    unlink(ncfname)
    ncout <- ncdf4::nc_create(ncfname, def_list, force_v4 = TRUE)
  }
  on.exit({
    ncdf4::nc_close(ncout)
  })

  # create netCDF file and put arrays
  index <- 1
  ncdf4::ncvar_put(ncout, def_list[[index]] ,x_efi[,1:length(depths), , 1])
  index <- index + 1
  ncdf4::ncvar_put(ncout, def_list[[index]] ,as.array(data_assimilation_flag))
  index <- index + 1
  ncdf4::ncvar_put(ncout, def_list[[index]] ,as.array(forecast_flag))
  index <- index + 1
  ncdf4::ncvar_put(ncout, def_list[[index]] ,as.array(da_qc_flag))
  index <- index + 1
  if(config$model_settings$model ==  "GLM") {
    # ncdf4::ncvar_put(ncout, def_list[[index]], snow_ice_thickness)
    # index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], snow_thickness)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], white_ice_thickness)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], blue_ice_thickness)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], lake_depth)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], avg_surf_temp)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], restart_variables)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], model_internal_depths)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], salt)
  } else if(config$model_settings$model ==  "GOTM") {
    # z vars
    ncdf4::ncvar_put(ncout, def_list[[index]], x_efi[,1:length(depths), , 2]) #Salt
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], u)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], uo)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], v)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], vo)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], xP)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], h)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], ho)
    index <- index + 1

    #zi vars
    ncdf4::ncvar_put(ncout, def_list[[index]], tke)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], tkeo)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], eps)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], num)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], nuh)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], nus)
  } else if(config$model_settings$model ==  "Simstrat") {

    ncdf4::ncvar_put(ncout, def_list[[index]], u)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], v)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], temp)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], x_efi[,1:length(depths), , 2]) #Salt
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], k)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], eps)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], num)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], nuh)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], seicheE)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], b_ice)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], w_ice)
    index <- index + 1
    ncdf4::ncvar_put(ncout, def_list[[index]], snow)
  }


  if(npars > 0){
    for(par in 1:npars){
      ncdf4::ncvar_put(ncout,def_list[[index + par]] ,pars[, par, ])
    }
  }

  #if(config$include_wq){
  #  for(s in 2:length(states_config$state_names)){
  #    ncdf4::ncvar_put(ncout,def_list[[index+npars+s-1]],x_efi[,states_config$wq_start[s]:states_config$wq_end[s], ])
  #  }
  #}

  if(length(config$output_settings$diagnostics_names) > 0){
    for(s in 1:length(config$output_settings$diagnostics_names)){
      ncdf4::ncvar_put(ncout, def_list[[index+npars+length(states_config$state_names) - 2 + s]],diagnostics_efi[s, , ,])
    }
  }

  #tmp_index <- index+npars+length(states_config$state_names)+length(config$output_settings$diagnostics_names)-1
  # for(s in 1:length(obs_config$state_names_obs)){
  #   if(!obs_config$state_names_obs[s] %in% states_config$state_names){
  #     tmp_index <- tmp_index + 1
  #     first_index <- 1
  #     for(ii in 1:length(states_config$state_names)){
  #       if(s %in% states_config$states_to_obs[[ii]]){
  #         temp_index <- which(states_config$states_to_obs[[ii]] == s)
  #         if(first_index == 1){
  #           temp_var <- x_efi[, states_config$wq_start[ii]:states_config$wq_end[ii], ] * states_config$states_to_obs_mapping[[ii]][temp_index]
  #           first_index <- 2
  #         }else{
  #           temp_var <- temp_var + x_efi[, states_config$wq_start[ii]:states_config$wq_end[ii], ] * states_config$states_to_obs_mapping[[ii]][temp_index]
  #         }
  #       }
  #     }
  #     ncdf4::ncvar_put(ncout,def_list[[tmp_index]] , temp_var)
  #
  #   }
  # }

  time_of_forecast <- lubridate::with_tz(da_forecast_output$time_of_forecast, tzone = "UTC")

  #Global file metadata
  ncdf4::ncatt_put(ncout,0,"title", config$metadata$forecast_title, prec =  "text")
  ncdf4::ncatt_put(ncout,0,"forecast_iteration_id" ,da_forecast_output$forecast_iteration_id, prec =  "text")
  ncdf4::ncatt_put(ncout,0,"forecast_project_id", config$metadata$forecast_project_id, prec =  "text")
  ncdf4::ncatt_put(ncout,0,"forecast_model_id", config$metadata$model_description$forecast_model_id, prec =  "text")
  ncdf4::ncatt_put(ncout,0,"forecast_model", config$model_settings$model, prec =  "text")
  ncdf4::ncatt_put(ncout,0,"time_zone_of_simulation","UTC", prec =  "text")

  invisible(ncfname)

}
