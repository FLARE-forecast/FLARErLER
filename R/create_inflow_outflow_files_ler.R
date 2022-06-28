#' @title Generating inflow and output files in the LakeEnsemblR format
#' @details Processes historical inflow data from inflow_obs and from files in the inflow_file_dir into the LakeEnsemblR format
#' @param inflow_file_dir string; full directory path that contains forecasted inflow and outflow files
#' @param inflow_obs string; full path to cleaned inflow observation in the specified format
#' @param working_directory string; full directory where FLARE executes
#' @param state_names vector; vector of state names that will be included in the inflow files
#' @return list with two vectors. One vector is the matrix of inflow_file_names and the other is the matrix of outflow_file_names
#' @export
#' @examples
create_inflow_outflow_files_ler <- function(inflow_file_dir = NULL,
                                            inflow_obs,
                                            working_directory,
                                            config,
                                            state_names)

{

  VARS <- c("time", "FLOW", "TEMP", "SALT")

  if(config$model_settings$model == "glm_aed"){
    VARS <- c(VARS,
              'OXY_oxy',
              'CAR_dic',
              'CAR_ch4',
              'SIL_rsi',
              'NIT_amm',
              'NIT_nit',
              'PHS_frp',
              'OGM_doc',
              'OGM_docr',
              'OGM_poc',
              'OGM_don',
              'OGM_donr',
              'OGM_pon',
              'OGM_dop',
              'OGM_dopr',
              'OGM_pop',
              'PHY_cyano',
              'PHY_green',
              'PHY_diatom')

  }

  start_datetime <- lubridate::as_datetime(config$run_config$start_datetime)
  if(is.na(config$run_config$forecast_start_datetime)){
    end_datetime <- lubridate::as_datetime(config$run_config$end_datetime)
    forecast_start_datetime <- end_datetime
  }else{
    forecast_start_datetime <- lubridate::as_datetime(config$run_config$forecast_start_datetime)
    end_datetime <- forecast_start_datetime + lubridate::days(config$run_config$forecast_horizon)
  }

  obs_inflow <- readr::read_csv(inflow_obs, col_types = readr::cols())
  hour_step <- lubridate::hour(start_datetime)

  if(config$inflow$use_forecasted_inflow){
    obs_inflow <- obs_inflow %>%
      dplyr::filter(time >= lubridate::as_date(start_datetime) & time <= lubridate::as_date(forecast_start_datetime)) %>%
      dplyr::mutate(inflow_num = 1)
  }else{
    obs_inflow <- obs_inflow %>%
      dplyr::filter(time >= lubridate::as_date(start_datetime) & time <= lubridate::as_date(end_datetime)) %>%
      dplyr::mutate(inflow_num = 1)
  }

  obs_outflow <- obs_inflow %>%
    dplyr::select(time, FLOW) %>%
    dplyr::mutate(outflow_num = 1)

  all_files <- NULL

  if(!is.null(inflow_file_dir)){
    if(dir.exists(inflow_file_dir)){
      all_files <- list.files(inflow_file_dir, full.names = TRUE)
    }
  }

  if(config$run_config$forecast_horizon > 16) {
    all_files <- all_files[!grepl("ens00", all_files)]
  }

  inflow_files <- all_files[stringr::str_detect(all_files,"INFLOW")]
  outflow_files <- all_files[stringr::str_detect(all_files,"OUTFLOW")]



  if(length(inflow_files) > 0){
    d <- readr::read_csv(inflow_files[1], col_types = readr::cols())
    num_inflows <- max(c(d$inflow_num,obs_inflow$inflow_num))
  }else{
    num_inflows <- max(obs_inflow$inflow_num)
  }

  inflow_file_names <- array(NA, dim = c(max(c(1, length(inflow_files))), num_inflows))

  for(j in 1:num_inflows) {

    if(length(inflow_files) == 0 | end_datetime == forecast_start_datetime){

      obs_inflow_tmp <- obs_inflow %>%
        dplyr::filter(inflow_num == j) %>%
        dplyr::select(dplyr::all_of(VARS))

      obs_inflow_tmp <- as.data.frame(obs_inflow_tmp)
      obs_inflow_tmp[, 1] <- format(obs_inflow_tmp[, 1], format="%Y-%m-%d %H:%M:%S")
      obs_inflow_tmp[, 1] <- lubridate::with_tz(obs_inflow_tmp[, 1]) + lubridate::hours(hour_step)
      obs_inflow_tmp[, 1] <- format(obs_inflow_tmp[, 1], format="%Y-%m-%d %H:%M:%S")
      colnames(obs_inflow_tmp) <- c("datetime", "Flow_metersCubedPerSecond", "Water_Temperature_celsius", "Salinity_practicalSalinityUnits")

      inflow_file_name <- file.path(working_directory, paste0("inflow",j,".csv"))

      readr::write_csv(x = obs_inflow_tmp,
                       file = inflow_file_name,
                       quote = "none")
      inflow_file_names[, j] <- inflow_file_name
    } else {

      for(i in 1:length(inflow_files)) {
        d <- readr::read_csv(inflow_files[i], col_types = readr::cols()) %>%
          dplyr::filter(inflow_num == j) %>%
          dplyr::select(dplyr::all_of(VARS)) %>%
          dplyr::mutate_at(dplyr::vars(dplyr::all_of(VARS)), list(~round(., 4)))

        obs_inflow_tmp <- obs_inflow %>%
          dplyr::filter(inflow_num == j,
                        time < lubridate::as_date(forecast_start_datetime)) %>%
          dplyr::select(dplyr::all_of(VARS))


        inflow <- rbind(obs_inflow_tmp, d)
        inflow <- as.data.frame(inflow)
        # inflow[, 1] <- as.POSIXct(inflow[, 1], tz = tz) + lubridate::hours(hour_step)
        inflow[, 1] <- format(inflow[, 1], format="%Y-%m-%d %H:%M:%S")
        inflow[, 1] <- lubridate::with_tz(inflow[, 1]) + lubridate::hours(hour_step)
        inflow[, 1] <- format(inflow[, 1], format="%Y-%m-%d %H:%M:%S")
        colnames(inflow) <- c("datetime", "Flow_metersCubedPerSecond", "Water_Temperature_celsius", "Salinity_practicalSalinityUnits")

        inflow_file_name <- file.path(working_directory, paste0("inflow",j,"_ens",i,".csv"))

        inflow_file_names[i, j] <- inflow_file_name

        if(config$inflow$use_forecasted_inflow){
          readr::write_csv(x = inflow,
                           inflow_file_name,
                           quote = "none")
        }else{
          readr::write_csv(x = obs_inflow_tmp,
                           inflow_file_name,
                           quote = "none")

        }
      }
    }
  }


  if(length(outflow_files) > 0){
    d <- readr::read_csv(outflow_files[1], col_types = readr::cols())
    num_outflows <- max(c(d$outflow_num,obs_outflow$outflow_num))
  }else{
    num_outflows <- max(obs_outflow$outflow_num)
  }

  outflow_file_names <- array(NA, dim = c(max(c(1, length(outflow_files))),num_outflows))

  for(j in 1:num_outflows){

    if(length(inflow_files) == 0 | end_datetime == forecast_start_datetime){

      obs_outflow_tmp <- obs_outflow %>%
        dplyr::filter(outflow_num == j) %>%
        dplyr::select(c("time", "FLOW"))

      outflow_file_name <- file.path(working_directory, paste0("outflow",j,".csv"))

      readr::write_csv(x = obs_outflow_tmp,
                       outflow_file_name,
                       quote = "none")
      outflow_file_names[, j] <- outflow_file_name
    } else {

      for(i in 1:length(outflow_files)){
        d <- readr::read_csv(outflow_files[i], col_types = readr::cols())%>%
          dplyr::filter(outflow_num == j) %>%
          dplyr::select(time, FLOW) %>%
          dplyr::mutate_at(dplyr::vars(c("FLOW")), list(~round(., 4)))

        obs_outflow_tmp <- obs_outflow %>%
          dplyr::filter(outflow_num == j,
                        time < lubridate::as_date(forecast_start_datetime)) %>%
          dplyr::select(-outflow_num)

        outflow <- rbind(obs_outflow_tmp, d)
        outflow <- as.data.frame(outflow)
        outflow[, 1] <- format(outflow[, 1], format="%Y-%m-%d %H:%M:%S")
        colnames(outflow) <- c("datetime", "Flow_metersCubedPerSecond")

        outflow_file_name <- file.path(working_directory, paste0("outflow",j,"_ens",i,".csv"))

        outflow_file_names[i, j]  <- outflow_file_name

        if(config$inflow$use_forecasted_inflow){
          readr::write_csv(x = outflow,
                           outflow_file_name,
                           quote = "none")
        }else{
          readr::write_csv(x = obs_outflow_tmp,
                           outflow_file_name,
                           quote = "none")
        }
      }
    }
  }

  return(list(inflow_file_names = as.character(inflow_file_names),
              outflow_file_names = as.character(outflow_file_names)))
}
