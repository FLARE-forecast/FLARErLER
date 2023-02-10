

set_up_model_ler <- function(model,
                             config,
                             working_directory,
                             state_names,
                             met_file_names,
                             inflow_file_names = NULL,
                             outflow_file_names = NULL,
                             member = "1",
                             start_datetime,
                             end_datetime){

  oldwd <- getwd()
  ens_working_directory <- file.path(working_directory, member)
  setwd(ens_working_directory)
  on.exit({
    setwd(oldwd)
  })

  file.copy(from = file.path(config$file_path$configuration_directory, config$model_settings$base_ler_yaml),
            to = ens_working_directory, overwrite = TRUE)
  yaml_file <- config$model_settings$base_ler_yaml

  ler_yaml <- yaml::read_yaml(yaml_file)
  if(is.null(inflow_file_names)) {
    ler_yaml$inflows$use <- FALSE
  }


  idx <- sample(1:length(met_file_names), 1)
  ler_yaml$input$meteo$file <- paste0("../", basename(met_file_names[idx]))
  if(!is.null(inflow_file_names)) {
    ler_yaml$inflows$file <- paste0("../", basename(inflow_file_names[idx]))
  }
  ler_yaml$time$start <- format(start_datetime, format = "%Y-%m-%d %H:%M:%S")
  ler_yaml$time$stop <- format(end_datetime, format = "%Y-%m-%d %H:%M:%S")

  # yml <- yaml::read_yaml(file.path(ens_working_directory, ler_yaml))
  ler_directory <- ens_working_directory

  if(model == "GLM") {

    dir.create(file.path(ens_working_directory, "GLM"), showWarnings = FALSE)
    file.copy(from = file.path(config$file_path$configuration_directory, config$model_settings$base_GLM_nml),
              to = file.path(ens_working_directory, "GLM", "glm3.nml"), overwrite = TRUE)

    inflow_var_names <- c("FLOW","TEMP","SALT")

    if(ler_yaml$inflows$use) {
      ler_yaml$model_parameters$GLM$`inflow/num_inflows` <- ncol(inflow_file_names)
      ler_yaml$model_parameters$GLM$`outflow/num_outlet` <- ncol(outflow_file_names)
    } else {
      ler_yaml$model_parameters$GLM$`inflow/num_inflows` <- as.integer(0)
      ler_yaml$model_parameters$GLM$`outflow/num_outlet` <- as.integer(0)
    }

    ler_yaml$model_parameters$GLM$`init_profiles/num_wq_vars` <- 0
    ler_yaml$model_parameters$GLM$`inflow/inflow_varnum` <- length(inflow_var_names)
    ler_yaml$model_parameters$GLM$`output/out_dir` <- "'output'"

    if(config$include_wq){

      file.copy(from =  file.path(config$run_config$forecast_location,config$base_AED_nml),
                to = paste0(ens_working_directory, "/", "aed2.nml"), overwrite = TRUE)

      file.copy(from =  file.path(config$run_config$forecast_location,config$base_AED_phyto_pars_nml),
                to = paste0(ens_working_directory, "/", "aed2_phyto_pars.nml"), overwrite = TRUE)

      file.copy(from =  file.path(config$run_config$forecast_location,config$base_AED_zoop_pars_nml),
                to = paste0(ens_working_directory, "/", "aed2_zoop_pars.nml"), overwrite = TRUE)

    }
  }

  gotmtools::write_yaml(ler_yaml, yaml_file)

  suppressMessages({
    LakeEnsemblR::export_config(config_file = basename(yaml_file), model = model, dirs = TRUE,
                              time = FALSE, location = TRUE, output_settings = TRUE,
                              meteo = TRUE, init_cond = FALSE, extinction = TRUE,
                              inflow = TRUE, model_parameters = TRUE,
                              folder = ens_working_directory)
    })

    if(file.exists(file.path(ens_working_directory,model,"restart.nc"))){
      unlink(file.path(ens_working_directory,model,"restart.nc"))
    }


}
