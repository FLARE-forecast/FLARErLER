#' @importFrom ncdf4 nc_open ncvar_get nc_close
#' @importFrom LakeEnsemblR get_output
#' @importFrom rLakeAnalyzer get.offsets

get_ler_var_all <- function(model,
                            working_dir,
                            z_out, vars_depth,
                            vars_no_depth,
                            diagnostic_vars,
                            ler_yaml,
                            run_success = run_success) {

  temp <- LakeEnsemblR::get_output(yaml = ler_yaml, model = model, vars = "temp", obs_depths = z_out, run_success = run_success)$temp
  salt <- LakeEnsemblR::get_output(yaml = ler_yaml, model = model, vars = "salt", obs_depths = z_out, run_success = run_success)$salt
  ice <- LakeEnsemblR::get_output(yaml = ler_yaml, model = model, vars = "ice_height", run_success = run_success)$ice_height
  deps <- rLakeAnalyzer::get.offsets(temp)
  # Subset to z_out
  idx <- which(deps %in% z_out) + 1
  temp <- temp[, c(1, idx)]
  salt <- salt[, c(1, idx)]
  deps <- deps[which(deps %in% z_out)]

  final_time_step <- nrow(temp)

  # No varying water level in Simstrat
  heights_surf <- max(deps)
  heights <- deps
  # heights_out <- rep()

  temps <- apply(temp[, -1], 2, mean)
  salt <- unlist(salt[final_time_step, -1])

  output <- array(NA, dim=c(length(temps), length(vars_depth)))
  for(v in 1:length(vars_depth)){
    output[,v] <- temps
  }

  depths_enkf = rev(heights_surf - heights)

  restart_vars <- LakeEnsemblR::read_restart(working_dir, model)

  if(model == "GLM") {

    glm_nc <- ncdf4::nc_open(file.path(working_dir, model, "output", "output.nc"))
    on.exit({
      ncdf4::nc_close(glm_nc)
    })
    tallest_layer <- ncdf4::ncvar_get(glm_nc, "NS")
    final_time_step <- length(tallest_layer)
    tallest_layer <- tallest_layer[final_time_step] # Edited
    heights <- matrix(ncdf4::ncvar_get(glm_nc, "z"), ncol = final_time_step)
    heights_surf <- heights[tallest_layer, final_time_step]
    heights <- heights[1:tallest_layer, final_time_step]

    output <- array(NA,dim=c(tallest_layer, length(vars_depth)))
    for(v in 1:length(vars_depth)){
      var_modeled <-  matrix(ncdf4::ncvar_get(glm_nc, vars_depth[v]), ncol = final_time_step)
      output[, v] <- var_modeled[1:tallest_layer,final_time_step]
    }

    depths_enkf = rev(heights_surf - heights)

    output_no_depth <- NA

    if(length(diagnostic_vars) > 0){
      diagnostics_output <- array(NA,dim=c(length(z_out), length(diagnostic_vars)))
      for(v in 1:length(diagnostic_vars)){
        var_modeled <- ncdf4::ncvar_get(glm_nc, diagnostic_vars[v])[1:tallest_layer, final_time_step]
        var_modeled <- approx(heights, var_modeled, xout = z_out, rule = 2)$y
        diagnostics_output[,v] <- var_modeled
      }
    }else{
      diagnostics_output <- NA
    }

  }

  # GOTM ----
  if( model == "GOTM") {

    output_no_depth <- NA

    # if(length(diagnostic_vars) > 0){
    #   diagnostics_output <- array(NA,dim=c(tallest_layer, length(diagnostic_vars)))
    #   for(v in 1:length(diagnostic_vars)){
    #     var_modeled <- ncdf4::ncvar_get(nc, diagnostic_vars[v])[, final_time_step]
    #     diagnostics_output[,v] <- var_modeled[1:tallest_layer]
    #   }
    # }else{
    diagnostics_output <- NULL
    # }
  }

  # Simstrat ----
  if( model == "Simstrat") {

    # snow <- read.delim(file.path(model, "output", "SnowH_out.dat"), sep = ",")[final_time_step, 2]
    # ice_white <- read.delim(file.path(model, "output", "WhiteIceH_out.dat"), sep = ",")[final_time_step, 2]

    # Extract variables for restarting initial conditions
    # U <- read.delim(file.path(model, "output", "U_out.dat"), sep = ",")[final_time_step, -1]
    # deps2 <- colnames(U) %>%
    #   regmatches(., gregexpr("[[:digit:]]+\\.*[[:digit:]]*", .)) %>%
    #   unlist() %>%
    #   as.numeric()
    # U <- approx(deps2, U, z_out, rule = 2)$y
    # V <- read.delim(file.path(model, "output", "V_out.dat"), sep = ",")[final_time_step, -1] %>%
    #   approx(deps2, ., z_out, rule = 2) %>%
    #   .[[2]]
    # k <- read.delim(file.path(model, "output", "k_out.dat"), sep = ",")[final_time_step, -1] %>%
    #   approx(deps2, ., z_out, rule = 2) %>%
    #   .[[2]]
    # eps <- read.delim(file.path(model, "output", "eps_out.dat"), sep = ",")[final_time_step, -1] %>%
    #   approx(deps2, ., z_out, rule = 2) %>%
    #   .[[2]]

    output_no_depth <- NA

    if(length(diagnostic_vars) > 0){
      diagnostics_output <- array(NA,dim=c(tallest_layer, length(diagnostic_vars)))
      for(v in 1:length(diagnostic_vars)){
        var_modeled <- ncdf4::ncvar_get(nc, diagnostic_vars[v])[, final_time_step]
        diagnostics_output[,v] <- var_modeled[1:tallest_layer]
      }
    }else{
      diagnostics_output <- NULL
    }
  }


  return(list(output = output,
              salt = salt,
              output_no_depth = output_no_depth,
              lake_depth = heights_surf,
              depths_enkf = depths_enkf,
              restart_vars = restart_vars,
              diagnostics_output = diagnostics_output))
}
