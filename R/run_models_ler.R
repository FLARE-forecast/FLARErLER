#' Run the lake model from LakeEnsemblR
#'
#'
#' @param model vector; model to export configuration file.
#'   Options include c("GOTM", "GLM", "Simstrat", "FLake", "MyLake")
#' @param folder folder
#' @param verbose Boolean; Should model output be shown in the console. Defaults to FALSE
#' @keywords methods
#'
#' @import GLM3r
#' @import GOTMr
#' @import SimstratR
#'
#' @export
run_models_ler <- function(model, folder, verbose, restart, member, the_temps, model_depths) {

  if(model == "GLM") {

    GLM3r::run_glm(sim_folder = file.path(folder, "GLM"), verbose = verbose)
    # glmtools::plot_temp(file.path(folder, "GLM", "output.nc"))

    if(verbose) message("GLM run is complete! ", paste0("[", Sys.time(), "]"))
  }

  # GOTM ----
  if(model == "GOTM") {

    # if(restart) {
    #   file.copy(from = file.path(folder, "GOTM", paste0("restart_", member, ".nc")),
    #             to = file.path(folder, "GOTM", "restart.nc"),
    #             overwrite = TRUE)
    #   nc <- ncdf4::nc_open(file.path(folder, "GOTM", "restart.nc"), write = TRUE)
    #   gotm_depths <- ncdf4::ncvar_get(nc, "z")
    #   input_temps <- approx(model_depths, the_temps, xout = abs(gotm_depths), rule = 2)$y
    #   ncdf4::ncvar_put(nc, "temp", input_temps)
    #   ncdf4::nc_close(nc)
    # }

    GOTMr::run_gotm(sim_folder = file.path(folder, "GOTM"), verbose = verbose)

    file.copy(from = file.path(folder, "GOTM", "restart.nc"),
              to = file.path(folder, "GOTM", paste0("restart_", member, ".nc")),
              overwrite = TRUE)

    if(verbose) message("GOTM run is complete! ", paste0("[", Sys.time(), "]"))
  }

  # Simstrat ----
  if(model == "Simstrat") {

    SimstratR::run_simstrat(sim_folder = file.path(folder, "Simstrat"), par_file = "simstrat.par", verbose = verbose)

    if(verbose) message("Simstrat run is complete! ", paste0("[", Sys.time(), "]"))
  }

}
