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
    if(verbose) message("GLM run is complete! ", paste0("[", Sys.time(), "]"))
  }

  # GOTM ----
  if(model == "GOTM") {
    GOTMr::run_gotm(sim_folder = file.path(folder, "GOTM"), verbose = verbose)
    if(verbose) message("GOTM run is complete! ", paste0("[", Sys.time(), "]"))
  }

  # Simstrat ----
  if(model == "Simstrat") {
    SimstratR::run_simstrat(sim_folder = file.path(folder, "Simstrat"), par_file = "simstrat.par", verbose = verbose)
    if(verbose) message("Simstrat run is complete! ", paste0("[", Sys.time(), "]"))
  }
}
