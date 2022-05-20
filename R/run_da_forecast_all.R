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

run_da_forecast_all <- function(states_init,
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
                                debug = FALSE) {
  if(config$model_settings$use_ler) {
    FLARErLER::run_da_forecast_ler(states_init = states_init,
                                pars_init = pars_init,
                                aux_states_init = aux_states_init,
                                obs = obs,
                                obs_sd = obs_sd,
                                model_sd = model_sd,
                                working_directory = working_directory,
                                met_file_names = met_file_names,
                                inflow_file_names = inflow_file_names,
                                outflow_file_names = outflow_file_names,
                                config = config,
                                pars_config = pars_config,
                                states_config = states_config,
                                obs_config = obs_config,
                                management = management,
                                da_method = da_method,
                                par_fit_method = par_fit_method,
                                debug = debug)
  } else {
    FLARErLER::run_da_forecast(states_init = states_init,
                            pars_init = pars_init,
                            aux_states_init = aux_states_init,
                            obs = obs,
                            obs_sd = obs_sd,
                            model_sd = model_sd,
                            working_directory = working_directory,
                            met_file_names = met_file_names,
                            inflow_file_names = inflow_file_names,
                            outflow_file_names = outflow_file_names,
                            config = config,
                            pars_config = pars_config,
                            states_config = states_config,
                            obs_config = obs_config,
                            management = management,
                            da_method = da_method,
                            par_fit_method = par_fit_method,
                            debug = debug)
  }
}
