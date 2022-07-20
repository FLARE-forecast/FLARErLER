#' @param folder filepath; to folder which contains the model
#' @param model character; name of the model
#' @noRd

check_model_output <- function(folder, model) {

  if(model == "GLM") {
    outfile <- file.path(folder, model, "output", "output.nc")
    file_check <- file.exists(outfile)

    if(file_check) {
      nc <- ncdf4::nc_open(outfile)
      on.exit({
        ncdf4::nc_close(nc)
      })
      tmp <- ncdf4::ncvar_get(nc, "temp")
      nan_check <- all(is.nan(unlist(tmp)))
      if(nan_check) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }
  } else if(model == "GOTM") {
    outfile <- file.path(folder, model, "output", "output.nc")
    file_check <- file.exists(outfile)

    if(file_check) {
      nc <- ncdf4::nc_open(outfile)
      on.exit({
        ncdf4::nc_close(nc)
      })
      tmp <- ncdf4::ncvar_get(nc, "temp")
      nan_check <- all(is.nan(unlist(tmp)))
      if(nan_check) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }
  }

  if(model == "Simstrat") {
    outfile <- file.path(folder, model, "output", "T_out.dat")
    file_check <- file.exists(outfile)

    if(file_check) {

      tmp <- tryCatch({
        read.csv(outfile)
      },
      error = function(e) return(NaN))
      nan_check <- any(is.nan(unlist(tmp)) | any(unlist(tmp) == "NaN"))
      if(nan_check) {
        return(FALSE)
      } else {
        if(nrow(tmp) == 0) {
          return(FALSE)
        } else {
          return(TRUE)
        }
      }
    } else {
      return(FALSE)
    }
  }
}
