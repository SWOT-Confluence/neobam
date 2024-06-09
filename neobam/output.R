#' Write neoBAM output to NetCDF file.
#'
#' discharge time series
#' posteriors (mean and sd): r, logn, logWb, logDb

# Libraries
library(RNetCDF,quietly=TRUE,warn.conflicts = FALSE)


# Constants
FILL = -999999999999

#' Insert NA values back into discharge list vectors to account for invalid
#' time steps.
#'
#' @param discharge list of discharge time series to insert NA at time steps
#' @param invalid_times list of invalid time indexes
#'
#' @return list of data with NA in place of invalid nodes
concatenate_invalid = function(discharge, invalid_times) {

  # Time-level data
  for (index in invalid_times) {
    discharge = append(discharge, NA, after=index-1)
  }
  return(discharge)
}

#' Write posteriors data to NetCDF file.
#'
#' @param chain integer number of neoBAM run
#' @param nc_out NetCDF file pointer to write to
#' @param posteriors list of posteriors
write_posteriors = function(nc_out, posteriors) {

  # Chain

  # Posteriors
  r = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "r"),
    grp.inq.nc(nc_out, "r")$self
  )
  print("here are posteriors---------")
  print(posteriors$r$mean)
  print('there were posterirs......')
  var.def.nc(r, "mean", "NC_DOUBLE", "nx")
  att.put.nc(r, "mean", "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(r, "mean", posteriors$r$mean)
  var.def.nc(r, "sd", "NC_DOUBLE", NA)
  att.put.nc(r, "sd", "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(r, "sd", posteriors$r$sd)

  logn = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logn"),
    grp.inq.nc(nc_out, "logn")$self
  )
  var.def.nc(logn, "mean", "NC_DOUBLE", "nx")
  att.put.nc(logn, "mean", "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logn, "mean", posteriors$logn$mean)
  var.def.nc(logn, "sd", "NC_DOUBLE", NA)
  att.put.nc(logn, "sd", "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logn, "sd", posteriors$logn$sd)

  logWb = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logWb"),
    grp.inq.nc(nc_out, "logWb")$self
  )
  var.def.nc(logWb, "mean", "NC_DOUBLE", "nx")
  att.put.nc(logWb, "mean", "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logWb, "mean", posteriors$logWb$mean)
  var.def.nc(logWb, "sd", "NC_DOUBLE", NA)
  att.put.nc(logWb, "sd", "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logWb, "sd", posteriors$logWb$sd)

  logDb = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "logDb"),
    grp.inq.nc(nc_out, "logDb")$self
  )
  var.def.nc(logDb, "mean", "NC_DOUBLE", "nx")
  att.put.nc(logDb, "mean", "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logDb, "mean", posteriors$logDb$mean)
  var.def.nc(logDb, "sd", "NC_DOUBLE", NA)
  att.put.nc(logDb, "sd", "_FillValue", "NC_DOUBLE", FILL)
  var.put.nc(logDb, "sd", posteriors$logDb$sd)

}

#' Write discharge data to NetCDF file.
#'
#' @param chain integer number that indicates neoBAM run
#' @param nc_out NetCDF file pointer to write to
#' @param discharge list of discharge values
write_discharge = function(chain, nc_out, discharge,discharge_sd, is_valid) {

  # # time
  # discharge_time_var = tryCatch(
  #   error = function(cond) grp.def.nc(discharge_time, "discharge_time"),
  #   grp.inq.nc(nc_out, "discharge_time")$self
  # )
  # var.def.nc(discharge_time_var, "discharge_time", "NC_DOUBLE", "nt")
  # att.put.nc(discharge_time_var, "discharge_time", "_FillValue", "NC_DOUBLE", FILL)
  # # discharge[is.nan(discharge)] = NA
  # # if (is_valid){
  # var.put.nc(q, "discharge_time", as.numeric(unlist(discharge)))

  # # } else{
  # #   var.put.nc(q,"discharge_time",FILL)
  # # }

  # Discharge
  q = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "q"),
    grp.inq.nc(nc_out, "q")$self
  )
  var.def.nc(q, "q", "NC_DOUBLE", "nt")
  att.put.nc(q, "q", "_FillValue", "NC_DOUBLE", FILL)
  # discharge[is.nan(discharge)] = NA
  if (is_valid){
      var.put.nc(q, "q", as.numeric(unlist(discharge)))

  } else{
    var.put.nc(q,"q",FILL)
  }


    # Discharge error
  q = tryCatch(
    error = function(cond) grp.def.nc(nc_out, "q"),
    grp.inq.nc(nc_out, "q")$self
  )
  var.def.nc(q, "q_sd", "NC_DOUBLE", NA)
  att.put.nc(q, "q_sd", "_FillValue", "NC_DOUBLE", FILL)
  # discharge[is.nan(discharge)] = NA
  if (is_valid){
      var.put.nc(q, "q_sd", discharge_sd)

  } else{
    var.put.nc(q,"q_sd",FILL)
  }

}

#' Write discharge and posteriors to NetCDF file.
#'
#' @param data named list of metadata
#' @param posteriors list of posterior "chains"
#' @param discharge list of discharge "chains"
#' @param out_dir string to output directory
#'
#' @export
write_output = function(in_data, data, posteriors, discharge, out_dir, is_valid) {

  # Concatenate invalid times back into discharge
  discharge = lapply(discharge, concatenate_invalid, invalid_times=data$invalid_times)

  for (i in 1:length(in_data$reach_id[[1]])){

    rid = in_data$reach_id[[1]][i]

    # File creation
    nc_file = paste(out_dir, paste0(rid, "_geobam.nc"), sep=.Platform$file.sep)
    nc_out = create.nc(nc_file, format="netcdf4")

    # Global attributes
    att.put.nc(nc_out, "NC_GLOBAL", "reach_id", "NC_INT64", rid)
    att.put.nc(nc_out, "NC_GLOBAL", "set_ids", "NC_INT64", unlist(in_data$reach_id[[1]]))
    att.put.nc(nc_out, "NC_GLOBAL", "discharge_time", "NC_STRING", unlist(in_data$time))

    # Dimensions
    # dim.def.nc(nc_out, "set_length", length(in_data$reach_id[[1]]))
    # var.def.nc(nc_out, "set_ids", "NC_INT", "set_length")
    # att.put.nc(nc_out, "set_ids", "units", "NC_STRING", "set_length")
    # var.put.nc(nc_out, "set_ids", unlist(in_data$reach_id[[1]]))


    dim.def.nc(nc_out, "nt", length(in_data$time))
    var.def.nc(nc_out, "nt", "NC_INT", "nt")
    att.put.nc(nc_out, "nt", "units", "NC_STRING", "time")
    var.put.nc(nc_out, "nt", seq(from = 0, by = 1, length.out = length(in_data$time)))
    # var.put.nc(nc_out, "nt", in_data$time)
  
    dim.def.nc(nc_out, "nx", length(posteriors$r$mean))
    var.def.nc(nc_out, "nx", "NC_INT", "nx")
    att.put.nc(nc_out, "nx", "units", "NC_STRING", "num_nodes")
    var.put.nc(nc_out, "nx", seq(from = 0, by = 1, length.out = length(posteriors$r$mean)))

    # # Write data
    # lapply(list(1,2,3), write_posteriors, nc_out=nc_out, posteriors=posteriors)
    write_posteriors(nc_out=nc_out, posteriors=posteriors)

    # # Discharge
    # lapply(list(1,2,3), write_discharge, nc_out=nc_out, discharge=discharge)
    # print('sd')
    # print(data$posterior_Q_sd)
    write_discharge(nc_out=nc_out, discharge=discharge, discharge_sd = data$posterior_Q_sd,is_valid=is_valid)

    # Close file
    close.nc(nc_out)
  }
}
