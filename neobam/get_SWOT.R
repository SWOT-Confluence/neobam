#' Retrieve SWOT data.
#'
#' @param swot_file string path to SWOT NetCDF file
#'
#' @return list of width, slope2, and time matrices
get_swot = function(swot_file) {
  swot = open.nc(swot_file)
  nx = var.get.nc(swot, "nx")
  nt = var.get.nc(swot, "nt")

  node_grp = grp.inq.nc(swot, "node")$self
  width = t(var.get.nc(node_grp, "width"))
  wse = t(var.get.nc(node_grp, "wse"))
  slope2 = t(var.get.nc(node_grp, "slope2"))
  time = t(var.get.nc(node_grp, "time"))

  reach_grp = grp.inq.nc(swot, "reach")$self
  reach_wse=  t(var.get.nc(reach_grp, "wse"))  
  obs_times = t(var.get.nc(reach_grp, "time_str"))   


  close.nc(swot)

  return(list(nx=nx, nt=nt, wse=wse, reach_wse=reach_sewidth=width,slope2=slope2, time=time, obs_times=obs_times))

}