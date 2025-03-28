source('/nas/cee-water/cjgleason/colin/neobam/neobam/get_invalid_nodes_times.R')
source('/nas/cee-water/cjgleason/colin/neobam/neobam/remove_invalid.R')
source('/nas/cee-water/cjgleason/colin/neobam/neobam/get_SWOT.R')
source('/nas/cee-water/cjgleason/colin/neobam/neobam/get_sos.R')

library(RNetCDF,quietly=TRUE,warn.conflicts = FALSE)
#' Title
#'
#' @param swot_file string path to SWOT NetCDF file
#' @param sos_file string path to SOS NetCDF file
#' @param reach_id integer unique reach identifier
#'
#' @return named list of data needed to run neoBAM
#' @export
get_input = function(swot_file, sos_file, reach_id) {


  # Get SWOT
  swot_data = get_swot(swot_file)

  # Get SOS
  sos_data = get_sos(sos_file, reach_id)

  # Check for valid number of observations
  cleaned_data = check_observations(swot_data, sos_data)

  # Return list of valid observations
  if (length(data) == 0) {
    # print('in get_input(), neobam has decided the data are invalid')
         # print(swot_data)
    return(list(valid=FALSE, time=swot_data$time,reach_id=reach_id, nx=swot_data$nx, nt=swot_data$nt, node_ids=sos_data$Q_priors$nids))
  } else {

      #drop bad nodes on the cumulative distance AFTER calculating it node by node
      if(!identical(cleaned_data$invalid_nodes, integer(0))){
      cleaned_data$sos_data$node_cum_dist=cleaned_data$sos_data$node_cum_dist[-cleaned_data$invalid_nodes]}
      
    

    # Create a list of data with reach identifier
    return(list(valid=TRUE, reach_id = reach_id,
                time=cleaned_data$swot_data$time,
                obs_times=cleaned_data$swot_data$obs_times,
                swot_data=cleaned_data$swot_data,
                sos_data=cleaned_data$sos_data,
                invalid_nodes=cleaned_data$invalid_nodes,
                invalid_times=cleaned_data$invalid_times,
                node_ids=sos_data$Q_priors$nids))
  }
}


check_observations = function(swot_data, sos_data) {

  # Q priors
  qhat = sos_data$Q_priors$logQ_hat
  qmax = sos_data$Q_priors$upperbound_logQ
  qmin = sos_data$Q_priors$lowerbound_logQ
  qsd = sos_data$Q_priors$logQ_sd

  if (all(is.na(qhat[[1]])) || is.na(qmax[[1]]) || is.na(qmin[[1]]) || is.na(qsd[[1]])) { return(vector(mode = "list")) }
    
  
  # SWOT data
  swot_data$width[swot_data$width < 0] = NA
  swot_data$wse[swot_data$wse < 0] = NA
  swot_data$slope2[swot_data$slope2 < 0] = NA
  # print(swot_data$width)
  # print(swot_data$slope2)

    
 
  invalid = get_invalid_nodes_times(swot_data$wse, swot_data$slope2, swot_data$time)
 
  # Return valid data (or empty list if invalid)
       #here, we pass the inpute data do the file, so this is where we rewrite missing qhats
  return(remove_invalid(swot_data, sos_data, invalid$invalid_nodes, invalid$invalid_times))
}




