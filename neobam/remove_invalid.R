remove_invalid = function(swot_data, sos_data, invalid_nodes, invalid_times){
  # print(invalid_nodes)
    #make a dataframe with a month keyf ield
    Q_hat_df= data.frame(logQ_hat=sos_data$Q_priors$logQ_hat, month=as.numeric(1:12))
    #change the date to a datetime object  using the !#$!@#$ origin of the data
    #time is nx by nt, we want an nt vector only
    obs_date= as.Date(as.Date(swot_data$time[1,]/86400,origin = '2000-01-01'),format='%Y%m%d')
  
    #make it a month
    obs_month=data.frame(month=as.numeric(format(obs_date,'%m')))
      # print(obs_month)
    #now join the Qhat_df by month
    Q_hat_month=left_join(obs_month,Q_hat_df,by='month')
    #some of these will be blank, so we need to drop the NAs in order to make the data work out
    # print(Q_hat_month)
     Q_hat_month$logQ_hat[is.na(Q_hat_month$logQ_hat)]=mean(Q_hat_month$logQ_hat,na.rm=TRUE)
#     print(Q_hat_month)
#                                 bonk
    
  # All valid
  if (identical(invalid_nodes, integer(0)) && identical(invalid_times, integer(0))) {
    # print('all valid')
      sos_data$Q_priors$logQ_hat=Q_hat_month$logQ_hat
    return(list(swot_data=swot_data, sos_data=sos_data,
                invalid_nodes=invalid_nodes, invalid_times=invalid_times))
    # Valid nodes
  } else if (identical(invalid_nodes, integer(0))) {
    # print('valid nodes')
    swot_data$width = swot_data$width[, -invalid_times]
    swot_data$wse= swot_data$wse[, -invalid_times]
    swot_data$slope2 = swot_data$slope2[, -invalid_times]
    swot_data$time = swot_data$time[, -invalid_times]
    swot_data$obs_times = swot_data$obs_times[, -invalid_times]
    sos_data$Q_priors$logQ_hat=Q_hat_month$logQ_hat[-invalid_times]

    # Valid time steps
  } else if (identical(invalid_times, integer(0))) {
    # print('valid time')
    swot_data$width = swot_data$width[-invalid_nodes,]
      swot_data$wse = swot_data$wse[-invalid_nodes,]
    swot_data$slope2 = swot_data$slope2[-invalid_nodes,]
     swot_data$time = swot_data$time[-invalid_nodes,]
       swot_data$obs_times = swot_data$obs_times[-invalid_nodes,]
    sos_data$window_params$logWb_hat = sos_data$window_params$logWb_hat[-invalid_nodes]
    sos_data$window_params$logWb_sd = sos_data$window_params$logWb_sd[-invalid_nodes]
    sos_data$window_params$logDb_hat = sos_data$window_params$logDb_hat[-invalid_nodes]
    sos_data$window_params$logDb_sd = sos_data$window_params$logDb_sd[-invalid_nodes]
    sos_data$window_params$r_hat = sos_data$window_params$r_hat[-invalid_nodes]
    sos_data$window_params$r_sd = sos_data$window_params$r_sd[-invalid_nodes]
    sos_data$window_params$logn_hat = sos_data$window_params$logn_hat[-invalid_nodes]
    sos_data$window_params$logn_sd = sos_data$window_params$logn_sd[-invalid_nodes]

    # Both invalid
  } else {
    # print('both invalid')
    swot_data$width = swot_data$width[-invalid_nodes, -invalid_times]
    swot_data$wse = swot_data$wse[-invalid_nodes, -invalid_times]
    swot_data$slope2 = swot_data$slope2[-invalid_nodes, -invalid_times]
    swot_data$time = swot_data$time[-invalid_nodes, -invalid_times]
    swot_data$obs_times = swot_data$obs_times[-invalid_nodes, -invalid_times]
    sos_data$Q_priors$logQ_hat=Q_hat_month$logQ_hat[-invalid_times]  
      
      
    sos_data$window_params$logWb_hat = sos_data$window_params$logWb_hat[-invalid_nodes]
    sos_data$window_params$logWb_sd = sos_data$window_params$logWb_sd[-invalid_nodes]
    sos_data$window_params$logDb_hat = sos_data$window_params$logDb_hat[-invalid_nodes]
    sos_data$window_params$logDb_sd = sos_data$window_params$logDb_sd[-invalid_nodes]
    sos_data$window_params$r_hat = sos_data$window_params$r_hat[-invalid_nodes]
    sos_data$window_params$r_sd = sos_data$window_params$r_sd[-invalid_nodes]
    sos_data$window_params$logn_hat = sos_data$window_params$logn_hat[-invalid_nodes]
    sos_data$window_params$logn_sd = sos_data$window_params$logn_sd[-invalid_nodes]

  }
    


  # Return list to indicate invalid data
  # if (is.null(dim(swot_data$width)) || nrow(swot_data$width) < 3 || ncol(swot_data$width) < 3 ) { return(vector(mode = "list")) }
  if (is.null(dim(swot_data$wse)) || nrow(swot_data$wse) < 3 || ncol(swot_data$wse) < 3 ) { return(vector(mode = "list")) }
  if (is.null(dim(swot_data$slope2)) || nrow(swot_data$slope2) < 3 || ncol(swot_data$slope2) < 3 ) { return(vector(mode = "list")) }
  if (is.null(dim(swot_data$time)) || nrow(swot_data$time) < 3 || ncol(swot_data$time) < 3 ) { return(vector(mode = "list")) }

      swot_data$nx=1:nrow(swot_data$wse)
    # print(sos_data$Q_priors$logQ_hat)
    # bonk

    
  # Return list of remaining valid observation data
  return(list(swot_data=swot_data, sos_data=sos_data,
              invalid_nodes=invalid_nodes, invalid_times=invalid_times))

}