get_invalid_nodes_times = function(wse, slope2, time) {
    get_invalid = function(obs) {
      invalid_nodes = rowSums(is.na(obs)) >= (ncol(obs) - 3)
      invalid_times = colSums(is.na(obs)) >= (nrow(obs) - 3)
      return(list(invalid_nodes=invalid_nodes, invalid_times=invalid_times))
    }

      invalid_wse = get_invalid(wse)
      # invalid_width = get_invalid(width)
      # print(invalid_width)
      invalid_slope2 = get_invalid(slope2)
      invalid_time = get_invalid(time)
      invalid_nodes = unique(c(which(invalid_wse$invalid_nodes == TRUE),
                               which(invalid_slope2$invalid_nodes == TRUE),
                               which(invalid_time$invalid_nodes == TRUE)))
      invalid_times = unique(c(which(invalid_wse$invalid_times == TRUE),
                               which(invalid_slope2$invalid_times == TRUE),
                               which(invalid_time$invalid_times == TRUE)))
      return(list(invalid_nodes=invalid_nodes, invalid_times=invalid_times))
}