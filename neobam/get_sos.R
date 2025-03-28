#' Retrieve SOS data.
#'
#' @param sos_file string path to SOS NetCDF file
#' @param reach_id integer unique reach identifier
#'
#' @return list of Q priors
#' @export
get_sos = function(sos_file, reach_id) {
reload_sos = function(sos_file, retry_number){
  tries = retry_number
  while (tries > 0){
    tryCatch (
      {
        sos = open.nc(sos_file)
        # print("sos loaded...")
        tries = 0

      },
      error=function(e) {
              message('An Error Occurred')
              # print(e)
              tries = tries - 1
              Sys.sleep(runif(1, min=100, max=500))

          },
      warning=function(w) {
            message('A Warning Occurred')
            # print(w)
            return(NA)
        
  }
    )

  }
  return(sos)

}
  Q_priors = list()
  sos = reload_sos(sos_file, 5)



  tries = 5
  while (tries >0){
    tryCatch ( {
    # print('trying...')
    reach_grp = grp.inq.nc(sos, "reaches")$self
    rids = var.get.nc(reach_grp, "reach_id")
    index = which(rids == reach_id, arr.ind=TRUE)

    node_grp = grp.inq.nc(sos, "nodes")$self
    nrids = var.get.nc(node_grp, "reach_id")
    indexes = which(nrids == reach_id, arr.ind=TRUE)
    nids = var.get.nc(node_grp, "node_id")
    Q_priors$nids = nids[indexes]

    model_grp = grp.inq.nc(sos, "model")$self
    # print("made it to the model group")
    # print(model_grp)
    Q_priors$logQ_hat = log(var.get.nc(model_grp, "monthly_q")[,index])
  
    Q_priors$upperbound_logQ = log(var.get.nc(model_grp, "max_q")[index])
    min_q = var.get.nc(model_grp, "min_q")[index]   # Check action taken
    # if(is.null(min_q)){Q_priors$lowerbound_logQ = NA}
    if(is.na(min_q)){Q_priors$lowerbound_logQ = NA} else {
        Q_priors$lowerbound_logQ = log(min_q)
        if(min_q ==0){
            min_q=0.01
         Q_priors$lowerbound_logQ = log(min_q)}
        
    }
        
      
    

    r_grp = grp.inq.nc(sos, "gbpriors/reach")$self
    Q_priors$logQ_sd = var.get.nc(r_grp, "logQ_sd")[index]

    window_params = list()
    n_grp = grp.inq.nc(sos, "gbpriors/node")$self
    window_params$logWb_hat = var.get.nc(n_grp, "logWb_hat")[indexes]
    window_params$logWb_sd = var.get.nc(n_grp, "logWb_sd")[indexes]
    window_params$lowerbound_logWb = min(var.get.nc(n_grp, "lowerbound_logWb")[index])
    window_params$upperbound_logWb = max(var.get.nc(n_grp, "upperbound_logWb")[index])

    window_params$logDb_hat = var.get.nc(n_grp, "logDb_hat")[indexes]
    window_params$logDb_sd = var.get.nc(n_grp, "logDb_sd")[indexes]
    window_params$lowerbound_logDb = min(var.get.nc(n_grp, "lowerbound_logDb")[index])
    window_params$upperbound_logDb = max(var.get.nc(n_grp, "upperbound_logDb")[index])

    window_params$r_hat = exp(var.get.nc(n_grp, "logr_hat")[indexes])
    window_params$r_sd = exp(var.get.nc(n_grp, "logr_sd")[indexes])
    window_params$lowerbound_r = min(exp(var.get.nc(n_grp, "lowerbound_logr")[index]))
    window_params$upperbound_r = max(exp(var.get.nc(n_grp, "upperbound_logr")[index]))

    window_params$logn_hat = var.get.nc(n_grp, "logn_hat")[indexes]
    window_params$logn_sd = var.get.nc(n_grp, "logn_sd")[indexes]
    window_params$lowerbound_logn = min(var.get.nc(n_grp, "lowerbound_logn")[index])
    window_params$upperbound_logn = max(var.get.nc(n_grp, "upperbound_logn")[index])
    close.nc(sos)
    tries = 0


  },
      error=function(e) {
              message('An Error Occurred, reloading sos')
              close.nc(sos)


              # print(e)
              tries = tries - 1
              Sys.sleep(runif(1, min=10, max=600))
              sos = reload_sos(sos_file, 5)
          },
      warning=function(w) {
            message('A Warning Occurred')
            # print(w)
            return(NA)
      }
  )
}

  # print("Read was successful...")

    sos_in=open.nc(sos_file)
    sos=read.nc(sos_in,recursive=TRUE)
    close.nc(heightfile)
          
    #from the top level function to disambiguate
    reach_id_in=reach_id

    reach_index=heights$nodes$reach_id==reach_id_in
    #calculate a station vector
        x=sos$nodes$x[reach_index]
        y=sos$nodes$y[reach_index]
        A=c(x,0)
        B=c(0,x)
        C=c(y,0)
        D=c(0,y)
        #drop the first row, but keep the last column
        #this will make a matrix where the 1st - 2nd columns gives you 2-1, 3-2, etc
        #for consecutive distances
        distmatrix=as.matrix(cbind(A,B,C,D),ncol=4)[2:length(x),]

        distfunc=function(matrix_in){
            distm(c(matrix_in[['A']], matrix_in[['C']]), c(matrix_in[['B']], matrix_in[['D']]), fun = distHaversine)
        }

        #pad with a leading zero
        node_dist=c(0,apply(distmatrix,1,distfunc))

        #take a cumulative to make a station vector
        node_cum_dist =cumsum(node_dist)
    # print(node_cum_dist)

       return(list(Q_priors=Q_priors, window_params=window_params,node_cum_dist=node_cum_dist))

 }