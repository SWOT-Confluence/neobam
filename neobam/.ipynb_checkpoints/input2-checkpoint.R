
get_input2 = function(swot_file, sos_file, reach_id_in) {
    suppressMessages({library(RNetCDF)})
    source('/nas/cee-water/cjgleason/colin/neobam/neobam/calculate_cum_dist.R')
    source('/nas/cee-water/cjgleason/colin/neobam/neobam/cleanup_data.R')
    source('/nas/cee-water/cjgleason/colin/neobam/neobam/rejection_sample.R')


  
  #getSWOT
    swot_in=open.nc(swot_file)
    swot_data=read.nc(swot_in,recursive=TRUE)
    close.nc(swot_in)

  # Get SOS
    sos_in=open.nc(sos_file)
    sos=read.nc(sos_in,recursive=TRUE)
    close.nc(sos_in)
          
  #calculate node cumulative distance
 
    # print('open file')

    node_cum_dist=calculate_cum_dist(sos,reach_id_in)

    # print('cumulative dist')

  #data cleanup.
    #1. filter the node data based on the reach data
    #2. check for remaining dates and times with enough data to invert
    #3. trim all the data to ensure that inputs are the same size
    cleaned_data = cleanup_data(sos=sos, 
                                swot_data=swot_data,
                                reach_id_in=reach_id_in)


 
  # print('cleaned data')
    # print(cleaned_data)

    if(cleaned_data$valid==FALSE){
            cleaned_data$Q_priors=NA
            cleaned_data$nodepriors=NA
            return(cleaned_data)
        }
    #now we've got cleaned data, we need to do some rejection sampling
    #and Q prior tighening. This makes the input ops longer but the actual
    #run of the sampler bits shorter
    nodepriors=make_priors(cleaned_data)
    
    if(nodepriors$validpriors==FALSE){
            cleaned_data$valid=FALSE
            cleaned_data$Q_priors=NA
            cleaned_data$nodepriors=NA
            return(cleaned_data)
        }
       
    Q_priors=cleaned_data$Q_priors%>%
    #two digit month vs one digit mmonth.....
        mutate(month=as.character(month))%>%
        mutate(month=ifelse(nchar(month)==1,paste0(0,month),month))

  
    obs_month=data.frame(month=format(as.Date(cleaned_data$obs_times), "%m"))%>%
      left_join(Q_priors,by='month')
    
    Q_hat=obs_month$monthly_q
    Q_sd=obs_month$qsd
    qmin=obs_month$qmin[1]
    qmax=obs_month$qmax[1]
    
    Q_priors=list('Q_hat'=Q_hat,
                  'Q_sd'=Q_sd,
                  'Q_min'=qmin,
                  'Q_max'=qmax)
    
    cleaned_data$Q_priors=Q_priors
    cleaned_data$nodepriors=nodepriors
   
       
      return( cleaned_data)
 }