cleanup_data=function(sos,swot_data,reach_id_in){
    source('/nas/cee-water/cjgleason/colin/neobam/neobam/nsync.R')
    library(dplyr)
    library(tidyr)

    #sos comes in by continent. Need to filter to this reach
    sos_reach_index=which(sos$reaches$reach_id == reach_id_in)
    sos_node_index=which(sos$nodes$reach_id == reach_id_in)
    

    # Q priors are based at the reach level
        #monthly Q is stored as 12 by nreach
            monthly_q= sos$model$monthly_q[,sos_reach_index]
        
        #we have a flow duration curve as 20 by n reach
        #in descending order (1st value is 96th percentile flow, 20th is 1st)
        #define Q
            probability=sos$model$probability
        #take the 66th percentile as the sd
            probability_index=which(probability == 66)
       
            qsd = sos$model$flow_duration_q[probability_index,sos_reach_index]

        #exception handling
        #if everything is NA or if any of the quanitiles are empty
         if (all(is.na(monthly_q)) || any(is.na(sos$model$flow_duration_q[,sos_reach_index])))
                                          {
             
              print('no flow duration curve')
               return(list('valid'=FALSE,
                   'Hobs'=NA,
                  'Wobs'=NA,
                  'Sobs'=NA,
                  'stationvec'=NA,
                  'nx'=NA,
                  'nt'=NA,
                  'obs_times'=NA,
                  'reach_id'=NA,
                  'Q_priors'=NA))
         }

        #20th position is 1st percentile
        qmin=sos$model$flow_duration_q[20,sos_reach_index]
        #1st position is 96th percentile
        qmax=sos$model$flow_duration_q[1,sos_reach_index]

        #we know that some monthly q exists
        monthly_q[is.na(monthly_q)]==mean(monthly_q,na.rm=TRUE)

        #build a little dataframe by month for later joining to swot data
        Q_prior_df=data.frame(month=1:12,monthly_q=monthly_q,qmin=qmin,qmax=qmax,qsd=qsd)
        
    #we need the station vector from the SoS as well
        stationvec=data.frame(stationvec=calculate_cum_dist(sos,reach_id_in),
                             xs= paste0('xs_',1:length(calculate_cum_dist(sos,reach_id_in))))
    
    if(is.null(stationvec)){
           return(list('valid'=FALSE,
                   'Hobs'=NA,
                  'Wobs'=NA,
                  'Sobs'=NA,
                  'stationvec'=NA,
                  'nx'=NA,
                  'nt'=NA,
                  'obs_times'=NA,
                  'reach_id'=NA,
                  'Q_priors'=NA))
        
        }

    #now get SWOT data
        obs_times=swot_data$reach$time_str
        reach_wse=data.frame(time=obs_times, reach_wse=swot_data$reach$wse)
        #remove outliers in time
        reach_wse_outliers=boxplot.stats(reach_wse$reach_wse)$out
        reach_wse$reach_wse[reach_wse$reach_wse %in% reach_wse_outliers]= NA

        reach_slope=data.frame(time=obs_times, reach_slope=swot_data$reach$slope2)

        #node data are nt by nx
        node_wse=data.frame(swot_data$node$wse)
        names(node_wse)=paste0('xs_',1:ncol(node_wse))

        #widths are expected to be super noisy. Do not filter
      
    node_width=data.frame(swot_data$node$width)
        names(node_width)=paste0('xs_',1:ncol(node_wse))

        node_wse_df=mutate(node_wse,time=obs_times)%>%
            gather(xs,node_wse,-time)
        node_width_df=mutate(node_width,time=obs_times)%>%
            gather(xs,node_width,-time)
 suppressWarnings({    
        #combine
    # print(nrow(node_wse_df))
    # print(nrow(node_width_df))
    # print(nrow(reach_wse))
    # print(nrow(reach_slope))
     
        swot_df=left_join(reach_wse,reach_slope,by='time')%>%
             filter(time!='no_data')%>%
             filter(!(is.na(reach_wse) & is.na(reach_slope)))%>%
   
     
            left_join(node_wse_df,by='time')%>%
            left_join(node_width_df,by=c('time','xs')) %>%
            group_by(xs)%>%
            #filter node wse outliers in time at each node
            mutate(node_wse= ifelse(node_wse %in% boxplot.stats(node_wse)$out,NA,node_wse))

})
    
 
        swot_df$node_wse[swot_df$node_wse=='NaN']=NA
        swot_df$node_width[swot_df$node_width=='NaN']=NA
        swot_df$reach_wse[swot_df$reach_wse=='NaN']=NA
    #first cleanup the node heights using the reach heights
        # 1. smooth each of the reach and node heights over time using
        #    the zoo package 'rollmean' function
        # 2. classify the hydrograph into rising, falling, peak, trough
        # 3. calculate the disagreement between those classifications of
        #     the smoothed hydrographs
        # the idea here is that the reach data are higher quality and reflect
        # the actual dynamics, so we only keep nodes that are 'in step' with
        # the reach dynamics
        xs_ids=unique(swot_df$xs)
   
    b=Sys.time()
   

    
    swot_df2=do.call(rbind,lapply(xs_ids,nsync,swot_df=swot_df))
    
    if(is.null(swot_df2)){
          return(list('valid'=FALSE,
                   'Hobs'=NA,
                  'Wobs'=NA,
                  'Sobs'=NA,
                  'stationvec'=NA,
                  'nx'=NA,
                  'nt'=NA,
                  'obs_times'=NA,
                  'reach_id'=NA,
                  'Q_priors'=NA))
        }
    swot_df2=swot_df2%>%
                   left_join(stationvec,by='xs')%>%
            filter(agreement>0.8)
    
    # print('agreement check')

        #sweet! some testing needs to be done with this 0.5 number
        #now, check how many xs are left
        nx=length(unique(swot_df2$xs))
        if (nx<3){
            
             print('not enough nodes')
            # print(unique(swot_df2$xs))
            
              return(list('valid'=FALSE,
                   'Hobs'=NA,
                  'Wobs'=NA,
                  'Sobs'=NA,
                  'stationvec'=NA,
                  'nx'=NA,
                  'nt'=NA,
                  'obs_times'=NA,
                  'reach_id'=NA,
                  'Q_priors'=NA))
                 }

       
        
        #count how many times we have repeating obs across the same nodes
        count_square_obs=swot_df2%>%
            group_by(time)%>%
            #this number is how many xs have data for a given time
            summarize(goodcount=sum(!is.na(node_wse)))

        #join the count to the original datagrame
        swot_df3=swot_df2%>%
            left_join(count_square_obs,by='time')%>%
        #filter by obs nx >5
            filter(goodcount>5)#%>%
        # #calculate a slope from the observed heights
        #     group_by(time,xs)%>%
        #     mutate(Sact=lm(node_wse~stationvec)$coefficients[2])%>%
        # #these slopes shoudl be NEGATIVE as they go downhill
        #     filter(Sact<0)
   


        #drop if empty
        if(nrow(swot_df3)==0){
            
              print('not enough time/node overlap')
               return(list('valid'=FALSE,
                   'Hobs'=NA,
                  'Wobs'=NA,
                  'Sobs'=NA,
                   'stationvec'=NA,
                  'nx'=NA,
                  'nt'=NA,
                  'obs_times'=NA,
                  'reach_id'=NA,
                  'Q_priors'=NA))
        }
        
        #nx and nt are now defined
        nx=length(unique(swot_df3$xs))
        nt=length(unique(swot_df3$time))
    

        
        #go back to obs matrices
        Hobs=select(swot_df3,time,xs,node_wse)%>%
            pivot_wider(names_from=time,values_from=node_wse)
        
        Wobs=select(swot_df3,time,xs,node_width)%>%
            pivot_wider(names_from=time,values_from=node_width)

        # Sobs=select(swot_df3,time,xs,reach_slope)%>%
        #     pivot_wider(names_from=time,values_from=reach_slope)
                              
        stationvec_df=select(swot_df3,stationvec,xs)%>%
            group_by(xs)%>%
            summarize(station_vec=first(stationvec))%>%
    #this is now out of order, so we need to substract the first value to make a zero based one for Zo in gureu
            mutate(station_vec=station_vec-min(station_vec))
    
    
        
        stationvec=stationvec_df

        timeHobs=names(Hobs)[2:ncol(Hobs)]
        timeWobs=names(Wobs)[2:ncol(Wobs)]
        # timeSobs=names(Sobs)[2:ncol(Sobs)]

        #check to make sure the ordering is correct. It shoudl be, but...
        # check1=all(timeHobs == timeSobs)
        check2=all(timeHobs == timeWobs)
        
        if( check2 ==FALSE){
             print('time mapping not consistent')
               return(list('valid'=FALSE,
                   'Hobs'=NA,
                  'Wobs'=NA,
                  'Sobs'=NA,
                   'stationvec'=NA,
                  'nx'=NA,
                  'nt'=NA,
                  'obs_times'=NA,
                  'reach_id'=NA,
                  'Q_priors'=NA))
            }
        
        
        obs_times=timeHobs
        #we have now defined Hobs, Wobs, Sobs, nx, nt, obs_times, stationvec, reach_id
        #and Q priors. return these
    
    Hobs=as.matrix(Hobs[,2:ncol(Hobs)])
    Wobs=as.matrix(Wobs[,2:ncol(Wobs)])
    # Sobs=as.matrix(Sobs[,2:ncol(Sobs)])
    
    #the system will fall apart if slope doesn't work with teh heights as measured
         slopefunc=function(H,stationvec){
                return(lm(H~stationvec)$coefficients[2])
            }

            #make Zo station vectored
            # Zo=Zo1 + S*stationvec
     Sact=apply(Hobs,2,slopefunc,stationvec=stationvec_df$station_vec)
    

    
    #repeat it over rows
    Sobs=matrix(rep(Sact,times=nrow(Hobs)),nrow=nrow(Hobs),ncol=ncol(Hobs),byrow=TRUE)

  print(Q_prior_df)
    bonk
    return(list('valid'=TRUE,
                  'Hobs'=Hobs,
                  'Wobs'=Wobs,
                  'Sobs'=Sobs,
                  'stationvec'=stationvec,
                  'nx'=nx,
                  'nt'=nt,
                  'obs_times'=obs_times,
                  'reach_id'=reach_id_in,
                  'Q_priors'=Q_prior_df))
}
    
    

