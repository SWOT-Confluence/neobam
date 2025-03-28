process_data = function(cleaned_data, stan_file) {


    source('/nas/cee-water/cjgleason/colin/neobam/neobam/neobam_functions.R')

  if(cleaned_data$valid==FALSE){
 output=list('posterior_Q'= NA,
                'neobam_time'= NA,
                'posterior_Q_sd'= NA,
                'posteriors'=NA,
                'diditrun'=FALSE)
     return(output)}
    
    

#this function does data prep and generate priors
  neobam_parameters=list(
    reachid=cleaned_data$reach_id,
    date=cleaned_data$obs_times,
    Sobs=cleaned_data$Sobs,
    # Wobs=data$swot_data$width,
    Hobs=cleaned_data$Hobs,

    #performative variables
    # Werr_sd=    0.3*mean(data$swot_data$width,na.rm=T), #30% width error, version D
    Herr_sd=    0.11, #11cm, version D
    Serr_sd=    0.000018, #changed december 2024 to mathc version D stats
    #logWerr_sd=  norm_to_lognorm(mean(data$swot_data$width,na.rm=T),20)$sigma,
    logSerr_sd=  norm_to_lognorm(abs(mean(cleaned_data$Sobs,na.rm=T)),0.000018)$sigma,
    logHerr_sd=  norm_to_lognorm(mean(cleaned_data$Hobs,na.rm=T),0.11)$sigma,
    iter=   4000
)

    # print(mean(neobam_parameters$Sobs,na.rm=TRUE))
    # print(neobam_parameters$logSerr_sd)
#assign priors

    #log Q prior taken directly from the priors
    #monthly Q prior as of December 2024
   neobam_parameters$logQ_hat=log(cleaned_data$Q_priors$Q_hat)
   neobam_parameters$lowerbound_logQ= log(cleaned_data$Q_priors$Q_min)
   neobam_parameters$upperbound_logQ=log(cleaned_data$Q_priors$Q_max)
   neobam_parameters$logQ_sd= log(cleaned_data$Q_priors$Q_sd)
     
   neobam_parameters$logQ_sd[neobam_parameters$logQ_sd <0] =1
 # neobam_parameters$logQ_sd=matrix(rep(neobam_parameters$logQ_sd,times=nrow(neobam_parameters$Hobs)),
 #                                  nrow=nrow(neobam_parameters$Hobs),ncol=ncol(neobam_parameters$Hobs),byrow=TRUE)
    
        neobam_parameters$r_hat= cleaned_data$nodepriors$nodepriors$r# new_priors$r_hat
        neobam_parameters$r_sd=  cleaned_data$nodepriors$nodepriors$rsd #new_priors$r_sd
        neobam_parameters$lowerbound_r=min(cleaned_data$nodepriors$nodepriors$r)
        neobam_parameters$upperbound_r=max(cleaned_data$nodepriors$nodepriors$r)

        neobam_parameters$logWb_hat= log(cleaned_data$nodepriors$nodepriors$wb)# new_priors$r_hat
        neobam_parameters$logWb_sd=  norm_to_lognorm(cleaned_data$nodepriors$nodepriors$wb,
                                                          cleaned_data$nodepriors$nodepriors$wbsd)$sigma #new_priors$r_sd
        neobam_parameters$lowerbound_logWb=min(log(cleaned_data$nodepriors$nodepriors$wb))
        neobam_parameters$upperbound_logWb=max(log(cleaned_data$nodepriors$nodepriors$wb))

        neobam_parameters$logDb_hat= log(cleaned_data$nodepriors$nodepriors$db)# new_priors$r_hat
        neobam_parameters$logDb_sd= norm_to_lognorm(cleaned_data$nodepriors$nodepriors$db,
                                                         cleaned_data$nodepriors$nodepriors$dbsd)$sigma #new_priors$r_sd
        neobam_parameters$lowerbound_logDb=min(log(cleaned_data$nodepriors$nodepriors$db))
        neobam_parameters$upperbound_logDb=max(log(cleaned_data$nodepriors$nodepriors$db))

        neobam_parameters$Zo1_hat= mean(cleaned_data$nodepriors$nodepriors$Zo)# new_priors$r_hat
        neobam_parameters$Zo1_sd=  mean(cleaned_data$nodepriors$nodepriors$Zosd) #new_priors$r_sd
        neobam_parameters$lowerbound_Zo1=min(cleaned_data$nodepriors$nodepriors$Zo)
        neobam_parameters$upperbound_Zo1=max(cleaned_data$nodepriors$nodepriors$Zo)

        neobam_parameters$logn_hat= rep(log(0.03),times=nrow(neobam_parameters$Hobs))# new_priors$r_hat
        neobam_parameters$logn_sd=rep(0.001,times=nrow(neobam_parameters$Hobs)) #new_priors$r_sd
        neobam_parameters$lowerbound_logn=log(0.028)
        neobam_parameters$upperbound_logn=log(0.032)
    


        #drop NA for stan
        # neobam_parameters$Wobs[is.na(neobam_parameters$Wobs)]=0
        neobam_parameters$Sobs[is.na(neobam_parameters$Sobs)]=0
        neobam_parameters$Hobs[is.na(neobam_parameters$Hobs)]=0
    
        #set nx and nt
        neobam_parameters$nx=cleaned_data$nx
        neobam_parameters$nt=cleaned_data$nt
        
    

        #set stan-required ntot
        neobam_parameters$ntot=sum( neobam_parameters$hasdat)

          neobam_parameters$sigma_man=matrix(0.26,
                         nrow=neobam_parameters$nx,ncol=neobam_parameters$nt)
    
        #add the cumulative sum for bed regularization
      stationvec=cleaned_data$stationvec$station_vec
        neobam_parameters$stationvec=matrix(rep(stationvec,
                                                times=neobam_parameters$nt),
                                            nrow=neobam_parameters$nx,
                                            ncol=neobam_parameters$nt)
    

#stan can't handle negative slopes
    #if we got here, we can try
    neobam_parameters$Sobs=abs(neobam_parameters$Sobs)

    
        #set stan-required hastdat
neobam_parameters$hasdat= unname(as.matrix((neobam_parameters$Hobs * neobam_parameters$Sobs)>0,nrow=neobam_parameters$nx,ncol=neobam_parameters$nt))
neobam_parameters$ntot=sum(neobam_parameters$hasdat)

 
   #run the MCMC 
    fit_out = run_neobam_stan(neobam_parameters,stan_file)
    
    # print(fit_out)

        #pull the posteriors
    posteriors=fit_out$posteriors
    hydrograph_posterior=fit_out$hydrograph_posterior
    hydrograph_posterior_sd= fit_out$hydrograph_post_sd

     hydrograph_recon= remake_discharge(H=neobam_parameters$Hobs,
                                        S=neobam_parameters$Sobs,
                                        stationvec=neobam_parameters$stationvec,
                                        posteriors=posteriors) #lin space

        #write to file
    output=list('posterior_Q'= hydrograph_posterior,
                'recon_Q'=hydrograph_recon$Q,
                'recon_sd'=hydrograph_recon$sd,
                'neobam_time'= neobam_parameters$date,
                'posterior_Q_sd'= hydrograph_posterior_sd,
                'posteriors'=posteriors,
                'diditrun'=TRUE)


    return(output)

  }


