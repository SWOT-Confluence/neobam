
run_neobam_stan = function(neobam_data_and_priors,sourcefile){
    
  library(BH,quietly=TRUE,warn.conflicts = FALSE)
  library(rstan,quietly=TRUE,warn.conflicts = FALSE)
  library(dplyr,quietly=TRUE,warn.conflicts = FALSE)
  library(tidyr,quietly=TRUE,warn.conflicts = FALSE)
  library(stringr,quietly=TRUE,warn.conflicts = FALSE)
  library(reshape2)
    
  rstan_options(auto_write = TRUE)


    
  return_posterior_mean_sd= function(parameter,stan_output){
      
      
    
    var_names=gsub("\\[.*\\]" ,"",   row.names(stan_output) )
    parameter_index=!is.na(str_match(var_names,parameter))
    string_length_index=unlist(lapply(var_names,nchar)) == nchar(parameter)
    
    combined_index=which(parameter_index*string_length_index==1)
      
      #print(stan_output[combined_index,])
    
    posterior_mean=stan_output[combined_index,]$mean
    posterior_sd=mean(stan_output[combined_index,]$sd^2)
    
    return(list('mean'=posterior_mean,'sd'=posterior_sd))
    
  }
    
   
    
  
iter=neobam_data_and_priors$iter
    

        
  fit1 <- stan(
    
    file = sourcefile,  # Stan program
    data = neobam_data_and_priors,    # named list of data
    chains = 3,             # number of Markov chains
    warmup = floor(0.1*iter),          # number of warmup iterations per chain
    iter = iter,   # total number of iterations per chain
    cores = 1,              # number of cores (could use one per chain)
    refresh = 0
    #        boost_lib='/home/cjgleason_umass_edu/R/x86_64-pc-linux-gnu-library/4.0',
   #   eigen_lib='/home/cjgleason_umass_edu/R/x86_64-pc-linux-gnu-library/4.0'
  )
  
  output= data.frame(rstan::summary(fit1)$summary)
    

  #OK! we now need a list of all the posteriors so we can return them to stuff back into the flow law
  # posterior_list=c('r','logn','betaslope','betaint','logbeta','logQ')
     posterior_list=c('r','logn','logWb','logDb','logQ','Zo1','Zo','logQ_transform')
    
  posteriors=lapply(posterior_list,return_posterior_mean_sd,output)
  names(posteriors)=posterior_list
    
    hydrograph_vector=posteriors$logQ_transform$mean

        replacements=which(as.vector(neobam_data_and_priors$hasdat==FALSE)) -1
        values=NA
        for (i in replacements){
            hydrograph_vector=append(hydrograph_vector, values, after =i)
         }

        hydromat=matrix(hydrograph_vector,
                          nrow=neobam_data_and_priors$nx,
                          ncol=neobam_data_and_priors$nt)
 

        hydrograph_posterior=apply(exp(hydromat),2,mean,na.rm=TRUE)
        Qsd=apply(exp(hydromat),2,sd,na.rm=TRUE)
    
    
  return(list('posteriors'=posteriors,
              'hydrograph_posterior'=hydrograph_posterior,
              'hydrograph_post_sd'=Qsd))

}

remake_discharge =function (H,S,stationvec,posteriors){
    

  

    Zo1=posteriors$Zo1$mean
    r=posteriors$r$mean
    wb=exp(posteriors$logWb$mean)
    db=exp(posteriors$logDb$mean)
    n=0.03
  
    r=matrix(rep(r,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
    wb=matrix(rep(wb,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
    db=matrix(rep(db,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))

    
    #make Zo station vectored
            Zo=Zo1 + S*stationvec
    #can use an abs on the slope, they are all nevative
    Q=abs(S)^(0.5) * (H-Zo)^(1.66 +(1/r))*db^(-1/r)*(wb)*(1/n)*(r/(r+1))^(1.66)
 
    

    Qmean=apply(Q,2,mean,na.rm=TRUE)
    Qsd=apply(Q,2,sd,na.rm=TRUE)

  
    return(list('Q'=Qmean,
                 'sd'=Qsd))

  
}

norm_to_lognorm=function(mu,sigma){
 
    lognorm_mu=2*log(mu)-0.5*log(mu^2+sigma^2)
    lognorm_sigma= log(mu^2 + sigma^2) - 2*(log(mu))
    return(list('mu'=lognorm_mu,'sigma'=lognorm_sigma))
    
}








