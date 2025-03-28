optimizeneobam=function(input_file){
source('/nas/cee-water/cjgleason/colin/neobam/neobam/make_prior_Q.R')
    suppressMessages({
library('optimx')
library('dplyr')})

    neo_input_object=readRDS(input_file)
    
 remake_discharge_opt =function (priors,H,S,stationvec,qmin,qmax){

        Zo1=priors[(3*nrow(H)+1):(4*nrow(H))]
        r=priors[1:nrow(H)]
        wb=priors[(2*nrow(H)+1):(3*nrow(H))]
        db=priors[(nrow(H)+1):(2*nrow(H))]

    #      print(r)
    #      print(db)
    #      print(wb)
    #      print(Zo1)
    #      bonk
        n=0.03

        r=matrix(rep(r,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
        wb=matrix(rep(wb,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
        db=matrix(rep(db,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
        Zo1=matrix(rep(Zo1,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))

        #make Zo station vectored
        Zo=Zo1 + S*stationvec
        #can use an abs on the slope, they are all nevative
        Q=abs(S)^(0.5) * (H-Zo)^(1.66 +(1/r))*db^(-1/r)*(wb)*(1/n)*(r/(r+1))^(1.66)
         Q[Q=='NaN']=NA

         # if(all(is.na(Q))){return(10000000)}

         # print(Q)
        Qmean=apply(Q,2,mean,na.rm=TRUE)
        # Qsd=apply(Q,2,sd,na.rm=TRUE)
        error=apply(Q,1,function(x){(x-Qmean)/Qmean})

           #don't let Q go too small or too big

        error[which(Q>qmax)]=5
        error[which(Q<qmin)]=5


         #Jeff style tulip function
         # error[abs(error)>2]=NA
         # print(error)
         loss=sum(error^2,na.rm=TRUE)

         percent_NA=sum(is.na(error))/(nrow(error)*ncol(error))


         return(loss)

  
 }

 remake_discharge =function (posteriors,H,S,stationvec){

        Zo1=posteriors[(3*nrow(H)+1):(4*nrow(H))]
        r=posteriors[1:nrow(H)]
        wb=posteriors[(2*nrow(H)+1):(3*nrow(H))]
        db=posteriors[(nrow(H)+1):(2*nrow(H))]


    #      print(r)
    #      print(db)
    #      print(wb)
    #      print(Zo1)
    #      bonk
        n=0.03

        r=as.numeric(matrix(rep(r,times=ncol(H)),nrow=nrow(H),ncol=ncol(H)))
        wb=as.numeric(matrix(rep(wb,times=ncol(H)),nrow=nrow(H),ncol=ncol(H)))
        db=as.numeric(matrix(rep(db,times=ncol(H)),nrow=nrow(H),ncol=ncol(H)))
        Zo1=as.numeric(matrix(rep(Zo1,times=ncol(H)),nrow=nrow(H),ncol=ncol(H)))


        #make Zo station vectored
        Zo=Zo1 + (S*stationvec)
        #can use an abs on the slope, they are all nevative
        Q=abs(S)^(0.5) * (H-Zo)^(1.66 +(1/r))*db^(-1/r)*(wb)*(1/n)*(r/(r+1))^(1.66)
         Q[Q=='NaN']=NA

         #don't let Q go too small or too big


        Qmean=apply(Q,2,mean,na.rm=TRUE)
        Qsd=apply(Q,2,sd,na.rm=TRUE)
        date=as.Date(names(Qmean))


         return(data.frame(posteriorsd=Qsd,
                                posteriorQ=Qmean,
                                date=date))


 }
    
   
#input ops-----------------
# names(neo_input_object$Q_priors)
priors=neo_input_object$nodepriors$nodepriors%>%
    select(r,db,wb,Zo)

# names(priors)
stationvec=matrix(rep(neo_input_object$stationvec$station_vec,
                                            times=neo_input_object$nt),
                                            nrow=neo_input_object$nx,
                                            ncol=neo_input_object$nt)
#-----------------------------
    
#----------------------------
 # remake_discharge_opt(transpriors,H=neo_input_object$Hobs,S=neo_input_object$Sobs,stationvec=stationvec,transpriors=c(priors$r,priors$db,priors$wb,priors$Zo)
#                       qmin=neo_input_object$Q_priors$Q_min,
#                       qmax=neo_input_object$Q_priors$Q_max)

#-----------------------------

#dfeinte lower bounds
lower=priors%>%
  summarize(r=0.01,
              db=min(db,na.rm=TRUE),#-
              # 2*mean(neo_input_object$nodepriors$nodepriors$dbsd),
              wb=min(wb,na.rm=TRUE),#-
              # 2*mean(neo_input_object$nodepriors$nodepriors$wbsd),
              Zo=min(Zo,na.rm=TRUE))%>%#-
              # 2*mean(neo_input_object$nodepriors$nodepriors$Zosd))%>%
as.double()

lower=c(rep(lower[1],nrow(neo_input_object$Hobs)),
        rep(lower[2],nrow(neo_input_object$Hobs)),
        rep(lower[3],nrow(neo_input_object$Hobs)),
        rep(lower[4],nrow(neo_input_object$Hobs)))
#------------------------------
    
#define upper bounds---------
upper=priors%>%
    summarize(r=max(r,na.rm=TRUE),#+
              # 2*mean(neo_input_object$nodepriors$nodepriors$rsd),
              db=max(db,na.rm=TRUE),#+
              # 2*mean(neo_input_object$nodepriors$nodepriors$dbsd),
              wb=max(wb,na.rm=TRUE),#+
              # 2*mean(neo_input_object$nodepriors$nodepriors$wbsd),
              Zo=max(Zo,na.rm=TRUE))%>%#+
              # 2*mean(neo_input_object$nodepriors$nodepriors$Zosd))%>%
as.double()

upper=c(rep(upper[1],nrow(neo_input_object$Hobs)),
        rep(upper[2],nrow(neo_input_object$Hobs)),
        rep(upper[3],nrow(neo_input_object$Hobs)),
        rep(upper[4],nrow(neo_input_object$Hobs)))
#------------------------------
 
    #transform the priors-----------------------------
transpriors=c(priors$r,priors$db,priors$wb,priors$Zo)
    #-----------------------------

#need a lot of parameters to account for space differences
    
    suppressWarnings({
optparams=optimx(transpriors, fn=remake_discharge_opt,
            H=neo_input_object$Hobs,
            S=neo_input_object$Sobs,
           stationvec=stationvec,
           qmin=neo_input_object$Q_priors$Q_min,
           qmax=neo_input_object$Q_priors$Q_max,
          method='L-BFGS-B',
           lower=lower,
            upper=upper)})
#---------------------------
    
#remake Q from params
H=neo_input_object$Hobs
newpar2=c(optparams[1:nrow(H)],
          optparams[(nrow(H)+1):(2*nrow(H))],
          optparams[(2*nrow(H)+1):(3*nrow(H))],
          optparams[(3*nrow(H)+1):(4*nrow(H))])

output= remake_discharge(posteriors=newpar2,
                         H=neo_input_object$Hobs,
                         S=neo_input_object$Sobs,
                         stationvec=stationvec)
    
    return(mutate(output,reach_id=neo_input_object$reach_id[1]))
}

