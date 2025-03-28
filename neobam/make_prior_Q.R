
make_prior_Q=function(input_file){


neobam_data_object=readRDS(input_file)
# neobam_data_object=readRDS('/nas/cee-water/cjgleason/colin/neobam20_files/input/23229000521SWOT_neobam_input.rds')

this_reach_id=neobam_data_object$reach_id[1]

    H=neobam_data_object$Hobs
    S=neobam_data_object$Sobs
    Zo1=neobam_data_object$nodepriors$nodepriors$Zo[1]
    r=neobam_data_object$nodepriors$nodepriors$r
    wb=neobam_data_object$nodepriors$nodepriors$wb
    db=neobam_data_object$nodepriors$nodepriors$db
    n=0.03
    stationvec1=neobam_data_object$stationvec$station_vec
    stationvec=matrix(rep(stationvec1,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
    r=matrix(rep(r,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
    wb=matrix(rep(wb,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
    db=matrix(rep(db,times=ncol(H)),nrow=nrow(H),ncol=ncol(H))
    
    #make Zo station vectored
            Zo=Zo1 + S*stationvec
    #can use an abs on the slope, they are all nevative
    Q=abs(S)^(0.5) * (H-Zo)^(1.66 +(1/r))*db^(-1/r)*(wb)*(1/n)*(r/(r+1))^(1.66)

    Qmean=apply(Q,2,mean,na.rm=TRUE)
    Qsd=apply(Q,2,sd,na.rm=TRUE)


return(data.frame(priorQ=Qmean,priorsd=Qsd,reach_id=this_reach_id,date=as.Date(names(Qmean))))
    
    }