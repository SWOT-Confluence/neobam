# H=neobam_data_object$swot_data$wse[,1]
# W=neobam_data_object$swot_data$width[,1]
# count=0
# save=list()
rejection_sample_single_node=function(seed,H,W){
    
 #all base R functions
# for(i in 1:1000){
    #sample an r, bankfull depth, and bankfull width
    r= abs(rnorm(1, mean = 1, sd = 1))
    db=abs(rnorm(1, mean = 1, sd = 1))
    #rescale to truncated normal
        x1=0 # not exactly
        x2=3 # not exactly
        y1=2
        y2=20
        m=(y2-y1)/(x2-x1)
        b=y1-m*x1
    db = db*m +b
    #wb won't go negative, but just in case
    wb=abs(rnorm(1, mean=max(W*1.2,na.rm=TRUE),sd=sd(W,na.rm=TRUE)))

    #by definition:
    #max depth is strictly less than bankfull.   
    #Iterate until you find this to be true
    max_depth=db*(W/wb)^r
    if(all(max_depth<db,na.rm=TRUE)){
       #now, we can calibrate zo to match ovserved heights
        Zo= H-max_depth
        
        #since this is a single cross section, then in theory the variation of Zo should be quite small.
        #this is different than saying the max depth is constant, since that needs to trade with H
        #test some meaningful value.
        #if we use SWOT uncertainties of W of about 30% and height of 11cm, we can say the level to which
        #this is resolved
        #the width error doens't necessarily dominate
        width_error=((mean(W,na.rm=TRUE)*0.3)^r)*(db/(wb^r))
        height_error=0.11
        if(is.na(sd(Zo,na.rm=TRUE))){return(data.frame(r=NA,
                          db=NA,
                          wb=NA,
                          Zo=NA))}
        if(sd(Zo,na.rm=TRUE) < width_error+height_error){

           
            return(data.frame(r=r,
                          db=db,
                          wb=wb,
                          Zo=mean(Zo,na.rm=TRUE)))
            # count=count+1
            # save[[count]]=list('r'=r,'db'=db,'wb'=wb,'Zo'=mean(Zo,na.rm=TRUE))
        }else{return(data.frame(r=NA,
                          db=NA,
                          wb=NA,
                          Zo=NA))}
    }else{return(data.frame(r=NA,
                          db=NA,
                          wb=NA,
                          Zo=NA))}

}

make_priors=function(cleaned_data){
    
    H=cleaned_data$Hobs
    W=cleaned_data$Wobs
    # library(parallel)
    minifunc=function(row_index,H,W){
        library(dplyr)

        output=do.call(rbind,lapply(1:500,rejection_sample_single_node,
                                    H=H[row_index,],W=W[row_index,]))%>%
        mutate(node_id=row_index)%>%
        filter(!is.na(r))

    }


    all_nodes=do.call(rbind,lapply(1:nrow(H),
                                      minifunc,
                                      H=H,
                                      W=W)) 

    #if too few nodes have systems, skip
    good_nodes=unique(all_nodes$node_id)
    #check how many
    physical_pass_rate=length(good_nodes)/cleaned_data$nx
   
    if(physical_pass_rate<0.5){
     return(list('validpriors'=FALSE,
            'nodepriors'=NA))
        }
    # stopCluster(clust)
    #what we have now is pretty dope! we have a node-by-node channel system that guarantees a few things:

    #1. bankfull depth is always greater than the instantaneous max depth
    #2. the bed elevation is always within error of the height and width obeservations
    #3. r, db, wb, and zo all exist in a physically plausible hydraulic system per noce

    #generate priors per node
    nodepriors=all_nodes%>%
        group_by(node_id)%>%
        summarize(rsd=sd(r,na.rm=TRUE),
                  dbsd=sd(db,na.rm=TRUE),
                  wbsd=sd(wb,na.rm=TRUE),
                  Zosd=sd(Zo,na.rm=TRUE),

                  r=mean(r,na.rm=TRUE),
                  db=mean(db,na.rm=TRUE),
                  wb=mean(wb,na.rm=TRUE),
                  Zo=mean(Zo,na.rm=TRUE),
                 
                  n=n())
     percent_samples_less_than_20= sum(nodepriors$n>20)/nrow(nodepriors)
        if(percent_samples_less_than_20<0.3){
             return(list('validpriors'=FALSE,
            'nodepriors'=NA))}
    
    #this throws some NA values sometimes. catch
    nodepriors$r[is.na(nodepriors$r)]=mean(nodepriors$r,na.rm=TRUE)
    nodepriors$db[is.na(nodepriors$db)]=mean(nodepriors$db,na.rm=TRUE)
    nodepriors$wb[is.na(nodepriors$wb)]=mean(nodepriors$wb,na.rm=TRUE)
    nodepriors$Zo[is.na(nodepriors$Zo)]=mean(nodepriors$Zo,na.rm=TRUE)
    nodepriors$rsd[is.na(nodepriors$rsd)]=mean(nodepriors$rsd,na.rm=TRUE)
    nodepriors$dbsd[is.na(nodepriors$dbsd)]=mean(nodepriors$dbsd,na.rm=TRUE)
    nodepriors$wbsd[is.na(nodepriors$wbsd)]=mean(nodepriors$wbsd,na.rm=TRUE)
    nodepriors$Zosd[is.na(nodepriors$Zosd)]=mean(nodepriors$Zosd,na.rm=TRUE)
    
        return(list('validpriors'=TRUE,
            'nodepriors'=nodepriors))
   

        
}
    
    
