#' neoBAM input operations
#'
#' SOS: max_q, mean_q, min_q
#' SWOT: node/time, node/width, node/slope2

# Libraries
library(RNetCDF,lib.loc='/home/cjgleason_umass_edu/.conda/pkgs/r-rnetcdf-2.6_2-r42h498a2f1_0/lib/R/library/',quietly=TRUE,warn.conflicts = FALSE)
#' Title
#'
#' @param swot_file string path to SWOT NetCDF file
#' @param sos_file string path to SOS NetCDF file
#' @param reach_id integer unique reach identifier
#'
#' @return named list of data needed to run neoBAM
#' @export
get_input = function( sos_file,set_json,set_index) {


  # Get SWOT
  swot_data = get_swot(set_json,set_index)
   
    if(typeof(swot_data) != "list"){
        return(valid=FALSE)}

  # Get Qpriors
  Q_priors_all = lapply(swot_data$width$reach_id,get_Qpriors,sos_file=sos_file)
  Q_priors= all_to_one_Q_priors(Q_priors_all)
    #these come out one per reach. average over the cycle


    width=swot_data$width%>%
    select(sort(names(swot_data$width)))

slope=swot_data$slope%>%
    select(sort(names(swot_data$slope)))
    
    
    
    
  return(list('reach_id'=Q_priors$reach_ids,'width'=width,'slope'=slope,'Q_priors'=Q_priors,'time'=names(width)[1:(ncol(width)-1)]))

 }

get_swot = function(set_json,set_index) {
    

  json_data = rjson::fromJSON(file=file.path(set_json))

count=0
reachid=0
setid=0
# json_data
# json_data[[1]][[1]]
for (i in 1:length(json_data)){
    this_set=json_data[[i]]
    # print(length(this_set))
    for (j in 1:length(this_set)){
        count=count+1
        # print(json_data[[i]][[j]])
        reachid[count]=this_set[[j]]$reach_id  
        setid[count]=i
    }
 }

        sets=data.frame('set_id'=setid,'reach_id'=reachid)

get_swot_from_set=function(setid,setdf){
    library(dplyr)
    library(tidyr)
    

    this_set=dplyr::filter(setdf,set_id==setid)
 
      thisset=data.frame('width'=NA,'slope'=NA,'time'=NA,
                       'reach_id'=NA, 'setid'=NA,'nt'=NA)
    
 
        
    for (reach_id in this_set$reach_id){
         swot_file=paste0('/nas/cee-water/cjgleason/SWOT_Q_UMASS/mnt/input/swot/',reach_id,'_SWOT.nc') 


            if(!file.exists(swot_file)){return(valid=FALSE)}  

          swot = open.nc(swot_file)        
          nx = var.get.nc(swot, "nx")
          nt = var.get.nc(swot, "nt")

          reach_grp = grp.inq.nc(swot, "reach")$self
          width = var.get.nc(reach_grp, "width")
          slope = var.get.nc(reach_grp, "slope2")
          time = var.get.nc(reach_grp, "time")


          thisset=rbind(thisset,data.frame('width'=width,'slope'=slope,'time'=time,
                               'reach_id'=reach_id, 'setid'=setid,'nt'=length(width)))

            close.nc(swot)
        
        }
    
    # print(nrow(thisset))

    
        if(nrow(thisset)<3){return(valid=FALSE)}


  width_matrix=select(thisset,width,time,reach_id)%>%
    filter(!is.na(width))%>%
    filter(!is.na(time))%>%
    mutate(time=as.integer(time))%>%
     pivot_wider(names_from = time, values_from = width)
   
  slope_matrix=select(thisset,slope,time,reach_id)%>%
    filter(!is.na(slope))%>%
    filter(!is.na(time))%>%
    mutate(time=as.integer(time))%>%
    mutate(slope=ifelse(slope<0,NA,slope))%>%
     pivot_wider(names_from = time, values_from = slope)
    
 
  
    return(list(width=width_matrix,slope=slope_matrix))
 
}
      

return(get_swot_from_set(setid=set_index,setdf=sets))


}



#' Retrieve SOS data.
#'
#' @param sos_file string path to SOS NetCDF file
#' @param reach_id integer unique reach identifier
#'
#' @return list of Q priors
#' @export
get_Qpriors = function(sos_file, reach_id) {

  Q_priors = list()

  sos=open.nc(sos_file)

    reach_grp = grp.inq.nc(sos, "reaches")$self
    rids = var.get.nc(reach_grp, "reach_id")
    index = which(rids == reach_id, arr.ind=TRUE)

    node_grp = grp.inq.nc(sos, "nodes")$self
    nrids = var.get.nc(node_grp, "reach_id")
    indexes = which(nrids == reach_id, arr.ind=TRUE)

    model_grp = grp.inq.nc(sos, "model")$self
    # print("made it to the model group")
    # print(model_grp)
    Q_priors$logQ_hat = log(var.get.nc(model_grp, "mean_q")[index])
    Q_priors$upperbound_logQ = log(var.get.nc(model_grp, "max_q")[index])
    min_q = var.get.nc(model_grp, "min_q")[index]   # Check action taken
    if ((min_q < 0) | (is.na(min_q))) {
      Q_priors$lowerbound_logQ = NA
    } else {
      Q_priors$lowerbound_logQ = log(min_q)
    }

    r_grp = grp.inq.nc(sos, "gbpriors/reach")$self
    Q_priors$logQ_sd = var.get.nc(r_grp, "logQ_sd")[index]

   
  close.nc(sos)


  # print("Read was successful...")

  return(list(Q_priors=Q_priors,reach_id=reach_id))

}

check_valid=function(data){
    date=data$time
    Sobs=data$slope
    Wobs=data$width
    
    if(nrow(Sobs) <3){return(valid=FALSE)}
    if(nrow(Wobs) <3){return(valid=FALSE)}
    if(ncol(Sobs) <3){return(valid=FALSE)}
    if(ncol(Wobs) <3){return(valid=FALSE)}
    if(sum(is.na(Sobs)) > ((ncol(Sobs)*nrow(Sobs))*0.4)){return(valid=FALSE)}
    if(sum(is.na(Wobs)) > ((ncol(Wobs)*nrow(Wobs))*0.4)){return(valid=FALSE)}
    
  
    if(all(names(Sobs) == names(Wobs)) != TRUE){return(valid=FALSE)}
    return(valid=TRUE)
    
    
    }
    
   all_to_one_Q_priors= function(Q_priors_all){
       #comes in as a list
       
       process_one_level=function(this_row){
    
           
             output=  data.frame('logQ_hat'=this_row$Q_priors$logQ_hat,
                  'logQ_sd'=this_row$Q_priors$logQ_sd,
                  'lowerbound_logQ'=this_row$Q_priors$lowerbound_logQ,
                  'upperbound_logQ'=this_row$Q_priors$upperbound_logQ,
                        'reach_id'=this_row$reach_id)
           }
       
       all_Q_df=do.call(rbind,lapply(Q_priors_all,process_one_level))
       
       Q_priors=all_Q_df%>%
       dplyr::select(-reach_id)%>%
       summarize(logQ_hat=mean(logQ_hat),
                 logQ_sd=mean(logQ_sd),
                 lowerbound_logQ=mean(lowerbound_logQ),
                 upperbound_logQ=mean(upperbound_logQ))
       

       
       return(list('reach_ids'=dplyr::select(all_Q_df,reach_id),'Q_priors'=Q_priors) )
       
       
       }
    
    
