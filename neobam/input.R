#' neoBAM input operations
#'
#' SOS: max_q, mean_q, min_q
#' SWOT: node/time, node/width, node/slope2

# Libraries
library(RNetCDF,quietly=TRUE,warn.conflicts = FALSE)
library(dplyr,quietly=TRUE,warn.conflicts = FALSE)
#' Title
#'
#' @param swot_file string path to SWOT NetCDF file
#' @param sos_file string path to SOS NetCDF file
#' @param reach_id integer unique reach identifier
#'
#' @return named list of data needed to run neoBAM
#' @export
get_input = function(set_json,set_index) {


  # Get SWOT
  swot_data = get_swot(set_json,set_index)
  print("here is the swot data we are looking for")
  print(swot_data$sos[4])
  print("There was the swot data")
  sos_file = paste0("/mnt/data/input/sos/",swot_data$sos, sep="")
  swot_data = get_swot_from_set(set_index, swot_data)
   
  if(typeof(swot_data) != "list"){
      return(valid=FALSE)}

  # Get Qpriors
  Q_priors_all = lapply(swot_data$width$reach_id,get_Qpriors,sos_file=sos_file)
  print("--------------------start------------------")
  print("here is q priors all")
  # print(Q_priors_all)
    print("--------------------end------------------")

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
  for (i in 1:length(json_data)){
      this_set=json_data[[i]]
      for (j in 1:length(this_set)){
          count=count+1
          reachid[count]=this_set[[j]]$reach_id  
          setid[count]=i
      }
  }
  print("set index")
  print(set_index)
  this_set = json_data[[set_index]]
  print("this set")
  print(this_set)
  a_reach = this_set[[1]]
  print("a reach")
  print(a_reach)
  sos_file = a_reach$sos
  print("Here is the new sos")
  print(sos_file)
  return(sets=data.frame('set_id'=setid,'reach_id'=reachid, 'sos_file'=rep(sos_file, length(reachid))))
        }

get_swot_from_set=function(setid,setdf){
    library(dplyr)
    library(tidyr)
    

    this_set=dplyr::filter(setdf,set_id==setid)
 
    thisset=data.frame('width'=NA,'slope'=NA,'time'=NA,
                       'reach_id'=NA, 'setid'=NA,'nt'=NA)
    
 
        
    for (reach_id in this_set$reach_id){
         swot_file=paste0('/mnt/data/input/swot/',reach_id,'_SWOT.nc') 


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
    print("here are the r ids")
    print(rids[4])
    print("there it was")
    print("here is the reach d")
    print(reach_id)
    print("yaaaa")
    index = which(rids == reach_id, arr.ind=TRUE)

    node_grp = grp.inq.nc(sos, "nodes")$self
    nrids = var.get.nc(node_grp, "reach_id")
    indexes = which(nrids == reach_id, arr.ind=TRUE)

    model_grp = grp.inq.nc(sos, "model")$self
    print("--------------log hat group------------------")
    print(index)
    print('thats the index----------')
    print("--------------end log hat--------------------")



    Q_priors$logQ_hat = log(var.get.nc(model_grp, "mean_q")[index])
    Q_priors$upperbound_logQ = log(var.get.nc(model_grp, "max_q")[index])
    min_q = var.get.nc(model_grp, "min_q")[index]   # Check action taken
    # print("here is whole group")
    # print(var.get.nc(model_grp, "min_q"))
    # print("here is minq")
    # print(min_q)
    # if ((min_q < 0) | (is.na(min_q))) {
    #   Q_priors$lowerbound_logQ = NA
    # } else {
    #   Q_priors$lowerbound_logQ = log(min_q)
    # }

    # Check for numeric(0) explicitly
    if (is.numeric(min_q) && length(min_q) == 0) {
      min_q = NA
    }

    if (is.na(min_q) | min_q < 0) {
      Q_priors$lowerbound_logQ = NA
    } else {
      Q_priors$lowerbound_logQ = log(min_q)
    }

    r_grp = grp.inq.nc(sos, "gbpriors/reach")$self
    Q_priors$logQ_sd = var.get.nc(r_grp, "logQ_sd")[index]

   
  close.nc(sos)



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
      #  dplyr::select(-reach_id)%>%
       summarize(logQ_hat=mean(logQ_hat),
                 logQ_sd=mean(logQ_sd),
                 lowerbound_logQ=mean(lowerbound_logQ),
                 upperbound_logQ=mean(upperbound_logQ))
       

       
       return(list('reach_ids'=dplyr::select(all_Q_df,reach_id),'Q_priors'=Q_priors) )
       
       
       }
    
    
