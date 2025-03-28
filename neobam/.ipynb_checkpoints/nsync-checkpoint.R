nsync=function(swot_df,xs_id){
     library(zoo)
    
classify_hydrograph=function(hydrograph){
             #first forward difference
        FFD_hydrograph_this= (c(hydrograph,0)-c(0,hydrograph))[2:length(hydrograph)]
        #get the next point in the series
        FFD_hydrograph_next= FFD_hydrograph_this[2:length(FFD_hydrograph_this)]
        #trim the first series to get a same-lenght vector
        FFD_hydrograph_this=FFD_hydrograph_this[1:(length(FFD_hydrograph_this)-1)]

        this_key=FFD_hydrograph_this
        tk1=which(this_key<0)
        tk2=which(this_key==0)
        tk3=which(this_key>0)

        this_key[tk1]=1 #neg
        this_key[tk2]=2 #zero
        this_key[tk3]=3 #pos

        next_key=FFD_hydrograph_next
        nk1=which(next_key<0)
        nk2=which(next_key==0)
        nk3=which(next_key>0)

        next_key[nk1]=10 #neg
        next_key[nk2]=100 #zero
        next_key[nk3]=1000 #pos

        classification=this_key+next_key
        classification[classification==11]='falling'
        classification[classification==101]='falling'
        classification[classification==1001]='trough'
        classification[classification==12]='falling'
        classification[classification==102]='flat'
        classification[classification==1002]='rising'
        classification[classification==13]='peak'
        classification[classification==103]='rising'
        classification[classification==1003]='rising'
        classification=c(NA,classification,NA)

#             for (i in 2:(length(FFD_hydrograph)-1)){
#                 this_point=FFD_hydrograph[i]
#                 prev_point=FFD_hydrograph[i-1]
#                 next_point=FFD_hydrograph[i+1]   

#                 #pos to negative = peak
#                 #pos to pos = rising
#                 #neg to pos = trough
#                 #neg to neg = falling

#                 if(is.na(next_point) | is.na(this_point)){class=NA
#                                                          next()}

#                 if(next_point <0 & this_point<0){class[i]='falling'}
#                 if(next_point >0 & this_point>0){class[i]='rising'}
#                 if(next_point <0 & this_point>0){class[i]='peak'}
#                 if(next_point >0 & this_point<0){class[i]='trough'}

#                 if(next_point ==0 & this_point >0){class[i]='rising'}
#                 if(next_point ==0 & this_point <0){class[i]='falling'}
#                 if(next_point <0 & this_point ==0){class[i]='falling'}
#                 if(next_point >0 & this_point ==0){class[i]='rising'}
   return(classification)
        }
         

    temp_df=filter(swot_df,xs==xs_id)
    

    smoothed_reachwse=rollapply(temp_df$reach_wse, width=10, mean,na.rm=TRUE)
    smoothed_nodewse=rollapply(temp_df$node_wse, width=10, mean,na.rm=TRUE)
  
    class_reach=classify_hydrograph(smoothed_reachwse)
    class_node=classify_hydrograph(smoothed_nodewse)
    
  
    
    # print(smoothed_reachwse)
    # print(smoothed_nodewse)
    # print(class_reach)
    # print(class_node)

    suppressWarnings({
    agreement=sum(class_reach ==class_node,na.rm=TRUE)/length(class_reach)
        })

    output_df=temp_df%>%
                mutate(agreement=agreement)

     
     }