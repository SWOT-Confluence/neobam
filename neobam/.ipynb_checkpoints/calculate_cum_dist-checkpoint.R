calculate_cum_dist=function(sos,reach_id_in){
    library(geosphere)
  
    reach_index=sos$nodes$reach_id==reach_id_in
    #calculate a station vector
        x=sos$nodes$x[reach_index]
        y=sos$nodes$y[reach_index]
  
        A=c(x,0)
        B=c(0,x)
        C=c(y,0)
        D=c(0,y)
        #drop the first row, but keep the last column
        #this will make a matrix where the 1st - 2nd columns gives you 2-1, 3-2, etc
        #for consecutive distances
        distmatrix=as.matrix(cbind(A,B,C,D),ncol=4)[2:length(x),]

        distfunc=function(matrix_in){
            distm(c(matrix_in[['A']], matrix_in[['C']]), c(matrix_in[['B']], matrix_in[['D']]), fun = distHaversine)
        }

        #pad with a leading zero
        node_dist=c(0,apply(distmatrix,1,distfunc))

        #take a cumulative to make a station vector
        node_cum_dist =cumsum(node_dist)
    # print(node_cum_dist)

       return(node_cum_dist)
    
    }