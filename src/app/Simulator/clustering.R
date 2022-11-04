clustering <- function(serie_cluster,nb_cluster,mon_iter.max,mon_nstart,ma_date_sat,elevation_point_sat,tol){

      a_prendre <- list()
      maxi <- 1
      ecarts <- 0.01
      while((sum((maxi/ecarts)>=tol)>=1) & (nb_cluster>=1)){

            a <- kmeans(serie_cluster,nb_cluster,iter.max=mon_iter.max,nstart=mon_nstart)
            a_prendre[[1]] <- a$cluster==1
            a_prendre[[2]] <- a$cluster==2
            a_prendre[[3]] <- a$cluster==3
            a_prendre[[4]] <- a$cluster==4
            a_prendre[[5]] <- a$cluster==5
            a_prendre[[6]] <- a$cluster==6
            a_prendre[[7]] <- a$cluster==7
            clusters <- a$centers
            
            ecarts <- (clusters)
            ecarts <- ecarts[order(ecarts)]
            ecarts <- abs(diff(ecarts))
            maxi <- max(ecarts)            
            nb_cluster <- nb_cluster - 1
      }  
      nb_cluster <- nb_cluster + 1
      a_prendre <- a_prendre[1:nb_cluster]

      return(a_prendre)
}