rm(list=ls())
graphics.off()
Sys.setlocale("LC_TIME", "English")

######################### To be defined by the user before any use 
nom_date_1 <- "date_h.txt"
nom_serie_1 <- "h.txt"
nom_date_2 <- "t_ADCP.txt"
nom_serie_2 <- "h_ADCP.txt"
chiffres_significatifs <- 2

analyse_dephasage <- TRUE
increment_dephasage <- 60
dephasage_max <- 3*3600
######################### End






reechantillonne_series <- function(t1,t2,s1,s2){
      T1 <- mean(diff(t1))
      T2 <- mean(diff(t2))
        
      T_final <- min(T1,T2)
      
      min_date <- max(min(t1),min(t2))
      max_date <- min(max(t1),max(t2))
      if(min_date>max_date){
            print('There is an error with the dates. Please check.'); PLOUF
      }
      t_finale <- seq(min_date,max_date,T_final)
      
      s1 <- approx(t1,s1,t_finale)$y
      s2 <- approx(t2,s2,t_finale)$y

      sortie <- list()
      sortie[[1]] <- t_finale
      sortie[[2]] <- s1
      sortie[[3]] <- s2
      return(sortie)      
}

tracer <- function(x,y1,y2){
      #par(mfrow=c(2,1))
      x11()
      plot(as.POSIXct(x,origin='1970-01-01',tz="UTC"),y1,type='l',ylim=c(min(serie_1_reech,serie_1_reech),max(serie_1_reech,serie_2_reech)),xlab='Date',ylab='Series',font.axis=2,cex.axis=2,cex.lab=1.5,lwd=2)
      lines(as.POSIXct(x,origin='1970-01-01',tz="UTC"),y2,col='red',lwd=2)
      #legend("topright",inset=0.05,legend=c('Series 1','Series 2'),col=c('black','red'),lty=c(1,1))
      model <- lm(y2 ~ y1)
      trend <- model$fitted.values
      coeff <- model$coefficients
      R2 <- summary(model)$r.squared  
      x11() 
      plot(y1,y2,main='Scatter plot',xlab=substr(nom_serie_1,1,nchar(nom_serie_1)-4),ylab=substr(nom_serie_2,1,nchar(nom_serie_2)-4),font.axis=2,cex.axis=2,cex.lab=1.5,pch=16,cex=2)
      mtext(paste('R2 =',round(R2,chiffres_significatifs),"; n =",length(y1)),side=1,line=-2,cex=2,col='blue')
      lines(y1,y1*coeff[2]+coeff[1],col='blue',lwd=3)   
}



date_1 <- read.table(nom_date_1)[[1]]
date_2 <- read.table(nom_date_2)[[1]]
serie_1 <- read.table(nom_serie_1)[[1]]
serie_2 <- read.table(nom_serie_2)[[1]]

if(length(date_1) != length(serie_1)){
      print(paste("Lengths of",nom_date_1,"and",nom_serie_1,"are different and should be the same. Please check")); PLOUF
}

if(length(date_2) != length(serie_2)){
      print(paste("Lengths of",nom_date_2,"and",nom_serie_2,"are different and should be the same. Please check")); PLOUF
}

result <- reechantillonne_series(date_1,date_2,serie_1,serie_2) 

date_finale <- result[[1]]
serie_1_reech <- result[[2]]
serie_2_reech <- result[[3]]

biais <- mean(serie_1_reech) - mean(serie_2_reech)
serie_2_reech <- serie_2_reech + biais
correlation <- cor(serie_1_reech,serie_2_reech)

print('****************************************',quote=FALSE)
print(paste('Comparison between', nom_serie_1,'and',nom_serie_2),quote=FALSE)
print('****************************************',quote=FALSE)
print('',quote=FALSE)

if(!analyse_dephasage){
      print(paste('Correlation =',round(correlation,chiffres_significatifs)),quote=FALSE) 
      print(paste('Bias =',round(biais,chiffres_significatifs),'m'),quote=FALSE)
      tracer(date_finale,serie_1_reech,serie_2_reech)
}


if(analyse_dephasage){

      seq_dephasage <- seq(-dephasage_max,dephasage_max,increment_dephasage)
      biais <- 0
      correlation <- 0
      dephasage <- 0
      incr <- 1
      for(incr_dephasage in seq_dephasage){
            date_2_deph <- date_2 + incr_dephasage
            if(min(max(date_2_deph),max(date_1))>max(min(date_2_deph),min(date_1))){
                  result <- reechantillonne_series(date_2_deph,date_1,serie_2,serie_1)  
                  date_finale <- result[[1]]
                  serie_2_reech <- result[[2]]
                  serie_1_reech <- result[[3]]
                  biais[incr] <- mean(serie_2_reech) - mean(serie_1_reech)
                  serie_2_reech <- serie_2_reech - biais[incr]
                  correlation[incr] <- cor(serie_2_reech,serie_1_reech) 
            }
            else{
                  biais[incr] <- 0
                  correlation[incr] <- 0
            } 
            incr <- incr + 1
      }
      
      ma_seq <- seq(1,length(seq_dephasage),1)
      index <- ma_seq[abs(correlation)==max(abs(correlation))][1] 
      dephasage <- seq_dephasage[index]
      biais <- biais[index]
      x11()
      plot(seq_dephasage/60,correlation,xlab='Shift (min)',ylab='Correlation',type='l',lwd=4)
      correlation <- correlation[index]
      
      serie_2 <- serie_2-biais
      date_2 <- date_2 + dephasage
      result <- reechantillonne_series(date_2,date_1,serie_2,serie_1)  
      date_finale <- result[[1]]
      serie_2_reech <- result[[2]]
      serie_1_reech <- result[[3]]
                  
      tracer(date_finale,serie_1_reech,serie_2_reech)
      print(paste('Correlation =',round(correlation,chiffres_significatifs)),quote=FALSE)
      print(paste('Bias =',round(biais,chiffres_significatifs),'m'),quote=FALSE)
      print(paste('Shift =',round(dephasage/60,chiffres_significatifs),'min'),quote=FALSE)
      if(abs(dephasage) == dephasage_max)
      {
            print('!!!! WARNING !!!!! The shift seems to be wrong, maybe check the border....', quote=FALSE,col='red')
      }
                         
}

      


