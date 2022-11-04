##################################################################
#
# Programme convertit les jours juliens en jour GPS
#
# Input : 
#   - annee, mois, jour en julien
# Output :                
#   - date_gps : la date GPS correspondante
#
#   juillet 2014 - Nicolas ROUSSEL
#
##################################################################


jul_2_gps <- function(annee,mois,jour){   

        annee <- as.numeric(annee)
        mois <- as.numeric(mois)
        jour <- as.numeric(jour)
      
        mamat <- NaN
      
        if(mois>12){
            annee <- annee + 1
            mois <- 1
            jour <- 1
        }
      
        if(jour==0){
           nb_jour <- julian(as.Date(paste(annee,'-',mois,'-',jour+1,sep='')))-1
        }
      
        else{try(mamat <- as.Date(paste(annee,'-',mois,'-',jour,sep='')),silent=TRUE)
          	  if(is.na(mamat)){
          		  mois <- mois + 1
          		  jour <- 1  
          	  } 
          
          	  if(mois>12){
          				annee <- annee + 1
            			mois <- 1
          				jour <- 1
          		  }
          
          	  nb_jour <- julian(as.Date(paste(annee,'-',mois,'-',jour,sep=''))) 
        }

        nb_jour_gps <- nb_jour - julian(as.Date('1980-01-06'))
        nb_weeks <- trunc(nb_jour_gps/7)
        nb_jours_gps <- nb_jour_gps-nb_weeks*7
        date_gps <- paste(nb_weeks,nb_jours_gps,sep='_')
        return(date_gps)
}

