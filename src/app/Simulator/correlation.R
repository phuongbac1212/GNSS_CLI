calcul_correlation <- function(x1_ori,y1_ori,x2_ori,y2_ori,corriger_shift){
      seq_shift <- seq(-700,700,10)*60
      correlation_2 <- 0
      RMSE_2 <- 0
      R2_2 <- 0
      bias_2 <- 0
      shift <- 0
      y2 <- approx(x2_ori,as.numeric(y2_ori),xout=x1_ori)$y;
      y2 <- (y2-min(y2))*(max(y1_ori)-min(y1_ori))/(max(y2)-min(y2))+min(y1_ori)
      x1 <- x1_ori
      y1 <- y1_ori
      bias <- mean(y2)-mean(y1)
      model <- lm(y2 ~ y1)
      trend <- model$fitted.values
      correlation <- cor(y1,y2)
  
      RMSE <- sqrt(sum((y2-y1)^2)/length(y2))
      R2 <- summary(model)$r.squared
      ## en enlevant le shift
      if(corriger_shift){
            incr <- 1
            correlation_avec_shift <- 0
            for(u in seq_shift){
                  x1 <- x1_ori+u
                  y2 <- approx(x2_ori,as.numeric(y2_ori),xout=x1)$y;
                  y2 <- (y2-min(y2))*(max(y1)-min(y1))/(max(y2)-min(y2))+min(y1)
                  correlation_avec_shift[incr] <- cor(y1,y2)
                  incr <- incr + 1
            }
            shift <- seq_shift[abs(correlation_avec_shift)==max(abs(correlation_avec_shift))]
            x1 <- x1_ori + shift
            y2 <- approx(x2_ori,as.numeric(y2_ori),xout=x1)$y;
            y2 <- (y2-min(y2))*(max(y1)-min(y1))/(max(y2)-min(y2))+min(y1)

            correlation_2 <- cor(y1,y2)
          # plot(seq_shift,correlation_avec_shift)
#           x11()
            bias_2 <- mean(y2)-mean(y1)
            #y2 <- y2-bias
            model <- lm(y2 ~ y1)
            trend <- model$fitted.values
            RMSE_2 <- sqrt(sum((y2-y1)^2)/length(y2))
            R2_2 <- summary(model)$r.squared
      }
      return(c(correlation,RMSE,R2,bias,correlation_2,RMSE_2,R2_2,bias_2,shift))
}