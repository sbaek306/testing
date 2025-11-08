#library(rootSolve)
library(stats)

###############
## functions ##
###############

# approximate observerd lilekelihood

lik.i <- function(par,hI,o.d,v.d){
  
  beta <- par[1:2]
  a <- par[3]
  tau <- par[4]
  
  #print(a)
  
  res <- NULL
  
  tmp.I <- 1
  
  #while (tmp.I<I){
  for (ii in 1:hI){
    
    ind.o <- which(o.d$I==ii)
    ind.v <- which(v.d$I==ii)
    o.data <- o.d[ind.o,]
    v.data <- v.d[ind.v,]  
    
    t.data <- rbind(o.data,v.data)
    
    Po <- plogis(as.matrix(o.data[,c(1,2)]) %*% beta)
    Pv <- plogis(as.matrix(v.data[,c(1,2)]) %*% beta)
    
    mi <- length(Pv)
    #const1 <- mean(1 - exp(a)/(1+exp(a))*Pv)
    #const2 <- mean(Pv)
    
    #c1 <- sum(o.data$s * exp(a)/(1+exp(a))^2)
    #c2 <- sum( (1 - o.data$s) * exp(a) * Po * (1 - exp(2*a) + exp(2*a)*Po) 
    #           / ( (1 + exp(a))^2 * (1 + exp(a) - exp(a) * Po)^2 ) )
    #c3 <- sum(exp(a)/(1+exp(a))^2 * v.data$y) 
    
    #print(c(c1,c2,c3,sum(c1,c2,c3)))
    #print(c(sum(o.data$s),sum(1-o.data$s),sum(v.data$y)))
    
    #if (cond > 0){
    l.i0F.2 <- ( - sum(o.data$s * exp(a)/(1+exp(a))^2)
                 - sum( (1 - o.data$s) * exp(a) * Po * (1 - exp(2*a) + exp(2*a)*Po) 
                        / ( (1 + exp(a))^2 * (1 + exp(a) - exp(a) * Po)^2 ) )
                 - sum(exp(a)/(1+exp(a))^2 * v.data$y) 
                # + mi*const2*(exp(a)/(1+exp(a))^2/const1 - 2*exp(2*a)/(1+exp(a))^3/const1 
                #              + exp(2*a)/(1+exp(a))^3/const1) 
    )
    #} else {
    #  l.i0F.2 <- ( - sum(o.data$s * exp(a)/(1+exp(a))^2)
    #               - sum( (1 - o.data$s) * exp(a) * (Po+0.1) * (1 - exp(2*a) + exp(2*a)*(Po+0.1)) 
    #                      / ( (1 + exp(a))^2 * (1 + exp(a) - exp(a) * (Po+0.1))^2 ) )
    #               - sum(exp(a)/(1+exp(a))^2 * v.data$y) 
    #  )
    #}
    
    #if (l.i0F.2 < 0){
    
    l.i0F <- ( sum(o.data$s * ( log(exp(a)/(1+exp(a))) + log(Po) )
                   + (1 - o.data$s) * log(1-exp(a)/(1+exp(a))*Po) )
               + sum( v.data$y * ( -log(1+exp(a)) + log(Pv) )
                      + (1 - v.data$y) * log(1 - Pv) )
               + sum(dnorm(t.data[,2],1,1,log=TRUE)) #+ sum(dnorm(v.data[,2],1,sqrt(1),log=TRUE))
              #  - mi*log(const1)
    )
    
    # l.i0F...normal density parameter -> need to be careful
    
    l.i0F.1 <- ( sum(o.data$s /(1+exp(a)) - (1 - o.data$s) * exp(a) * Po /
                       ((1 + exp(a)) * (1 + exp(a) - exp(a) * Po))) 
                 - sum(v.data$y * exp(a)/(1+exp(a)))   
               #  + mi*exp(a)*const2/(1+exp(a))^2/const1
    )
    #res[tmp.I] <- 0.5 * log(1 - tau*l.i0F.2) - l.i0F - 0.5 * tau*(l.i0F.1^2) / (1-tau*l.i0F.2) # minus logL
    res[ii] <- 0.5 * log(1 - tau*l.i0F.2) - l.i0F - 0.5 * tau*(l.i0F.1^2) / (1-tau*l.i0F.2) # minus logL
    #  tmp.I <- tmp.I + 1
    
    # } 
    #print(o.d[1,])
  }
  #print(dat[[1]][1,])
  #ind <- which(is.na(res)==1)
  tmp <- res[complete.cases(res)]
  #h5 <- quantile(tmp,0.95);l5 <- quantile(tmp,0.05)
  #tmp2 <- tmp[which(tmp <= h5 & tmp >= l5)]
  #print(tmp)
  return(sum(tmp))
}

lik.i.H0 <- function(par,hI,o.d,v.d){
  
  beta <- par[1:2]
  a <- par[3]
  tau <- 0
  
  res <- NULL
  
  tmp.I <- 1
  
  #while (tmp.I<I){
  for (ii in 1:hI){
    
    ind.o <- which(o.d$I==ii)
    ind.v <- which(v.d$I==ii)
    o.data <- o.d[ind.o,]
    v.data <- v.d[ind.v,]  
    
    t.data <- rbind(o.data,v.data)
    
    Po <- plogis(as.matrix(o.data[,c(1,2)]) %*% beta)
    Pv <- plogis(as.matrix(v.data[,c(1,2)]) %*% beta)
    
    #mi <- length(Pv)
    
    #l.i0F.2 <- ( - sum(o.data$s * exp(a)/(1+exp(a))^2)
    #             - sum( (1 - o.data$s) * exp(a) * Po * (1 - exp(2*a) + exp(2*a)*Po) 
    #                    / ( (1 + exp(a))^2 * (1 + exp(a) - exp(a) * Po)^2 ) )
    #             - sum(exp(a)/(1+exp(a))^2 * v.data$y) 
    #)
    
    l.i0F <- ( sum(o.data$s * ( log(exp(a)/(1+exp(a))) + log(Po) )
                   + (1 - o.data$s) * log(1-exp(a)/(1+exp(a))*Po) )
               + sum( v.data$y * ( -log(1+exp(a)) + log(Pv) )
                      + (1 - v.data$y) * log(1 - Pv) )
               + sum(dnorm(t.data[,2],1,1,log=TRUE)) #+ sum(dnorm(v.data[,2],1,sqrt(1),log=TRUE))
              #  - mi*log(const1)
    )
        
    #l.i0F.1 <- ( sum(o.data$s /(1+exp(a)) - (1 - o.data$s) * exp(a) * Po /
    #                   ((1 + exp(a)) * (1 + exp(a) - exp(a) * Po))) 
    #             - sum(v.data$y * exp(a)/(1+exp(a)))   
               #  + mi*exp(a)*const2/(1+exp(a))^2/const1
    #)
    
    #res[ii] <- 0.5 * log(1 - tau*l.i0F.2) - l.i0F - 0.5 * tau*(l.i0F.1^2) / (1-tau*l.i0F.2) # minus logL
    res[ii] <- - l.i0F 
    #  tmp.I <- tmp.I + 1
    
    # } 
    #print(o.d[1,])
  }
  #print(dat[[1]][1,])
  #ind <- which(is.na(res)==1)
  tmp <- res[complete.cases(res)]
  #h5 <- quantile(tmp,0.95);l5 <- quantile(tmp,0.05)
  #tmp2 <- tmp[which(tmp <= h5 & tmp >= l5)]
  #print(tmp)
  return(sum(tmp))
}

lrt.stat <- function(par,hI,o.d,v.d){
  
  beta <- par[1:2]
  a <- par[3]
  tau <- par[4]
  
  res <- NULL
  
  tmp.I <- 1
  
  #while (tmp.I<I){
  for (ii in 1:hI){
    
    ind.o <- which(o.d$I==ii)
    ind.v <- which(v.d$I==ii)
    o.data <- o.d[ind.o,]
    v.data <- v.d[ind.v,]  
    
    t.data <- rbind(o.data,v.data)
    
    Po <- plogis(as.matrix(o.data[,c(1,2)]) %*% beta)
    Pv <- plogis(as.matrix(v.data[,c(1,2)]) %*% beta)
    
    mi <- length(Pv)
    
    l.i0F.2 <- ( - sum(o.data$s * exp(a)/(1+exp(a))^2)
                 - sum( (1 - o.data$s) * exp(a) * Po * (1 - exp(2*a) + exp(2*a)*Po) 
                        / ( (1 + exp(a))^2 * (1 + exp(a) - exp(a) * Po)^2 ) )
                 - sum(exp(a)/(1+exp(a))^2 * v.data$y) 
    )
    
    l.i0F <- ( sum(o.data$s * ( log(exp(a)/(1+exp(a))) + log(Po) )
                   + (1 - o.data$s) * log(1-exp(a)/(1+exp(a))*Po) )
               + sum( v.data$y * ( -log(1+exp(a)) + log(Pv) )
                      + (1 - v.data$y) * log(1 - Pv) )
               + sum(dnorm(t.data[,2],1,1,log=TRUE)) #+ sum(dnorm(v.data[,2],1,sqrt(1),log=TRUE))
              #  - mi*log(const1)
    )
    
    
    l.i0F.1 <- ( sum(o.data$s /(1+exp(a)) - (1 - o.data$s) * exp(a) * Po /
                       ((1 + exp(a)) * (1 + exp(a) - exp(a) * Po))) 
                 - sum(v.data$y * exp(a)/(1+exp(a)))   
    )

    res[ii] <- 0.5 * log(1 - tau*l.i0F.2) - l.i0F - 0.5 * tau*(l.i0F.1^2) / (1-tau*l.i0F.2) # minus logL

  }
  tmp <- res[complete.cases(res)]
  #ind <- which(tmp<quantile(tmp,0.04))
  #tmp <- tmp[-ind]
  return(sum(-tmp))
}

score.stat <- function(par,hI,o.d,v.d){
  
  beta <- par[1:2]
  a <- par[3]
  tau <- par[4]
  
  res = fisher <- NULL
  
  tmp.I <- 1
  
  #while (tmp.I<I){
  for (ii in 1:hI){
    
    ind.o <- which(o.d$I==ii)
    ind.v <- which(v.d$I==ii)
    o.data <- o.d[ind.o,]
    v.data <- v.d[ind.v,]  
    
    t.data <- rbind(o.data,v.data)
    
    Po <- plogis(as.matrix(o.data[,c(1,2)]) %*% beta)
    Pv <- plogis(as.matrix(v.data[,c(1,2)]) %*% beta)
    
    mi <- length(Pv)
    
    l.i0F.2 <- ( - sum(o.data$s * exp(a)/(1+exp(a))^2)
                 - sum( (1 - o.data$s) * exp(a) * Po * (1 - exp(2*a) + exp(2*a)*Po) 
                        / ( (1 + exp(a))^2 * (1 + exp(a) - exp(a) * Po)^2 ) )
                 - sum(exp(a)/(1+exp(a))^2 * v.data$y) 
    )
    
    l.i0F <- ( sum(o.data$s * ( log(exp(a)/(1+exp(a))) + log(Po) )
                   + (1 - o.data$s) * log(1-exp(a)/(1+exp(a))*Po) )
               + sum( v.data$y * ( -log(1+exp(a)) + log(Pv) )
                      + (1 - v.data$y) * log(1 - Pv) )
               + sum(dnorm(t.data[,2],1,1,log=TRUE)) #+ sum(dnorm(v.data[,2],1,sqrt(1),log=TRUE))
              #  - mi*log(const1)
    )
        
    l.i0F.1 <- ( sum(o.data$s /(1+exp(a)) - (1 - o.data$s) * exp(a) * Po /
                       ((1 + exp(a)) * (1 + exp(a) - exp(a) * Po))) 
                 - sum(v.data$y * exp(a)/(1+exp(a)))   
    )

    fisher[ii] <- l.i0F.2^2 + 2*l.i0F.1^2*l.i0F.2
    ## adjustments 
    #if (abs(l.i0F.1)>4) { 
    #  l.i0F.1 <- sign(l.i0F.1)*4
    #} 
    ## 
    res[ii] <- 0.5 * (l.i0F.2 + l.i0F.1^2)
    #print(c(l.i0F.1,l.i0F.2))

  }
  #ind <- which(res>quantile(res,0.90) | res<quantile(res,0.1))
  ind <- which(res>quantile(res,0.95))
  res <- res[-ind]
  fisher <- fisher[-ind]
  #asy.var <- -2*hI/sum(fisher)
  asy.var <- -2*(hI-length(ind))/sum(fisher)
  tmp <- sum(res) / sqrt(asy.var)
  #tmp <- sum(res) * (-2*hI/sum(fisher))^(-0.5)
  #print(c(sum(res),sqrt(asy.var)))
  #return(tmp^2) #two-sides
  #print(sum(res))
  return(tmp)
}

#################
## simulations ##
#################

simu <- function(simuno,seedn,toI,n,true.beta,true.a,tt,rr,r.sel){
  
  #simuno <- 1
  est <- matrix(0,simuno,4)
  est.0 <- matrix(0,simuno,3)
  score <- rep(NA,simuno)
  lrt <- matrix(0,simuno,2)
  err <- NULL; LRT <- NULL
  #pv <- rep(0,simuno)
  # true beta
  #beta <- c(1,-5) # prevalence of Y: 10%
  #I <- 1000
  #n <- 100
  #a <- 3
  s.tt <- sqrt(tt)
  beta <- true.beta 
  a <- true.a 
  
  set.seed(seedn)
  
  tmp.iter <- 1
  iter <- 1
  
  while (simuno >= tmp.iter){
  #for (jj in 1:simuno){
    
    br <- rnorm(toI,0,s.tt) # random bi's
    
    v.d = o.d <- NULL
    
    for (ii in 1:toI){
      # generate X
      x1 <- rnorm(n,1,1)
      x <- as.matrix(cbind(rep(1,n),x1))
      
      # generate Y
      p <- plogis(x %*% beta)
      y <- rbinom(n,1,p)
      
      # generate S
      b <- br[ii]
      s <- rep(0,n)
      sp <- exp(a+b)/(1+exp(a+b))
      for (i in 1:n) {
        if (y[i]==1){
          s[i] <- rbinom(1,1,sp)} else {s[i]==0}
      }
      
      # data
      full.data <- data.frame(x,s,y)
      #r <- rr # the proportion of validation subset
      #tmp <- rbinom(sum(full.data$s==0),1,r) # available Y 
      #R <- rep(0,n)
      #R[which(full.data$s==0)] <- tmp
      
      #v.ind <- which(R==1)
      tmp.n <- ceiling(r.sel/toI)
      v.ind <- sample(which(full.data$s==0),tmp.n)
      v.data <- data.frame(full.data[v.ind,]) # validation data whose Y will be used.
      o.data <- data.frame(full.data[-v.ind,]) # y will not be used.
      
      v.data$I <- ii # add indices for i
      o.data$I <- ii # add indices for i
      
      v.d <- rbind(v.d,v.data)
      o.d <- rbind(o.d,o.data)
    }
    
    #pert <- abs(c(true.beta,true.a,tt) * 0.1)
    #par0 <- c(true.beta,true.a,tt) + c(runif(3,-pert,pert),runif(1,0,0.001)) ################ for simulation
    par0 <- jitter(c(true.beta,true.a,tt))#;par0[3] <- runif(1,-0.01,0.01)
    A <- matrix(c(0,0,0,1),1,4,byrow=TRUE)
    optim.est <- try(constrOptim(par0,lik.i,grad=NULL,ui=A,ci=c(0),method="Nelder-Mead",control=list(maxit=200),hI=toI,o.d=o.d,v.d=v.d),silent=TRUE)
    
    par0.H0 <- jitter(c(true.beta,true.a))
    optim.est.H0 <- try(optim(par0.H0,lik.i.H0,gr=NULL,method="Nelder-Mead",control=list(maxit=200),hI=toI,o.d=o.d,v.d=v.d),silent=TRUE)

    if (is.character(optim.est)){ 
      est[tmp.iter,] <- NA 
      #est[tmp.iter,5] <- optim.est$convergence
    } else if (optim.est$convergence==0 & optim.est.H0$convergence==0 ){
      #est[tmp.iter,5] <- optim.est$convergence
      est[tmp.iter,] <- optim.est$par
      #st.0[tmp.iter,4] <- optim.est.H0$convergence
      est.0[tmp.iter,] <- optim.est.H0$par

      lrt[tmp.iter,1] <- lrt.stat(optim.est$par,hI=toI,o.d=o.d,v.d=v.d)
      lrt[tmp.iter,2] <- lrt.stat(c(optim.est.H0$par,0),hI=toI,o.d=o.d,v.d=v.d)
      score[tmp.iter] <- score.stat(c(optim.est.H0$par,0),hI=toI,o.d=o.d,v.d=v.d)
      
      cat('tmp.iter = ',tmp.iter,'\n')
      
      tmp.iter <- tmp.iter + 1
    
    }
    #cat('iter = ',tmp.iter,'\n',est[tmp.iter,],'\n')
    cat('iter = ',iter,'\n')
    iter <- iter + 1
    #print(err[jj])
  }
  
  # Likelihood ration test 
  #tt=lrt 
  #ts <- -2*tt[,2] + 2*tt[,1] ## Tl
  #pv <- (1-pchisq(ts,1))/2 ## p-value 
  #LRT <- mean(pv<0.05) # the mean number of rejections
  
  return(list(est,est.0,lrt,score))
}


est <- simu(simuno=1000,seedn=432,toI=50,n=100,true.beta=c(1,-5),true.a=0.1,tt=0.1,rr=0.5,r.sel=200)

save(est, file = "p2-1d.RData")


