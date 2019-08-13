### 29/10/15: SIStrain

siStrainNT<- function(beta=2, gamma=1, theta=0.01, N=25, t.max=100, M=100, t.step=5, t.burn=0, t.delay=1, lambda=0.5, 
                      csv=F,v=NULL, w=NULL){
##### siStrainNT simulates QSD for SIS epidemic with evolving strains by obtaining multiple samples after burn-in
    ## OUT
    # list consisting of
    # s = samples           samples of QSD particles
    # w = samples.w         weights of particles in s
	# iw = iw				weighted average of I, S, K
    # rho = rho             average proportion of individual-based immunity
##### initialisation ##### 
    if(is.null(v)){
    v <- array(0, dim=c(2,N,M)) # 1 initial infective of strain 0
        v[1,1,] <- 1 # 1 infective strain 0
        v[2,2,] <- N-1 # N-1 susceptibles strain -1
    }
    if(is.null(w)){
        w <- rep(1/M,M) #particle weights
    }
    #tock <- 0
    
    s <- 1 # event count
    t <- 0 # time
    t.res <- 0 #time since last resampling step
    frames <- floor((t.max-t.burn)/t.delay)+1 #number of samples to take
    f <- 0 # current frame number
    
    samples <- array(NA,dim=c(2,N,M*frames))
    samples.w <- rep(NA,M*frames)

    while(t < t.max && sum(v[1,,])>0 ){ # while not hit end and infectives remain
        #if(t > tock+1){print(tock <- tock+1)}
##### resampling step #####
        if(sum(rep(1,N)%*%(v[1,,]!=0)!=0) < lambda*M || t.res > t.step){ 
            t.res <- 0
            for(m in which(w!=0)){
                if(is.na(v[1,1,m])){
                    next
                }
                tr <- which( rep(1,dim(v)[2])%*%(abs(v[1,,]-rep(v[1,,m], times=M))+abs(v[2,,]-rep(v[2,,m], times=M)))==0)
                if(length(tr)>1){
                    tr <- tr[-1]
                    w[m] <- w[m] + sum(w[tr])
                    w[tr] <- 0
                    v[,,tr] <- NA
                }
            }
            new.p <- sample.int(M, sum(w==0), replace=T, prob=w)
            new.n <- c(which(w!=0), new.p)
            v <- v[,,new.n]
            w <- w[new.n]
            for(m in 1:M){
                w[m] <- w[m]/sum(new.n==new.n[m])
            }
        }
##### simulation step #####
        rates.S <- (rep(beta/N,N)%*%v[1,,])*(rep(1,N)%*%v[2,,]) #infection rates
        rates.I <- rep(gamma,N)%*%v[1,,] # recovery rates
        rate <- sum(rates.S,rates.I)
        rate.M <- rates.S + rates.I # particle rates
        
        t.next <- rexp(1,rate)
        if(t + t.next >= t.max){  #if too long to next event
            s <- s+1
            t <- t.max
            next
        }else{
            m <- max(which(cumsum(c(0, rate.M)/sum(rate.M)) < runif(1))) #which particle event occurs on
            if(runif(1) < sum(rates.I[m])/rate.M[m]){
                ### Recovery Event Ik -> Sk
                j <- sample.int(N,1,prob=v[1,,m]) #pick infected strain
                v[,j,m] <- v[,j,m] + c(-1,1)
                #i.count[m] <- i.count[m] - 1
                if(v[1,j,m] == 0){ #if strain becomes empty
                    #k.active[m] <- k.active[m] - 1
                    if(j < N){
                        v[1,j:N,m] <- c(v[1,(j+1):N,m],0)
                        if(j < N-1){
                            v[2,j:N,m] <- c(v[2,j,m]+v[2,j+1,m], v[2,(j+2):N,m], 0)
                        }else{
                            v[2,(N-1):N,m] <- c(v[2,j,m]+v[2,j+1,m], 0)
                        }
                    }
                }
                t <- t + t.next # next event time
                t.res <- t.res + t.next # add to resample clock
                s <- s + 1
            }else{
                js <- sample.int(N,1,prob=v[2,,m]) #pick susceptible strain
                if(runif(1) < theta){ #if mutation event
                    #k.count[m] <- k.count[m] + 1
                    #k.active[m] <- k.active[m] + 1
                    #i.count[m] <- i.count[m] + 1
                    v[2,js,m] <- v[2,js,m] - 1
                    v[,,m] <- c(1,0,v[,-N,m])
                    t <- t + t.next # next event time
                    t.res <- t.res + t.next # add to resample clock
                    s <- s + 1
                }else{
                    ji <- sample.int(N,1,prob=v[1,,m]) # pick infective to infect
                    if(ji < js){
                        v[1,ji,m] <- v[1,ji,m] + 1
                        v[2,js,m] <- v[2,js,m] - 1
                        t <- t + t.next # next event time
                        t.res <- t.res + t.next # add to resample clock
                        s <- s + 1
                        #i.count[m] <- i.count[m] + 1
                    }else{
                        next
                    }
                }
            }
        }
        
############# reweighting #####
        w[m] <- w[m]*(sum(v[1,,m])!=0)
        w <- w/sum(w)
        if(is.nan(w[1])){
            print(w[1])
        }
        
############# taking sample step #####
        if(t > t.burn + f*t.delay){ 
            f <- f+1
            samples.w[((f-1)*M + 1):(f*M)] <- w
            samples[,,((f-1)*M + 1):(f*M)] <- v
        }
    }
    #### Final Sampling Step
    f <- f+1
    samples.w[((f-1)*M + 1):(f*M)] <- w
    samples[,,((f-1)*M + 1):(f*M)] <- v
    
    samples.w[samples.w==0] <- NA
    samples <- samples[,,!is.na(samples.w)]
    samples.w <- samples.w[!is.na(samples.w)]
    samples.w <- samples.w/sum(samples.w)
    
    v.summary <- matrix(0,ncol=N,nrow=2)
    v.summary[1,] <- samples[1,,]%*%samples.w
    v.summary[2,] <- samples[2,,]%*%samples.w
    
    I <- samples[1,,]
    S <- samples[2,,]
    
    A <- matrix(0,N,N)
    for(i in 1:N){
        for(j in 1:N){
            A[i,j] <- 1*(j >= i)
        }
    }
    rho <- (1/N) * { (rep(1,N) %*% I) + 
                     ((rep(1,N) %*% (S * (A%*%I)))/(rep(1,N)%*%I))  }%*%samples.w
    #B <- sapply(k.count, function(x){y <- rep(0,N); y[1:x] <- x:1;y})
    iw <- rep(NA,3)
    iw[1] <- (rep(1,N)%*%samples[1,,])%*%samples.w #I
    iw[2] <- (rep(1,N)%*%(samples[1,,]>0))%*%samples.w #K
    iw[3] <- rho
    #iw[4] <- rep(1,N)%*%((samples[1,,]+samples[2,,])*B)%*%t(samples.w/(N*k.count))
    
    
    ## see uI(k) = \sum_m=1^M w^{(m)} I_k^(m)/|I|^(m)
    Irep <- matrix(rep(rep(1,N)%*%samples[1,,], each=N), nrow=N)
    uI <- (samples[1,,]/Irep)%*%samples.w
    
    ### uL(j) = \sum_m=1^M w^(m) \sum_{j<k} (I_k^(m) + S_k^(m) + R_k^(m))/N
    Kback <- samples[1,(N:1),] + samples[2,(N:1),]
    revcumsum <- function(x){c(rev(cumsum(x[1:(length(x)-1)])), 0)}
    uL <- (1/N) * (apply(Kback, 2, revcumsum) %*% samples.w)
    
    ### R_Q
    rq <- beta*(theta + (1-theta)*(t(uL) %*% uI))/gamma
    
    iw[4] <- rq
    #write samples table
    if(csv){
    write.csv(cbind(matrix(aperm(samples,c(3,1,2)),nrow=length(samples.w), ncol=2*N),matrix(samples.w, ncol=1)), 
              file=paste0("./samples_b",beta,"_t",theta,"_n",N,".csv"), row.names=F)
        #### How to read samples into V and W again.
        # D <- read.csv("samples_b2_t0.05_n25.csv")
        # WW <- D[,dim(D)[2]]
        # VV <- array(t(D[,1:50]),dim=c(2,25,dim(D)[1]))
        ####
    }

    return(list(v=v.summary[,N:1], iw=iw, samples=samples, sw = samples.w, rq=rq))
}