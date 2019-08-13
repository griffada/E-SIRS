### 29/10/15: SIRStrain

sirStrainNT <- function(beta=2, gamma=1, delta=0.5, theta=0.01, N=25, t.max=100, M=100, t.step=5, t.burn=0, 
                        t.delay=1, lambda=0.5, csv=F){
##### SirStrainNT simulates QSD for SIRS epidemic with evolving strains by obtaining multiple samples after burn-in. #####
    ## IN
    # beta, gamma, delta    infection, recovery, loss of immunity rates > 0
    # theta                 mutation rate >= 0
    # N                     population size >= 1
    # t.max                 maximum time to run for > 0
    # t.step                maximum time between resampling steps
    # t.burn            `   time before samples are taken < t.max
    # t.delay               time between samples
    # lambda                proportion of particles dead to trigger resampling 
    # M                     number of particles to run
    # plot                  if TRUE, plot trace of epidemic
    ## OUT
    # list consisting of
    # s = samples           samples of QSD particles
    # w = samples.w         weights of particles in s
    # v = summary           
    
######### initialisation ##### 
    v <- array(0, dim=c(3,N,M)) # 1 initial infective of strain 0 rows (I,R,S)
    v[1,1,] <- 1
    v[3,2,] <- N-1
    w <- rep(1/M,M) #particle weights
    
    s <- 1 # event count
    t <- 0 # time
    t.res <- 0 #time since last resampling step
    frames <- floor((t.max-t.burn)/t.delay)+1 #number of samples to take
    f <- 0 # current frame number
    
    samples <- array(NA,dim=c(3,N,M*frames))
    samples.w <- rep(NA,M*frames)
    
    while(t < t.max && sum(v[1,,])>0 ){ # while not hit end and infectives remain
        
############# resampling step #####
        if(sum((rep(1,N)%*%v[1,,])!=0) < lambda*M || t.res > t.step){ 
            t.res <- 0
            for(m in which(w!=0)){
                if(is.na(v[1,1,m])){
                    next
                }
                tr <- which( rep(1,dim(v)[2])%*%(abs(v[1,,]-rep(v[1,,m], times=M))
                                                 +abs(v[2,,]-rep(v[2,,m], times=M))
                                                 +abs(v[3,,]-rep(v[3,,m], times=M)))==0)
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
        
############# simulation step #####
        rates.S <- (rep(beta/N,N)%*%v[1,,])*(rep(1,N)%*%v[3,,]) #infection rates
        rates.I <- rep(gamma,N)%*%v[1,,] # recovery rates
        rates.R <- rep(delta,N)%*%v[2,,] # loss of immunity rates
        rate.M <- rates.S + rates.I + rates.R # particle rates
        
        t.next <- rexp(1,sum(rates.S,rates.I, rates.R))
        if(t + t.next >= t.max){
            s <- s+1
            t <- t.max
            next
        }else{
            m <- max(which(cumsum(c(0,rate.M)/sum(rate.M))<runif(1))) # which particle event occurs on
            u1 <- runif(1)
            if(u1 < rates.R[m]/rate.M[m]){
                ### Loss of Immunity Event Rk -> Sk
                jr <- sample.int(N,1,prob=v[2,,m]) # pick immune strain
                v[2:3,jr,m] <- v[2:3,jr,m] + c(-1,1) # +1 S, -1 R
                t <- t + t.next
                t.res <- t.res + t.next
                s <- s + 1
                
            }else if( u1 < (rates.R[m] + rates.I[m])/rate.M[m]){
                ### Recovery Event Ik -> Rk
                j <- sample.int(N,1,prob=v[1,,m]) # pick infected strain
                v[1:2,j,m] <- v[1:2,j,m] + c(-1,1) #+1 R, -1 I
                if(v[1,j,m]==0){
                    if(j < N){
                        v[1,j:N,m] <- c(v[1,(j+1):N,m],0)
                        if(j < N-1){
                            v[2,j:N,m] <- c(v[2,j,m]+v[2,j+1,m], v[2,(j+2):N,m], 0)
                            v[3,j:N,m] <- c(v[3,j,m]+v[3,j+1,m], v[3,(j+2):N,m], 0)
                        }else{
                            v[2,(N-1):N,m] <- c(v[2,N-1,m]+v[2,N,m],0)
                            v[3,(N-1):N,m] <- c(v[3,N-1,m]+v[3,N,m],0)
                        }
                    }
                }
                t <- t + t.next
                t.res <- t.res + t.next
                s <- s + 1
                
            }else{
                js <- sample.int(N,1,prob=v[3,,m]) #pick susceptible
                if(runif(1) < theta){
                ### Mutation Event Sk -> I_K+1
                    v[3,js,m] <- v[3,js,m]-1
                    v[,,m] <- c(1,0,0,v[,-N,m])
                    t <- t + t.next
                    t.res <- t.res + t.next
                    s <- s + 1
                }else{
                ### Infection Event Sk -> Ij
                    ji <- sample.int(N,1,prob=v[1,,m]) # pick infective to infect
                    if(js > ji){
                        v[1,ji,m] <- v[1,ji,m]+1
                        v[3,js,m] <- v[3,js,m]-1
                        t <- t + t.next
                        t.res <- t.res + t.next
                        s <- s + 1
                    }else{
                        next
                    }
                }
            }
        }
        
############# reweighting #####
        w[m] <- w[m]*(sum(v[1,,m])!=0)
        w <- w/sum(w)
        
############# taking sample step #####
        if(t > t.burn + f*t.delay){ 
            f <- f+1
            samples.w[((f-1)*M + 1):(f*M)] <- w
            samples[,,((f-1)*M + 1):(f*M)] <- v
        }
    }
    
    
    
##### Final Sampling Step #####
    f <- f+1
    samples.w[((f-1)*M + 1):(f*M)] <- w
    samples[,,((f-1)*M + 1):(f*M)] <- v
    
    samples.w[samples.w==0] <- NA
    samples <- samples[,,!is.na(samples.w)]
    samples.w <- samples.w[!is.na(samples.w)]
    samples.w <- samples.w/sum(samples.w)
    
    v.summary <- matrix(0,ncol=N,nrow=3)
    for(j in 1:3){
        v.summary[j,] <- samples[j,,]%*%samples.w
    }
    
    K <- rep(1,N)%*%(samples[1,,]>0)
    I <- rep(1,N)%*%samples[1,,]
    R <- rep(1,N)%*%samples[2,,]
    A <- matrix(0,N,N)
    for(i in 1:N){
        for(j in 1:N){
            A[i,j] <- 1*(j >= i)
        }
    }
    B <- matrix(0,N,dim(samples)[3])
    for(i in 1:N){
        for(j in 1:dim(samples)[3]){
            B[i,j] <- max(K[j]+i-1,0)
        }
    }
    
    z <- v.summary[1,]/sum(v.summary[1,])
    suma <- c(0,unlist(cumsum((v.summary[1,] + v.summary[2,] + v.summary[3,])/N)[-N]))
    sumb <- sum(z*suma)
    R0 <- (sumb/N*(1-theta) + theta)*beta*(N-1)/N
    
    
    iw <- rep(NA,3)
    I <- samples[1,,] #I
    S <- samples[3,,] #S
    R <- samples[2,,] #R
    
    iw[1] <- (rep(1,N)%*%samples[1,,])%*%samples.w
    
    iw[2] <- K %*% samples.w #K mean
    
    rho <- (1/N) * { (rep(1,N) %*% (I+R)) + 
                     ((rep(1,N) %*% (S * (A%*%I)))/(rep(1,N)%*%I))  }%*%samples.w
    
    iw[3] <- rho #rhoI
    #iw[6] <- {(1 / (N * K)) * {(rep(1, N)%*%(samples[3,,] * B)) + (I * (K - 1)) + (R * K)}} %*% samples.w #rhoK
    #iw[7] <- R0 #r0
    
    Irep <- matrix(rep(rep(1,N)%*%samples[1,,], each=N), nrow=N)
    uI <- (samples[1,,]/Irep)%*%samples.w
    
    ### uL(j) = \sum_m=1^M w^(m) \sum_{j<k} (I_k^(m) + S_k^(m) + R_k^(m))/N
    Kback <- samples[1,(N:1),] + samples[2,(N:1),] + samples[3,(N:1),]
    revcumsum <- function(x){c(rev(cumsum(x[1:(length(x)-1)])), 0)}
    uL <- (1/N) * (apply(Kback, 2, revcumsum) %*% samples.w)
    
    ### R_Q
    rq <- beta*(theta + (1-theta)*(t(uL) %*% uI))/gamma
    iw[4] <- rq
    
    if(csv){
        write.csv(cbind(matrix(aperm(samples,c(3,1,2)),nrow=length(samples.w), ncol=3*N),matrix(samples.w, ncol=1)), 
              file=paste0("./samples_b",beta,"_d",delta,"_t",theta,"_n",N,".csv"), row.names=F)
        #### How to read samples into V and W again.
        # D <- read.csv("samples_b2_t0.05_n25.csv")
        # WW <- D[,dim(D)[2]]
        # VV <- array(t(D[,1:(2*N)]),dim=c(2,N,dim(D)[1]))
        ####
    }
        
        return(list(v=v.summary[c(1,3,2),N:1], iw=iw))
}