#### New compilation for running on one PC

library("stats")

#### EDIT THIS LINE TO SHOW THE RIGHT FOLDER
#setwd("./EvolvingStrains/EvolvePaper")
#source('SiStrain_NoTrace2.R')
#source('sirStrain_NoTrace2.R')


#####  Expected Proportion of Infectives [Figure 5] #####

betaRange <- 1:40
NRange <- seq(10, 200, by=10)
thetaRange <- c(0.05, 0.5)
deltaRange <- c(NA, 0.5, 2.5)
betaRange2 <- c(0.5, 2.5)
gamma <- 1
theta <- 0.4
beta <- 2
N <- 100

figA <- expand.grid(beta=betaRange,
                   theta=thetaRange,
                   delta=deltaRange,
                   gamma=gamma,
                   N=N)

figB <- expand.grid(beta=betaRange2,
                   theta=theta,
                   delta=deltaRange,
                   gamma=gamma,
                   N=NRange)

RQ <- matrix(NA, nrow=2, ncol=length(betaRange)) #matrix to include data for graphing
#file.remove("./GraphImmThetaDataC.csv") #resets the file
st<-Sys.time()
for(i in 1:nrow(figA)){
    if(is.na(figA$delta[i])){
      x <- siStrainNT(beta = figA$beta[i], gamma=figA$gamma[i],
                      theta=figA$theta[i], N=figA$N[i], 
                      M = 100, t.max = 100,  
                      t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
    }else{
      x <- sirStrainNT(beta = figA$beta[i], gamma=figA$gamma[i],
                       theta=figA$theta[i], delta=figA$delta[i], N=figA$N[i], 
                       M = 100, t.max = 100,  
                       t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
    }
    figA$I[i] <- x$iw[1] / figA$N[i]
    write.table(matrix(c(figA[i,], x$iw),nrow=1), 
                file="./Fig5_NewData.csv", 
                sep=",", row.names=F, col.names=F,append=T)
  #print(i,"Estimated duration:",round(length(thetaRange)*difftime(Sys.time(),st,units="mins")/i,0),"minutes.\n")
  cat(i,"Estimated duration:",
      round(length(thetaRange)*difftime(Sys.time(),st,units="mins")/i,0),"minutes.\n")
}

st<-Sys.time()
for(i in 1:nrow(figB)){
  if(is.na(figB$delta[i])){
    x <- siStrainNT(beta = figB$beta[i], gamma=figB$gamma[i],
                    theta=figB$theta[i], N=figB$N[i], 
                    M = 100, t.max = 100,  
                    t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
  }else{
    x <- sirStrainNT(beta = figB$beta[i], gamma=figB$gamma[i],
                     theta=figB$theta[i], delta=figB$delta[i], N=figB$N[i], 
                     M = 100, t.max = 100,  
                     t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
  }
  figB$I[i] <- x$iw[1] / figB$N[i]
  write.table(matrix(c(figB[i,], x$iw),nrow=1), 
              file="./Fig5_NewData.csv", 
              sep=",", row.names=F, col.names=F,append=T)
  #print(i,"Estimated duration:",round(length(thetaRange)*difftime(Sys.time(),st,units="mins")/i,0),"minutes.\n")
  cat(i,"Estimated duration:",
      round(length(thetaRange)*difftime(Sys.time(),st,units="mins")/i,0),"minutes.\n")
}

library(reshape2)
castA <- dcast(figA, delta+theta ~ beta, value.var="I")
castB <- dcast(figB, delta+beta ~ N, value.var="I")


pdf("./Fig5a_Newdata.pdf", width=4, height=5, pointsize=10, title="beta vary, theta=0.05")
par(mar=c(3,3,0.4,0.4), mgp=c(2,1,0))
plot(betaRange, castA[1,], type='l', lwd=2, xlim=c(0,40), ylim=c(0,35),
     xlab=expression(beta), ylab="Expected Proportion of Infectives")
lines(betaRange, castA[2,], lty=2, lwd=2)
lines(betaRange, castA[3,], col="grey70", lwd=2)
lines(betaRange, castA[4,], col="grey40", lwd=2)
lines(betaRange, castA[5,], col="grey70", lwd=2, lty=2)
lines(betaRange, castA[6,], col="grey40", lwd=2, lty=2)
legend("topleft", legend=c(expression(list("SIS", theta==0.05)), expression(list("SIS",theta==0.5)),
                           expression(list(delta==0.5, theta==0.05)), expression(list(delta==0.5,theta==0.5)),
                           expression(list(delta==2.5, theta==0.05)), expression(list(delta==2.5,theta==0.5))),
       col=c(1,1,"grey70","grey70","grey40","grey40"), lty=c(1,2,1,2,1,2))

dev.off()


pdf("./Fig5b_Newdata.pdf", width=4, height=5, pointsize=10, title="beta vary, theta=0.05")
par(mar=c(3,3,0.4,0.4), mgp=c(2,1,0))
plot(NRange, castB[1,], type='l', lwd=2, xlim=c(0,40), ylim=c(0,35),
     xlab=expression(beta), ylab="Expected Proportion of Infectives")
lines(NRange, castB[2,], lty=2, lwd=2)
lines(NRange, castB[3,], col="grey70", lwd=2)
lines(NRange, castB[4,], col="grey40", lwd=2)
lines(NRange, castB[5,], col="grey70", lwd=2, lty=2)
lines(NRange, castB[6,], col="grey40", lwd=2, lty=2)
legend("topleft", legend=c(expression(list("SIS",      beta==0.5)), expression(list("SIS",      beta==2.5)),
                           expression(list(delta==0.5, beta==0.5)), expression(list(delta==0.5, beta==2.5)),
                           expression(list(delta==2.5, beta==0.5)), expression(list(delta==2.5, beta==2.5))),
       col=c(1,1,"grey70","grey70","grey40","grey40"), lty=c(1,2,1,2,1,2))

dev.off()