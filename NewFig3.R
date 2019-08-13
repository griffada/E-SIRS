#### New compilation for running on one PC

library("stats")

#### EDIT THIS LINE TO SHOW THE RIGHT FOLDER
#setwd("./EvolvingStrains/EvolvePaper")
#source('SiStrain_NoTrace2.R')
#source('sirStrain_NoTrace2.R')


#####  RQ figure [Figure 3] #####
thetaRange <- seq(0,1,by=0.025)[-1]
betaRange <- c(0.5,2)
gamma <- 1
N <- 100

RQ <- matrix(NA, nrow=2, ncol=length(thetaRange)) #matrix to include data for graphing
#file.remove("./GraphImmThetaDataC.csv") #resets the file
st<-Sys.time()
for(i in 1:length(thetaRange)){
  for(j in 1:2){
    x <- siStrainNT(beta = betaRange[j], gamma=1, theta=thetaRange[i], N=N, 
                    M = 100, t.max = 100,  
                    t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
    write.table(matrix(c(betaRange[j], gamma, NA, thetaRange[i], N, x$iw),nrow=1), 
                file="./Fig3_NewData.csv", 
                sep=",", row.names=F, col.names=F,append=T)
    rhoTheta[j,i] <- x$iw[4]
  }
  #print(i,"Estimated duration:",round(length(thetaRange)*difftime(Sys.time(),st,units="mins")/i,0),"minutes.\n")
  cat(i,"Estimated duration:",
      round(length(thetaRange)*difftime(Sys.time(),st,units="mins")/i,0),"minutes.\n")
}

R0 <- c(0.5,2)
Rstar <- outer(betaRange, thetaRange,
               FUN=function(x,y){
                 ifelse(x*(1-y) < 1, (x*y)/(1 - x*(1-y)) ,Inf)})

pdf("Fig3_New.pdf", width=4, height=4, pointsize=9,
    title="Reproduction Numbers")
  par(mfrow=c(1,2), mgp=c(2,1,0), mar=c(3,3,0.4,0.4))
  plot(c(0,1), c(R0[1], R0[1]),
       ylim=c(0,1), xlim=c(0,1),
       xlab=expression(theta), ylab ="Reproduction Number",
       type='l', lwd=2)
  lines(c(0, thetaRange), c(0,Rstar[1,]), lwd=2, lty=2)
  lines(c(0, thetaRange), c(0,RQ[1,]),    lwd=2, lty=3)
  legend("topleft",
         legend=c(expression(R[0]),expression(R["*"]),expression(R[Q])),
         lty=1:3, lwd=2)
  
  plot(c(0,1), c(R0[2], R0[2]),
       ylim=c(0,20), xlim=c(0,1),
       xlab=expression(theta), ylab ="Reproduction Number",
       type='l', lwd=2)
  lines(c(0, thetaRange), c(0,Rstar[2,]), lwd=2, lty=2)
  lines(c(0, thetaRange), c(0,RQ[2,]),    lwd=2, lty=3)
  legend("topleft",
         legend=c(expression(R[0]),expression(R["*"]),expression(R[Q])),
         lty=1:3, lwd=2)
dev.off()
