#### New compilation for running on one PC

library("stats")
library("foreach")
library("doParallel")
registerDoParallel(cores=7)

#### EDIT THIS LINE TO SHOW THE RIGHT FOLDER
#setwd("./EvolvePaper")
source('SiStrain_NoTrace2.R')
source('sirStrain_NoTrace2.R')


#####  Strain Index Variability [Figure 7] #####

betaRange <- c(0.5, 1, 2)
br_labels <- c("05", "1", "2")
thetaRange <- c(0.05,0.5, 0.95)
tr_labels <- c("05", "50", "95")
deltaRange <- c(2.5,0.5)
dr_labels <- c("SIS", "SIRS")
Nrange <- c(50, 100, 200)
beta <- 2
theta <- 0.5
gamma <- 1
delta <- 2
N <- 100


if(!dir.exists("./VObjects")){ # make a VObjects folder
  dir.create("./VObjects")
}

### Vary beta  
st<-Sys.time()
output<-foreach(i=1:3) %dopar% {
  
  # Run through each choice of beta for each selection of theta and delta.
  x1 <- siStrainNT(beta = betaRange[i], gamma = gamma, theta = theta, N = N, 
                   M = 100, t.max = 100,  
                   t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
  write.table(matrix(c(betaRange[i], gamma, NA, theta, N, x1$iw),nrow=1),
              file="./Fig7a_NewData.csv", 
              sep=",", row.names=F, col.names=F,append=T)
  write.csv(x1$v, paste0("./VObjects/Fig7_beta", br_labels[i],"_v.csv"))
  
}


### Vary theta
st<-Sys.time()
output<-foreach(i=1:3) %dopar% {
  
  # Run through each choice of beta for each selection of theta and delta.
  x1 <- siStrainNT(beta = beta, gamma = gamma, theta = thetaRange[i], N = N, 
                   M = 100, t.max = 100,  
                   t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
  write.table(matrix(c(beta, gamma, NA, thetaRange[i], N, x1$iw),nrow=1),
              file="./Fig7b_NewData.csv", 
              sep=",", row.names=F, col.names=F,append=T)
  write.csv(x1$v, paste0("./VObjects/Fig7_theta", tr_labels[i],"_v.csv"))
}


### Vary N
st<-Sys.time()
output<-foreach(i=1:3) %dopar% {
  
  # Run through each choice of N for each selection of theta and delta.
  x1 <- siStrainNT(beta = betaRange[i], gamma = gamma, theta = theta, N = Nrange[i], 
                   M = 100, t.max = 100,  
                   t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
  write.table(matrix(c(beta, gamma, NA, theta, Nrange[i], x1$iw),nrow=1),
              file="./Fig7c_NewData.csv", 
              sep=",", row.names=F, col.names=F,append=T)
  write.csv(x1$v, paste0("./VObjects/Fig7_N", Nrange[i],"_v.csv"))
}


### Vary delta

  # Run through each choice of beta for each selection of theta and delta.
  x1 <- siStrainNT(beta = betaRange[i], gamma = gamma, theta = theta, N = N, 
                   M = 100, t.max = 100,  
                   t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
  write.table(matrix(c(betaRange[i], gamma, NA, thetaRange[1], N, x1$iw),nrow=1),
              file="./Fig7d_NewData.csv", 
              sep=",", row.names=F, col.names=F,append=T)
  write.csv(x1$v, paste0("./VObjects/Fig7_delta", dr_labels[i],"_v.csv"))
  
  x2 <- sirStrainNT(beta = betaRange[i], gamma = gamma, theta = theta, delta=delta, N = N, 
                   M = 100, t.max = 100,  
                   t.step = 10, t.delay = 1, lambda = 0.4, csv = T)
  write.table(matrix(c(beta, gamma, delta, theta, N, x1$iw),nrow=1),
              file="./Fig7d_NewData.csv", 
              sep=",", row.names=F, col.names=F,append=T)
  write.csv(x2$v, paste0("./VObjects/Fig7_delta", dr_labels[i],"_v.csv"))

  
  

  
  
  
  
  vdata1 <- read.csv(paste0("./VObjects/Fig7_beta", br_labels[1],"_v.csv"), header=F) #beta=.5
  vdata2 <- read.csv(paste0("./VObjects/Fig7_beta", br_labels[2],"_v.csv"), header=F) #beta=1
  vdata3 <- read.csv(paste0("./VObjects/Fig7_beta", br_labels[3],"_v.csv"), header=F) #beta=2
 
  pdf("./Figures/Fig7a.pdf", width=3.5, height=5, pointsize=10, title="Strain Index, Vary Beta")
  par(mar=c(4,4,0.5,0.5))
  plot((1-N):0, vdata1[1,]/N, xlim=c(-20,0), ylim=c(0,.6), type='l', lwd=2, xlab="Strain Index", ylab="Mean Proportion")
  lines((1-N):0, (vdata1[1,]+vdata1[2,])/N, lwd=2, lty=2)
  lines((1-N):0, vdata2[1,]/N, col="grey40", lwd=2)
  lines((1-N):0, (vdata2[1,]+vdata2[2,])/N, lwd=2, col="grey40", lty=2)
  lines((1-N):0, vdata3[1,]/N, col="grey70", lwd=2)
  lines((1-N):0, (vdata3[1,]+vdata3[2,])/N, lwd=2, col="grey70", lty=2)
  legend("topleft", legend=c(expression(paste(beta==0.5,   " - Infectives")),
                             expression(paste(beta==0.5, " - Total Indivs")),
                             expression(paste(beta==1,     " - Infectives")),
                             expression(paste(beta==1,   " - Total Indivs")),
                             expression(paste(beta==2,     " - Infectives")),
                             expression(paste(beta==2,   " - Total Indivs"))),
         col=c(1,1,"grey40","grey40","grey70","grey70"), lwd=2, lty=c(1,2,1,2,1,2))
  dev.off()
  
  
  
  vdata1 <- read.csv(paste0("./VObjects/Fig7_theta", tr_labels[1],"_v.csv"), header=F) 
  vdata2 <- read.csv(paste0("./VObjects/Fig7_theta", tr_labels[2],"_v.csv"), header=F) 
  vdata3 <- read.csv(paste0("./VObjects/Fig7_theta", tr_labels[3],"_v.csv"), header=F) 
  
  pdf("./Figures/Fig7b.pdf", width=3.5, height=5, pointsize=10, title="Strain Index, Vary Theta")
  par(mar=c(4,4,0.6,0.5))
  plot((1-N):0, vdata1[1,]/N, xlim=c(-20,0), ylim=c(0,0.5), type='l', lwd=2, xlab="Strain Index", ylab="Mean Proportion")
  lines((1-N):0, (vdata1[1,]+vdata1[2,])/N, lwd=2, lty=2)
  lines((1-N):0, vdata2[1,]/N, col="grey40", lwd=2)
  lines((1-N):0, (vdata2[1,]+vdata2[2,])/N, lwd=2, col="grey40", lty=2)
  lines((1-N):0, vdata3[1,]/N, col="grey70", lwd=2)
  lines((1-N):0, (vdata3[1,]+vdata3[2,])/N, lwd=2, col="grey70", lty=2)
  legend("topleft", legend=c(expression(paste(theta==0.05,   " - Infectives")),
                             expression(paste(theta==0.05, " - Total Indivs")),
                             expression(paste(theta==0.5,     " - Infectives")),
                             expression(paste(theta==0.5,   " - Total Indivs")),
                             expression(paste(theta==0.95,     " - Infectives")),
                             expression(paste(theta==0.95,   " - Total Indivs"))),
         col=c(1,1,"grey40","grey40","grey70","grey70"), lwd=2, lty=c(1,2,1,2,1,2))
  dev.off()
  
  
  
  
  vdata1 <- read.csv(paste0("./VObjects/Fig7_N", Nrange[1],"_v.csv"), header=F) 
  vdata2 <- read.csv(paste0("./VObjects/Fig7_N", Nrange[2],"_v.csv"), header=F) 
  vdata3 <- read.csv(paste0("./VObjects/Fig7_N", Nrange[3],"_v.csv"), header=F) 
  
  pdf("./Figures/Fig7c.pdf", width=3.5, height=5, pointsize=10, title="Strain Index, Vary N")
  par(mar=c(4,4,0.5,0.5))
  plot(-49:0, vdata1[1,]/Nrange[1], xlim=c(-30,0), ylim=c(0,0.2), type='l', lwd=2, xlab="Strain Index", ylab="Mean Proportion")
  lines(-49:0, (vdata1[1,]+vdata1[2,])/20, lwd=2, lty=2)
  lines(-99:0, vdata2[1,]/Nrange[2], col="grey40", lwd=2)
  lines(-99:0, (vdata2[1,]+vdata2[2,])/Nrange[2], lwd=2, col="grey40", lty=2)
  lines(-199:0, vdata3[1,]/Nrange[3], col="grey70", lwd=2)
  lines(-199:0, (vdata3[1,]+vdata3[2,])/Nrange[3], lwd=2, col="grey70", lty=2)
  legend("topleft", legend=c(expression(paste(N==50,   " - Infectives")),
                             expression(paste(N==50, " - Total Indivs")),
                             expression(paste(N==100,     " - Infectives")),
                             expression(paste(N==100,   " - Total Indivs")),
                             expression(paste(N==200,     " - Infectives")),
                             expression(paste(N==200,   " - Total Indivs"))),
         col=c(1,1,"grey40","grey40","grey70","grey70"), lwd=2, lty=c(1,2,1,2,1,2))
  dev.off()
  
  
  
  
  vdata1 <- read.csv(paste0("./VObjects/Fig7_delta", dr_labels[1],"_v.csv"), header=F) 
  vdata2 <- read.csv(paste0("./VObjects/Fig7_delta", dr_labels[2],"_v.csv"), header=F) 
  
  pdf("./Figures/Fig7d.pdf", width=3.5, height=5, pointsize=10, title="Strain Index, Vary Delta")
  
  par(mar=c(4,4,0.5,0.5))
  plot((1-N):0, vdata1[1,]/N, type='l', lwd=2, xlab="Number of Infectives", ylab="Mean Proportion", ylim=c(0,0.3))
  lines((1-N):0, vdata1[2,]/N + vdata1[1,]/N, col=1, lty=2, lwd=2)
  
  lines((1-N):0, vdata2[1,]/N, col="grey50", lty=1, lwd=2)
  lines((1-N):0, vdata2[2,]/N + vdata2[1,]/N, col="grey50", lty=2, lwd=2)
  lines((1-N):0, vdata2[2,]/N + vdata2[1,]/N + vdata2[3,]/N, col="grey50", lty=4, lwd=2)

  legend("topleft", legend=c(expression(paste("E-SIS",   " - Infectives")),
                             expression(paste("E-SIS",  " - Total Indivs")),
                             expression(paste("E-SIRS, ", delta==2,    " - Infectives")),
                             expression(paste("E-SIRS, ", delta==2,  " - Infects + Immunes")),
                             expression(paste("E-SIRS, ", delta==2,     " - Total"))),
         col=c(1,1,"grey50", "grey50", "grey50"), lwd=2, lty=c(1,2,1,4,2))
  dev.off()
