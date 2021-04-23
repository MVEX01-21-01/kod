library(spatstat)
library(GET)
moderatePoints<-readRDS(file.choose())
normalPoints<-readRDS(file.choose())

moderateBranch<-readRDS(file.choose())
normalBranch<-readRDS(file.choose())
#----------------count number of points in each pattern


moderatePointsNumber<-c()
moderateBranchNumber<-c()
normalPointsNumber<- c()
normalBranchNumber<-c()

for(i in 1:7){
  
  moderatePointsNumber[i]<-moderatePoints[[i]]$n
  moderateBranchNumber[i]<-moderateBranch[[i]]$n
}
for(i in 1:8){
  normalPointsNumber[i]<-normalPoints[[i]]$n
  normalBranchNumber[i]<-normalBranch[[i]]$n
}
#----------------create multitype ppp object
moderateMulti<-list()
normalMulti<-list()
for (i in 1:7)
  moderateMulti[[i]]<-superimpose(point=moderatePoints[[i]],branch=moderateBranch[[i]])
for (i in 1:8)
  normalMulti[[i]]<-superimpose(point=normalPoints[[i]],branch=normalBranch[[i]])
#---------------calculate Kcross for multi and Kest for points data
moderateKcross<-list()
normalKcross<-list()
normalKest<-list()
moderateKest<-list()

for(i in 1:8){
  normalKcross[[i]]<-Kcross(normalMulti[[i]])
  normalKest[[i]]<-Kest(normalPoints[[i]])
}
for(i in 1:7){
  moderateKcross[[i]]<-Kcross(moderateMulti[[i]])
  moderateKest[[i]]<-Kest(moderatePoints[[i]])
}

#----------------compare Kest with Kcross
par(mfrow=c(2,8),mar=c(1,1,1,1))

for(i in 1:8){
  plot(normalKcross[[i]]$r,normalKcross[[i]]$iso,type="l",col="red",ylim=c(0,40000))
  lines(normalKcross[[i]]$r,normalKcross[[i]]$theo,type="l",col="green")
  lines(normalKest[[i]]$r,normalKest[[i]]$iso,type="l",col="black")
}
for(i in 1:7){
  plot(moderateKcross[[i]]$r,moderateKcross[[i]]$iso,type="l",col="red",ylim=c(0,40000))
  lines(moderateKest[[i]]$r,moderateKest[[i]]$theo,col="green")
  lines(moderateKest[[i]]$r,moderateKest[[i]]$iso,col="black")
}

#----------------plot points with branching points
par(mfrow=c(2,4),mar=c(1,1,1,1))
for(i in 1:7){ 
  plot(moderatePoints[[i]])
  points(moderateBranch[[i]],pch=4,col="red")
}
par(mfrow=c(2,4),mar=c(1,1,1,1))
for(i in 1:8){ 
  plot(normalPoints[[i]])
  points(normalBranch[[i]],pch=4,col="red")
}

#--------------estimate groupwise Kcross function
moderateKcrossGroup<-do.call(pool,c(moderateKcross,moderateKest))
normalKcrossGroup<-do.call(pool,c(normalKcross,normalKest))
KbarM<-do.call(pool,moderateKest)
KbarN<-do.call(pool,normalKest)





mod1<-do.call(pool,moderateKest)


#--------------estimate params for amtern and thomas process for each group using kcross
#moderateMultiMatclust<-matclust.estfun(moderateKcrossGroup)
#moderateMultiThomas<-thomas.estfun(moderateKcrossGroup)

#normalMuliMatclust<-matclust.estfun(normalKcrossGroup)
#normalMuliThomas<-thomas.estfun(normalKcrossGroup)

thomson<-list()
lambda<-grouped(intensitybar,data)
moderateMultiMatclust<-matclust.estK(moderateKcrossGroup,lambda=lambda$MODERATE)
moderateMultiThomas<-thomas.estK(moderateKcrossGroup,lambda=lambda$MODERATE)

normalMultiMatclust<-matclust.estK(normalKcrossGroup,lambda=lambda$NORMAL)
normalMultiThomas<-thomas.estK(normalKcrossGroup,lambda=lambda$NORMAL)

moderateThomas<-thomas.estK(KbarM,lambda=lambda$MODERATE)
#------------envelope test for multiKcross
handlers(handler_rstudio())
with_progress({
  p <- progressor(2*length(data$ppp))
  rp <- progress_aggregator(p)
  rp(envs.thomas   <- grouped(multiGET.composite, data, thomson, c(Gest), alpha=0.1, type='erl'))
  rp(envs.matclust <- grouped(multiGET.composite, data, morten, c(Gest), alpha=0.1, type='erl'))
})

#------------plot random models

par(mfrow=c(2,4),mar=c(1,1,1,1))
for(i in 1:4)
  plot(moderatePoints[[i]])
for(i in 1:4)
  plot(rMatClust(fit.matclust.K$MODERATE$modelpar[1],fit.matclust.K$MODERATE$modelpar[2],fit.matclust.K$MODERATE$modelpar[3],win=moderatePoints[[1]]$window ))

par(mfrow=c(2,4),mar=c(1,1,1,1))
for(i in 1:4)
  plot(rMatClust(moderateMultiMatclust$modelpar[1],moderateMultiMatclust$modelpar[2],3.5,win=moderatePoints[[1]]$window ))
for(i in 1:4)
  plot(moderatePoints[[i]])

#-------------simult.env test

SIMS<-list()
for(i in 1:7) SIMS[[i]] <- as.solist(rMatClust(moderateMultiMatclust$modelpar[1],moderateMultiMatclust$modelpar[2],
                                               moderateMultiMatclust$modelpar[3],saveparents = TRUE,win=ModerateData[[i]]$window,nsim=998))

for(i in 8:15) SIMS[[i]] <- as.solist(rMatClust(normalMultiMatclust$modelpar[1],normalMultiMatclust$modelpar[2],
                                                normalMultiMatclust$modelpar[3],saveparents = TRUE,win=NormalData[[i-7]]$window,nsim=998))
a<-cbind(data, hyperframe(Sims=SIMS))
EE <- with(a, envelope(ppp, Gest, nsim=499,simulate=Sims,global=TRUE))
plot(EE)


#----------------
#fit matern and thomas parameters for branching point data
Km<-list()
Kn<-list()
for(i in 1:7){
  Km[[i]]<-Kest(moderateBranch[[i]])
}
for(i in 1:8){
  Kn[[i]]<-Kest(normalBranch[[i]])
}
KbarMB<-do.call(pool,Km)
KbarNB<-do.call(pool,Kn)
lambdaBM<-intensitybar(moderateBranch)
lambdaBN<-intensitybar(normalBranch)

fitBranchMM<-matclust.estK(KbarMB,lambda=lambdaBM)
fitBranchMN<-matclust.estK(KbarNB,lambda=lambdaBN)
fitBranchTM<-thomas.estK(KbarMB,lambda=lambdaBM)
fitBranchTN<-thomas.estK(KbarNB,lambda=lambdaBN)
#-----------------try fits, test different parameters
#-------------
c=0.5
k1<-fitBranchMM$modelpar[[1]]
mu1<-fitBranchMM$modelpar[[3]]
mu2<-avgEndMod/avgBrMod
R1<-fitBranchMM$modelpar[[2]]
R2<-20

listMod<-list()
listMod2<-list()
listMod3<-list()
#par(mfrow=c(3,7))

for(i in 1:100){
  
  listMod[[i]]<-rMatClust(k1,R1,mu1,win=moderateBranch[[1]]$window)
  #rPMatClus...matern around parents in listMod[[i]]
  listMod2[[i]]<-rPMatClust(listMod[[i]],scale=R2,mu=mu2,win = listMod[[1]]$window)
  #plot(listMod2[[i]])
}

for(i in 1:7){
  #plot(ModerateData[[i]])
}

for(i in 1:100){
  listMod3[[i]]<-rMatClust(fit.matclust$MODERATE$modelpar[[1]],fit.matclust$MODERATE$modelpar[[2]],fit.matclust$MODERATE$modelpar[[3]],win=ModerateData[[1]]$window)
  #plot(listMod3[[i]])
  
}
par(mfrow=c(1,1))

plot(do.call(pool,lapply(listMod2,Kest))$pooliso,type="l",xlim=c(0,100),ylim=c(0,8000))
lines(do.call(pool,lapply(listMod3,Kest))$pooliso,col="red")
lines(do.call(pool,lapply(ModerateData,Kest))$pooliso,col="green")
print(getAvgSim(rPMatClust(listMod[[i]],scale=R2,mu=mu2,win = listMod[[1]]$window,nsim=1000)))
#-----------------
#---------------------- plot simulations from rPmatClust,data and rmatClust
par(mfrow=c(3,7))
for(i in 1:7){
  plot(listMod2[[i]])
}
for(i in 1:7){
  plot(ModerateData[[i]])
}
for(i in 1:7){
  plot(listMod3[[i]])
}
#---------------
#

#-----------------
#how does change of parameters affect k function
cols=c("gray8","gray20","gray39","gray48","gray60")

k1<-fitBranchMM$modelpar[[1]]
mu1<-fitBranchMM$modelpar[[3]]
mu2<-avgEndMod/avgBrMod
R1<-fitBranchMM$modelpar[[2]]
R2<-25

listMod<-list()
listMod2<-list()
listMod3<-list()
#par(mfrow=c(3,7))
for(i in 1:10){
  
  listMod[[i]]<-rMatClust(k1,R1,mu1,win=moderateBranch[[1]]$window)
  listMod2[[i]]<-rPMatClust(listMod[[i]],scale=R2,mu=mu2,win = listMod[[1]]$window)
  #plot(listMod2[[i]])
}
par(mfrow=c(1,1))

plot(do.call(pool,lapply(listMod2,Kest))$pooliso,type="l",xlim=c(0,100),ylim=c(0,8000))
lines(do.call(pool,lapply(listMod3,Kest))$pooliso,col="red")
lines(do.call(pool,lapply(ModerateData,Kest))$pooliso,col="green")
print(getAvgSim(rPMatClust(listMod[[i]],scale=R2,mu=mu2,win = listMod[[1]]$window,nsim=1000)))
#------------------
















