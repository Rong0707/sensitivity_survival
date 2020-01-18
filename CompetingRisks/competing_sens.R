#two methods to get cumulative baseline hazard function
#newdata = data.frame(X1 = 0, X2 = 0, U=0, Z=0)
#ps = survfit(t1.fit, newdata=newdata)
#plot(ps$time, -log(ps$surv), type="l") #cumulative baseline hazard

#bh = basehaz(t1.fit, centered=F)
#plot(bh$time, bh$hazard, type="l") #cumulative baseline hazard
source("SimulateU.R")
source("emU.R")

competing_sens <- function(data, zetat, zetaz, B = 40, offset){
  t = data$t
  d1 = as.numeric(data$d==1)
  d2 = as.numeric(data$d==2)
  Z = data$Z
  X = data$X
  nx = dim(X)[2]
  U = data$U
  
  Z.fit = glm(Z ~ X, offset = zetaz * U, family = binomial(link="probit"))
  Coefz.trU = Z.fit$coefficients
  t1.fit = coxph(Surv(t, d1) ~ X + Z + offset(zetat[1] * U))
  Coeft1.trU = t1.fit$coefficients
  t2.fit = coxph(Surv(t, d2) ~ X + Z + offset(zetat[2] * U))
  Coeft2.trU = t2.fit$coefficients
  
  Z.fit = glm(Z ~ X, family = binomial(link="probit"))
  Coefz.noU = Z.fit$coefficients
  t1.fit = coxph(Surv(t, d1) ~ X + Z)
  Coeft1.noU = t1.fit$coefficients
  t2.fit = coxph(Surv(t, d2) ~ X + Z)
  Coeft2.noU = t2.fit$coefficients
  
  Uem = emU(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz)
  
  if(offset){
    #coefficients with simulated U
    Coefz = matrix(0, nrow = B, ncol = nx+1)
    Coeft1 = matrix(0, nrow = B, ncol = nx+1)
    Coeft2 = matrix(0, nrow = B, ncol = nx+1)
    #Usim = SimulateU(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, iter = B-1, offset=TRUE)
    
    for (j in 1:B){
      Usim = SimulateU(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, offset=TRUE)
      
      Z.fit = glm(Z ~ X, offset = zetaz * Usim$U, family=binomial(link="probit"))
      Coefz[j,] = Z.fit$coefficients
      
      t1.fit = coxph(Surv(t, d1) ~ X + Z + offset(zetat[1] * Usim$U))
      Coeft1[j,] = t1.fit$coefficients
      
      t2.fit = coxph(Surv(t, d2) ~ X + Z + offset(zetat[2] * Usim$U))
      Coeft2[j,] = t2.fit$coefficients
    }
    colnames(Coefz) = names(Z.fit$coefficients)
    colnames(Coeft1) = colnames(Coeft2)  = names(t1.fit$coefficients)
  }
  else{
    #coefficients with simulated U
    Coefz = matrix(0, nrow = B, ncol = nx+2)
    Coeft1 = matrix(0, nrow = B, ncol = nx+2)
    Coeft2 = matrix(0, nrow = B, ncol = nx+2)
    #Usim = SimulateU(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, offset=FALSE, iter=B-1)
    
    for (j in 1:B){
      Usim = SimulateU(t, cbind(d1,d2), Z, X, zetat = zetat, zetaz = zetaz, offset=FALSE)

      Z.fit = glm(Z ~ X + Usim$U, family=binomial(link="probit"))
      Coefz[j,] = Z.fit$coefficients
      
      t1.fit = coxph(Surv(t, d1) ~ X + Z + Usim$U)
      Coeft1[j,] = t1.fit$coefficients
      
      t2.fit = coxph(Surv(t, d2) ~ X + Z + Usim$U)
      Coeft2[j,] = t2.fit$coefficients
    }
    colnames(Coefz) = names(Z.fit$coefficients)
    colnames(Coeft1) = colnames(Coeft2)  = names(t1.fit$coefficients)
  }
  
  return (list(Coefz.simU = colMeans(Coefz), Coefz.trU = Coefz.trU, Coefz.noU = Coefz.noU, Coefz.emU = Uem$z.coef,
               Coeft1.simU = colMeans(Coeft1), Coeft1.trU = Coeft1.trU, Coeft1.noU = Coeft1.noU, Coeft1.emU = Uem$t.coef1,
               Coeft2.simU = colMeans(Coeft2), Coeft2.trU = Coeft2.trU, Coeft2.noU = Coeft2.noU, Coeft2.emU = Uem$t.coef2))
}