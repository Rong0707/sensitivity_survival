library(survival)
library(ggplot2)
library(directlabels)
source("competing_sens_var.R")
source("emU.R")
source("competing_final.R")
source("competing_sens_ipw.R")

#out1:clinical remission
a2cd_out1 <- read.csv("/CD_PS/a2cd_out1.csv", header = T)
t = a2cd_out1$clinremt
#clinremt is time to event.
d = a2cd_out1$clinrem2
#clinrem2 is indicator for events, if clinrem2=0 then censored, if clinrem2=1 then event, if clinrem2=2 then competing risk happened.
Z = a2cd_out1$gr2
#gr2 is the indicator for treatment, if gr2=1, then the patients receive the treatment, Vedolizumab.
X = qnorm(a2cd_out1$ps_ate)
ipw = a2cd_out1$ipw_ate

data = list(t=t, d=d, Z=Z, X= X)
data1 = data.frame(data)

#regression ignoring U
Z.fit = glm(Z ~ X, family = binomial(link="probit"), data=data)
Coefz.noU = Z.fit$coefficients
t1.fit = coxph(Surv(t, d==1) ~ X + Z, data = data)
Coeft1.noU = t1.fit$coefficients
t2.fit = coxph(Surv(t, d==2) ~ X + Z, data = data)
Coeft2.noU = t2.fit$coefficients

#ipw
t1.ipw = coxph(Surv(t, d==1) ~ Z, weights = ipw, robust = TRUE)

tau1.stoEM = tau1.em = tau1.ipw = data.frame()
#tau2.stoEM = tau2.em = data.frame()

zetaZ = seq(-2,2,by=0.5)
zetaT = seq(-2,2,by=0.5)

set.seed(2019)
zetat2 = 0

for(i in 1:length(zetaZ)){
  for(j in 1:length(zetaT)){
    #stochastic EM/regression
    tau.sto = competing_sens_var(data, zetat = c(zetaT[j],zetat2), zetaz = zetaZ[i], theta = 0.3, B=100)
    tau1.stoEM = rbind(tau1.stoEM, data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2, tau1 = tau.sto$tau1, tau1.se = tau.sto$tau1.se))
    #tau2.stoEM = rbind(tau2.stoEM, data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2, tau2 = tau.sto$tau2, tau2.se = tau.sto$tau2.se))
    
    #stochastic EM/IPW
    tau.ipw = competing_sens_ipw(data, zetat = c(zetaT[j],zetat2), zetaz = zetaZ[i], theta = 0.3, B=100)
    tau1.ipw = rbind(tau1.ipw, data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2, tau1 = tau.ipw$tau1, tau1.se = tau.ipw$tau1.se))
   
    #EM/regression
    tau.em = emU(t, cbind(d==1, d==2), Z, X, zetat = c(zetaT[j],zetat2), zetaz = zetaZ[i], theta = 0.3)
    data1$p = tau.em$p
    tau.em.final = competing_final(data1, zetat = c(zetaT[j],zetat2), zetaz = zetaZ[i], z.coef = tau.em$z.coef[1:2])
    tau1.em = rbind(tau1.em, data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2, tau1 = tau.em.final$t1.coef[2], tau1.se = tau.em.final$t1.coef.se[2]))
    #tau2.em = rbind(tau2.em, data.frame(zetaz = zetaZ[i], zetat1 = zetaT[j], zetat2 = zetat2, tau2 = tau.em.final$t2.coef[2], tau2.se = tau.em.final$t2.coef.se[2]))
  }
}

tau1.stoEM$t = tau1.stoEM$tau1/tau1.stoEM$tau1.se
tau1.em$t = tau1.em$tau1/tau1.em$tau1.se
tau1.ipw$t = tau1.ipw$tau1/tau1.ipw$tau1.se

#tau2.stoEM$t = tau2.stoEM$tau2/tau2.stoEM$tau2.se
#tau2.em$t = tau2.em$tau2/tau2.em$tau2.se

pdf("CD_tau1_0_EM.pdf", width=6, height=5)
g = ggplot(tau1.em, aes(zetaz, zetat1)) +
  stat_contour(aes(z = tau1, colour = ..level..)) +
  stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
  xlab("Coef. on U in model for treatment") +
  ylab("Coef. on U in model for type 1 event") +
  annotate("text",x=0,y=0, label = "0.0605") + 
  theme_bw()
direct.label(g, method="bottom.pieces")
dev.off()

pdf("CD_tau1_0_stoEM.pdf", width=6, height=5)
g = ggplot(tau1.stoEM, aes(zetaz, zetat1)) +
  stat_contour(aes(z = tau1, colour = ..level..)) +
  stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
  xlab("Coef. on U in model for treatment") +
  ylab("Coef. on U in model for type 1 event") +
  annotate("text",x=0,y=0, label = "0.0605") + 
  theme_bw()
direct.label(g, method="bottom.pieces")
dev.off()


pdf("CD_tau1_0_ipw.pdf", width=6, height=5)
g = ggplot(tau1.ipw, aes(zetaz, zetat1)) +
  stat_contour(aes(z = tau1, colour = ..level..)) +
  stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
  xlab("Coef. on U in model for treatment") +
  ylab("Coef. on U in model for type 1 event") +
  annotate("text",x=0,y=0, label = "-0.0563") + 
  theme_bw()
direct.label(g, method="bottom.pieces")
dev.off()

tau1.em$value = paste(round(tau1.em$tau1,4), " (", round(tau1.em$tau1.se,4), ")", sep="")
tau1.stoEM$value = paste(round(tau1.stoEM$tau1,4), " (", round(tau1.stoEM$tau1.se,4), ")", sep="")
tau1.ipw$value = paste(round(tau1.ipw$tau1,4), " (", round(tau1.ipw$tau1.se,4), ")", sep="")

write.csv(tau1.stoEM, "CD_tau1_0_stoEM.csv")
write.csv(tau1.em, "CD_tau1_0_EM.csv")
write.csv(tau1.ipw, "CD_tau1_0_ipw.csv")