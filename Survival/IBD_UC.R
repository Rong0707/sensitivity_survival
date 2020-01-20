library(survival)
library(ggplot2)
library(directlabels)
source("surv_sens_var.R")
source("emU_surv.R")
source("surv_final.R")
source("surv_sens_ipw.R")

#out1:clinical remission
a2uc_out1 <- read.csv("/UC Effective/a2uc_out1.csv", header = T)
t = a2uc_out1$clinremt
#clinremt is time to event.
d = a2uc_out1$clinrem2==1
#clinrem2 is indicator for events, if clinrem2=0 then censored, if clinrem2=1 then event, if clinrem2=2 then competing risk happened.
Z = a2uc_out1$gr2
#gr2 is the indicator for treatment, if gr2=1, then the patients receive the treatment, Vedolizumab.
X = qnorm(a2uc_out1$ps_ate)
ps = a2uc_out1$ps_ate
#ipw = (sum(Z)/719) * Z/ps + (1-(sum(Z)/719))*(1-Z)/(1-ps)
ipw = a2uc_out1$ipw_ate

data = list(t=t, d=d, Z=Z, X= X)
data1 = data.frame(data)

#regression ignoring U
Z.fit = glm(Z ~ X, family = binomial(link="probit"))
Coefz.noU = Z.fit$coefficients
t1.fit = coxph(Surv(t, d) ~ X + Z)
Coeft1.noU = t1.fit$coefficients

#ipw
t1.ipw = coxph(Surv(t, d) ~ Z, weights = ipw, robust = TRUE)

tau1.sto = data.frame()
tau1.em = data.frame()
tau1.ipw = data.frame()
zetaZ = seq(-2,2,by=0.5)
zetaT = seq(-2,2,by=0.5)

set.seed(2019)
for(i in 1:length(zetaZ)){
  for(j in 1:length(zetaT)){
    #stochastive EM/regression
    tau.sto = surv_sens_var(data, zetat = zetaT[j], zetaz = zetaZ[i], theta = 0.3, B=100)
    tau1.sto = rbind(tau1.sto, data.frame(zetaz = zetaZ[i], zetat = zetaT[j], tau1 = tau.sto$tau1, tau1.se = tau.sto$tau1.se))
    
    #stochastive EM/IPW
    tau.ipw = surv_sens_ipw(data, zetat = zetaT[j], zetaz = zetaZ[i], theta = 0.3, B=100)
    tau1.ipw = rbind(tau1.ipw, data.frame(zetaz = zetaZ[i], zetat = zetaT[j], tau1 = tau.ipw$tau1, tau1.se = tau.ipw$tau1.se))
    
    #EM/regression
    tau.em = emU_surv(t, d, Z, X, zetat = zetaT[j], zetaz = zetaZ[i], theta = 0.3)
    data1$p = tau.em$p
    tau.em.final = surv_final(data1, zetat = zetaT[j], zetaz = zetaZ[i], z.coef = tau.em$z.coef[1:2])
    tau1.em = rbind(tau1.em, data.frame(zetaz = zetaZ[i], zetat = zetaT[j],tau1 = tau.em.final$coef[2], tau1.se = tau.em.final$coef.se[2]))
  }
}

tau1.sto$t = tau1.sto$tau1/tau1.sto$tau1.se
tau1.em$t = tau1.em$tau1/tau1.em$tau1.se
tau1.ipw$t = tau1.ipw$tau1/tau1.ipw$tau1.se

pdf("out1_stoEM.pdf", width=6, height=5)
g = ggplot(tau1.sto, aes(zetaz, zetat)) +
  stat_contour(aes(z = tau1, colour = ..level..)) +
  stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
  xlab("Coef. on U in model for treatment") +
  ylab("Coef. on U in model for response") +
  theme_bw() +
  annotate("text",x=0,y=0, label = "0.5756")
direct.label(g, method="bottom.pieces")
dev.off()

pdf("out1_EM.pdf", width=6, height=5)
g = ggplot(tau1.em, aes(zetaz, zetat)) +
  stat_contour(aes(z = tau1, colour = ..level..)) +
  stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
  xlab("Coef. on U in model for treatment") +
  ylab("Coef. on U in model for response") +
  theme_bw() + 
  annotate("text",x=0,y=0, label = "0.5756")
direct.label(g, method="bottom.pieces")
dev.off()

pdf("out1_ipw.pdf", width=6, height=5)
g = ggplot(tau1.ipw, aes(zetaz, zetat)) +
  stat_contour(aes(z = tau1, colour = ..level..)) +
  stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
  xlab("Coef. on U in model for treatment") +
  ylab("Coef. on U in model for response") +
  theme_bw()+
  annotate("text",x=0,y=0, label = "0.5241")
direct.label(g, method="bottom.pieces")
dev.off()

tau1.sto$value = paste(round(tau1.sto$tau1,4), " (", round(tau1.sto$tau1.se,4), ")", sep="")
write.csv(tau1.sto, "out1_stoEM.csv")
tau1.em$value = paste(round(tau1.em$tau1,4), " (", round(tau1.em$tau1.se,4), ")", sep="")
write.csv(tau1.em, "out1_EM.csv")
tau1.ipw$value = paste(round(tau1.ipw$tau1,4), " (", round(tau1.ipw$tau1.se,4), ")", sep="")
write.csv(tau1.ipw, "out1_stoEM_ipw.csv")