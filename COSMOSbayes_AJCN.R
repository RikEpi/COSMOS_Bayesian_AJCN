#Project: Re-examining the COSMOS trial results using a Bayesian approach
#Made by Rikuta Hamaya

library(tidyverse)
library(haven)
library(survHE)
library(rstan)
library(bayesplot)
library(tidybayes)
library(scales)


# priors ------------------------------------------------------------------
#These calculations are based on previous literatures
#Priors were defined by three mechanisms: CF-BP-CVD, CF-FMD-CVD, and CF-CRP-CVD
#Please see the method section for each numerics

#BP-CVD
HR = 0.91
L = 0.88
U = 0.95
logHR = log(HR)
logL = log(L)
logU = log(U)
SE = ((logHR - logL)/1.96 + (logU - logHR)/1.96)/2

SEc = ((2.37 - 0.44)/1.96 + (4.3 - 2.37)/1.96)/2
BP1 = rnorm(100000,2.37,SEc)
BP2 = rnorm(100000,logHR,SE)
BPall=BP1*BP2/5

Est_BP=mean(BPall)
SE_BP=sd(BPall)

exp(Est_BP)
exp(Est_BP-1.96*SE_BP)
exp(Est_BP+1.96*SE_BP)


#FMD-CVD
HR = 0.88
L = 0.84
U = 0.91
logHR = log(HR)
logL = log(L)
logU = log(U)
SE = ((logHR - logL)/1.96 + (logU - logHR)/1.96)/2

SEc = ((1.34 - 1)/1.96 + (1.68 - 1.34)/1.96)/2
FMD1 = rnorm(100000,1.34,SEc)
FMD2 = rnorm(100000,logHR,SE)
FMDall=FMD1*FMD2

Est_FMD=mean(FMDall)
SE_FMD=sd(FMDall)

exp(Est_FMD)
exp(Est_FMD-1.96*SE_FMD)
exp(Est_FMD+1.96*SE_FMD)

##CRP-CVD
##300mg arm: per mg/l
HR = 0.86
L = 0.75
U = 0.99
logHR = log(HR)
logL = log(L)
logU = log(U)
SE = ((logHR - logL)/1.96 + (logU - logHR)/1.96)/2

SEc = ((0.978 - 0.269)/1.96 + (1.687 - 0.978)/1.96)/2
CRP1 = rnorm(100000,0.978,SEc)
CRP2 = rnorm(100000,logHR,SE)
CRPall=CRP1*CRP2/1.85

Est_CRP=mean(CRPall)
SE_CRP=sd(CRPall)

exp(Est_CRP)
exp(Est_CRP+1.96*SE_CRP)
exp(Est_CRP-1.96*SE_CRP)



##LDL-CVD
HR = 0.77
L = 0.75
U = 0.79
logHR = log(HR)
logL = log(L)
logU = log(U)
SE = ((logHR - logL)/1.96 + (logU - logHR)/1.96)/2

SEc = (11.13 - 0.61)/3.92
LDL1 = rnorm(100000,5.87,SEc)
LDL2 = rnorm(100000,logHR,SE)
LDLall=LDL1*LDL2/38.67

Est_LDL=mean(LDLall)
SE_LDL=sd(LDLall)

exp(Est_LDL)
exp(Est_LDL+1.96*SE_LDL)
exp(Est_LDL-1.96*SE_LDL)

#combine all three pathways: 14 senarios
CRP = rnorm(100000,Est_CRP,SE_CRP)
FMD = rnorm(100000,Est_FMD,SE_FMD)
BP = rnorm(100000,Est_BP,SE_BP)
LDL = rnorm(100000,Est_LDL,SE_LDL)
CRP2 = rnorm(100000,Est_CRP,SE_CRP*2)
FMD2 = rnorm(100000,Est_FMD,SE_FMD*2)
LDL2 = rnorm(100000,Est_LDL,SE_LDL*2)

#original scenarios 1 and 2 were excluded for the too strong assumptions
#scenarios in papers coresponds to [scenario - 2] in this code (e.g. scenario 1 id based on dis3)
# dis1=CRP+FMD+BP
# dis2=CRP2+FMD2+BP2
dis3=FMD
dis4=rnorm(100000,Est_FMD,SE_FMD*2)
dis5=CRP+BP
dis6=CRP2+BP2
dis7=CRP*0.75+BP*0.75
dis8=CRP2*0.75+BP2*0.75
dis9=CRP*0.75+FMD2*0.75+BP*0.75
dis10=CRP2*0.75+FMD2*0.75+BP2*0.75
dis11=CRP+FMD2*0.5+BP*0.5
dis12=CRP2+FMD2*0.5+BP2*0.5
dis13=CRP*0.5+FMD2*0.5+BP*0.5
dis14=CRP2*0.5+FMD2*0.5+BP2*0.5
dis15=CRP*0.5+FMD2*0.5+BP*0.5+LDL2*0.5
dis16=CRP2*0.5+FMD2*0.5+BP2*0.5+LDL2*0.5

est1=mean(dis1)
se1=sd(dis1)
est2=mean(dis2)
se2=sd(dis2)
est3=mean(dis3)
se3=sd(dis3)
est4=mean(dis4)
se4=sd(dis4)
est5=mean(dis5)
se5=sd(dis5)
est6=mean(dis6)
se6=sd(dis6)
est7=mean(dis7)
se7=sd(dis7)
est8=mean(dis8)
se8=sd(dis8)
est9=mean(dis9)
se9=sd(dis9)
est10=mean(dis10)
se10=sd(dis10)
est11=mean(dis11)
se11=sd(dis11)
est12=mean(dis12)
se12=sd(dis12)
est13=mean(dis13)
se13=sd(dis13)
est14=mean(dis14)
se14=sd(dis14)
est15=mean(dis15)
se15=sd(dis15)
est16=mean(dis16)
se16=sd(dis16)

#Table 1: Output of prior distributions
PriorHR=function(dis){
  paste0(round(exp(mean(dis)),2)," (",
         round(exp(quantile(dis,0.025)),2),", ",
         round(exp(quantile(dis,0.975)),2),")")
}
PriorHR(dis1)
PriorHR(dis2)
PriorHR(dis3)
PriorHR(dis4)
PriorHR(dis5)
PriorHR(dis6)
PriorHR(dis7)
PriorHR(dis8)
PriorHR(dis9)
PriorHR(dis10)
PriorHR(dis11)
PriorHR(dis12)
PriorHR(dis13)
PriorHR(dis14)
PriorHR(dis15)
PriorHR(dis16)

#SFig 2
dat=data.frame(A=exp(c(dis3,dis4,dis5,dis6,dis7,dis8,dis9,dis10,dis11,dis12,dis13,dis14,dis15,dis16)),
               name=c(rep("Senario 1",1000000),
                      rep("Senario 2",1000000),
                      rep("Senario 3",1000000),
                      rep("Senario 4",1000000),
                      rep("Senario 5",1000000),
                      rep("Senario 6",1000000),
                      rep("Senario 7",1000000),
                      rep("Senario 8",1000000),
                      rep("Senario 9",1000000),
                      rep("Senario 10",1000000),
                      rep("Senario 11",1000000),
                      rep("Senario 12",1000000),
                      rep("Senario 13",1000000),
                      rep("Senario 14",1000000)))
ggplot(dat, aes(x=A,color=name)) + geom_density(size=1)+ 
  scale_x_continuous(trans = log_trans(),limits=c(0.5,1.5)) +
  ylab("Prior density") + xlab("Hazard ratio")



# Read original data ------------------------------------------------------
#define main_use dataframe


# Frequentist analysis ----------------------------------------------------
#Cox
coxph(Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi, data = main_use)%>% summary()
coxph(Surv(t2mistrrevcvddth,mistrrevcvddth)~interventionC, data = main) %>% summary()

#Weibull models
survreg(Surv(t2mistrrevcvddthoth, mistrrevcvddthoth) ~ interventionC+agerand+men+interventionM+whi, dist="weibull", data=main_use) %>% summary()



# Bayesian analysis -------------------------------------------------------
#CVD results
CVD0 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                  data = main_use, distr =c("weibullPH"),method = "hmc",
                  priors = list(
                    exp = list(mu_beta = c(0,0,0,0,0,0), sigma_beta= c(10,10,10,10,10,10)),
                    weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                     mu_beta = c(0,0,0,0,0,0), sigma_beta= c(10,10,10,10,10,10) )),
                  chains=3,iter=4000,seed=1233)

# original scenarios 1 and 2 are not used
# CVD1 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
#                      data = main_use, distr =c("weibullPH"),method = "hmc",
#                      priors = list(
#                        exp = list(mu_beta = c(0,est1,0,0,0,0), sigma_beta= c(10,se1,10,10,10,10)),
#                        weibullPH = list(a_alpha=0.1, b_alpha=0.1,
#                                         mu_beta = c(0,est1,0,0,0,0), sigma_beta= c(10,se1,10,10,10,10) )),
#                      chains=3,iter=4000,seed=1233)
# 
# CVD2 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
#                   data = main_use, distr =c("weibullPH"),method = "hmc",
#                   priors = list(
#                     exp = list(mu_beta = c(0,est2,0,0,0,0), sigma_beta= c(10,se2,10,10,10,10)),
#                     weibullPH = list(a_alpha=0.1, b_alpha=0.1,
#                                      mu_beta = c(0,est2,0,0,0,0), sigma_beta= c(10,se2,10,10,10,10) )),
#                   chains=3,iter=4000,seed=1233)

CVD3 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                  data = main_use, distr =c("weibullPH"),method = "hmc",
                  priors = list(
                    exp = list(mu_beta = c(0,est3,0,0,0,0), sigma_beta= c(10,se3,10,10,10,10)),
                    weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                     mu_beta = c(0,est3,0,0,0,0), sigma_beta= c(10,se3,10,10,10,10) )),
                  chains=3,iter=4000,seed=1233)

CVD4 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                  data = main_use, distr =c("weibullPH"),method = "hmc",
                  priors = list(
                    exp = list(mu_beta = c(0,est4,0,0,0,0), sigma_beta= c(10,se4,10,10,10,10)),
                    weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                     mu_beta = c(0,est4,0,0,0,0), sigma_beta= c(10,se4,10,10,10,10) )),
                  chains=3,iter=4000,seed=1233)

CVD5 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                  data = main_use, distr =c("weibullPH"),method = "hmc",
                  priors = list(
                    exp = list(mu_beta = c(0,est5,0,0,0,0), sigma_beta= c(10,se5,10,10,10,10)),
                    weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                     mu_beta = c(0,est5,0,0,0,0), sigma_beta= c(10,se5,10,10,10,10) )),
                  chains=3,iter=4000,seed=1233)

CVD6 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                  data = main_use, distr =c("weibullPH"),method = "hmc",
                  priors = list(
                    exp = list(mu_beta = c(0,est6,0,0,0,0), sigma_beta= c(10,se6,10,10,10,10)),
                    weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                     mu_beta = c(0,est6,0,0,0,0), sigma_beta= c(10,se6,10,10,10,10) )),
                  chains=3,iter=4000,seed=1233)

CVD7 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                  data = main_use, distr =c("weibullPH"),method = "hmc",
                  priors = list(
                    exp = list(mu_beta = c(0,est7,0,0,0,0), sigma_beta= c(10,se7,10,10,10,10)),
                    weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                     mu_beta = c(0,est7,0,0,0,0), sigma_beta= c(10,se7,10,10,10,10) )),
                  chains=3,iter=4000,seed=1233)

CVD8 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                  data = main_use, distr =c("weibullPH"),method = "hmc",
                  priors = list(
                    exp = list(mu_beta = c(0,est8,0,0,0,0), sigma_beta= c(10,se8,10,10,10,10)),
                    weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                     mu_beta = c(0,est8,0,0,0,0), sigma_beta= c(10,se8,10,10,10,10) )),
                  chains=3,iter=4000,seed=1233)

CVD9 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                  data = main_use, distr =c("weibullPH"),method = "hmc",
                  priors = list(
                    exp = list(mu_beta = c(0,est9,0,0,0,0), sigma_beta= c(10,se9,10,10,10,10)),
                    weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                     mu_beta = c(0,est9,0,0,0,0), sigma_beta= c(10,se9,10,10,10,10) )),
                  chains=3,iter=4000,seed=1233)

CVD10 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                   data = main_use, distr =c("weibullPH"),method = "hmc",
                   priors = list(
                     exp = list(mu_beta = c(0,est10,0,0,0,0), sigma_beta= c(10,se10,10,10,10,10)),
                     weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                      mu_beta = c(0,est10,0,0,0,0), sigma_beta= c(10,se10,10,10,10,10) )),
                   chains=3,iter=4000,seed=1233)

CVD11 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                   data = main_use, distr =c("weibullPH"),method = "hmc",
                   priors = list(
                     weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                      mu_beta = c(0,est11,0,0,0,0), sigma_beta= c(10,se11,10,10,10,10) )),
                   chains=3,iter=4000,seed=1233)

CVD12 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                   data = main_use, distr =c("weibullPH"),method = "hmc",
                   priors = list(
                     weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                      mu_beta = c(0,est12,0,0,0,0), sigma_beta= c(10,se12,10,10,10,10) )),
                   chains=3,iter=4000,seed=1233)

CVD13 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                   data = main_use, distr =c("weibullPH"),method = "hmc",
                   priors = list(
                     exp = list(mu_beta = c(0,est13,0,0,0,0), sigma_beta= c(10,se13,10,10,10,10)),
                     weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                      mu_beta = c(0,est13,0,0,0,0), sigma_beta= c(10,se13,10,10,10,10) )),
                   chains=3,iter=4000,seed=1233)

CVD14 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                   data = main_use, distr =c("weibullPH"),method = "hmc",
                   priors = list(
                     exp = list(mu_beta = c(0,est14,0,0,0,0), sigma_beta= c(10,se14,10,10,10,10)),
                     weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                      mu_beta = c(0,est14,0,0,0,0), sigma_beta= c(10,se14,10,10,10,10) )),
                   chains=3,iter=4000,seed=1233)

CVD15 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                   data = main_use, distr =c("weibullPH"),method = "hmc",
                   priors = list(
                     weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                      mu_beta = c(0,est15,0,0,0,0), sigma_beta= c(10,se15,10,10,10,10) )),
                   chains=3,iter=4000,seed=1233)

CVD16 = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                   data = main_use, distr =c("weibullPH"),method = "hmc",
                   priors = list(
                     weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                      mu_beta = c(0,est16,0,0,0,0), sigma_beta= c(10,se16,10,10,10,10) )),
                   chains=3,iter=4000,seed=1233)


CVD<-grep("CVD",names(.GlobalEnv),value=TRUE)
CVD_list<-do.call("list",mget(CVD))






# Main Results: Table 3 -----------------------------------------------------------------

#results output function
comb = function(res){
  temp = res$models$`Weibull (PH)`
  n4 = temp@.MISC$summary$quan[2,5]
  n5 = temp@.MISC$summary$quan[2,1]
  n6 = temp@.MISC$summary$quan[2,9]
  
  print(paste0('Weibull: ',
               round(exp(n4),2)," (",
               round(exp(n5),2),", ",
               round(exp(n6),2),")"))
}

#Results of fits
for (d in 1:17) {
  d2=CVD_list[[d]]
  print(d2$models$`Weibull (PH)`)
}

#Table 2: Results of estimates::: not use CVD1 and CVD2
for (i in 1:17){
  print(names(CVD_list)[[i]])
  comb(CVD_list[[i]])
}



# risk difference: Table 3 (response to reviewer's comment) ---------------------------------------------------------

t_eval <- 365*4

RDcalc = function(cvdfit){
  temp=cvdfit$models$`Weibull (PH)`
  post <- rstan::extract(temp)
  alpha <- post$alpha
  beta <- post$beta
  
  # Set covariates
  X_control <- c(1, 0, mean(main_use$agerand), mean(main_use$men), 0, 0)  # intercept + covariates
  X_treat   <- c(1, 1, mean(main_use$agerand), mean(main_use$men), 0, 0)
  
  # Calculate linear predictors
  lp_control <- as.vector(beta %*% X_control)
  lp_treat   <-  as.vector(beta %*% X_treat)
  
  S_control <- exp(- (t_eval ^ alpha) * exp(lp_control))
  S_treat   <- exp(- (t_eval ^ alpha) * exp(lp_treat))
  
  risk_control <- 1 - S_control
  risk_treat <- 1 - S_treat
  RD <- risk_treat - risk_control
  
  # Summary
  return(c(mean(RD),quantile(RD, c(0.025, 0.975)) ))
}


model_list <- list(CVD0, CVD1, CVD2, CVD3, CVD4,
                   CVD5, CVD6, CVD7, CVD8, CVD9,
                   CVD10, CVD11, CVD12, CVD13, CVD14, CVD15, CVD16)

# Apply RDcalc to each model
rd_results <- lapply(model_list, RDcalc)

# Convert to a data frame
rd_df <- do.call(rbind, rd_results)%>% as.data.frame()
colnames(rd_df) <- c("RD_mean", "RD_lower", "RD_upper") 

rd_df$est = paste0(round(rd_df$RD_mean*100,2)," (",round(rd_df$RD_lower*100,2),", ",round(rd_df$RD_upper*100,2),")")
write.csv(rd_df,"rd_df.csv")



# Posterior likelihood: Figure 2 ---------------------------------------------------------

#prior likelihood
set.seed(123)
dis0=rnorm(100000,0,10)
priors<-grep("dis",names(.GlobalEnv),value=TRUE)
prior_list<-do.call("list",mget(priors))


pl=replicate(7,rep(NA,17)) %>% data.frame()

for (i in 1:17){
  pl[i,1]=names(prior_list)[[i]]
  AA=prior_list[[i]]
  vec=AA %>% exp() 
  pl[i,2:4] = c(mean(vec<1),mean(vec<0.9),mean(vec<0.8))*100 
}

#posterior likelihood
desired_order <- paste0("CVD", 0:16)
CVD_list <- CVD_list[desired_order]

for (i in 1:17){
  AA=CVD_list[[i]]
  vec=as.matrix(AA$models$`Weibull (PH)`)[,2] %>% exp() 
  pl[i,5:7] =  c(mean(vec<1),mean(vec<0.9),mean(vec<0.8),mean(vec<0.7))*100 
}

colnames(pl)=c("olds","pre1","pre0.9","pre0.8","post1","post0.9","post0.8")

plnew = pl[c(1,4:17),]
plnew$scenario = 0:14


plnew_long2 <- plnew %>%
  pivot_longer(
    cols = c("pre1", "pre0.9", "pre0.8", "post1", "post0.9", "post0.8"),
    names_to = "prob_type",
    values_to = "prob_val"
  ) %>%
  mutate(
    # Classify as prior/posterior
    Group = ifelse(grepl("^pre", prob_type), "A: Prior probability", "B: Posterior probability"),
    # Classify threshold
    Threshold = case_when(
      prob_type %in% c("pre1", "post1") ~ "HR<1.0",
      prob_type %in% c("pre0.9", "post0.9") ~ "HR<0.9",
      TRUE ~ "HR<0.8"
    )
  )

plnew_long2$scenario = factor(plnew_long2$scenario, 
                              levels = sort(unique(plnew_long2$scenario),decreasing = TRUE))
plnew_long2$Threshold = factor(plnew_long2$Threshold,
                               levels = c("HR<1.0", "HR<0.9", "HR<0.8"))

ggplot(plnew_long2, aes(x = Threshold, y = scenario, fill = prob_val)) +
  geom_tile() +
  geom_text(aes(label = ifelse(prob_val >99.9, ">99.9%", sprintf("%.1f", prob_val))), size = 3) +
  facet_wrap(~ Group) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = NULL, y = "Scenario", fill = "Probability") +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  )



# K-M: Figure 3 ---------------------------------------------------------------------

cvd9km = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~as.factor(interventionC),
                    data = main_use, distr =c("weibullPH"),method = "hmc",
                    priors = list(
                      weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                       mu_beta = c(0,est9), sigma_beta= c(10,se9) )),
                    chains=3,iter=4000,seed=1233)


plot(cvd9km,add.km=F) + theme_classic() +ylim(c(0.90,1))+xlab("Follow-up (days)")

main_use$interventionC=as.factor(main_use$interventionC)
km = survfit(Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC, data = main_use)
autoplot(km, conf.int = FALSE, censor = FALSE) +  # Remove confidence intervals
  ylim(0.90, 1) +                    # Set y-axis limits
  xlab("Follow-up (days)") +         # Label for x-axis
  ylab("Survival") +                  # Label for y-axis
  theme_classic() +                  # Apply a classic theme
  scale_color_manual(values = c("blue", "blue")) +      # Set both lines to red
  scale_linetype_manual(values = c("dashed", "solid")) + # Set line types
  theme(legend.title = element_blank())   

ggsurvplot(
  km,
  data = main_use,
  linetype = "interventionC",
  palette = c("blue", "blue"),
  censor = FALSE,
  ylim = c(0.90, 1),
  xlab = "Follow-up (days)",
  ylab = "Survival",
  legend.title = "InterventionC",
  legend.labs = c("Group 0", "Group 1"),
  ggtheme = theme_minimal(),
  risk.table = FALSE
) +
  theme(
    legend.position = "bottom",         # Position the legend at the bottom
    legend.key = element_blank(),       # Remove the legend key background
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 10)    # Increase legend text size
  )



km_plot <- ggsurvplot(
  km,
  data = main_use,
  linetype = "interventionC",
  palette = c("blue", "blue"),
  censor = FALSE,
  ylim = c(0.90, 1),
  xlab = "Follow-up (days)",
  ylab = "Survival",
  legend.title = "InterventionC",
  legend.labs = c("Group 0", "Group 1"),
  ggtheme = theme_classic(),
  risk.table = FALSE
)

# Reverse the line types (solid <-> dashed)
km_plot$plot <- km_plot$plot + 
  scale_linetype_manual(values = c("dashed", "solid"))

km_plot$plot <- km_plot$plot +
  scale_x_continuous(breaks = c(500, 1000, 1500))

# Adjust the legend and other theme elements for better clarity
km_plot$plot <- km_plot$plot + 
  theme(
    legend.position = "bottom",         # Position the legend at the bottom
    legend.key = element_blank(),       # Remove the legend key background
    legend.title = element_text(size = 12),  # Increase legend title size
    legend.text = element_text(size = 10)    # Increase legend text size
  )



# diagnosis: SFig 3 ---------------------------------------------------------------
#only representative case is shown
rstan::traceplot(CVD9$models$`Weibull (PH)`)
bayesplot::mcmc_acf(as.matrix(CVD9$models$`Weibull (PH)`))





# sensitivity analysis responding to reviewer's comment -------------------

#skeptical prior [estimate of HR was 1, and SEs were defined such that a 
#probability of observing the treatment effect shown in scenario 7 less than 5%]). 
est_skep = 0
skeptical_sd = function(est){
  for (i in seq(0.001,1,by=0.001)){
    if (pnorm(est9,mean=est,sd=i) > 0.05){
      return(i)
      break}
  }
}
se_skep = skeptical_sd(0)

skep_res = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                      data = main_use, distr =c("weibullPH"),method = "hmc",
                      priors = list(
                        weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                         mu_beta = c(0,est_skep,0,0,0,0), sigma_beta= c(10,se_skep,10,10,10,10) )),
                      chains=3,iter=4000,seed=1233)


#pessimistic  prior [estimate of HR was 1.1, and SEs were defined such that a 
#probability of observing the treatment effect shown in scenario 7 less than 5%]). 
est_neg = log(1.1)
se_neg = skeptical_sd(est_neg)

neg_res = fit.models(formula = Surv(t2mistrrevcvddthoth,mistrrevcvddthoth)~interventionC+agerand+men+interventionM+whi,
                     data = main_use, distr =c("weibullPH"),method = "hmc",
                     priors = list(
                       weibullPH = list(a_alpha=0.1, b_alpha=0.1,
                                        mu_beta = c(0,est_neg,0,0,0,0), sigma_beta= c(10,se_neg,10,10,10,10) )),
                     chains=3,iter=4000,seed=1233)


#results of sensitivity analysis
temp=skep_res$models$`Weibull (PH)`
exp(temp@.MISC$summary$quan[2,c(1,5,9)])

vec=as.matrix(temp)[,2] %>% exp() 
c(mean(vec<1),mean(vec<0.9),mean(vec<0.8),mean(vec<0.7))*100 

