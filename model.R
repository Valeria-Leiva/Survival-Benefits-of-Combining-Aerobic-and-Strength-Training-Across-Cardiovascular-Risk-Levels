rm(list=ls())
library(haven)

metadata <- read_sas('....../aeandre2.sas7bdat')
metadata <- as.data.frame(metadata)
names(metadata)

metadata$gender <- factor(metadata$female, levels = c(0,1), labels = c('Hombre','Mujer'))
metadata$FollowMortalYr <- ifelse(metadata$FollowMortalYr==0, 
                                  metadata$FollowMortalYr + 0.001,
                                  metadata$FollowMortalYr)
metadata$AGE <- as.numeric(metadata$AGE)
metadata$CHOLSTRL <- as.numeric(metadata$CHOLSTRL)
metadata$HDL <- as.numeric(metadata$HDL)
metadata$RSTSYSBP <- as.numeric(metadata$RSTSYSBP)


metadata$treated <- ifelse(metadata$abnHTN == 1, 1, 0)

calculate_framingham <- function(age, sex, chol, hdl, sbp, smoker, diabetes, treated) {
  if (sex == 1) { # Mujeres
    age_coef <- 2.32888
    chol_coef <- 1.20904
    hdl_coef <- -0.70833
    sbp_coef <- ifelse(treated == 1, 2.82263, 2.76157)
    smoker_coef <- 0.52873
    diabetes_coef <- 0.69154
    s10 <- 0.95012
    mean_lp <- 26.1931
  } else { # Hombres
    age_coef <- 3.06117
    chol_coef <- 1.12370
    hdl_coef <- -0.93263
    sbp_coef <- ifelse(treated == 1, 1.99881, 1.93303)
    smoker_coef <- 0.65451
    diabetes_coef <- 0.57367
    s10 <- 0.88936
    mean_lp <- 23.9802
  }
  
  lp <- age_coef * log(age) + 
    chol_coef * log(chol) + 
    hdl_coef * log(hdl) + 
    sbp_coef * log(sbp) + 
    smoker_coef * smoker + 
    diabetes_coef * diabetes
  
  risk <- 1 - (s10 ^ exp(lp - mean_lp))
  return(risk * 100)  # riesgo en %
}
metadata$fram_risk <- mapply(
  calculate_framingham,
  age = metadata$AGE,
  sex = metadata$female,
  chol = metadata$CHOLSTRL,
  hdl = metadata$HDL,
  sbp = metadata$RSTSYSBP,
  smoker = metadata$smokenow1,
  diabetes = metadata$abnDM,
  treated = metadata$treated
)
metadata$FS_cat_new <- cut(
  metadata$fram_risk,
  breaks = c(-Inf, 10, 20, Inf),
  labels = c("Low", "Intermediate", "High")
)

metadata$SportCat <- 0
metadata$SportCat[metadata$RecAE==1 & metadata$RecRE==0] <- 1
metadata$SportCat[metadata$RecAE==0 & metadata$RecRE==1] <- 2
metadata$SportCat[metadata$RecAE==1 & metadata$RecRE==1] <- 3
metadata$SportCat <- factor(metadata$SportCat , levels=c(0,1,2,3), 
                            labels=c('None','Aerobic', 'Strength','Both'))

metadata$SB_cat <- 0
metadata$SB_cat[metadata$RSTSYSBP<130 & metadata$RSTDIABP <85] <- 1
metadata$SB_cat[(metadata$RSTSYSBP>=130 & metadata$RSTSYSBP<140) | (metadata$RSTDIABP>=85 & metadata$RSTDIABP <90)] <- 2
metadata$SB_cat[(metadata$RSTSYSBP>=140 & metadata$RSTSYSBP<160) | (metadata$RSTDIABP>=90 & metadata$RSTDIABP <100)] <- 3
metadata$SB_cat[metadata$RSTSYSBP>=160 | metadata$RSTDIABP >=100] <- 4
metadata$SB_cat <- factor(metadata$SB_cat , levels=c(0,1,2,3,4), 
                          labels=c('NA','Normal', 'Prehypertension','Stage 1 hypertension','Stage 2 hypertension'))


mydata <- metadata[-which(metadata$SB_cat=='NA'),]
mydata$SB_cat <- droplevels(mydata$SB_cat)


mydata$SportCat_reduced <- mydata$SportCat
levels(mydata$SportCat_reduced)[levels(mydata$SportCat_reduced) %in% c("Strength", "Both")] <- "Strength_or_Both"

mydata2 <- mydata[mydata$SportCat!="Strength",]
mydata2$SportCat <- droplevels(mydata2$SportCat)

mymetada <- mydata2[mydata2$female==0 & (mydata2$AGE>=30 & mydata2$AGE<=74), 
                    c('FollowMortalYr','deceased','BMI', 'AGE', 'TRIG', 'abnEcg', 'abnDM',
                      'FS_cat_new', 'SportCat_reduced', 'SB_cat','smokenow1','SportCat')]
mymetada <- na.omit(mymetada)

K <- length(table(mymetada$SB_cat))# number of stratum (SB)
a <- seq(0,max(mymetada$FollowMortalYr)+0.001, length.out = K+1)

# int.obs: vector that tells us at which interval each observation is
int.obs <- matrix(data = NA, nrow = nrow(mymetada), ncol = length(a)-1)
d <- matrix(data = NA, nrow = nrow(mymetada), ncol = length(a)-1)
for(i in 1:nrow(mymetada)){
  for(k in 1:(length(a)-1)){
    d[i,k] <- ifelse(mymetada$FollowMortalYr[i]-a[k] > 0,1,0)*ifelse(a[k+1]-mymetada$FollowMortalYr[i] > 0,1,0)
    int.obs[i,k] <- d[i,k]*k
  }
}
int.obs <- rowSums(int.obs)

# X: design matrix
X <- model.matrix(~BMI + AGE+ TRIG +abnEcg +FS_cat_new * SportCat, 
                  data = mymetada)
dim(X)

library(rjags)
library(R2jags)

model.bayes <- 'model{
	for(i in 1:n){
		for(k in 1:int.obs[i]){
			cond[i,k] <- step(time[i]-a[k+1])
			HH[i,k] <- cond[i,k]*(a[k+1]-a[k])*lambda[k] +(1-cond[i,k])*(time[i]-a[k])*lambda[k]
		}
		H[i] <- sum(HH[i,1:int.obs[i]])
	}
	for(i in 1:n){
		elinpred[i] <- exp(inprod(beta[],X[i,]))
		logHaz[i] <- log(lambda[int.obs[i]]*elinpred[i])
		logSurv[i] <- -H[i]*elinpred[i]
		phi[i] <- 100000- delta[i]*logHaz[i] -logSurv[i]
	    zeros[i] ~ dpois(phi[i])
	}
	
  for(l in 1:Nbetas){
  	beta[l] ~ dnorm(0,0.001)
  }
  for(k in 1:m){
  	lambda[k] ~ dgamma(0.001,0.001)
  }
	
}'

d.jags <- list(n = nrow(mymetada), m = length(a)-1, delta = mymetada$deceased, time = mymetada$FollowMortalYr,
               X = X, a = a, int.obs = int.obs, Nbetas = ncol(X), zeros = rep(0,nrow(mymetada)))
p.jags <- c('beta','lambda')
out <- jags(data=d.jags,
            parameters=p.jags,
            model = textConnection(model.bayes),
            n.chains=2,
            n.iter = 10000,
            n.thin=10,
            n.burnin=5000)

print(out, intervals=c(0.025,0.975),digits = 4)
