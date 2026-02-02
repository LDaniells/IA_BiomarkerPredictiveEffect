rm(list=setdiff(ls(), c()))

#### IMPORT LIBRARIES ######################################################
library(dplyr)
library(mgcv)
library(detectseparation)   # to detect infinite estimates of the coeff in the glm model 
# (https://cran.r-project.org/web/packages/brglm2/brglm2.pdf)
library(doParallel)
library(foreach)
library(purrr)
library(doRNG) #for fixing seeds
library(tidyr)

Linear = function(data, markdiff,seed=-1) {
  
  if(seed > 0 )
    set.seed(seed)
  
  #check separation and convergence of models with given dataset
  
  separationglm <- glm(y ~ x + trt +  trt*x, data = data, family = binomial, 
                       method = "detect_separation")$outcome #Can the linear model be fit to the data? (can't have all responses on experimental <- seperation)
  
  convergence <- FALSE
  
  if(separationglm==0){ 
    
    convergence <- TRUE
    
    ###################################################
    # Fit a linear model
    
    model = glm(y ~ x +  trt +  trt*x, data = data, family = binomial)
    BIC_lin = BIC(model)
    
    # Prob that slope is negative (warning: parameter is for trt - control)
    #P_linear = pnorm(-markdiff, mean=-summary(model)$coefficients["x:trt1", "Estimate"], sd=summary(model)$coefficients["x:trt1", "Std. Error"])
    P_linear_inter = summary(model)$coefficients["x:trt", "Pr(>|z|)"] #Extracting the p-value of the interaction term (interaction test method)
    
    # Obtain predicted values and their confidence interval 
    predicted  = predict(model, se.fit = TRUE) #logit of the predictions of each observations
    fit =  predicted$fit
    se = predicted$se.fit
    lwr = fit - qt(1 - 0.025, df = model$df.residual) * se
    upr = fit + qt(1 - 0.025, df = model$df.residual) * se
    data_gam=cbind(data, fit, lwr, upr)
    data_gam = data_gam[order(data_gam$x), ]
    
    x = seq(min(data$x), max(data$x), length.out = 100) #sequence of biomarker values
    pred1 = predict(model, newdata = data.frame(x, trt = 0), se.fit = TRUE) #predict for values of biomarker under placebo
    pred2 = predict(model, newdata = data.frame(x, trt = 1), se.fit = TRUE)#predict for values of biomarker under experimental
    
    # Calculate the contrast and its confidence interval
    contrast = pred1$fit - pred2$fit #Placebo-Experimental (hence why we do <0)
    se.constrat = sqrt(pred1$se.fit^2 + pred2$se.fit^2) #standard error of the difference 
    lwr.contrast = contrast - qnorm(1 - 0.025) * se.constrat
    upr.contrast = contrast + qnorm(1 - 0.025) * se.constrat
    
    max_contrast = sample(1:length(x), 10000, replace=T) #Random samples of the biomarker
    min_contrast = sample(1:length(x), 10000, replace=T) 
    diff = contrast[pmax(min_contrast, max_contrast)] - contrast[pmin(min_contrast, max_contrast)]
    se.diff = sqrt(se.constrat[min_contrast]^2 + se.constrat[max_contrast]^2)
    # Prob that the Diff bw the contrasts is negative
    y = rnorm(10000, mean=diff, sd=se.diff)
    P_linear =  mean(y < -markdiff)
    
  }
  if(convergence){
    return(list("P(delta<-markdiff)_lin"=P_linear, #Probability Delta is positive (placebo-experimental)
                "P-value_inter" = P_linear_inter, #Interaction p-value
                "BIC_lin"= BIC_lin
                
    ))
  }else {
    return(list("P(delta<-markdiff)_lin"=NA,
                "P-value_inter" = NA,
                "BIC_lin"= NA
                
    ))
  }
  
  
}


# Build table of scenarios

names_scenarios <- c("NPNPT","NPNPNT", 
                     "HPNP1_50","HPNP1_40", "MPNP1_30", "LPNP1_20",
                     "HPNP2_50", "HPNP2_40", "MPNP2_30", "LPNP2_20",
                     "NPPT","NPPNT", 
                     "HPMP1_50","HPMP1_40", "MPMP1_30", "LPMP1_20",
                     "HPMP2_50", "HPMP2_40", "MPMP2_30", "LPMP2_20",
                     "HPHP1_50","HPHP1_40", "MPHP1_30", "LPHP1_20",
                     "HPHP2_50", "HPHP2_40", "MPHP2_30", "LPHP2_20")

#coefficients for logistic response-biomarker relationships for each scenario
coefficients_b0123 <- rbind(c(-0.405, 0.811, 0, 0),
                            c(-0.405, 0, 0, 0),
                            
                            c(-0.405, -1.428, 0.049, 0),
                            c(-0.405, -0.868, 0.035,0),
                            c(-0.405, -0.45, 0.026,0),
                            c(-0.405, -0.02, 0.017,0),
                            c(-0.405, -1.863, 0.045, 0),
                            c(-0.405, -1.36, 0.034,0),
                            c(-0.405, -0.853, 0.024,0),
                            c(-0.405, -0.465, 0.017,0),
                            
                            c(-1.655, 0.8, 0, 0.026),
                            c(-1.655, 0, 0, 0.026))
colnames(coefficients_b0123) <- c("b0", "b1", "b2", "b3")

coefficients_b0123 <- rbind(coefficients_b0123, coefficients_b0123[3:10, ])
coefficients_b0123[13:20, "b3"] = 0.2*coefficients_b0123[3:10, "b2"]

coefficients_b0123 <- rbind(coefficients_b0123, coefficients_b0123[3:10, ])
coefficients_b0123[21:28, "b3"] = 0.5*coefficients_b0123[3:10, "b2"]

distributionv <- as.character(c(rep("gamma", length(names_scenarios))))
param_distributionv <- c(rep(0.0827,length(names_scenarios))) #rate for gamma distribution and shape =  650.5434 *rate^2; max value range for uniform

#Data frame with all possible scenarios
df_all <- data.frame(distributionv, param_distributionv,
                     coefficients_b0123,
                     (rep(names_scenarios,1)))
colnames(df_all) <- c("distributionbmk", "param_distribution", "b0", "b1", "b2", "b3", "names_scenarios")

##select only NULL, Medium and High prognostic cases

df <- df_all %>% filter(names_scenarios %in%  c("NPNPT","NPNPNT", 
                                                "NPPT","NPPNT", 
                                                "HPNP1_50","HPNP1_40", "MPNP1_30", "LPNP1_20",
                                                "HPMP1_50","HPMP1_40", "MPMP1_30", "LPMP1_20",
                                                "HPHP1_50","HPHP1_40", "MPHP1_30", "LPHP1_20"))
# etaSTOP <- 0.1
# etaACC <- 0.3


run <- 5000

#Interim 
#Number of patients on experimental arm 
N.exp = 40
#Number of patients on placebo arm 
N.pbo = 20
#Total number of patients
Ntot = N.exp+N.pbo
for(j in 1:16){
  task_id <- j
  #Distribution of the biomarker
  distribution = df$distributionbmk[task_id]
  #Set parameters of the BMK-distribution
  if(distribution == "gamma"){
    rate = df$param_distribution[task_id]
    shape = 650.5434*rate^2
  }
  
  if(distribution == "uniform"){
    max = df$param_distribution[task_id]
    min = 0
  }
  
  #Coefficients of logistic function for the given scenario
  b0 = df$b0[task_id]
  b1 = df$b1[task_id]
  b2 = df$b2[task_id]
  b3 = df$b3[task_id]
  name_scen = df$names_scenarios[task_id]
  
  Decision.Store <- rep(NA,run)
  interim.prob <- rep(NA,run)
  final.prob <- rep(NA,run)
  for(i in 1:run){
    # Function for generating biomarker data, with a logistic relationship with the response
    SimDataLogistic = function(N, b0, b1, b2, b3, trt, distribution){ 
      
      if(distribution=="gamma"){
        x = matrix(rgamma(N, shape=shape, rate=rate), N, 1)
        p = exp(b0 + b1*(trt) + b2*x*(trt) + b3*x) / (1 + exp(b0 + b1*(trt) + b2*x*(trt) + b3*x)) 
        y = rbinom(N, 1, p)
        mydata  = cbind.data.frame(x,y,p)
        names(mydata) = c("x","y", "p")
        return(list("data"=mydata))
      }
      
      if(distribution=="uniform"){
        x = matrix(runif(N, min = min, max = max), N, 1)
        p = exp(b0 + b1*(trt) + b2*x*(trt) + b3*x) / (1 + exp(b0 + b1*(trt) + b2*x*(trt) + b3*x)) 
        y = rbinom(N, 1, p)
        mydata  = cbind.data.frame(x,y,p)
        names(mydata) = c("x","y", "p")
        return(list("data"=mydata))
      }
    }
    N.exp = 40
    #Number of patients on placebo arm 
    N.pbo = 20
    myData.exp = SimDataLogistic(N.exp, b0=b0, b1=b1, b2=b2,b3=b3,trt=1, distribution = distribution)
    myData.pbo = SimDataLogistic(N.pbo, b0=b0, b1=b1, b2=b2,b3=b3,trt=0, distribution = distribution)
    data=cbind(rbind(myData.exp$data, myData.pbo$data), "trt"=c(rep(1, length(myData.exp$data$x)), rep(0, length(myData.pbo$data$x))))
    
    #GaelleMK
    res_GaelleMK_lin <- Linear(data, 0) #remove seed 
    res_GaelleMK_lin_v <- unlist(res_GaelleMK_lin)
    output <- res_GaelleMK_lin_v
    
    interim.prob[i] <- as.numeric(output[1])
    
#Expansion Cohort
    N.exp = 15
    #Number of patients on placebo arm 
    N.pbo = 15
    #Total number of patients
    Ntot = N.exp+N.pbo
    myData.exp2 = SimDataLogistic(N.exp, b0=b0, b1=b1, b2=b2,b3=b3,trt=1, distribution = distribution)
    myData.pbo2 = SimDataLogistic(N.pbo, b0=b0, b1=b1, b2=b2,b3=b3,trt=0, distribution = distribution)
    data2=cbind(rbind(rbind(myData.exp$data,myData.exp2$data), rbind(myData.pbo$data,myData.pbo2$data)), "trt"=c(rep(1, length(myData.exp2$data$x)+length(myData.exp$data$x)), rep(0, length(myData.pbo2$data$x)+length(myData.pbo$data$x))))
    
    res_GaelleMK_lin <- Linear(data2, 0) #remove seed 
    res_GaelleMK_lin_v <- unlist(res_GaelleMK_lin)
    output <- res_GaelleMK_lin_v
    
    final.prob[i] <- as.numeric(output[1])
    
  }
  
  Results <- data.frame('Interim Prob'=interim.prob,'Final Prob'=final.prob)
  write.csv(Results, paste0("RawData", name_scen,"SS4020Interim15_15.csv"))
  
}


SS <- c(10,15,20,25,30)
Bounds.List <- list()
InterimStop <- c(0.80,0.75,0.70)
for(z in 1:3){
  Bounds <- matrix(NA,nrow=3,ncol=5)
  for(p in 1:5){
    #Interim------------
    filenamesdelta_all <- list.files(pattern="RawData*.", 
                                     full.names=TRUE)
    filenamesdelta <- filenamesdelta_all[grep(paste0("SS", 4020), filenamesdelta_all)]
    filenamesdelta <- filenamesdelta[grep(paste0("Interim", SS[p]), filenamesdelta)]
    
    Q <- delta <- c()
    for(k in 13:16){
      scenario <- read.csv(filenamesdelta[k])
      Q[k] <- quantile(scenario$Interim.Prob, prob=InterimStop[z],na.rm=T)
      delta[k] <- deltaL <- c(mean(scenario$Interim.Prob<Q[k], na.rm = TRUE))
    }
    eta0 <- max(Q,na.rm=T)
    rbind(Q,delta)
    eta0
    
    
    #Comment out if 2-bound configuration
    Q2 <- delta2 <- c()
    for(k in 13:16){
      scenario <- read.csv(filenamesdelta[k])
      Q2[k] <- quantile(scenario$Interim.Prob, prob=1-0.075,na.rm=T)
      delta2[k] <- deltaL <- c(mean(scenario$Interim.Prob>Q2[k], na.rm = TRUE))
    }
    eta1 <- max(Q2,na.rm=T)
    
    eta2 <- seq(0.65,0.8,0.0005)
    Reject <- Stop <- Acc <- Average.N <- Pos <-  matrix(NA,nrow=length(eta2),ncol=16)
    for(k in 13:16){
      scenario <- read.csv(filenamesdelta[k])
      
      for(i in 1:length(eta2)){
        Decision <- c()
        for(j in 1:5000){
          if(is.na(scenario$Interim.Prob[j])==T|is.na(scenario$Final.Prob[j])==T){
            Decision[j] <- 'NA'
          }else if(scenario$Interim.Prob[j]<eta0){
            Decision[j] <- 'STOP'
          }else if(scenario$Interim.Prob[j]>eta1){ #For two bound configuration set eta2==eta1
            Decision[j] <- 'ACC'
          }else if(scenario$Final.Prob[j]<=eta2[i]){
            Decision[j] <- 'Futile'
          }else{
            Decision[j] <- 'Efficacy'
          }
        }
        Reject[i,k] <- length(which((Decision=='ACC')|(Decision=='Efficacy')))/5000
        Stop[i,k] <- length(which(Decision=='STOP'))/5000
        Acc[i,k] <- length(which(Decision=='ACC'))/5000
        N <- c()
        for(t in 1:5000){
          Dec <- Decision[t]
          if(Dec=='STOP'|Dec=='ACC'){
            N[t] <- 20+40
          }else{
            N[t] <- 20+40+15+15
          }
        }
        Average.N[i,k] <- mean(N)
        Pos[i,k] <- length(which(Decision=='Efficacy'))/(length(which(Decision=='Efficacy'))+length(which(Decision=='Futile'))) 
      }
    }
    id <- which.min(abs(apply(Reject[,13:16],1,max)-0.15))
    eta2 <- eta2[id]  
    
    Bounds[1,p] <- eta0
    Bounds[2,p] <- eta1
    Bounds[3,p] <- eta2
  }
  Bounds.List[[z]] <- Bounds
}



SS <- c(10,15,20,25,30)
Results.80 <- Results.70 <- Results.75 <- list()
for(z in 1:3){
  for(p in 1:5){
    
    filenamesdelta_all <- list.files(pattern="RawData*.", 
                                     full.names=TRUE)
    filenamesdelta <- filenamesdelta_all[grep(paste0("Interim", SS[p]), filenamesdelta_all)]
    
    
    Reject <- Stop <- Acc <- Average.N <- Pos <- name <-   c()
    Q <- Bounds.List[[z]][1,p]
    eta1 <- Bounds.List[[z]][2,p]
    eta2 <- Bounds.List[[z]][3,p]
    
    for(k in 1:16){
      scenario <- read.csv(filenamesdelta[k])
      Decision <- c()
      for(j in 1:5000){
        if(is.na(scenario$Interim.Prob[j])==T|is.na(scenario$Final.Prob[j])==T){
          Decision[j] <- 'NA'
        }else if(scenario$Interim.Prob[j]<Q){
          Decision[j] <- 'STOP'
        }else if(scenario$Interim.Prob[j]>eta1){
          Decision[j] <- 'ACC'
        }else if(scenario$Final.Prob[j]<=eta2){
          Decision[j] <- 'Futile'
        }else{
          Decision[j] <- 'Efficacy'
        }
      }
      Reject[k] <- length(which((Decision=='ACC')|(Decision=='Efficacy')))/5000
      Stop[k] <- length(which(Decision=='STOP'))/5000
      Acc[k] <- length(which(Decision=='ACC'))/5000
      N <- c()
      for(t in 1:5000){
        Dec <- Decision[t]
        if(Dec=='STOP'|Dec=='ACC'){
          N[t] <- 20+40
        }else{
          N[t] <- 20+40+SS[p]+SS[p]
        }
      }
      Average.N[k] <- mean(N)
      Pos[k] <- length(which(Decision=='Efficacy'))/(length(which(Decision=='Efficacy'))+length(which(Decision=='Futile'))) 
      xx <- str_extract_part(filenamesdelta[k], before = TRUE, pattern = 'SS')
      scen_name_nsim_1 <- str_extract_part(xx, before = FALSE, pattern = "RawData")
      name[k] <- scen_name_nsim_1
    }
    if(z==1){
      FB <- rep(80,16)
    }else if(z==2){
      FB <- rep(75,16)
    }else{
      FB <- rep(70,16)
    }
    data <- data.frame('Scenario'=name,'SS'=rep(SS[p],16),'FutileBound'=FB,'Reject'=round(Reject,2),'Stop'=round(Stop,2),'Acc'=round(Acc,2),'N'=round(Average.N,2),'Pos'=round(Pos,2))
    
    
    names_scenarios <- c("NPNPT","NPPT","NPNPNT", "NPPNT",
                         "HPNP1_50","HPMP1_50",'HPHP1_50',
                         "HPNP1_40","HPMP1_40",'HPHP1_40',
                         "MPNP1_30","MPMP1_30",'MPHP1_30',
                         "LPNP1_20","LPMP1_20",'LPHP1_20')
    
    data <- data %>% arrange(factor(name,levels=names_scenarios))
    if(z==1){
      Results.80[[p]] <- data}else if(z==2){
        Results.75[[p]] <- data}else if(z==3){
          Results.70[[p]] <- data}
  }
}

