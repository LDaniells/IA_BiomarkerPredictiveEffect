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

b0 <- -0.405
b1 <- -1.428
b2 <- 0.049
b3 <- 0

distribution <- 'gamma'

rate = 0.0827
shape = 650.5434*rate^2

myData.exp = SimDataLogistic(N.exp, b0=b0, b1=b1, b2=b2,b3=b3,trt=1, distribution = distribution)
myData.pbo = SimDataLogistic(N.pbo, b0=b0, b1=b1, b2=b2,b3=b3,trt=0, distribution = distribution)
data=cbind(rbind(myData.exp$data, myData.pbo$data), "trt"=c(rep(1, length(myData.exp$data$x)), rep(0, length(myData.pbo$data$x))))



#Number of patients on placebo arm 
N.pbo = 10
#Total number of patients
myData.pbo2 = SimDataLogistic(N.pbo, b0=b0, b1=b1, b2=b2,b3=b3,trt=0, distribution = distribution)
data2=cbind(rbind(rbind(myData.exp$data), rbind(myData.pbo$data,myData.pbo2$data)), "trt"=c(rep(1, length(myData.exp$data$x)), rep(0, length(myData.pbo2$data$x)+length(myData.pbo$data$x))))


#Number of patients on placebo arm 
N.pbo = 15
#Total number of patients
myData.pbo2 = SimDataLogistic(N.pbo, b0=b0, b1=b1, b2=b2,b3=b3,trt=0, distribution = distribution)
data3=cbind(rbind(rbind(myData.exp$data), rbind(myData.pbo$data,myData.pbo2$data)), "trt"=c(rep(1, length(myData.exp$data$x)), rep(0, length(myData.pbo2$data$x)+length(myData.pbo$data$x))))


#Number of patients on placebo arm 
N.pbo = 20
#Total number of patients
myData.pbo2 = SimDataLogistic(N.pbo, b0=b0, b1=b1, b2=b2,b3=b3,trt=0, distribution = distribution)
data4=cbind(rbind(rbind(myData.exp$data), rbind(myData.pbo$data,myData.pbo2$data)), "trt"=c(rep(1, length(myData.exp$data$x)), rep(0, length(myData.pbo2$data$x)+length(myData.pbo$data$x))))

#Number of patients on placebo arm 
N.pbo = 30
#Total number of patients
myData.pbo2 = SimDataLogistic(N.pbo, b0=b0, b1=b1, b2=b2,b3=b3,trt=0, distribution = distribution)
data5=cbind(rbind(rbind(myData.exp$data), rbind(myData.pbo$data,myData.pbo2$data)), "trt"=c(rep(1, length(myData.exp$data$x)), rep(0, length(myData.pbo2$data$x)+length(myData.pbo$data$x))))



model1 = glm(y ~ x +  trt +  trt*x, data = data, family = binomial)
model2 = glm(y ~ x +  trt +  trt*x, data = data2, family = binomial)
model3 = glm(y ~ x +  trt +  trt*x, data = data3, family = binomial)
model4 = glm(y ~ x +  trt +  trt*x, data = data4, family = binomial)
model5 = glm(y ~ x +  trt +  trt*x, data = data5, family = binomial)


x = seq(min(data$x), max(data$x), length.out = 100) #sequence of biomarker values
pred1 = predict(model1, newdata = data.frame(x, trt = 0), se.fit = TRUE) #predict for values of biomarker under placebo
pred2 = predict(model1, newdata = data.frame(x, trt = 1), se.fit = TRUE)#predict for values of biomarker under experimental

# Calculate the contrast and its confidence interval
contrast = -pred1$fit + pred2$fit #Placebo-Experimental (hence why we do <0)
se.constrat = sqrt(pred1$se.fit^2 + pred2$se.fit^2) #standard error of the difference 
lwr.contrast = contrast - qnorm(1 - 0.025) * se.constrat
upr.contrast = contrast + qnorm(1 - 0.025) * se.constrat
data_contrast_1 <- data.frame(x,contrast,lwr.contrast,upr.contrast)


x2 = seq(min(data2$x), max(data2$x), length.out = 100) #sequence of biomarker values
pred1 = predict(model2, newdata = data.frame(x2, trt = 0), se.fit = TRUE) #predict for values of biomarker under placebo
pred2 = predict(model2, newdata = data.frame(x2, trt = 1), se.fit = TRUE)#predict for values of biomarker under experimental

# Calculate the contrast and its confidence interval
contrast2 = -pred1$fit + pred2$fit #Placebo-Experimental (hence why we do <0)
se.constrat2 = sqrt(pred1$se.fit^2 + pred2$se.fit^2) #standard error of the difference 
lwr.contrast2 = contrast2 - qnorm(1 - 0.025) * se.constrat2
upr.contrast2 = contrast2 + qnorm(1 - 0.025) * se.constrat2
data_contrast_2 <- data.frame(x,contrast2,lwr.contrast2,upr.contrast2)

max_contrast = sample(1:length(x2), 10000, replace=T) #Random samples of the biomarker
min_contrast = sample(1:length(x2), 10000, replace=T) 
diff = contrast2[pmax(min_contrast, max_contrast)] - contrast2[pmin(min_contrast, max_contrast)]
se.diff = sqrt(se.constrat2[min_contrast]^2 + se.constrat2[max_contrast]^2)
# Prob that the Diff bw the contrasts is negative
y = rnorm(10000, mean=diff, sd=se.diff)
P_linear2 =  mean(y > 0)

x3 = seq(min(data3$x), max(data3$x), length.out = 100) #sequence of biomarker values
pred1 = predict(model3, newdata = data.frame(x3, trt = 0), se.fit = TRUE) #predict for values of biomarker under placebo
pred2 = predict(model3, newdata = data.frame(x3, trt = 1), se.fit = TRUE)#predict for values of biomarker under experimental

# Calculate the contrast and its confidence interval
contrast3 = -pred1$fit + pred2$fit #Placebo-Experimental (hence why we do <0)
se.constrat3 = sqrt(pred1$se.fit^2 + pred2$se.fit^2) #standard error of the difference 
lwr.contrast3 = contrast3 - qnorm(1 - 0.025) * se.constrat3
upr.contrast3 = contrast3 + qnorm(1 - 0.025) * se.constrat3
data_contrast_3 <- data.frame(x,contrast3,lwr.contrast3,upr.contrast3)

max_contrast = sample(1:length(x3), 10000, replace=T) #Random samples of the biomarker
min_contrast = sample(1:length(x3), 10000, replace=T) 
diff = contrast3[pmax(min_contrast, max_contrast)] - contrast3[pmin(min_contrast, max_contrast)]
se.diff = sqrt(se.constrat3[min_contrast]^2 + se.constrat3[max_contrast]^2)
# Prob that the Diff bw the contrasts is negative
y = rnorm(10000, mean=diff, sd=se.diff)
P_linear3 =  mean(y > 0)


x4 = seq(min(data4$x), max(data4$x), length.out = 100) #sequence of biomarker values
pred1 = predict(model4, newdata = data.frame(x4, trt = 0), se.fit = TRUE) #predict for values of biomarker under placebo
pred2 = predict(model4, newdata = data.frame(x4, trt = 1), se.fit = TRUE)#predict for values of biomarker under experimental

# Calculate the contrast and its confidence interval
contrast4 = -pred1$fit + pred2$fit #Placebo-Experimental (hence why we do <0)
se.constrat4 = sqrt(pred1$se.fit^2 + pred2$se.fit^2) #standard error of the difference 
lwr.contrast4 = contrast4 - qnorm(1 - 0.025) * se.constrat4
upr.contrast4 = contrast4 + qnorm(1 - 0.025) * se.constrat4
data_contrast_4 <- data.frame(x,contrast4,lwr.contrast4,upr.contrast4)


x5 = seq(min(data5$x), max(data5$x), length.out = 100) #sequence of biomarker values
pred1 = predict(model5, newdata = data.frame(x5, trt = 0), se.fit = TRUE) #predict for values of biomarker under placebo
pred2 = predict(model5, newdata = data.frame(x5, trt = 1), se.fit = TRUE)#predict for values of biomarker under experimental

# Calculate the contrast and its confidence interval
contrast5 = -pred1$fit + pred2$fit #Placebo-Experimental (hence why we do <0)
se.constrat5 = sqrt(pred1$se.fit^2 + pred2$se.fit^2) #standard error of the difference 
lwr.contrast5 = contrast5 - qnorm(1 - 0.025) * se.constrat5
upr.contrast5 = contrast5 + qnorm(1 - 0.025) * se.constrat5
data_contrast_5 <- data.frame(x,contrast5,lwr.contrast5,upr.contrast5)




# ================================
# Overlay plot: interim + final
# ================================
# Rename columns of data_contrast_2
names(data_contrast_2) <- c("x", "contrast", "lwr.contrast", "upr.contrast")
names(data_contrast_3) <- c("x", "contrast", "lwr.contrast", "upr.contrast")
names(data_contrast_4) <- c("x", "contrast", "lwr.contrast", "upr.contrast")
names(data_contrast_5) <- c("x", "contrast", "lwr.contrast", "upr.contrast")


# Add stage labels for plotting
data_contrast_1$stage <- "20"
data_contrast_2$stage <- "30"
data_contrast_3$stage <- "35"
data_contrast_4$stage <- "40"
data_contrast_5$stage <- "50"

df_plot <- rbind(data_contrast_1, data_contrast_2,data_contrast_3,data_contrast_4,data_contrast_5)

library(ggplot2)
library(patchwork)

# Interim plot (left)
p_interim <- ggplot(subset(df_plot, stage=="20"), aes(x=x, y=contrast, color=stage, fill=stage)) +
  geom_ribbon(aes(ymin=lwr.contrast, ymax=upr.contrast), alpha=1, color=NA) +
  geom_line(size=1.5) +
  geom_line(aes(y=lwr.contrast), linetype="dashed", size=1) +
  geom_line(aes(y=upr.contrast), linetype="dashed", size=1) +
  geom_hline(yintercept=0, linetype=2, color="black") +
  xlab("Biomarker value") +
  ylab("Contrast (Experimental - Placebo)") +
  scale_color_manual(values=c("20"="blue")) +
  scale_fill_manual(values=c("20"="lightblue")) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.position="none") +
  ggtitle("Interim Analysis")+ylim(-5,30)

# Final + Final 2 plot (right, overlaid)
p_final_overlay <- ggplot(subset(df_plot, stage %in% c("30", "35","40","50")), 
                          aes(x=x, y=contrast, color=stage, fill=stage)) +
  geom_ribbon(aes(ymin=lwr.contrast, ymax=upr.contrast), alpha=1, color=NA) +
  geom_line(size=1.5) +
  geom_line(aes(y=lwr.contrast), linetype="dashed", size=1) +
  geom_line(aes(y=upr.contrast), linetype="dashed", size=1) +
  geom_hline(yintercept=0, linetype=2, color="black") +
  xlab("Biomarker value") +
  ylab("Contrast (Experimental - Placebo)") +
  scale_color_manual(
    name="Post-Interim Allocation", # Legend title
    values=c("30"="purple", "Final 2"="orange"),
    labels=c("30", "35","40","50") # Legend labels
  ) +
  scale_fill_manual(
    name="Post-Interim Allocation", # Legend title
    values=c("30"="lightpink", "Final 2"="lightyellow"),
    labels=c("30", "35","40","50")
  ) +
  theme_bw()  +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  ggtitle("Final Analyses")+ylim(-5,30)

# Combine side by side
combined_plot <- p_interim | p_final_overlay
print(combined_plot)



# Interim plot (left)
p_interim <- ggplot(df_plot, aes(x=x, y=contrast, color=stage, fill=stage)) +
  geom_ribbon(aes(ymin=lwr.contrast, ymax=upr.contrast), alpha=1, color=NA) +
  geom_line(size=1.5) +
  geom_line(aes(y=lwr.contrast), linetype="dashed", size=1) +
  geom_line(aes(y=upr.contrast), linetype="dashed", size=1) +
  geom_hline(yintercept=0, linetype=2, color="black") +
  xlab("Biomarker value") +
  ylab("Contrast (Experimental - Placebo)") +
  scale_color_manual(values=c("20"="blue")) +
  scale_fill_manual(values=c("20"="lightblue")) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.position="none") +
  ggtitle("Interim Analysis")+ylim(-5,30)
