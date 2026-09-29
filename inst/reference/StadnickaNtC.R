# code for Julita in R, which will later become part of the shiny app. 

# the libraries below need to be installed
# load libraries for plotting
library(ggplot2)
library(cowplot)
library(nlmrt)
library(nls2)

## Read data from a csv file. 
# the first column in the file are the concentrations, while the other columns are replicates
# !!!!!!!! CHANGE PATH TO FILE !!!!!! ONLY WORKS WITH CSV FILES WITH SEMICOLON BELOW !!!!!!!!!!!
data = read.csv(file = "PATH_TO_FILE", sep = ";")
conc = data[2:nrow(data),1]
logconc = log10(data[2:nrow(data),1])
effect_rep = data[2:nrow(data),2:ncol(data)]
effect_mean = rowMeans(effect_rep, na.rm = TRUE)

# prepare a data frame for plotting replicates as points
data_for_point_plot = cbind(logconc,effect_rep)
data_for_point_plot = as.data.frame(data_for_point_plot)

## Make a new dataframe, ready for fitting (remove controls)
data_for_fit = cbind(logconc,rowMeans(effect_rep, na.rm = TRUE))
data_for_fit = as.data.frame(data_for_fit)
names(data_for_fit) = c("concentration", "effect")

## Formula for the sigmoidal dose response curve
drc_formula = effect ~ 100/(1+10^((logEC50-concentration)*slope))

## fit sigmoidal with nlmrt
fit1 = nlxb(drc_formula,
                start = c(logEC50 = log10(conc[ceiling(length(conc)/2)]), slope = -2),
                trace = FALSE,
                data = data_for_fit)
curvefit <- nls2(drc_formula, data = data_for_fit, start = fit1$coefficients,
                algorithm = "brute-force")

# get parameter values of the fit
coefficient = coef(curvefit)
logEC50_value = coefficient[1]
slope_value = coefficient[2]


## get the confidence intervals
# confidence intervals are calculated beyond the concentrations, at which measurements were taken
# these are the concentration values at which to plot the curve
if (logconc[1]>0){
  x_start = 0
} else {
  x_start = logconc[1] + 2*logconc[1]
}
x_values = seq(x_start,logconc[length(logconc)]*2,0.001)

# get the fitted values for all x_values
curve.predict = predict(curvefit, newdata = data.frame(concentration = c(x_values)))

# get the covariance matrix for the parameters
covmatrix= vcov(curvefit)

# get the jacobian for all x_values
pp = c(logEC50_value, slope_value)
myjacfun = model2jacfun(drc_formula,pp)
myjacobian = myjacfun(pp, effect = curve.predict, concentration = x_values)

# get the confidence interval for all x_values
error = myjacobian %*% covmatrix %*% t(myjacobian)
error = diag(error)
error = sqrt(error)
# t for the inverse studen't distribution depends on the number of measurements
t_student = qt(.975, df = nrow(data)-3)
CI = error*t_student
upperCI = curve.predict + CI
lowerCI = curve.predict - CI

## NtC upper 
# find concentration at which the upper confidence interval falls below
# 100 
inds = max(which(upperCI >= 99.999999999))
value_atupper = 100 - curve.predict[inds]
if (value_atupper > 0){
  NtC_upper = 10^(x_values[inds])
} else {
  NtC_upper = 0}

## NtC lower
# find concentration at which lower confidence interval is 90
inds = max(which(lowerCI > 90))
value_atlower = 100 - curve.predict[inds]
if (value_atlower > 0){
  NtC_lower = 10^(x_values[inds])
} else {
  NtC_lower = 0}

## Choose NtC based on upper and lower CI criteria
if (NtC_upper > NtC_lower) {
  NtC = NtC_lower
  Effect_NtC = value_atlower
} else {
  if (value_atlower/value_atupper > 10) {
    Effect_NtC = value_atlower/10
    NtC = 10^(logEC50_value-log10(100/(100-Effect_NtC)-1)/slope_value)
  } else {
    NtC = NtC_upper
    Effect_NtC = value_atupper
  }
} 

## Include measured data if the effect of any replicate (at concentration below
# the concentration, which has been thus far chosen as the NtC) is > 10%
if (length(which(conc < NtC)) > 0) {
  inds = which(effect_rep[which(conc < NtC),] <= 90)
  if (length(inds) > 0) {
    remained = inds%%nrow(effect_rep[which(conc < NtC),])
    remained = replace(remained, remained == 0,nrow(effect_rep[which(conc < NtC),]))
    NtC = conc[min(remained)]
  }
}
  
## plotting
# prepare data in data.frame for plotting
data_for_plot = data.frame(cbind(x_values,curve.predict, upperCI, lowerCI, logconc, c(as.matrix(effect_rep))))
# find the correct value of the predicted curve for plotting (value at NtC)
temp = abs(x_values -log10(NtC))
index = which(temp == min(temp))
# plot itself
ggplot(data_for_plot, aes(x_values, curve.predict)) + geom_line(aes(y = curve.predict, colour = "red")) +
  geom_line(aes(y = upperCI), colour = "red", linetype = 2) + geom_line(aes(y = lowerCI), colour = "red", linetype = 2) +
  geom_segment(aes(x = log10(NtC), xend = log10(NtC), y = -15, yend = curve.predict[index]), colour = "blue", linetype = 2) + 
  geom_point(aes(x = log10(NtC), y = curve.predict[index]), shape = 4, size = 8, colour = "blue") + 
  geom_point(aes(x = logconc, y = V6)) +
  xlim(min(log10(NtC)/2,log10(conc[1])), log10(conc[length(conc)])*1.3) + ylim(-15,115) + 
  scale_y_continuous(breaks = seq(0,100, by = 20),"100-Effect, %") +
  xlab("log(conc), unit") + ylab("100-Effect, %") + theme(legend.position = "none") +ggtitle("your title") 


  



  


