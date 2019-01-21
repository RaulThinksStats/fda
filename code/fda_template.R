#**************************************************************************
# FDA - fda_template.R
#
# Author: Raul J.T.A.
#
# Summary: Cursory overview of analysis involving a functional object.
#***************************************************************************

#preliminary Settings----------------------------------------------------


library(fda)
library(fda.usc)
library(tidyverse)
library(VIM)

setwd("data/")


#setting up Data--------------------------------------------------------

temp = list.files(pattern="*.csv")
rfm = map(temp, read_csv)

#chop off data where mask is leaking, i.e. -1
rfm <- map(rfm, function(x) x[(which(x$Leak != -1)[1]):length(x$Leak),])

#matrix of replicates needed for fda and easier to plot
flow <- as.matrix(map_dfc(rfm, function(x) x[1:600,"Flow"]))
vte <- as.matrix(map_dfc(rfm, function(x) x[1:600,"Vte"]))
vti <- as.matrix(map_dfc(rfm, function(x) x[1:600,"Vti"]))



#plotting raw data--------------------------------------------------------

par(mfrow = c(4,3))

for (i in 1:12){
  plot(flow[,i], type = "l", xlab = "Time (every .5s sec)",
       ylab = "Flow Metric", main = i)
}

for (i in 1:12){
  plot(vte[,i], type = "l", xlab = "Time (every .5s sec)",
       ylab = "Vte Metric", main = i)
}
for (i in 1:12){
  plot(vti[,i], type = "l", xlab = "Time (every .5s sec)",
       ylab = "Vti Metric", main = i)
}
par(mfrow = c(1,1))

#missing data + outliers ---------------------------------------------------------

#missingness simulation
miss_prob <- runif(13, 0, .3)
miss_ind <- vti_miss <- matrix(nrow = 600, ncol = 13)

for(i in 1:13){miss_ind[,i] <- rbinom(600, 12, miss_prob[i])
vti_miss[,i] <- miss_ind[,i] * vti[,i]}
vti_miss <- apply(vti_miss, 2, function(x) ifelse(x == 0, NA, x))

#visualizing missingness to explore potential patterns
aggr_plot <- aggr(vti_miss, col=c('navyblue','red'), numbers=TRUE, sortVars=FALSE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))


#creating functional object--------------------------------------------------------


fdnames = vector("list", 3)
patient = vector("list", 13)
patient[[1]]  = "Patient 1"
patient[[2]]  = "Patient 2"
patient[[3]]  = "Patient 3"
patient[[4]]  = "Patient 4"
patient[[5]]  = "Patient 5"
patient[[6]]  = "Patient 6"
patient[[7]]  = "Patient 7"
patient[[8]]  = "Patient 8"
patient[[9]]  = "Patient 9"
patient[[10]] = "Patient 10"
patient[[11]] = "Patient 11"
patient[[12]] = "Patient 12"
patient[[13]] = "Patient 13"
fdnames = list("Time (sec)",
               "Patient ID" = patient,
               "Flow Metric")


time = (1:600)
flowbasis = create.bspline.basis(c(0, 600), 65)
vtebasis = create.bspline.basis(c(0, 600), 65)
vtibasis = create.bspline.basis(c(0, 600), 65)

flowfd = smooth.basis(time, flow, flowbasis)$fd
flowfd$fdnames = list("Time (recorded every 0.5 seconds)",
                      "Patient ID",
                      "Flow Metric")
vtefd = smooth.basis(time, vte, vtebasis)$fd
vtefd$fdnames = list("Time (recorded every 0.5 seconds)",
                      "Patient ID",
                      "VTE Metric")
vtifd = smooth.basis(time, vti, vtibasis)$fd
vtifd$fdnames = list("Time (recorded every 0.5 seconds)",
                      "Patient ID",
                      "VTI Metric")

# visualizing functional object--------------------------------------------------------------


summary(flow)
par(mfrow = c(3,1))
plot(flowfd, lty=1, main = "Flow Curves for 13 Patients")
lines(mean(flowfd), lty = 1, lwd = 4, col = "black")

plot(vtefd, lty=1, main = "VTE Curves for 13 Patients", ylim = c(0,22.5))
lines(mean(vtefd), lty = 1, lwd = 4, col = "black")

plot(vtifd, lty=1, main = "VTI Curves for 13 Patients", ylim = c(0,85))
lines(mean(vtifd), lty = 1, lwd = 4, col = "black")
ar(mfrow = c(1,1))

#commented lines yield naive CIs

# principal components analysis---------------------------------------------------

par(mfrow = c(2,3))
flow.pcalist = pca.fd(flowfd, 3)
print(flow.pcalist$values)
plot.pca.fd(flow.pcalist, main = "PC for Flow")

flow.rotpcalist = varmx.pca.fd(flow.pcalist)
plot.pca.fd(flow.rotpcalist, main = "Rotated PC for Flow")

par(mfrow = c(2,3))
vte.pcalist = pca.fd(vtefd, 3)
print(vte.pcalist$values)
plot.pca.fd(vte.pcalist, main = "PC for VTE")

vte.rotpcalist = varmx.pca.fd(vte.pcalist)
plot.pca.fd(vte.rotpcalist, main = "Rotated PC for VTE")

par(mfrow = c(2,2))
vti.pcalist = pca.fd(vtifd, 2)
print(vti.pcalist$values)
plot.pca.fd(vti.pcalist, main = "PC for VTI")

vti.rotpcalist = varmx.pca.fd(vti.pcalist)
plot.pca.fd(vti.rotpcalist, main = "Rotated PC for VTI")

