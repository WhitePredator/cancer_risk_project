t=c(32,37,42,47,52,57,62,67)


par(mfrow=c(3,2))

for (i in 1:6) {
    plot(log10(t),log10(bladder_list[[2*i-1]]),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(bladder_list[2*i-1])),14,50))
    abline(lm(log10(bladder_list[[2*i-1]])~log10(t)),col="blue")
    points(log10(t),log10(bladder_list[[2*i]]),col="red")
    abline(lm(log10(bladder_list[[2*i]])~log10(t)),col="red")
    
}


for (i in 1:6) {
  plot(log10(t),log10(colon_list[[2*i-1]]),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(colon_list[2*i-1])),12,50))
  abline(lm(log10(colon_list[[2*i-1]])~log10(t)),col="blue")
  points(log10(t),log10(colon_list[[2*i]]),col="red")
  abline(lm(log10(colon_list[[2*i]])~log10(t)),col="red")
  
}


for (i in 1:6) {
  plot(log10(t),log10(melanoma_list[[2*i-1]]),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(melanoma_list[2*i-1])),15,50))
  abline(lm(log10(melanoma_list[[2*i-1]])~log10(t)),col="blue")
  points(log10(t),log10(melanoma_list[[2*i]]),col="red")
  abline(lm(log10(melanoma_list[[2*i]])~log10(t)),col="red")
  
}

##test   colon  melanoma
t2=log10(rep(t,2))
sex=c(rep(1,8),rep(0,8))
colon_test=log10(rbind(colon_list[[1]],colon_list[[2]]))
data_colon=data.frame(test=colon_test,t2=t2,sex=c(rep(1,8),rep(0,8)),sex_t=sex*t2)
colon_fit=summary(lm(colon_test~t2+sex+sex_t,data=data_colon))

melanoma_test=log10(rbind(melanoma_list[[1]],melanoma_list[[2]]))
data_melanoma=data.frame(test=melanoma_test,t2=t2,sex=c(rep(1,8),rep(0,8)),sex_t=sex*t2)
melanoma_fit=summary(lm(melanoma_test~t2+sex+sex_t,data=data_melanoma))