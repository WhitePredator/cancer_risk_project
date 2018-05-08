t=c(32,37,42,47,52,57,62,67)


par(mfrow=c(3,3))

for (i in 1:9) {
    plot(log10(t),log10(bladder_list[[2*i-1]]),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(bladder_list[2*i-1])),14,50))
    abline(lm(log10(bladder_list[[2*i-1]])~log10(t)),col="blue")
    points(log10(t),log10(bladder_list[[2*i]]),col="red")
    abline(lm(log10(bladder_list[[2*i]])~log10(t)),col="red")
    
}


for (i in 1:9) {
  plot(log10(t),log10(colon_list[[2*i-1]]),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(colon_list[2*i-1])),12,50))
  abline(lm(log10(colon_list[[2*i-1]])~log10(t)),col="blue")
  points(log10(t),log10(colon_list[[2*i]]),col="red")
  abline(lm(log10(colon_list[[2*i]])~log10(t)),col="red")
  
}


for (i in 1:9) {
  plot(log10(t),log10(melanoma_list[[2*i-1]]),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(melanoma_list[2*i-1])),15,50))
  abline(lm(log10(melanoma_list[[2*i-1]])~log10(t)),col="blue")
  points(log10(t),log10(melanoma_list[[2*i]]),col="red")
  abline(lm(log10(melanoma_list[[2*i]])~log10(t)),col="red")
  
}

for (i in 1:9) {
  plot(log10(t),log10(lusc_list[[2*i-1]]+10^-100),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(lusc_list[2*i-1])),11,50))
  abline(lm(log10(lusc_list[[2*i-1]]+10^-100)~log10(t)),col="blue")
  points(log10(t),log10(lusc_list[[2*i]]+10^-100),col="red")
  abline(lm(log10(lusc_list[[2*i]]+10^-100)~log10(t)),col="red")
  
}

for (i in 1:9) {
  plot(log10(t),log10(stomach_list[[2*i-1]]+10^-100),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(stomach_list[2*i-1])),14,50))
  abline(lm(log10(stomach_list[[2*i-1]]+10^-100)~log10(t)),col="blue")
  points(log10(t),log10(stomach_list[[2*i]]+10^-100),col="red")
  abline(lm(log10(stomach_list[[2*i]]+10^-100)~log10(t)),col="red")
  
}


for (i in 1:9) {
  plot(log10(t),log10(liver_list[[2*i-1]]+10^-100),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(liver_list[2*i-1])),12,50))
  abline(lm(log10(liver_list[[2*i-1]]+10^-100)~log10(t)),col="blue")
  points(log10(t),log10(liver_list[[2*i]]+10^-100),col="red")
  abline(lm(log10(liver_list[[2*i]]+10^-100)~log10(t)),col="red")
  
}

for (i in 1:9) {
  plot(log10(t),log10(luad_list[[2*i-1]]+10^-100),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(luad_list[2*i-1])),11,50))
  abline(lm(log10(luad_list[[2*i-1]]+10^-100)~log10(t)),col="blue")
  points(log10(t),log10(luad_list[[2*i]]+10^-100),col="red")
  abline(lm(log10(luad_list[[2*i]]+10^-100)~log10(t)),col="red")
  
}


for (i in 1:9) {
  plot(log10(t),log10(lusc_list[[2*i-1]]+10^-100),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(lusc_list[2*i-1])),11,50))
  abline(lm(log10(lusc_list[[2*i-1]]+10^-100)~log10(t)),col="blue")
  points(log10(t),log10(lusc_list[[2*i]]+10^-100),col="red")
  abline(lm(log10(lusc_list[[2*i]]+10^-100)~log10(t)),col="red")
  
}