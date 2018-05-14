#load("D:/my_study/Project/cancer_risk/data_git/new_result.Rdata")
log_t=log10(t)


country_id=11
cancer_id=9

par(pin=c(3,2.5))
plot(log_t,log10(country_cancer_list[[country_id]][2*cancer_id-1,]),col="blue",xlab = "age",ylab = "incidence",main = paste(rownames(cancer_dict)[cancer_id],'_',names(country_cancer_list)[country_id],sep=''))
abline(lm(log10(country_cancer_list[[country_id]][2*cancer_id-1,])~log_t),col="blue")
points(log_t,log10(country_cancer_list[[country_id]][2*cancer_id,]),col="red")
abline(lm(log10(country_cancer_list[[country_id]][2*cancer_id,])~log_t),col="red")
    






#for (i in 1:9) {
  #plot(log10(t),log10(luad_list[[2*i-1]]+10^-100),col="blue",xlab = "age",ylab = "incidence",main = substr(paste(names(luad_list[2*i-1])),11,50))
  #abline(lm(log10(luad_list[[2*i-1]]+10^-100)~log10(t)),col="blue")
  #points(log10(t),log10(luad_list[[2*i]]+10^-100),col="red")
  #abline(lm(log10(luad_list[[2*i]]+10^-100)~log10(t)),col="red")
  
#}


