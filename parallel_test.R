#### non parallel  test
t=c(32,37,42,47,52,57,62,67)
sex=c(rep(1,8),rep(0,8))
t2=log10(rep(t,2))

load("D:\\my_study\\Project\\cancer_risk\\data\\lm.Rdata")

delta=0.6  #等效（平行）检验边界
beta3_threshold=0.3



###colon 
p_value_colon=matrix(data=0,nrow = length(colon_list)/2,ncol = 6)


for (i in 1:(length(colon_list)/2)) {
  p_value_colon[i,1]=substr(paste(names(colon_list[2*i-1])),12,50)
  colon_test=log10(rbind(colon_list[[2*i-1]],colon_list[[2*i]])+10^-100)
  data_colon=data.frame(test=colon_test,t2=t2,sex=sex,sex_t=sex*t2)
  colon_fit=summary(lm(colon_test~t2+sex+sex_t,data=data_colon))
  p_value_colon[i,2]=colon_fit[[4]][4,1] 
  p_value_colon[i,3]=colon_fit[[4]][4,2]
  p_value_colon[i,4]=colon_fit[[4]][4,4] 
  p_value_colon[i,5]=colon_fit[[8]]
  if ((colon_fit[[4]][4,1] + 1.96*colon_fit[[4]][4,2] <= delta) & (colon_fit[[4]][4,1] - 1.96*colon_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_colon[i,6]=1  #  1代表平行
  }
  else{
    p_value_colon[i,6]=0
  }
}


###melanoma
p_value_melanoma=matrix(data=0,nrow = length(melanoma_list)/2,ncol = 6)


for (i in 1:(length(melanoma_list)/2)) {
  p_value_melanoma[i,1]=substr(paste(names(melanoma_list[2*i-1])),15,50)
  melanoma_test=log10(rbind(melanoma_list[[2*i-1]],melanoma_list[[2*i]])+10^-100)
  data_melanoma=data.frame(test=melanoma_test,t2=t2,sex=sex,sex_t=sex*t2)
  melanoma_fit=summary(lm(melanoma_test~t2+sex+sex_t,data=data_melanoma))
  p_value_melanoma[i,2]=melanoma_fit[[4]][4,1]
  p_value_melanoma[i,3]=melanoma_fit[[4]][4,2]
  p_value_melanoma[i,4]=melanoma_fit[[4]][4,4]
  p_value_melanoma[i,5]=melanoma_fit[[8]]
  if ((melanoma_fit[[4]][4,1] + 1.96*melanoma_fit[[4]][4,2] <= delta) & (melanoma_fit[[4]][4,1] - 1.96*melanoma_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_melanoma[i,6]=1  #  1代表平行
  }
  else{
    p_value_melanoma[i,6]=0
  }
}

###oesophagus
p_value_oesophagus=matrix(data=0,nrow = length(oesophagus_list)/2,ncol = 6)


for (i in 1:(length(oesophagus_list)/2)) {
  p_value_oesophagus[i,1]=substr(paste(names(oesophagus_list[2*i-1])),17,50)
  oesophagus_test=log10(rbind(oesophagus_list[[2*i-1]],oesophagus_list[[2*i]])+10^-100) #log(0)
  data_oesophagus=data.frame(test=oesophagus_test,t2=t2,sex=sex,sex_t=sex*t2)
  oesophagus_fit=summary(lm(oesophagus_test~t2+sex+sex_t,data=data_oesophagus))
  p_value_oesophagus[i,2]=oesophagus_fit[[4]][4,1]
  p_value_oesophagus[i,3]=oesophagus_fit[[4]][4,2]
  p_value_oesophagus[i,4]=oesophagus_fit[[4]][4,4]
  p_value_oesophagus[i,5]=oesophagus_fit[[8]]
  if ((oesophagus_fit[[4]][4,1] + 1.96*oesophagus_fit[[4]][4,2] <= delta) & (oesophagus_fit[[4]][4,1] - 1.96*oesophagus_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_oesophagus[i,6]=1  #  1代表平行
  }
  else{
    p_value_oesophagus[i,6]=0
  }
  
}


###stomach
p_value_stomach=matrix(data=0,nrow = length(stomach_list)/2,ncol = 6)

for (i in 1:(length(stomach_list)/2)) {
  p_value_stomach[i,1]=substr(paste(names(stomach_list[2*i-1])),14,50)
  stomach_test=log10(rbind(stomach_list[[2*i-1]],stomach_list[[2*i]])+10^-100)
  data_stomach=data.frame(test=stomach_test,t2=t2,sex=sex,sex_t=sex*t2)
  stomach_fit=summary(lm(stomach_test~t2+sex+sex_t,data=data_stomach))
  p_value_stomach[i,2]=stomach_fit[[4]][4,1]
  p_value_stomach[i,3]=stomach_fit[[4]][4,2]
  p_value_stomach[i,4]=stomach_fit[[4]][4,4]
  p_value_stomach[i,5]=stomach_fit[[8]]
  if ((stomach_fit[[4]][4,1] + 1.96*stomach_fit[[4]][4,2] <= delta) & (stomach_fit[[4]][4,1] - 1.96*stomach_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_stomach[i,6]=1  #  1代表平行
  }
  else{
    p_value_stomach[i,6]=0
  }
}



###pancreas
p_value_pancreas=matrix(data=0,nrow = length(pancreas_list)/2,ncol = 6)

for (i in 1:(length(pancreas_list)/2)) {
  p_value_pancreas[i,1]=substr(paste(names(pancreas_list[2*i-1])),15,50)
  pancreas_test=log10(rbind(pancreas_list[[2*i-1]],pancreas_list[[2*i]])+10^-100)
  data_pancreas=data.frame(test=pancreas_test,t2=t2,sex=sex,sex_t=sex*t2)
  pancreas_fit=summary(lm(pancreas_test~t2+sex+sex_t,data=data_pancreas))
  p_value_pancreas[i,2]=pancreas_fit[[4]][4,1]
  p_value_pancreas[i,3]=pancreas_fit[[4]][4,2]
  p_value_pancreas[i,4]=pancreas_fit[[4]][4,4]
  p_value_pancreas[i,5]=pancreas_fit[[8]]
  if ((pancreas_fit[[4]][4,1] + 1.96*pancreas_fit[[4]][4,2] <= delta) & (pancreas_fit[[4]][4,1] - 1.96*pancreas_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_pancreas[i,6]=1  #  1代表平行
  }
  else{
    p_value_pancreas[i,6]=0
  }
}

###lusc
p_value_lusc=matrix(data=0,nrow = length(lusc_list)/2,ncol = 6)


for (i in 1:(length(lusc_list)/2)) {
  p_value_lusc[i,1]=substr(paste(names(lusc_list[2*i-1])),11,50)
  lusc_test=log10(rbind(lusc_list[[2*i-1]],lusc_list[[2*i]])+10^-100)
  data_lusc=data.frame(test=lusc_test,t2=t2,sex=sex,sex_t=sex*t2)
  lusc_fit=summary(lm(lusc_test~t2+sex+sex_t,data=data_lusc))
  p_value_lusc[i,2]=lusc_fit[[4]][4,1]
  p_value_lusc[i,3]=lusc_fit[[4]][4,2]
  p_value_lusc[i,4]=lusc_fit[[4]][4,4]
  p_value_lusc[i,5]=lusc_fit[[8]]
  if ((lusc_fit[[4]][4,1] + 1.96*lusc_fit[[4]][4,2] <= delta) & (lusc_fit[[4]][4,1] - 1.96*lusc_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_lusc[i,6]=1  #  1代表平行
  }
  else{
    p_value_lusc[i,6]=0
  }
}


###luad
p_value_luad=matrix(data=0,nrow = length(luad_list)/2,ncol = 6)


for (i in 1:(length(luad_list)/2)) {
  p_value_luad[i,1]=substr(paste(names(luad_list[2*i-1])),11,50)
  luad_test=log10(rbind(luad_list[[2*i-1]],luad_list[[2*i]])+10^-100)
  data_luad=data.frame(test=luad_test,t2=t2,sex=sex,sex_t=sex*t2)
  luad_fit=summary(lm(luad_test~t2+sex+sex_t,data=data_luad))
  p_value_luad[i,2]=luad_fit[[4]][4,1]
  p_value_luad[i,3]=luad_fit[[4]][4,2]
  p_value_luad[i,4]=luad_fit[[4]][4,4]
  p_value_luad[i,5]=luad_fit[[8]]
  if ((luad_fit[[4]][4,1] + 1.96*luad_fit[[4]][4,2] <= delta) & (luad_fit[[4]][4,1] - 1.96*luad_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_luad[i,6]=1  #  1代表平行
  }
  else{
    p_value_luad[i,6]=0
  }
}


###liver
p_value_liver=matrix(data=0,nrow = length(liver_list)/2,ncol = 6)


for (i in 1:(length(liver_list)/2)) {
  p_value_liver[i,1]=substr(paste(names(liver_list[2*i-1])),12,50)
  liver_test=log10(rbind(liver_list[[2*i-1]],liver_list[[2*i]])+10^-100)
  data_liver=data.frame(test=liver_test,t2=t2,sex=sex,sex_t=sex*t2)
  liver_fit=summary(lm(liver_test~t2+sex+sex_t,data=data_liver))
  p_value_liver[i,2]=liver_fit[[4]][4,1]
  p_value_liver[i,3]=liver_fit[[4]][4,2]
  p_value_liver[i,4]=liver_fit[[4]][4,4]
  p_value_liver[i,5]=liver_fit[[8]]
  if ((liver_fit[[4]][4,1] + 1.96*liver_fit[[4]][4,2] <= delta) & (liver_fit[[4]][4,1] - 1.96*liver_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_liver[i,6]=1  #  1代表平行
  }
  else{
    p_value_liver[i,6]=0
  }
}



###kindey
p_value_kindey=matrix(data=0,nrow = length(kindey_list)/2,ncol = 6)


for (i in 1:(length(kindey_list)/2)) {
  p_value_kindey[i,1]=substr(paste(names(kindey_list[2*i-1])),13,50)
  kindey_test=log10(rbind(kindey_list[[2*i-1]],kindey_list[[2*i]])+10^-100)
  data_kindey=data.frame(test=kindey_test,t2=t2,sex=sex,sex_t=sex*t2)
  kindey_fit=summary(lm(kindey_test~t2+sex+sex_t,data=data_kindey))
  p_value_kindey[i,2]=kindey_fit[[4]][4,1]
  p_value_kindey[i,3]=kindey_fit[[4]][4,2]
  p_value_kindey[i,4]=kindey_fit[[4]][4,4]
  p_value_kindey[i,5]=kindey_fit[[8]]
  if ((kindey_fit[[4]][4,1] + 1.96*kindey_fit[[4]][4,2] <= delta) & (kindey_fit[[4]][4,1] - 1.96*kindey_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_kindey[i,6]=1  #  1代表平行
  }
  else{
    p_value_kindey[i,6]=0
  }
}



###bladder
p_value_bladder=matrix(data=0,nrow = length(bladder_list)/2,ncol = 6)


for (i in 1:(length(bladder_list)/2)) {
  p_value_bladder[i,1]=substr(paste(names(bladder_list[2*i-1])),14,50)
  bladder_test=log10(rbind(bladder_list[[2*i-1]],bladder_list[[2*i]])+10^-100)
  data_bladder=data.frame(test=bladder_test,t2=t2,sex=sex,sex_t=sex*t2)
  bladder_fit=summary(lm(bladder_test~t2+sex+sex_t,data=data_bladder))
  p_value_bladder[i,2]=bladder_fit[[4]][4,1]
  p_value_bladder[i,3]=bladder_fit[[4]][4,2]
  p_value_bladder[i,4]=bladder_fit[[4]][4,4]
  p_value_bladder[i,5]=bladder_fit[[8]]
  if ((bladder_fit[[4]][4,1] + 1.96*bladder_fit[[4]][4,2] <= delta) & (bladder_fit[[4]][4,1] - 1.96*bladder_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_bladder[i,6]=1  #  1代表平行
  }
  else{
    p_value_bladder[i,6]=0
  }
}


###rectum
p_value_rectum=matrix(data=0,nrow = length(rectum_list)/2,ncol = 6)


for (i in 1:(length(rectum_list)/2)) {
  p_value_rectum[i,1]=substr(paste(names(rectum_list[2*i-1])),13,50)
  rectum_test=log10(rbind(rectum_list[[2*i-1]],rectum_list[[2*i]])+10^-100)
  data_rectum=data.frame(test=rectum_test,t2=t2,sex=sex,sex_t=sex*t2)
  rectum_fit=summary(lm(rectum_test~t2+sex+sex_t,data=data_rectum))
  p_value_rectum[i,2]=rectum_fit[[4]][4,1]
  p_value_rectum[i,3]=rectum_fit[[4]][4,2]
  p_value_rectum[i,4]=rectum_fit[[4]][4,4]
  p_value_rectum[i,5]=rectum_fit[[8]]
  if ((rectum_fit[[4]][4,1] + 1.96*rectum_fit[[4]][4,2] <= delta) & (rectum_fit[[4]][4,1] - 1.96*rectum_fit[[4]][4,2] >= -delta) & (abs(colon_fit[[4]][4,1]) <= beta3_threshold)) {
    p_value_rectum[i,6]=1  #  1代表平行
  }
  else{
    p_value_rectum[i,6]=0
  }
}

colnames(p_value_colon)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_bladder)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_kindey)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_liver)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_luad)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_lusc)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_melanoma)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_oesophagus)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_pancreas)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_rectum)=c("country", "beta", "SE","p_value","R2","Parallel")
colnames(p_value_stomach)=c("country", "beta", "SE","p_value","R2","Parallel")


save.image(file = "D:\\my_study\\Project\\cancer_risk\\data\\non_parallel_result.Rdata")




