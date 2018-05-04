load("D:\\my_study\\Project\\cancer_risk\\data\\country.Rdata")

t_life=c(2,7,12,17,22,27,32,37,42,47,52,57,62,67,72,77,82,87)
t=c(32,37,42,47,52,57,62,67)


###oesophagus
oesophagus_num=c(24:27)   #TCGA上癌种为Carcinoma，所以选取oesophagus中属于Carcinoma的部分
oesophagus_list=list()
for (j in 1:length(country_list)) {
    oesophagus=matrix(data = 0,ncol = 1,nrow = length(t))
    population=matrix(data = 0,ncol = 1,nrow = length(t))
    for (i in oesophagus_num) {
        a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
        b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
        oesophagus=oesophagus+a
        population=population+b
    }
    #assign(paste('oesophagus_male_',names(country_list[j]),sep=''),oesophagus/population)
    oesophagus_list[[paste('oesophagus_male_',names(country_list[j]),sep='')]]= oesophagus/population
    for (i in oesophagus_num) {
      a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
      b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
      oesophagus=oesophagus+a
      population=population+b
    }
    #assign(paste('oesophagus_female_',names(country_list[j]),sep=''),oesophagus/population)
    oesophagus_list[[paste('oesophagus_female_',names(country_list[j]),sep='')]]= oesophagus/population
}
oesophagus_lm_list=list()
for (i in 1:length(oesophagus_list)) {
  #assign(names(oesophagus_list[i]),summary(lm(log(oesophagus_list[[i]])~log(t))))  
  oesophagus_lm_list[[paste(names(oesophagus_list[i]))]]=summary(lm(log(oesophagus_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}


###stomach
stomach_num=35   
stomach_list=list()
for (j in 1:length(country_list)) {
  stomach=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in stomach_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    stomach=stomach+a
    population=population+b
  }
  stomach_list[[paste('stomach_male_',names(country_list[j]),sep='')]]= stomach/population
  for (i in stomach_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    stomach=stomach+a
    population=population+b
  }
  stomach_list[[paste('stomach_female_',names(country_list[j]),sep='')]]= stomach/population
}
stomach_lm_list=list()
for (i in 1:length(stomach_list)) {
  stomach_lm_list[[paste(names(stomach_list[i]))]]=summary(lm(log(stomach_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}



###pancreas
pancreas_num=70   
pancreas_list=list()
for (j in 1:length(country_list)) {
  pancreas=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in pancreas_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    pancreas=pancreas+a
    population=population+b
  }
  pancreas_list[[paste('pancreas_male_',names(country_list[j]),sep='')]]= pancreas/population
  for (i in pancreas_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    pancreas=pancreas+a
    population=population+b
  }
  pancreas_list[[paste('pancreas_female_',names(country_list[j]),sep='')]]= pancreas/population
}
pancreas_lm_list=list()
for (i in 1:length(pancreas_list)) {
  pancreas_lm_list[[paste(names(pancreas_list[i]))]]=summary(lm(log(pancreas_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}



###melanoma
melanoma_num=100   
melanoma_list=list()
for (j in 1:length(country_list)) {
  melanoma=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in melanoma_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    melanoma=melanoma+a
    population=population+b
  }
  melanoma_list[[paste('melanoma_male_',names(country_list[j]),sep='')]]= melanoma/population
  for (i in melanoma_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    melanoma=melanoma+a
    population=population+b
  }
  melanoma_list[[paste('melanoma_female_',names(country_list[j]),sep='')]]= melanoma/population
}
melanoma_lm_list=list()
for (i in 1:length(melanoma_list)) {
  melanoma_lm_list[[paste(names(melanoma_list[i]))]]=summary(lm(log(melanoma_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}



###lusc
lusc_num=78   
lusc_list=list()
for (j in 1:length(country_list)) {
  lusc=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in lusc_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    lusc=lusc+a
    population=population+b
  }
  lusc_list[[paste('lusc_male_',names(country_list[j]),sep='')]]= lusc/population
  for (i in lusc_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    lusc=lusc+a
    population=population+b
  }
  lusc_list[[paste('lusc_female_',names(country_list[j]),sep='')]]= lusc/population
}
lusc_lm_list=list()
for (i in 1:length(lusc_list)) {
  lusc_lm_list[[paste(names(lusc_list[i]))]]=summary(lm(log(lusc_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}


###luad
luad_num=79   
luad_list=list()
for (j in 1:length(country_list)) {
  luad=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in luad_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    luad=luad+a
    population=population+b
  }
  luad_list[[paste('luad_male_',names(country_list[j]),sep='')]]= luad/population
  for (i in luad_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    luad=luad+a
    population=population+b
  }
  luad_list[[paste('luad_female_',names(country_list[j]),sep='')]]= luad/population
}
luad_lm_list=list()
for (i in 1:length(luad_list)) {
  luad_lm_list[[paste(names(luad_list[i]))]]=summary(lm(log(luad_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}

###liver  Liver Hepatocellular Carcinoma
liver_num=59   
liver_list=list()
for (j in 1:length(country_list)) {
  liver=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in liver_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    liver=liver+a
    population=population+b
  }
  liver_list[[paste('liver_male_',names(country_list[j]),sep='')]]= liver/population
  for (i in liver_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    liver=liver+a
    population=population+b
  }
  liver_list[[paste('liver_female_',names(country_list[j]),sep='')]]= liver/population
}
liver_lm_list=list()
for (i in 1:length(liver_list)) {
  liver_lm_list[[paste(names(liver_list[i]))]]=summary(lm(log(liver_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}



###kindey  
kindey_num=160   
kindey_list=list()
for (j in 1:length(country_list)) {
  kindey=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in kindey_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    kindey=kindey+a
    population=population+b
  }
  kindey_list[[paste('kindey_male_',names(country_list[j]),sep='')]]= kindey/population
  for (i in kindey_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    kindey=kindey+a
    population=population+b
  }
  kindey_list[[paste('kindey_female_',names(country_list[j]),sep='')]]= kindey/population
}
kindey_lm_list=list()
for (i in 1:length(kindey_list)) {
  kindey_lm_list[[paste(names(kindey_list[i]))]]=summary(lm(log(kindey_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}


###bladder    transitional cell carcinoma  Urothelial cells are transitional cells
bladder_num=165   
bladder_list=list()
for (j in 1:length(country_list)) {
  bladder=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in bladder_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    bladder=bladder+a
    population=population+b
  }
  bladder_list[[paste('bladder_male_',names(country_list[j]),sep='')]]= bladder/population
  for (i in bladder_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    bladder=bladder+a
    population=population+b
  }
  bladder_list[[paste('bladder_female_',names(country_list[j]),sep='')]]= bladder/population
}
bladder_lm_list=list()
for (i in 1:length(bladder_list)) {
  bladder_lm_list[[paste(names(bladder_list[i]))]]=summary(lm(log(bladder_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}



###colon    
colon_num=42   
colon_list=list()
for (j in 1:length(country_list)) {
  colon=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in colon_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    colon=colon+a
    population=population+b
  }
  colon_list[[paste('colon_male_',names(country_list[j]),sep='')]]= colon/population
  for (i in colon_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    colon=colon+a
    population=population+b
  }
  colon_list[[paste('colon_female_',names(country_list[j]),sep='')]]= colon/population
}
colon_lm_list=list()
for (i in 1:length(colon_list)) {
  colon_lm_list[[paste(names(colon_list[i]))]]=summary(lm(log(colon_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}


###rectum    
rectum_num=49   
rectum_list=list()
for (j in 1:length(country_list)) {
  rectum=matrix(data = 0,ncol = 1,nrow = length(t))
  population=matrix(data = 0,ncol = 1,nrow = length(t))
  for (i in rectum_num) {
    a=country_list[[j]][((i-1)*38+7):((i-1)*38+14),4]
    b=country_list[[j]][((i-1)*38+7):((i-1)*38+14),5]
    rectum=rectum+a
    population=population+b
  }
  rectum_list[[paste('rectum_male_',names(country_list[j]),sep='')]]= rectum/population
  for (i in rectum_num) {
    a=country_list[[j]][((i-1)*39+7):((i-1)*39+14),4]
    b=country_list[[j]][((i-1)*39+7):((i-1)*39+14),5]
    rectum=rectum+a
    population=population+b
  }
  rectum_list[[paste('rectum_female_',names(country_list[j]),sep='')]]= rectum/population
}
rectum_lm_list=list()
for (i in 1:length(rectum_list)) {
  rectum_lm_list[[paste(names(rectum_list[i]))]]=summary(lm(log(rectum_list[[i]]+exp(-100))~log(t)))  ###加exp(-100)，防止出现0取log出错
}


save.image(file = "D:\\my_study\\Project\\cancer_risk\\data\\lm.Rdata")