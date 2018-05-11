load("D:\\my_study\\Project\\cancer_risk\\data_git\\country.Rdata")
t=c(32,37,42,47,52,57,62,67)
country_cancer_list=list()
cancer_num=11   ###  关注的癌种数量
country_num=length(country_dict)
cancer_dict=matrix(data = 0, nrow = cancer_num, ncol = 3  )    
cancer_dict[,1]=1:cancer_num
row.names(cancer_dict)=c("bladder", "colon", "kindey", "liver", "luad", "lusc", "melanoma" , "oesophagus", "pancreas", "rectum", "stomach")
cancer_dict[,2]=c(165,42,160,59,79,78,100,24,70,49,35)
cancer_dict[,3]=c(165,42,160,59,79,78,100,27,70,49,35)
colnames(cancer_dict)=c("dict","ci5_start","ci5_end")   
#####"ci5_start","ci5_end"  指癌种在ci5 data里的起始序号。   注：目前只有oesophagus24:27有合并
###bladder    transitional cell carcinoma  Urothelial cells are transitional cells
#####Liver 只包括	Hepatocellular carcinoma   因为TCGA只有Hepatocellular carcinoma

for (j in 1:country_num) {
  eval(parse(text=paste(paste(country_dict[j],'_matrix',sep=''), '= matrix(data=0, ncol=length(t), nrow=cancer_num*2)')))
  eval(parse(text=paste('rownames(',paste(country_dict[j],'_matrix)',sep=''), '=c(rep("a",cancer_num*2))')))
  eval(parse(text=paste('colnames(',paste(country_dict[j],'_matrix)',sep=''), '=t')))
  for (i in 1:cancer_num) {
    #####male
    incidence_count=matrix(data = 0,ncol = 1,nrow = length(t))
    population_count=matrix(data = 0,ncol = 1,nrow = length(t))
    for (k in (cancer_dict[i,2]:cancer_dict[i,3])) {
      a=country_list[[j]][((k-1)*38+7):((k-1)*38+14),4]
      b=country_list[[j]][((k-1)*38+7):((k-1)*38+14),5]
      incidence_count=incidence_count+a
      population_count=population_count+b
    }
    eval(parse(text=paste(paste(country_dict[j],'_matrix',sep=''), '[2*i-1,]= incidence_count/population_count')))
    
    eval(parse(text=paste('rownames(',paste(country_dict[j],'_matrix)',sep=''), '[2*i-1]="', paste(rownames(cancer_dict)[i],'_male',sep=''),'"')))
    #####female
    incidence_count=matrix(data = 0,ncol = 1,nrow = length(t))
    population_count=matrix(data = 0,ncol = 1,nrow = length(t))
    for (k in (cancer_dict[i,2]:cancer_dict[i,3])) {
      a=country_list[[j]][((k-1)*38+19+7):((k-1)*38+19+14),4]
      b=country_list[[j]][((k-1)*38+19+7):((k-1)*38+19+14),5]
      incidence_count=incidence_count+a
      population_count=population_count+b
    }
    eval(parse(text=paste(paste(country_dict[j],'_matrix',sep=''), '[2*i,]= incidence_count/population_count')))
    eval(parse(text=paste('rownames(',paste(country_dict[j],'_matrix)',sep=''), '[2*i]="', paste(rownames(cancer_dict)[i],'_female',sep=''),'"')))
  }
eval(parse(text=paste(paste('country_cancer_list[["',country_dict[j],'"]]=',sep=''), paste(country_dict[j],'_matrix',sep=''))))
}


#### parallel  test

sex=c(rep(1,8),rep(0,8))
t_log=log10(rep(t,2))
delta=1  #等效（平行）检验边界
beta3_threshold=0.5
R2__threshold=0.9
country_cancer_test_list=list()


for (j in 1:country_num) {
  test_matrix=matrix(data=0,nrow = cancer_num ,ncol = 13)
  colnames(test_matrix)=c("beta3","SE_beta3","p_value","R2_beta3","unparallel_test","Equivalence_test","beta1","SE_beta1","R2_beta1","beta2","SE_beta2","R2_beta2","Gao_test")
  rownames(test_matrix)=rownames(cancer_dict)
  for (i in 1:cancer_num) {
    cancer_log=log10(rbind(country_cancer_list[[j]][2*i-1,],country_cancer_list[[j]][2*i,])+10^-100)   #+10^-100防止出现log(0)
    cancer_log=as.vector(t(cancer_log))
    fit=summary(lm(cancer_log~t_log+sex+sex*t_log))
    test_matrix[i,1]=fit[[4]][4,1] 
    test_matrix[i,2]=fit[[4]][4,2]
    test_matrix[i,3]=fit[[4]][4,4] 
    test_matrix[i,4]=fit[[8]]
    ######Equivalence_test
    if ((fit[[4]][4,1] + 1.96*fit[[4]][4,2] <= delta) & (fit[[4]][4,1] - 1.96*fit[[4]][4,2] >= -delta) & (abs(fit[[4]][4,1]) <= beta3_threshold) & fit[[8]] >= R2__threshold) {
      test_matrix[i,6]=1  #  1代表平行
    }
    else{
      test_matrix[i,6]=0
    }
    ######gao_test
    cancer_log_male=log10(country_cancer_list[[j]][2*i-1,]+10^-100)   #+10^-100防止出现log(0)
    beta1_test=summary(lm(cancer_log_male~log10(t)))
    cancer_log_female=log10(country_cancer_list[[j]][2*i,]+10^-100)   #+10^-100防止出现log(0)
    beta2_test=summary(lm(cancer_log_female~log10(t)))
    test_matrix[i,7]=beta1_test[[4]][2,1]
    test_matrix[i,8]=beta1_test[[4]][2,2]
    test_matrix[i,9]=beta1_test[[8]]
    test_matrix[i,10]=beta2_test[[4]][2,1]
    test_matrix[i,11]=beta2_test[[4]][2,2]
    test_matrix[i,12]=beta2_test[[8]]
    if ((test_matrix[i,7] >= test_matrix[i,10] - 1.96*test_matrix[i,11]) & (test_matrix[i,7] <= test_matrix[i,10] + 1.96*test_matrix[i,11]) 
        & (test_matrix[i,10] >= test_matrix[i,7] - 1.96*test_matrix[i,8]) & (test_matrix[i,10] <= test_matrix[i,7] + 1.96*test_matrix[i,8]) 
        & beta1_test[[8]]>=R2__threshold & beta2_test[[8]]>=R2__threshold) {
      test_matrix[i,13]=1
    }
  }
  test_matrix[which(test_matrix[,3]<=0.05),5]=1
  eval(parse(text=paste(paste('country_cancer_test_list[["',country_dict[j],'"]]=',sep=''), 'test_matrix')))
}

cancer_country_test_matrix=matrix(data=0,nrow =cancer_num*country_num ,ncol = 13)

for (k in 1:(cancer_num*country_num)) {
  if((k%%country_num)!=0){
    cancer_country_test_matrix[k,]=country_cancer_test_list[[k%%country_num]][ceiling(k/(country_num)),]
  }
  else{
    cancer_country_test_matrix[k,]=country_cancer_test_list[[country_num]][ceiling(k/(country_num)),]
  }
}
rownames(cancer_country_test_matrix)=paste(rep(rownames(cancer_dict),each=country_num),'_',rep(country_dict,cancer_num),sep='')
colnames(cancer_country_test_matrix)=colnames(test_matrix)
cancer_country_test_matrix[,3]=p.adjust(cancer_country_test_matrix[,3],method = "fdr")
cancer_country_test_matrix[which(cancer_country_test_matrix[,3]<=0.05),5]=1
cancer_country_test_matrix[which(cancer_country_test_matrix[,3]>0.05),5]=0
result_simple=cancer_country_test_matrix[,c(3,5,6,13)]
save(result_simple, country_cancer_list,country_cancer_test_list,cancer_country_test_matrix,country_dict,cancer_dict,t,t_log,delta,sex,beta3_threshold,cancer_num,country_num,
     file = "D:\\my_study\\Project\\cancer_risk\\data_git\\new_result.Rdata")
