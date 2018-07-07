library(dplyr)
library(MatchIt)

cancer_num=12   ###  关注的癌种数量
cancer_dict=matrix(data = 0, nrow = cancer_num, ncol = 2  )    
cancer_dict[,1]=1:cancer_num
cancer_dict[,2]=c("blca", "coad", "esca","kirc", "kirp","lihc", "luad", "lusc", "paad", "read","skcm","stad")

clinical_all=read.delim("D:\\my_study\\Project\\cancer_risk\\data_git\\clinical\\clinical_PANCAN_patient_with_followup.tsv", comment.char="#")
clinical_1=clinical_all[,c(2,4,10,20,29,32,33,34,35,39,56,57,78,79,80)]
clinical_2=clinical_1[,c(1,5,13,14)]
tumor_purity = read_excel("D:/my_study/Project/cancer_risk/data_git/clinical/tumor_purity.xlsx")
#去除同一样本a b 重复，取平均
purity_simple=tumor_purity[,c(1,7)]
purity_simple_substr=purity_simple
for (i in 1:nrow(purity_simple)) {
  purity_simple_substr[i,1]=substr(purity_simple[i,1],1,12)
}
purity_simple_substr$CPE=as.numeric(purity_simple_substr$CPE)

index=duplicated(purity_simple_substr[,1])
for (i in 1:(length(index)-1)) {
  if (!is.na( purity_simple_substr[i+1,2])) {
    purity_simple_substr[i,2]=(purity_simple_substr[i,2]+purity_simple_substr[i+1,2]*index[i+1])/(1+1*index[i+1])
  }
}
purity_simple_substr=purity_simple_substr[!duplicated(purity_simple_substr[,1]),]
colnames(purity_simple_substr)[2]="purity"
colnames(purity_simple_substr)[1]="id"


maf_pvalue_list=list()
maf_adjust_pvalue_list=list()
white_count_list=list()
gene_count_list=list()

for (j in 1:cancer_num) {
  print(cancer_dict[j,2])
  
  eval(parse(text=paste('tcga=read.delim("D:/my_study/Project/cancer_risk/data_git/TCGA_maf/',cancer_dict[j,2],'/maf/muse.maf", comment.char="#")',sep='')))
  eval(parse(text=paste('clinical=read.delim("D:/my_study/Project/cancer_risk/data_git/clinical/',cancer_dict[j,2],'.tsv", comment.char="#")',sep='')))
  clinical=clinical[,c(2,4,6,12,13,14)]
  sample_num=dim(clinical)[1]
  maf=tcga[,c(1,9,16)]
  
  ###maf_full=maf_full[-which(maf_full[,3]==""),]    如果id有空值
  
  maf_no_slient=maf[-c(which(maf[,2]=="Silent")),]     #去掉slient mutation
  #maf_no_slient=maf_no_slient[-c(which(maf_no_silient[,2]=="RNA")),] 
  #maf_no_slient=maf_no_slient[-c(which(maf_no_slient[,2]=="Intron")),]
  #maf_no_slient=maf_no_slient[-c(which(maf_no_slient[,2]=="IGR")),]
  #maf_join=maf_join[complete.cases(join),]
  
  sub_str=maf_no_slient
  sub_str=as.matrix(sub_str)    #加快下一步substr速度
  for(i in 1:nrow(maf_no_slient)){
    sub_str[i,3]=substr(maf_no_slient[i,3],1,12)
    }

  sub_str=data.frame(sub_str)
  clinical_copy=data.frame(clinical)
  colnames(clinical_2)[1]="id"
  colnames(sub_str)[3]="id"
  colnames(clinical_copy)[1]="id"
  clinical_join=left_join(clinical_copy,clinical_2,by="id")
  clinical_join=left_join(clinical_join,purity_simple_substr,by="id")
  
  
  clinical_num=dim(clinical)[1]
  ##去掉突变频率不足5%的基因
  sub_str_dup=sub_str[!duplicated(sub_str[,c(1,3)]),]  
  gene_count=count(sub_str_dup,Hugo_Symbol)
  gene_count_005=gene_count[which(gene_count[,2]>=clinical_num*0.05),]
  maf_005=sub_str_dup[which(sub_str_dup[,1]  %in% gene_count_005$Hugo_Symbol),]
  
  ##去掉突变数大于1000的样本
  sample_mutation_count=count(sub_str,id)
  sample_mutaion_over=sample_mutation_count[which(sample_mutation_count[,2]>=1000),]
  if (dim(sample_mutaion_over)[1]!=0) {
    clinical_1000=clinical_join[-which(clinical_join[,1]  %in% sample_mutaion_over$id),] 
  }else{clinical_1000=clinical_join}
  
  
  
  
  ### psm
  psm=clinical_1000
  psm[is.na(psm[,10]),10]=0   #na 用0表示
  psm=psm[which(psm$race=="white"),]  #只保留白人
  psm$Group <- as.logical(psm$gender == 'male')
  if (length(unique(psm[,7]))==1 & length(unique(psm[,8]))==1 & length(unique(psm[,9]))==1) {
    match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status+purity,data = psm, method="nearest", ratio=1)
  }else if(length(unique(psm[,7]))==1 & length(unique(psm[,8]))!=1 & length(unique(psm[,9]))!=1){
    match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status+tobacco_smoking_history+number_pack_years_smoked+purity,
                        data = psm, method="nearest", ratio=1)
  }else if(length(unique(psm[,7]))==1 & length(unique(psm[,8]))==1 & length(unique(psm[,9]))!=1){
    match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status+number_pack_years_smoked+purity,
                        data = psm, method="nearest", ratio=1)
  }else if(length(unique(psm[,7]))==1 & length(unique(psm[,8]))!=1 & length(unique(psm[,9]))==1){
    match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status+tobacco_smoking_history+purity,
                        data = psm, method="nearest", ratio=1)
  }else if(length(unique(psm[,7]))!=1 & length(unique(psm[,8]))==1 & length(unique(psm[,9]))==1){
    match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status+histological_type+purity,
                        data = psm, method="nearest", ratio=1)
  }else if(length(unique(psm[,7]))!=1 & length(unique(psm[,8]))!=1 & length(unique(psm[,9]))!=1){
    match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status+histological_type+tobacco_smoking_history+number_pack_years_smoked+purity,
                        data = psm, method="nearest", ratio=1)
  }else if(length(unique(psm[,7]))!=1 & length(unique(psm[,8]))==1 & length(unique(psm[,9]))!=1){
    match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status+histological_type+number_pack_years_smoked+purity,
                        data = psm, method="nearest", ratio=1)
  }else if(length(unique(psm[,7]))!=1 & length(unique(psm[,8]))!=1 & length(unique(psm[,9]))==1){
    match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status+histological_type+tobacco_smoking_history+purity,
                        data = psm, method="nearest", ratio=1)
    }
   
  
  df.match <- match.data(match.it)[1:ncol(psm)]
  sample_num=dim(df.match)[1]
    
  ##clinical和maf join
  maf_join=left_join(df.match,maf_005,by="id")
  maf_join_duplicated=maf_join[!duplicated(maf_join[,c(1,12)]),]   #去冗余
  
  
  
  
  ######################分组####################
  male_white=maf_join_duplicated[which(maf_join_duplicated$gender=="male" ),]
  female_white=maf_join_duplicated[which(maf_join_duplicated$gender=="female"),]
  
  female_white_num=length(unique(female_white[,1]))
  male_white_num=length(unique(male_white[,1]))
  
  ##############white count#################
  female_white_count=count(female_white,Hugo_Symbol)
  male_white_count=count(male_white,Hugo_Symbol)
  female_white_count=cbind(female_white_count,female_white_num-female_white_count[,2])   #补全未突变case数量
  male_white_count=cbind(male_white_count,male_white_num-male_white_count[,2])
  colnames(female_white_count)=c("gene","female_1","female_0")
  colnames(male_white_count)=c("gene","male_1","male_0")
  white_count=full_join(female_white_count,male_white_count,by="gene")
  white_count=data.frame(white_count)
  white_count[is.na(white_count[,2]),2]=0                                         #na赋值
  white_count[is.na(white_count[,3]),3]=female_white_num
  white_count[is.na(white_count[,4]),4]=0
  white_count[is.na(white_count[,5]),5]=male_white_num
  
  ####################white test#######################
  white_count=as.matrix(white_count)
  pvalue_matrix=matrix(0,nrow =nrow(white_count) ,ncol = 2)
  colnames(pvalue_matrix)=c("gene","p_value")
  pvalue_matrix[,1]=white_count[,1]
  adjust_pvalue_matrix=pvalue_matrix
  
  for(i in 1:nrow(white_count)){
    a=matrix(c(white_count[i,2:5]),nrow = 2,ncol = 2,byrow = 1)
    a=apply(a, 2, as.numeric)
    b=chisq.test(a)
    pvalue_matrix[i,2]=b[[3]]
  }
  
  adjust_pvalue_matrix[,2]=as.matrix(p.adjust(pvalue_matrix[,2]),method = "fdr")

  eval(parse(text=paste('maf_pvalue_list[["',cancer_dict[j,2],'"]]=pvalue_matrix',sep=''))) 
  eval(parse(text=paste('maf_adjust_pvalue_list[["',cancer_dict[j,2],'"]]=adjust_pvalue_matrix',sep='')))
  eval(parse(text=paste('white_count_list[["',cancer_dict[j,2],'"]]=white_count',sep='')))
  eval(parse(text=paste('gene_count_list[["',cancer_dict[j,2],'"]]=gene_count',sep='')))
   
  
}
for (i in 1:dim(cancer_dict)[1]) {
  a=maf_adjust_pvalue_list [[i]]
  eval(parse(text=paste(cancer_dict[i,2],'=a',sep='')))
}

#######################save########################
save(maf_pvalue_list,maf_adjust_pvalue_list,white_count_list,gene_count_list,cancer_dict,blca, coad, esca,kirc, kirp,lihc, luad, lusc, paad, read,skcm,stad,
     file="D:/my_study/Project/cancer_risk/data_git/TCGA_maf/maf_psm_result.Rdata")
#save.image(file = "D:\\my_study\\Project\\gao_group\\TCGA\\colorectal\\rectal\\result.RData")


