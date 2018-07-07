library(biomaRt)
library(GenomicRanges)
library(dplyr)
library(MatchIt)

#Set up an gene annotation template to use
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes = getBM(attributes=c("hgnc_symbol","chromosome_name","start_position","end_position"), mart=mart)
genes = genes[genes[,1]!="" & genes[,2] %in% c(1:22,"X","Y"),]
xidx = which(genes[,2]=="X")
yidx = which(genes[,2]=="Y")
genes[xidx, 2] = 23
genes[yidx, 2] = 24
genes[,2] = sapply(genes[,2],as.integer)
genes = genes[order(genes[,3]),]
genes = genes[order(genes[,2]),]
colnames(genes) = c("GeneSymbol","Chr","Start","End")
genes_GR = makeGRangesFromDataFrame(genes,keep.extra.columns = TRUE)

threshold_up=0.848      ###segment mean threshold
threshold_low=-0.737   ##doi:  10.1038/nature06358
cancer_num=12   ###  关注的癌种数量
cancer_dict=matrix(data = 0, nrow = cancer_num, ncol = 2  )    
cancer_dict[,1]=1:cancer_num
cancer_dict[,2]=c("blca", "coad", "esca","kirc", "kirp","lihc", "luad", "lusc", "paad", "read","skcm","stad")

cancer_pvalue_list=list()
cancer_count=list()


for (j in 1:cancer_num) {
  print(j)
  
  eval(parse(text=paste('tcga=read.delim("D:/my_study/Project/cancer_risk/data_git/TCGA_cnv/seg_clinic/',cancer_dict[j,2],'.txt")',sep='')))
  eval(parse(text=paste('clinical=read.delim("D:/my_study/Project/cancer_risk/data_git/TCGA_cnv/seg_clinic/',cancer_dict[j,2],'.tsv")',sep='')))
  
  #psm
  psm=clinical[,c(2,4,6,12,13,14)]
  psm=psm[which(psm$race=="white"),]  #只保留白人
  psm$Group <- as.logical(psm$gender == 'male')
  match.it <- matchit(Group ~ age_at_diagnosis+tumor_stage+vital_status, data = psm, method="nearest", ratio=1)
  df.match <- match.data(match.it)[1:ncol(psm)]
  
  #filter
  tcga[,1]=substr(tcga[,1],1,12)
  tcga_psm=tcga[tcga$Sample %in% df.match$submitter_id,]
  tcga_psm_filtered=tcga_psm[which(tcga_psm[,6]< threshold_low | tcga_psm[,6]> threshold_up),]
  tcga_psm_filtered_GR = makeGRangesFromDataFrame(tcga_psm_filtered, keep.extra.columns = TRUE)
  #Overlap your regions with the reference dataset that you created
  hits = findOverlaps(genes_GR, tcga_psm_filtered_GR, type="within")
  result = cbind(tcga_psm_filtered[subjectHits(hits),],genes[queryHits(hits),])
  
  #合并clinical
  result_simple=result[,c(1,7,6)]
  clinical_simple=clinical[,c(2,4,6)]
  colnames(clinical_simple)[1]=colnames(result_simple)[1]
  join=left_join(result_simple,clinical_simple,by=colnames(clinical_simple)[1])
  
  #删除赋值0.扩增赋值1
  join[which(join[,3]<0),3]=0
  join[which(join[,3]>0),3]=1
  #去重复
  join=join[!duplicated(join),]   #去完全重复冗余
  #检查部分重复项是否有矛盾（即删除又扩增）
  rep=join[duplicated(join[,c(1,2)]),]
  rep_test=paste(rep[,1],'_',rep[,2],sep = '')
  test_rep=paste(join[,1],'_',join[,2],sep = '')
  rep_id=which( test_rep %in% rep_test)
  if(!isEmpty(rep_id)){
    join=join[-rep_id,]#去除矛盾项
  }
  
  #男女分组
  male=join[which(join$gender=="male" ),]
  female=join[which(join$gender=="female" ),]
  female_num=length(which(df.match$gender =="female"))
  male_num=length(which(df.match$gender=="male"))
  
  ####不合并del和amp
  female_count=count(female,GeneSymbol,Segment_Mean)
  male_count=count(male,GeneSymbol,Segment_Mean)
  female_count=cbind(female_count,0)   #补全未突变case数量
  for (i in 1:dim(female_count)[1]) {
    female_count[i,4]=female_num-sum(female_count[which(female_count[,1]==female_count[i,1]),3])
  }
  female_count[,1]=paste(female_count[,1],'_',female_count[,2],sep = '')
  female_count=female_count[,-2]
  
  male_count=cbind(male_count,0)
  for (i in 1:dim(male_count)[1]) {
    male_count[i,4]=male_num-sum(male_count[which(male_count[,1]==male_count[i,1]),3])
  }
  male_count[,1]=paste(male_count[,1],'_',male_count[,2],sep = '')
  male_count=male_count[,-2]
  
  white_count=full_join(female_count,male_count,by= colnames(male_count)[1])
  white_count=data.frame(white_count)
  white_count=white_count[!is.na(white_count[,1]),]
  white_count[is.na(white_count[,2]),2]=0                                         #na赋值
  white_count[is.na(white_count[,3]),3]=female_num
  white_count[is.na(white_count[,4]),4]=0
  white_count[is.na(white_count[,5]),5]=male_num
  
  colnames(white_count)=c("gene","female_1","female_0","male_1","male_0")
  ####################white test#######################
  p_value_white=matrix(nrow = dim(white_count)[1],ncol = 1)
  for(i in 1:nrow(white_count)){
    a=matrix(c(white_count[i,2:5]),nrow = 2,ncol = 2,byrow = 1)
    a=apply(a, 2, as.numeric)
    b=chisq.test(a)
    p_value_white[i,1]=b[[3]]
  }
  rownames(p_value_white)=white_count[,1]
  colnames(p_value_white)="p_value"
  adjsut_p_value_white=as.matrix(p.adjust(p_value_white),method = "fdr")
  rownames(adjsut_p_value_white)=rownames(p_value_white)   
  
  ####合并del和amp
  female_count_bind=count(female,GeneSymbol)
  male_count_bind=count(male,GeneSymbol)
  female_count_bind=cbind(female_count_bind,female_num-female_count_bind[,2])   #补全未突变case数量
  male_count_bind=cbind(male_count_bind,male_num-male_count_bind[,2])
  colnames(female_count_bind)=c("gene","female_1","female_0")
  colnames(male_count_bind)=c("gene","male_1","male_0")
  white_count_bind=full_join(female_count_bind,male_count_bind,by=colnames(male_count_bind)[1])
  white_count_bind=data.frame(white_count_bind)
  white_count_bind=white_count_bind[!is.na(white_count_bind[,1]),]
  white_count_bind[is.na(white_count_bind[,2]),2]=0                                         #na赋值
  white_count_bind[is.na(white_count_bind[,3]),3]=female_num
  white_count_bind[is.na(white_count_bind[,4]),4]=0
  white_count_bind[is.na(white_count_bind[,5]),5]=male_num
  
  p_value_white_bind=matrix(nrow = dim(white_count_bind)[1],ncol = 1)
  for(i in 1:nrow(white_count_bind)){
    a=matrix(c(white_count_bind[i,2:5]),nrow = 2,ncol = 2,byrow = 1)
    a=apply(a, 2, as.numeric)
    b=chisq.test(a)
    p_value_white_bind[i,1]=b[[3]]
  }
  rownames(p_value_white_bind)=white_count_bind[,1]
  colnames(p_value_white_bind)="p_value"
  adjsut_p_value_white_bind=as.matrix(p.adjust(p_value_white_bind),method = "fdr")
  rownames(adjsut_p_value_white_bind)=rownames(p_value_white_bind)
  
  eval(parse(text=paste('cancer_pvalue_list[["',cancer_dict[j,2],'_p_value_white"]]=p_value_white',sep=''))) 
  eval(parse(text=paste('cancer_pvalue_list[["',cancer_dict[j,2],'_adjsut_p_value_white"]]=adjsut_p_value_white',sep='')))
  eval(parse(text=paste('cancer_pvalue_list[["',cancer_dict[j,2],'_p_value_white_bind"]]=p_value_white_bind',sep=''))) 
  eval(parse(text=paste('cancer_pvalue_list[["',cancer_dict[j,2],'_adjsut_p_value_white_bind"]]=adjsut_p_value_white_bind',sep='')))
  eval(parse(text=paste('cancer_count[["',cancer_dict[j,2],'"]]=white_count',sep='')))
  eval(parse(text=paste('cancer_count[["',cancer_dict[j,2],'_bind"]]=white_count_bind',sep='')))
  
}

save(cancer_pvalue_list,cancer_count,threshold_low,threshold_up,cancer_dict,file="D:/my_study/Project/cancer_risk/data_git/TCGA_cnv/cnv_result.Rdata")
