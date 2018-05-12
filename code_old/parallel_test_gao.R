#parallel_test
#bladder

parallel_test_bladder=matrix(data=0,nrow = length(bladder_lm_list)/2,ncol = 2)

for (i in 1:(length(bladder_lm_list)/2)) {
  beta1=bladder_lm_list[[2*i-1]][4][[1]][2,1]
  beta1_se=bladder_lm_list[[2*i-1]][4][[1]][2,2]
  beta2=bladder_lm_list[[2*i]][4][[1]][2,1]
  beta2_se=bladder_lm_list[[2*i]][4][[1]][2,2]
  if ((beta1 >= beta2 - 1.96*beta2_se) & (beta1 <= beta2 + 1.96*beta2_se) & (beta2 >= beta1 - 1.96*beta1_se) & (beta2 <= beta1 + 1.96*beta1_se)) {
    parallel_test=1
  }
  else{
    parallel_test=0
  }
  parallel_test_bladder[i,1]=substr(paste(names(bladder_lm_list[2*i-1])),14,50)
  parallel_test_bladder[i,2]=parallel_test
}


parallel_test_melanoma=matrix(data=0,nrow = length(melanoma_lm_list)/2,ncol = 2)

for (i in 1:(length(melanoma_lm_list)/2)) {
  beta1=melanoma_lm_list[[2*i-1]][4][[1]][2,1]
  beta1_se=melanoma_lm_list[[2*i-1]][4][[1]][2,2]
  beta2=melanoma_lm_list[[2*i]][4][[1]][2,1]
  beta2_se=melanoma_lm_list[[2*i]][4][[1]][2,2]
  if ((beta1 >= beta2 - 1.96*beta2_se) & (beta1 <= beta2 + 1.96*beta2_se) & (beta2 >= beta1 - 1.96*beta1_se) & (beta2 <= beta1 + 1.96*beta1_se)) {
    parallel_test=1
  }
  else{
    parallel_test=0
  }
  parallel_test_melanoma[i,1]=substr(paste(names(melanoma_lm_list[2*i-1])),14,50)
  parallel_test_melanoma[i,2]=parallel_test
}
