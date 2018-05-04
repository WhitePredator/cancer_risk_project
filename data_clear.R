setwd("D:\\my_study\\Project\\cancer_risk\\data")

name = list.files("world_age_2003-2007")
dir = paste("D:\\my_study\\Project\\cancer_risk\\data\\world_age_2003-2007\\",name,sep="")
n = length(dir)
csv_list=list()
for (i in 1:n) {
  csv_list[[i]]=read.csv(file=dir[i],header = F)
  #eval(parse(text=paste(paste('a',i,sep=''), '= read.csv(file=dir[i],header = F)')))
}

registry <- read.delim2("D:/my_study/Project/cancer_risk/data/world_age_2003-2007/registry.txt", header=FALSE)

country_list=list()

canada_num=34
canada=csv_list[[canada_num]]
country_list[["canada"]]=canada

usa_white_num=209
usa_white=csv_list[[usa_white_num]]
country_list[["usa_white"]]=usa_white

usa_black_num=211
usa_black=csv_list[[usa_black_num]]
country_list[["usa_black"]]=usa_black

south_korea_num=253
south_korea=csv_list[[south_korea_num]]
country_list[["south_korea"]]=south_korea

england_num=392
england=csv_list[[england_num]]
country_list[["england"]]=england

scotland_num=401
scotland=csv_list[[scotland_num]]
country_list[["scotland"]]=scotland


three_col=csv_list[[1]][,1:3]

#argentina
argentina_num=9:12
argentina=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in argentina_num) {
  a=csv_list[[i]][,4:5]
  argentina=argentina+a
}
argentina=cbind(three_col,argentina)
country_list[["argentina"]]=argentina

#brazil
brazil_num=13:18
brazil=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in brazil_num) {
  a=csv_list[[i]][,4:5]
  brazil=brazil+a
}
brazil=cbind(three_col,brazil)
country_list[["brazil"]]=brazil

#china
china_num=215:228
china=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in china_num) {
  a=csv_list[[i]][,4:5]
  china=china+a
}
china=cbind(three_col,china)
country_list[["china"]]=china

#japan
japan_num=245:252
japan=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in japan_num) {
  a=csv_list[[i]][,4:5]
  japan=japan+a
}
japan=cbind(three_col,japan)
country_list[["japan"]]=japan


#india
india_num=229:240
india=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in india_num) {
  a=csv_list[[i]][,4:5]
  india=india+a
}
india=cbind(three_col,india)
country_list[["india"]]=india

#france
france_num=299:309
france=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in france_num) {
  a=csv_list[[i]][,4:5]
  france=france+a
}
france=cbind(three_col,france)
country_list[["france"]]=france

#germany
germany_num=310:318
germany=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in germany_num) {
  a=csv_list[[i]][,4:5]
  germany=germany+a
}
germany=cbind(three_col,germany)
country_list[["germany"]]=germany

#italy
italy_num=321:353
italy=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in italy_num) {
  a=csv_list[[i]][,4:5]
  italy=italy+a
}
italy=cbind(three_col,italy)
country_list[["italy"]]=italy

#spain
spain_num=368:380
spain=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in spain_num) {
  a=csv_list[[i]][,4:5]
  spain=spain+a
}
spain=cbind(three_col,spain)
country_list[["spain"]]=spain

#switzerland
switzerland_num=382:390
switzerland=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in switzerland_num) {
  a=csv_list[[i]][,4:5]
  switzerland=switzerland+a
}
switzerland=cbind(three_col,switzerland)
country_list[["switzerland"]]=switzerland

#australia
australia_num=405:413
australia=matrix(0,nrow =2*244*19 ,ncol = 2)
for (i in australia_num) {
  a=csv_list[[i]][,4:5]
  australia=australia+a
}
australia=cbind(three_col,australia)
country_list[["australia"]]=australia

country_dict=matrix(data = NA,nrow =length(country_list) ,ncol = 1)
for (i in 1:length(country_list)) {
  country_dict[i,1]=names(country_list[i])
}

save(country_dict,country_list,file = "D:\\my_study\\Project\\cancer_risk\\data\\country.Rdata")
