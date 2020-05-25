library(DMwR)
library(pls)
library(stats)
library(MASS)
pls_part1 = function(x,y){
  #x,y定義
  x<-as.data.frame(x)
  
  y<-as.data.frame(y)
  
  #normalize x,y
  x_scale<-as.matrix(scale(x))
  
  y_scale<-as.matrix(scale(y))
  
  #ncomp pls
  pls_overview<-plsr(y_scale~x_scale)
  
  #取出 R2 跟 X variance
  # r2<-matrix(
  #   R2(pls_overview, estimate="train",
  #      intercept = FALSE)$val,nrow=max(col(x)),byrow = T)
  r2<-matrix(
    R2(pls_overview, estimate="train",
       intercept = FALSE)$val,nrow=max(pls_overview$ncomp),byrow = T)
  
  Xv<-as.data.frame(cumsum(explvar(pls_overview)/100))
  
  #修改名稱
  
  row.names(r2)<-row.names(Xv)
  colnames(r2)=colnames(y)
  # rsq<-`colnames<-`(r2,colnames(y))
  
  names(Xv)<-"X variance"
  
  R2_df<-data.frame(Xv,r2)
  
  outlist<-list("R2_df"=R2_df,"pls_overview"=pls_overview)
  #製作總覽表格
  return(outlist)
  
}
pls_part2 =function(ncomp,nruns){
  plsncomp_res<-nruns
  #把plstep1結果拿進來
  pls_res=list()
  
  pls_res$loadings<-as.data.frame(plsncomp_res$loadings[,1:ncomp])
  if(length(colnames(plsncomp_res$model$y_scale))>1){
    pls_res$Yloadings<-t(plsncomp_res$Yloadings[,1:ncomp])
  }else{
    pls_res$Yloadings<-plsncomp_res$Yloadings[,1:ncomp]
  }
  pls_res$model<-as.data.frame(plsncomp_res$model)
  return(pls_res)
}
pls_weighting = function(pls){
  res<-pls
  #x,y loading 取絕對值
  xloadings_abs<-as.matrix(abs(res$loadings[]))
  yloadings_abs<-as.matrix(abs(res$Yloadings[]))
  colnames(yloadings_abs)<-colnames(res$model$y_scale)
  
  x_los<-t(scale(t(xloadings_abs)))
  
  y_los<-scale(yloadings_abs)
  x_los[is.nan(x_los)]<-0
  y_los[is.nan(y_los)]<-0
  weight<-x_los%*%y_los/(ncol(x_los)-1)
  #weight<-xloadings_abs%*%t(yloadings_abs)
  wwe<-as.data.frame(weight)
  weight_t<-data.frame("names"=rownames(wwe),wwe)
  
  weight_sort<-weight_t[order(weight_t[,2],decreasing = T),]
  
  
  return(weight_sort)
}
Tsqr = function(dat,matrixcov){
  apply(dat, 1, function(i) { i %*% matrixcov %*% i })
}

df =read.csv("FHS_strk-4.csv")
x_cols = c('AGE','SEX','SYSBP','DIABP','BPMEDS','CURSMOKE','CIGPDAY','TOTCHOL','HDLC','LDLC','BMI','GLUCOSE','HEARTRTE')
x_df = df[,x_cols]
x_df = knnImputation(x_df)
y_cols = c('DIABETES')
y_df = df[,y_cols]
n =length(x_cols)
row = list()
p_array=c()
count = 1
for(i in 1:(n-2)){
  for(j in (i+1):(n-1)){
    for(k in (j+1):n){
      data = df[,x_cols[c(i,j,k)]]
      na_df = is.na(data)
      temp = (na_df[,1]+na_df[,2]+na_df[,3])>0
      data=data[!temp,]
      data=scale(data)
      clus = kmeans(data,3)$cluster
      anova_df = data.frame(x = clus,DIABETES = df$DIABETES[!temp])
      p_value = anova(lm(DIABETES~x,data = anova_df))$`Pr(>F)`
      row[[count]] = x_cols[c(i,j,k)]
      p_array=c(p_array,p_value[1])   
      count=count+1
    }
  }
}


result =table(unlist( row[which(p_array<0.05)]))
result=sort(result,decreasing = T)
result

#A組
data_scale<-scale(x_df[,names(result)[1:3]])
corr<-cor(data_scale,use = "complete.obs")
E<-eigen(corr)
score<-data.frame(data_scale%*%E$vectors)
evs<-round(E$values,digits = 5)
eigentable=data.frame("a"=evs,"b"=round(cumsum(evs)/sum(evs)*100,5))
names(eigentable)=c("特徵值","累積解釋量(%)")
row.names(eigentable)=paste("PC",1:length(eigentable[,1]),sep="")
eigentable

#A組
dat=scale(as.matrix(score[,1,drop=F]),T,F)
varcov=cov(dat)
matrixcov=ginv(varcov)
tscore=Tsqr(dat = dat,matrixcov = matrixcov)
mean_tscore = mean(tscore)
anova_df_A = data.frame(range = as.numeric(tscore >mean_tscore)+1,DIABETES = df$DIABETES)
anova(lm(DIABETES~range,data = anova_df_A))

#B組
data_scale<-scale(x_df[,names(result)[11:13]])
corr<-cor(data_scale,use = "complete.obs")
E<-eigen(corr)
score<-data.frame(data_scale%*%E$vectors)
evs<-round(E$values,digits = 5)
eigentable=data.frame("a"=evs,"b"=round(cumsum(evs)/sum(evs)*100,5))
names(eigentable)=c("特徵值","累積解釋量(%)")
row.names(eigentable)=paste("PC",1:length(eigentable[,1]),sep="")
eigentable
#B組
dat=scale(as.matrix(score[,1:2,drop=F]),T,F)
varcov=cov(dat)
matrixcov=ginv(varcov)
tscore=Tsqr(dat = dat,matrixcov = matrixcov)
mean_tscore = mean(tscore)
anova_df_B = data.frame(range = as.numeric(tscore >mean_tscore)+1,DIABETES = df$DIABETES)
anova(lm(DIABETES~range,data = anova_df_B))

# 參數說明  name 為欄位名稱  df 為資料集 type 連續選con 類別選class
single_fac_info = function(name,df,type = c("con","class")){
  data = df[,name]
  data = data[!is.na(data)]
  par(mfcol=c(2,1))
  switch(type,
         "con" = {
           temp = table(data)
           info_data = list(
           ma = mean(data), #平均
           med = median(data), #中位
           mode = names(temp[temp==max(temp)]), #眾數
           sd_data = sd(data), #標準差
           var_data = var(data), #變異數
           hi = hist(data,xlab = name), #直方
           scatter = plot(data,ylab = name) #散佈
           )
         },
         "class" = {
           temp = table(data)
           info_data = list(
           med = median(data), #中位
           iqr = IQR(data), #4分位距
           mode = names(temp[temp==max(temp)]), #眾數
           max_dis = max(data)-min(data), #全距
           box = boxplot(data,ylab = name), #盒鬚圖
           bar = hist(data,ylab = name)) #長條圖
         })
return(info_data)
}

single_fac_info("SYSBP",df,type="con") #高血壓
single_fac_info("DEATH",df,type="class") #死亡與否
