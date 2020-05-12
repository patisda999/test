library(DMwR)
library(pls)
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
df =read.csv("FHS_strk-4.csv")
x_cols = c('AGE','SEX','SYSBP','DIABP','BPMEDS','CURSMOKE','CIGPDAY','TOTCHOL','HDLC','LDLC','BMI','GLUCOSE','HEARTRTE')
x_df = df[,x_cols]
x_df = knnImputation(x_df)
y_cols = c('DIABETES','PREVAP','PREVCHD','PREVMI','PREVHYP')
y_df = df[,y_cols]



result = pls_part1(x_df,y_df)

pls<-pls_part2(9,result$pls_overview)
weight=plstep2_weighting(pls)
xloading=pls$loadings[]
yloading=pls$Yloadings[]
key=rep(F,length(pls_res$weight[,1]))



weight_hyp=weight[order(weight[,"PREVHYP"],decreasing = T),]
weight_dia=weight[order(weight[,"DIABETES"],decreasing = T),]

