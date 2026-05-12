library(tidyverse)
library(deSolve)
library(rootSolve)
library(geomtextpath)
library(metR)

# Functions
malemimic<-function(t,state,params){
  
  tm<-params["tm"]
  tf<-params["tf"]
  O<-params["O"]
  M<-params["M"]
  K<-params["K"]
  muM<-params["muM"]
  muF<-params["muF"]
  Pm<-params["Pm"]
  Pf<-params["Pf"]
  
  Cm<-1-abs(tf-muM)
  Cf<-1-abs(tm-muF)
  
  Nm<-state["Nm"]
  Nf<-state["Nf"]
  
  N<-Nm+Nf
  
  Bm<-Nm*Nf*Cm*Cf*O*M*(1-N/K)
  Bf<-Nm*Nf*Cm*Cf*O*(1-M)*(1-N/K)
  
  Dm<-Nm*(Pm*Nm+Pf*Nf*(1-abs(tm-tf)))/(Nm+Nf*(1-abs(tm-tf)))
  Df<-Nf*(Pf*Nf+Pm*Nm*(1-abs(tm-tf)))/(Nf+Nm*(1-abs(tm-tf)))
  
  dNm<-Bm-Dm
  dNf<-Bf-Df
  
  list(c(dNm,dNf))
}

# If reproduction is limited by # of females (e.g., when females can only mate infrequently)
# FLR stands for female-limited reproduction
malemimic_FLR<-function(t,state,params){
  
  tm<-params["tm"]
  tf<-params["tf"]
  O<-params["O"]
  M<-params["M"]
  K<-params["K"]
  muM<-params["muM"]
  muF<-params["muF"]
  Pm<-params["Pm"]
  Pf<-params["Pf"]
  
  Cm<-1-abs(tf-muM)
  Cf<-1-abs(tm-muF)
  
  Nm<-state["Nm"]
  Nf<-state["Nf"]
  
  N<-Nm+Nf
  
  Bm<-Nf*Cm*Cf*O*M*(1-N/K)
  Bf<-Nf*Cm*Cf*O*(1-M)*(1-N/K)
  
  Dm<-Nm*(Pm*Nm+Pf*Nf*(1-abs(tm-tf)))/(Nm+Nf*(1-abs(tm-tf)))
  Df<-Nf*(Pf*Nf+Pm*Nm*(1-abs(tm-tf)))/(Nf+Nm*(1-abs(tm-tf)))
  
  dNm<-Bm-Dm
  dNf<-Bf-Df
  
  list(c(dNm,dNf))
}


# value<-0.8
# lag<-0.05
# 
# time<-seq(0,1000,by = 0.1)
# params<-c(tm=0,tf=1,O=2,M=0.9,K=10000,muM=1,muF=0,Pm=0.0578,Pf=0.289)
# test<-rootSolve::runsteady(y = state,func = malemimic_FLR,parms = params,times = c(0,1000))
# test$y
# round(test$y[1]+test$y[2],0)

state<-c(Nm=10,Nf=10)
  
tf<-seq(0,1,length.out=50)
# lag<-seq(0,0.1,length.out=10)
M<-seq(0.1,0.9,length.out=10)
Pm<-seq(0.01,0.4,length.out=50)
Pf<-Pm

all_combo<-expand.grid(tf,M,Pm,Pf)
parameter.space<-as.data.frame(matrix(unlist(all_combo[all_combo$Var4>all_combo$Var3,]),ncol = 4))

data<-matrix(-99,nrow = nrow(parameter.space),ncol = 3)

# simulations when females mate only once
for (i in 1:nrow(parameter.space)){
  
  params<-c(tm=0,tf=parameter.space[i,1],O=2,M=parameter.space[i,2],K=10000,muM=parameter.space[i,1],muF=0,Pm=parameter.space[i,3],Pf=parameter.space[i,4])
  test<-rootSolve::runsteady(y = state,func = malemimic_FLR,parms = params,times = c(0,1000))
  data[i,]<-c(test$y,sum(test$y))
  
}

# simulations when females can mate multiple times
for (i in 1:nrow(parameter.space)){
  
  params<-c(tm=0,tf=parameter.space[i,1],O=2,M=parameter.space[i,2],K=10000,muM=parameter.space[i,1],muF=0,Pm=parameter.space[i,3],Pf=parameter.space[i,4])
  test<-rootSolve::runsteady(y = state,func = malemimic_FLR,parms = params,times = c(0,1000))
  data[i,]<-c(test$y,sum(test$y))
  
}

df<-cbind(parameter.space,data)
colnames(df)<-c("tf","M","Pm","Pf","Nm","Nf","N")

write_csv(df,"df_notFLR.csv")

# df_lag<-read_csv("df_lag.csv")
# df_lag0<-subset(df_lag,subset = df_lag$lag==0)
# 
# df_lag0$M<-0.5
# 
# df2<-merge(df,df_lag0,all = TRUE)

df2<-read_csv("df_notFLR.csv")

df2$MPmPf<-paste0(df2$M,df2$Pm,df2$Pf)

combo<-unique(df2$MPmPf)

df2$w<-100

for (i in 1:length(combo)){
  
  
  df2[df2$MPmPf==combo[i],]$w<-df2[df2$MPmPf==combo[i],]$N/max(df2[df2$MPmPf==combo[i],]$N)
  
}

write_csv(df2,"df_notFLR.csv")
# Select results from representative Pm values
df2<-read_csv("df_FLR.csv")

Pm_level<-unique(df2$Pm)[c(1:10,20,30)]
df_plot<-subset(df2,subset = df2$Pm%in%Pm_level)

df_plot[df_plot$N<2,]$w<-0

# Select different M values
M_level<-unique(df_plot$M)
df_plot2<-subset(df_plot,subset = df_plot$M==M_level[10])

Lab<-paste0("Pm = ",round(Pm_level,digits = 3))

df_plot2$Lab<-paste0("Pm = ",round(df_plot2$Pm,digits = 3))
df_plot2$sexratio<-df_plot2$Nf/df_plot2$Nm
df_plot2$W<-df_plot2$N/10000


breaks1<-round(seq(min(df_plot2$sexratio),max(df_plot2$sexratio),by = 0.1*(max(df_plot2$sexratio)-min(df_plot2$sexratio))),2)
# Heat map of relative fitness calculated for each Pm-Pf combo
ggplot(data = df_plot2,aes(x = Pf,y = tf))+
  geom_raster(aes(fill = w))+
  scale_fill_gradientn(aes(fill = w),colors = c(viridisLite::magma(8),"white"),values = c(0,0.8,0.9,0.95,0.99,0.995,0.999,1))+
  geom_contour(aes(z = sexratio),color = "blue",linetype=2,linewidth = 0.2,breaks = breaks1)+
  geom_text_contour(aes(z = sexratio),colour = "blue",stroke=0.2,size=3,alpha=1,breaks = breaks1,label.placer = label_placer_fraction(frac=0.5))+
  facet_wrap(~Lab,ncol=4)+
  theme_bw()

# breaks2<-round(seq(min(df_plot2$W),max(df_plot2$W),by = 0.1*(max(df_plot2$W)-min(df_plot2$W))),2)
breaks2<-round(seq(min(df_plot2$W),max(df_plot2$W),by = 0.1*(max(df_plot2$W)-min(df_plot2$W))),6)
# Heat map of absolute fitness expressed as ratio to K across all conditions
ggplot(data = df_plot2,aes(x = Pf,y = tf))+
  geom_raster(aes(fill = W))+
  scale_fill_viridis_c(aes(fill = W),option = "magma")+
  # geom_contour(aes(z = w),color = "blue",linetype=2,linewidth = 0.2,breaks = c(0,0.7,0.8,0.9,0.95,0.99,0.995,0.999))+
  # geom_contour(aes(z = W),color = "black",linetype=1,linewidth = 0.2)+
  geom_contour(aes(z = W),color = "blue",linetype = 2,linewidth=0.2,breaks=breaks2)+
  # geom_text_contour(aes(z = W),color = "blue",size=2.5,stroke=0.2,alpha=1,breaks=breaks2,label.placer = label_placer_fraction(frac=0.5))+
  facet_wrap(~Lab,ncol=4)+
  theme_bw()

# Supplemental fig showing the effect of lag (mismatch between males preference and female trait)

df_lag<-read_csv("df_lag.csv")

df_lag$lagPmPf<-paste0(df_lag$lag,df_lag$Pm,df_lag$Pf)

combo<-unique(df_lag$lagPmPf)

df_lag$w<-100

for (i in 1:length(combo)){
  
  
  df_lag[df_lag$lagPmPf==combo[i],]$w<-df_lag[df_lag$lagPmPf==combo[i],]$N/max(df_lag[df_lag$lagPmPf==combo[i],]$N)
  
}

write_csv(df_lag,"df_lag.csv")

lag_level<-unique(df_lag$lag)[c(2,4,8,10)]

df_lag_plot<-subset(df_lag,subset = df_lag$lag%in%lag_level&df_lag$Pm%in%unique(df_lag$Pm)[c(1,4,10)])
df_lag_plot$W<-df_lag_plot$N/10000
df_lag_plot$sexratio<-df_lag_plot$Nf/df_lag_plot$Nm

breaks3<-round(seq(min(df_lag_plot$sexratio),max(df_lag_plot$sexratio),by = 0.1*(max(df_lag_plot$sexratio)-min(df_lag_plot$sexratio))),2)

# Heat map of absolute fitness expressed as ratio to K across all conditions
ggplot(data = df_lag_plot,aes(x = Pf,y = tf))+
  geom_raster(aes(fill = w))+
  scale_fill_gradientn(aes(fill = w),colors = c(viridisLite::magma(8),"white"),values = c(0,0.8,0.9,0.95,0.99,0.995,0.999,1))+
  geom_contour(aes(z = sexratio),color = "blue",linetype=2,linewidth = 0.2,breaks = breaks3)+
  geom_text_contour(aes(z = sexratio),colour = "blue",stroke=0.2,size=3,alpha=1,breaks = breaks3,label.placer = label_placer_fraction(frac=0.5))+
  facet_wrap(~lag*Pm,ncol=3)+
  theme_bw()

breaks4<-round(seq(min(df_lag_plot$W),max(df_lag_plot$W),by = 0.1*(max(df_lag_plot$W)-min(df_lag_plot$W))),6)

ggplot(data = df_lag_plot,aes(x = Pf,y = tf))+
  geom_raster(aes(fill = W))+
  scale_fill_viridis_c(aes(fill = W),option = "magma")+
  # geom_contour(aes(z = w),color = "blue",linetype=2,linewidth = 0.2,breaks = c(0,0.7,0.8,0.9,0.95,0.99,0.995,0.999))+
  # geom_contour(aes(z = W),color = "black",linetype=1,linewidth = 0.2)+
  geom_contour(aes(z = W),color = "blue",linetype = 2,linewidth=0.2,breaks=breaks2)+
  # geom_text_contour(aes(z = W),color = "blue",size=2.5,stroke=0.2,alpha=1,breaks=breaks2,label.placer = label_placer_fraction(frac=0.5))+
  facet_wrap(~lag*Pm,ncol=3)+
  theme_bw()
