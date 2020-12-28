library(sqldf)
library(ggplot2)
library(splines)
library(plot3D)
library(rgl)
library(RColorBrewer)
library(data.table)
library(plm)
library(astsa)
library("gridExtra")
library("ggpubr")
library(plyr)
library(RColorBrewer)
library(grDevices)
library(gridExtra)


# B sline basis function
BB26=bs(data$tij/max(data$tij),degree=3,df=80,intercept = F)
tg4=seq(0,1,length.out = 200)
bb26=bs(tg4,degree=3,df=80,intercept = F)
amplitude = array(NA, dim = 80)
for (i in 1:80){
  amplitude[i] = length(which(bb26[,i] != 0)) / 200
}

# Import regression coefficients
Beta=read.table("ExternalBeta.txt",header=F)
as.data.frame(Beta)
Beta$iter=seq(1:1599)
Beta$V1

# Trace plots for regression coefficients
x = {ggplot()+geom_line(data=Beta,aes(x = iter, y = V1),alpha=1,size=0.3)+
    ggtitle(substitute("Regression Coeficient Sex"))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
    theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
    ylab("Value")+ theme(legend.position = "none")+
    geom_hline(aes(yintercept=mean(Beta$V1),alpha=0.2,col = "red"))+
    scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")}

y = {ggplot()+geom_line(data=Beta,aes(x = iter, y = V2),alpha=1,size=0.3)+
    ggtitle(substitute("Regression Coeficient Age"))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
    theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
    ylab("Value")+ theme(legend.position = "none")+
    geom_hline(aes(yintercept=mean(Beta$V2),alpha=0.2,col = "red"))+
    scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")}

z = {ggplot()+geom_line(data=Beta,aes(x = iter, y = V3),alpha=1,size=0.3)+
    ggtitle(substitute("Regression Coeficient Environment"))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
    theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
    ylab("Value")+ theme(legend.position = "none")+
    geom_hline(aes(yintercept=mean(Beta$V3),alpha=0.2,col = "red"))+
    scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")}

# Point estimates for Regression coefficients
mean(Beta$V1)
sd(Beta$V1)
x1 = Beta$V1[1:1598]
x2 = Beta$V1[2:1599]
sd(Beta$V1) / sqrt(length(Beta$V1))
sd(Beta$V1) / sqrt(length(Beta$V1)) * sqrt((1+cor(x1,x2))/(1-cor(x1,x2)))
quantile(Beta$V1, 0.025)
quantile(Beta$V2, 0.975)

mean(Beta$V2)
sd(Beta$V2)
x1 = Beta$V2[1:1598]
x2 = Beta$V2[2:1599]
sd(Beta$V2) / sqrt(length(Beta$V2))
sd(Beta$V2) / sqrt(length(Beta$V2)) * sqrt((1+cor(x1,x2))/(1-cor(x1,x2)))
quantile(Beta$V2, 0.025)
quantile(Beta$V2, 0.975)

mean(Beta$V3)
sd(Beta$V3)
x1 = Beta$V3[1:1598]
x2 = Beta$V3[2:1599]
sd(Beta$V3) / sqrt(length(Beta$V3))
sd(Beta$V3) / sqrt(length(Beta$V3)) * sqrt((1+cor(x1,x2))/(1-cor(x1,x2)))
quantile(Beta$V3, 0.025)
quantile(Beta$V3, 0.975)


# Import of beta GARCH coefficient
beta = read.table("beta.txt", header = FALSE)

beta$V1
beta$iter=seq(1:1599)

# Trace plot for beta GARCH coefficient
w = {ggplot()+geom_line(data=beta,aes(x = iter, y = V1),alpha=1,size=0.3)+
    ggtitle(substitute("GARCH coeffiecient"))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
    theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
    ylab("Value")+ theme(legend.position = "none")+
    geom_hline(aes(yintercept=mean(beta$V1),alpha=0.2,col = "red"))+
    scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")}

# Import of alpha GARCH coefficient
alpha = read.table("alpha.txt", header = FALSE)
alpha$iter=seq(1:1599)

# Trace plots for alpha GARCH coefficient
k = {ggplot()+geom_line(data=alpha,aes(x = iter, y = V1),alpha=1,size=0.3)+
    ggtitle(substitute("First ARCH coeffiecient"))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
    theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
    ylab("Value")+ theme(legend.position = "none")+
    geom_hline(aes(yintercept=mean(alpha$V1),alpha=0.2,col = "red"))+
    scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")}

ggplot()+geom_line(data=alpha,aes(x = iter, y = V2),alpha=1,size=0.3)+
    ggtitle(substitute("Second ARCH coeffiecient"))+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
    theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
    ylab("Value")+ theme(legend.position = "none")+
    geom_hline(aes(yintercept=mean(alpha$V2),alpha=0.2,col = "red"))+
    scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")

# Import of overall mean coefficient
m = read.table("m.txt", header=F) 
m$iter=seq(1:1599)

# Trace plot for overall mean coefficient
ggplot()+geom_line(data=m,aes(x = iter, y = V1),alpha=1,size=0.3)+
  ggtitle(substitute("Second ARCH coeffiecient"))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
  theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
  ylab("Value")+ theme(legend.position = "none")+
  geom_hline(aes(yintercept=mean(m$V1),alpha=0.2,col = "red"))+
  scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")

# Import seasonal random intercepts
mu=read.table("Mu.txt", header=F) 
mu$iter=seq(1:1599)
ggplot()+geom_line(data=mu,aes(x = iter, y = V1),alpha=1,size=0.3)+
  ggtitle(substitute("mu"))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
  theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
  ylab("Value")+ theme(legend.position = "none")+
  geom_hline(aes(yintercept=mean(mu$V1),alpha=0.2,col = "red"))+
  scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")

# Import GARCH variance 
h = read.table("h.txt", header = F)
plot(t(h)[1,],type='l')
hist(t(h)[1,],80)
mean(t(h)[1,])
which(h[1,]==0)

# Import GARCH error 
Zeta=read.table("z.txt",header=F)
plot(t(Zeta)[1,],type='l')
hist(t(Zeta)[1,],80)
mean(t(Zeta)[1,])
which(Zeta[1,]==0)


# Import number of non-local bases
kk = read.table("k.txt", header = F)
max(kk)

# Import B-sline basis coefficients
Theta=read.table("Theta.txt",header=F)
Theta$iter=seq(1:1599)
ggplot()+geom_line(data=Theta,aes(x = iter, y = V1),alpha=1,size=0.3)+
  ggtitle(substitute("Second ARCH coeffiecient"))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
  theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
  ylab("Value")+ theme(legend.position = "none")+
  geom_hline(aes(yintercept=mean(Theta$V1),alpha=0.2,col = "red"))+
  scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")

# Import shrinkage measure for non-local bases
Lambda_norm=read.table("Lambda_Norm.txt",header=F, fill = T, col.names = paste0("V", seq_len(max(kk))))
Lambda_norm[is.na(Lambda_norm)] <- 0
Lambda=read.table("Lambda.txt", header = F, fill=T, col.names = paste0("V", seq_len(max(kk)*80)))
Lambda[is.na(Lambda)] <- 0
# exponential decay of non-local bases norms
boxplot(Lambda_norm)

# Import error term
psi=read.table("Psi.txt",header=F)
psi$iter=seq(1:1599)
ggplot()+geom_line(data=psi,aes(x = iter, y = V1),alpha=1,size=0.3)+
  ggtitle(substitute("Second ARCH coeffiecient"))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
  theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
  ylab("Value")+ theme(legend.position = "none")+
  geom_hline(aes(yintercept=mean(psi$V1),alpha=0.2,col = "red"))+
  scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")

mean(psi$V1)

abline(h = mean(t(psi)) , col = "red")
hist(t(psi)[1,],80)

grid.arrange(x, y, z, w, j, k,
             ncol = 3, nrow = 2)

# Error analysis

var_m = matrix(0, nrow=dim(h)[1], ncol=1)
for(i in 1:653){
  for(s in 1:Si[i]){
    var_m[,1] = var_m[,1] + 1./h[, 19*(i-1)+s] 
  }
}
var_m[,1] = var_m[,1] + 1/Hyperparameters$s0
var_m[,1] = 1/var_m[,1]

varm = read.table("VARm.txt", header = F)

meanm = read.table("MEANm.txt", header = F)
mean_m = matrix(0, nrow=dim(h)[1], ncol=1)
for(i in 1:653){
  for(s in 1:Si[i]){
    mean_m[,1] = mean_m[,1] + mu[, 19*(i-1)+s]/h[, 19*(i-1)+s]
  }
}
mean_m[,1] = mean_m[,1] * var_m[,1]

# Residuals Analysis

Residuals_all=fread("Residuals.txt", sep="auto", dec=".", header=F, blank.lines.skip=T)
Residuals_all=as.data.frame(Residuals_all)
Residuals_all$iter=seq(1:1599)


tt = cbind(rep(0,653),tt)
sum_tt = t(apply(tt, 1, cumsum))


idx=cumsum(nni[unique(data$ID)])
idx=c(0,idx)  # incremental index of observations for each athlete 

#season changes for all athletes 

time.span = difftime(as.Date(max(data$Date[which(data$ID == 26)]), origin = "1899-12-30"), 
                     as.Date(sprintf("%d-01-01", min(data$Season[which(ID == 26)])), origin = "1899-12-30"))

time.span = as.numeric(time.span)

start_year = min(data$Season[which(data$ID == 26)])
seasonchange = array(NA, dim = max(data$season_number[which(ID == 26)])+1) 
for (i in 0:max(data$season_number[which(ID == 26)])){
  seasonchange[i+1] = difftime(as.Date(sprintf("%d-01-01", start_year+i)), as.Date(sprintf("%d-01-01", start_year)), units = 'days') / time.span 
}

seasonchange

sqldf("select ID, count(tij) from data
      group by ID
      having count(tij) > 200")
sqldf("select distinct(ID), max(season_number) from data
       group by season_number
       having max(season_number) > 17")

its = 1:((Hyperparameters$nit-Hyperparameters$burn)/Hyperparameters$thin - 1)


# evolution of error per season with respect to iterations averaged over time (in season) and athletes

attach(data)

for(i in 1:max(Si)){
  res_t = apply(Residuals_all[,which(season_number == i)], 1, function(x){mean(x[x != 0])})
  windows()
  plot(res_t, main=substitute(paste('Residuals season', k), list(k = i)), type = 'l')
  abline(h = 0.0, col = 'red')
}

graphics.off()
# error per season averaged over iterations, athletes and time (per season)

mean_err_over_it = apply(Residuals_all, 2 ,mean)
mean_err_over_seas = matrix(0, nrow = max(Si), ncol = 653)
for(i in 1:max(Si)){
  for(j in 1:653){
    ind = intersect(which(season_number == i), (idx[j]:idx[j+1]))
    if(length(ind) == 0){
      mean_err_over_seas[i,j] = 0
    } else {
      mean_err_over_seas[i,j] = mean(mean_err_over_it[ind])
    }
  }
}

errors=matrix(NA,nrow=19,ncol=2)
errors[,1]=c(1:19)
errors[,2]=apply(mean_err_over_seas, 1, function(x){mean(x[x != 0])})
errors=as.data.frame(errors)
colnames(errors)=c("t","epsilon")

ggplot()+geom_point(data=errors,aes(x=t, y=epsilon),alpha=1,size=1.8)

# errors per time per athlete averaged over iterations... for a selection of athletes 

athl_list = c(16,101,149,189,226,259,504)
get_plette = colorRampPalette(brewer.pal(12, "Paired"))
athl_col = get_plette(10)

pp = ggplot() + theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold")) + 
  theme(plot.title = element_text(size=11)) + 
  theme(legend.title=element_text(size=9),legend.text=element_text(size=8)) +
  xlab("Time") + ylab("Errors") + xlim(0,1) + ylim(-1,1)

ii = 1
for(athl in athl_list){
  err = mean_err_over_it[(idx[athl]:idx[athl+1])]
  err = cbind(data$tij[(idx[athl]:idx[athl+1])],err)
  err = as.data.frame(err)
  colnames(err)=c("t","epsilon")
  pp = pp + geom_line(data=err,aes(x=t, y=epsilon),alpha=1,size=0.3,col = athl_col[ii])
  ii = ii+1
}
pp 

# error per miniseason averaged over iterations, athletes and time (per season)

miniseasonlength = (seasonchange[2]-seasonchange[1])/4
miniseasonchange = array(NA, dim = max(data$season_number[which(ID == 26)])*4+1)  
for(i in 0:(max(data$season_number[which(ID == 26)])-1)){
  miniseasonchange[4*i+1] = seasonchange[i+1]
  miniseasonchange[4*i+2] = seasonchange[i+1] + miniseasonlength
  miniseasonchange[4*i+3] = seasonchange[i+1] + 2*miniseasonlength
  miniseasonchange[4*i+4] = seasonchange[i+1] + 3*miniseasonlength
}
miniseasonchange[77] = seasonchange[19]

miniseason = array(NA, dim = dim(data)[1])  
for(i in 1:dim(data)[1]){
  miniseason[i] = which.max(miniseasonchange > data$tij[i])-1
}
data$miniseason = miniseason

mean_err_over_miniseas = matrix(0, nrow = max(miniseason)+1, ncol = 653)
for(i in 1:max(miniseason)){
  for(j in 1:653){
    ind = intersect(which(data$miniseason == i), (idx[j]:idx[j+1]))
    if(length(ind) == 0){
      mean_err_over_miniseas[i,j] = 0
    } else {
      mean_err_over_miniseas[i,j] = mean(mean_err_over_it[ind])
    }
  }
}

errors2=matrix(NA,nrow=76,ncol=2)
errors2[,1]=c(1:76)
errors2[,2]=apply(mean_err_over_miniseas, 1, function(x){mean(x[x != 0])})
errors2=as.data.frame(errors2)
colnames(errors2)=c("t","epsilon")
col = c("blue","black","black","red")
col = rep(col, 19)
ggplot()+geom_point(data=errors2,aes(x=t, y=epsilon),alpha=1,size=1.8,col = col)


# covariance function matrix
idx=cumsum(nni[unique(data$ID)]) #if use b
idx=c(0,idx)  # index usefull  

#estimated regression coefficients
Xb=matrix(NA,nrow = length(tg4),ncol=dim(Beta)[2]-1)   
beta_pr=matrix(NA,nrow=dim(Beta)[2]-1,ncol=3)
beta_pr[,1]=apply(Beta[,1:3],2,mean)
beta_pr[,2:3]=t(apply(Beta[,1:3],2,quantile,probs=c(0.05,0.95)))
beta_pr 

# Trajectory estimation 
start_trajectory = function(i){
  atl=unique(data$ID)[i]
  data_single_atl=data[which((data$ID==atl)),]
  
  # covariate contribution
  Xb[,1]=seq(data$Gender[idx[i]+1], data$Gender[idx[i+1]], length.out = length(tg4))
  Xb[,2]=rep(data$age[idx[i]+1], length.out = length(tg4))
  Xb[,3]=rep(NA,length(tg4))
  index = 1
  for(t in tg4){
    Xb[index,3] = data_single_atl$Environment[which.min(abs(data_single_atl$tij-t))]
    index = index+1
  }
  beta_pred=Xb%*%beta_pr
  
  #seasonal contribution
  bet_som=as.matrix(Beta[,1:3])%*%t(Xb)
  dim(bet_som)
  which((round(apply(bet_som, 2, mean), digits = 4) == round(beta_pred[,1], digits = 4)) == FALSE)
  
  mu_som=matrix(0,nrow=dim(mu)[1],ncol=length(tg4))
  dim(mu_som)
  
  # season change for athlete dependent seasons
  # season_change <- c(0)
  # for (j in 2:dim(data_single_atl)[1]) {
  #   if (data_single_atl$season_number[j]==(data_single_atl$season_number[j-1]+1))
  #     season_change<-c(season_change,(data_single_atl$tij[j]+data_single_atl$tij[j-1])/(2*max(data$tij)))
  # }
  # 
  # season_change<-c(season_change,(data_single_atl$tij[dim(data_single_atl)[1]]+0.05)/max(data$tij))
  # 
  # seaso=season_change[2:length(season_change)]
  
  season.i = max(data_single_atl$season_number)
  
  seaso2 = seasonchange[1:season.i+1]
  
  seaso = seaso2
  
  for (it in 1:dim(mu)[1]) {
    mu_pass=mu[it,((i-1)*max(data$season_number)+1):((i-1)*max(data$season_number)+length(seaso))]
    g=1
    for (k in 1:length(tg4)) {
      if (tg4[k]>seaso[g]) 
        g=g+1
      if (tg4[k]>seaso[length(seaso)])
        break
      mu_som[it,k]=mu_pass[,g]
    }
  }
  dim(mu_som)
  
  
  #functional contribution 
  q=ncol(BB26)
  
  its=1:(dim(Theta)[1])
  thti = Theta[its,((i-1)*q+1):(i*q)];
  Eyi=t(bb26%*%t(thti))
  dim(Eyi)
  
  #additive model
  tot_som=matrix(NA,nrow=dim(mu)[1],ncol=length(tg4))
  for (it in its) {
    for (t in 1:length(tg4)){
      tot_som[it,t]=Eyi[it,t]+mu_som[it,t]+bet_som[it,t]
    }
  }
  
  
  # Seasonal component estimate 
  est0 = matrix(0,nrow=length(tg4),ncol=3)
  est0[,1]=apply(mu_som,2,mean)
  est0[,2:3] = t(apply(mu_som,2,quantile,probs=c(0.05 ,0.95)));   
  est0[]=est0[]#+mean(data_single_atl$Result)
  
  # Functional component estimate 
  est1 = matrix(0,nrow=length(tg4),ncol=3)
  est1[,1]=apply(Eyi,2,mean)
  est1[,2:3] = t(apply(Eyi,2,quantile,probs=c(0.05 ,0.95)));   
  est1[]=est1[]#+mean(data_single_atl$Result)
  
  # Regression component estimate 
  est2 = matrix(0,nrow=length(tg4),ncol=3)
  est2[,1]=apply(bet_som,2,mean)
  est2[,2:3] = t(apply(bet_som,2,quantile,probs=c(0.05 ,0.95)));   
  est2[]=est2[]#+mean(data_single_atl$Result)
  
  # overall trajectory estimate 
  est3 = matrix(0,nrow=length(tg4),ncol=3)
  est3[,1]=apply(tot_som,2,mean)
  est3[,2:3] = t(apply(tot_som,2,quantile,probs=c(0.05 ,0.95)));   
  est3[]=est3[]#+mean(data_single_atl$Result)
  
  
  df=matrix(NA,nrow=length(data_single_atl$ID),ncol=2)
  df[,1]=t(data_single_atl$tij/max(data$tij))*6796
  df[,2]=t(data_single_atl$Result) - mean(data_single_atl$Result)
  df=as.data.frame(df)
  colnames(df)=c("t","y")
  df2=as.data.frame(est3)
  df2$t=tg4*6796
  df3=as.data.frame(est0)
  df3$t=tg4*6796
  df4=as.data.frame(est1)
  df4$t=tg4*6796
  df5=as.data.frame(est2)
  df5$t=tg4*6796
  colnames(df2)=c("tr.est","lb","ub","t")
  colnames(df3)=c("tr.est","lb","ub","t")
  colnames(df4)=c("tr.est","lb","ub","t")
  colnames(df5)=c("tr.est","lb","ub","t")
  colnames(df2)
  drop = which(tg4<min(data_single_atl$tij))
  drop2 =  which(tg4>(max(data_single_atl$tij)))
  drop = c(drop[1:length(drop)-1], drop2[-1])
  if (length(drop) > 1){
    df2=df2[-drop,]
    df3=df3[-drop,]
    df4=df4[-drop,]
    df5=df5[-drop,]
  }
  d=""
  if (data_single_atl$dop[1]==1)
    d=", Doped"
  
  print(i)
  
  ggplot()+geom_point(data=df,aes(x=t, y=y),col = "black" ,size=1)+
    theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
    theme(legend.title=element_text(size=9),legend.text=element_text(size=8))+xlab("Time in Days")+
    ylab("Lenght of throw (in m)")+
    geom_line(data=df2,aes(x=t,y=tr.est),color='black', size=1.2)+
    #geom_line(data=df3,aes(x=t,y=tr.est),color='black', size=1)+
    #geom_line(data=df4,aes(x=t,y=tr.est),color='black', size=1)+
    #geom_line(data=df5,aes(x=t,y=tr.est),color='black', size=1)+
    geom_vline(aes(xintercept=c(0,seasonchange[1:(length(seasonchange))]*6796 )),col = "black" ,size=0.001)+
    geom_line(data=df2,aes(x=t,y=lb),color='black', size=0.6)+
    geom_line(data=df2,aes(x=t,y=ub),color='black', size=0.6)+
    #geom_ribbon(aes(ymin=df2$lb, ymax=df2$ub, x=df2$t), alpha = 0.2, size = 0.1, color='black')+
    #geom_ribbon(aes(ymin=df3$lb, ymax=df3$ub, x=df2$t), alpha = 0.2,color='blue')+
    #geom_line(data=df3,aes(x=t,y=lb),color='black', size=0.6)+
    #geom_line(data=df3,aes(x=t,y=ub),color='black', size=0.6)+
    #geom_ribbon(aes(ymin=df4$lb, ymax=df4$ub, x=df4$t), alpha = 0.2,color='green')+
    #geom_line(data=df4,aes(x=t,y=lb),color='black', size=0.6)+
    #geom_line(data=df4,aes(x=t,y=ub),color='black', size=0.6)+
    #geom_ribbon(aes(ymin=df5$lb, ymax=df5$ub, x=df2$t), alpha = 0.2,color='yellow')+
    #geom_line(data=df5,aes(x=t,y=lb),color='black', size=0.6)+
    #geom_line(data=df5,aes(x=t,y=ub),color='black', size=0.6)+
    ggtitle(substitute(paste('Athlete ', k), list(k = i)))+
    scale_colour_manual("",values="blue")+scale_fill_manual("",values="grey12")+
    scale_x_continuous(breaks=seq(0,max(tij)*6796+365 ,365))+
    theme(axis.text=element_text(size=6),axis.title=element_text(size=9,face="bold"))+
    theme(plot.title = element_text(size=13,face="bold"))
}
  
  
start_trajectory(44)
start_trajectory(101)
start_trajectory(149)
start_trajectory(189)
start_trajectory(226)
start_trajectory(259)
start_trajectory(504)

a = start_trajectory(144)
b = start_trajectory(226)
c = start_trajectory(303)
d = start_trajectory(504)

setEPS()
postscript("datatraj.eps")
grid.arrange(a, b, c, d,
             ncol = 1, nrow = 4)
dev.off()

x = start_trajectory(226)
y = start_trajectory(226)
z = start_trajectory(226)
w = start_trajectory(226)

setEPS()
postscript("components.eps")
grid.arrange(x,y,z,w,
             ncol = 2, nrow = 2)
dev.off()

# Alternaive error analysis

ath = 122
obs = 25
index = idx[ath]+obs
estimation = as.data.frame(mu[, 19*(ath-1) + data$season_number[index]] - 
  as.matrix(Beta[,1:3])%*%(Reg_Matrix[index,]) - 
  as.matrix(Theta[,((ath-1)*q+1):(ath*q)]) %*% as.matrix(BB26[index,])) + Alt_esiduals[,index]
estimation$iter = seq(1:1599)


ggplot()+geom_line(data=estimation,aes(x = iter, y = V1),alpha=1,size=0.3)+
  ggtitle(substitute("Second ARCH coeffiecient"))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=9,face="bold"))+ 
  theme(plot.title = element_text(size=13,face="bold"))+xlab("Iteration")+
  ylab("Value")+ theme(legend.position = "none")+
  geom_hline(aes(yintercept=data$Result[index] - mean(data$Result[which((data$ID==ath))]),alpha=0.2,col = "red"))+
  scale_colour_manual("",values="red")+scale_fill_manual("",values="grey12")

Alt_esiduals = matrix(NA, nrow = 1599, ncol = 41033)
for(n in 1:653){
  mean = mean(data$Result[which((data$ID==n))])
  for(j in ((idx[n]+1):(idx[n+1]))){
    Alt_esiduals[,j] = data$Result[j] - mu[, 19*(n-1) + data$season_number[j]] - 
      as.matrix(Beta[,1:3])%*%(Reg_Matrix[j,]) - 
      as.matrix(Theta[its,((n-1)*q+1):(n*q)]) %*% as.matrix(BB26[j,]) - mean
  }
}
Alt_esiduals[,index]

Alt_mean_err_over_it = apply(Alt_esiduals, 2 ,mean)
Alt_mean_err_over_miniseas = matrix(0, nrow = max(miniseason)+1, ncol = 653)
for(i in 1:max(miniseason)){
  for(j in 1:653){
    ind = intersect(which(data$miniseason == i), (idx[j]:idx[j+1]))
    if(length(ind) == 0){
      Alt_mean_err_over_miniseas[i,j] = 0
    } else {
      Alt_mean_err_over_miniseas[i,j] = mean(Alt_mean_err_over_it[ind])
    }
  }
}

Alt_err2=matrix(NA,nrow=76,ncol=2)
Alt_err2[,1]=c(1:76)
Alt_err2[,2]=apply(Alt_mean_err_over_miniseas, 1, function(x){mean(x[x != 0])})
Alt_err2=as.data.frame(Alt_err2)
colnames(Alt_err2)=c("t","epsilon")
col = c("blue","black","black","red")
col = rep(col, 19)
ggplot()+geom_point(data=Alt_err2,aes(x=t, y=epsilon),alpha=1,size=1.8,col = col)

Alt_err504 = apply(Alt_esiduals[,which(data$ID == 504)], 2, mean)

# seasonal errors and patterns
q = 80
Seas_residuals = matrix(NA, nrow = 1599, ncol = 41033)
for(n in 1:653){
  mean = mean(data$Result[which((data$ID==n))])
  for(j in ((idx[n]+1):(idx[n+1]))){
    Seas_residuals[,j] = data$Result[j] - as.matrix(Beta[,1:3])%*%(as.matrix(Reg_Matrix)[j,]) - 
      as.matrix(Theta[its,((n-1)*q+1):(n*q)]) %*% as.matrix(BB26[j,]) - mean 
    #- mu[, 19*(n-1) + data$season_number[j]]
  }
}
Mean_seas_res_over_it = apply(Seas_residuals, 2 ,mean)
Mean_seas_res = matrix(0, nrow = max(Si), ncol = 653)
for(i in 1:max(Si)){
  for(j in 1:653){
    ind = intersect(which(data$season_number == i), (idx[j]:idx[j+1]))
    if(length(ind) == 0){
      Mean_seas_res[i,j] = 0
    } else {
      Mean_seas_res[i,j] = mean(Mean_seas_res_over_it[ind])
    }
  }
}

Alt_err2=matrix(NA,nrow=19,ncol=2)
Alt_err2[,1]=c(1:19)
Alt_err2[,2]=apply(Mean_seas_res, 1, function(x){mean(x[x != 0])})
Alt_err2=as.data.frame(Alt_err2)
colnames(Alt_err2)=c("t","epsilon")
ggplot()+geom_point(data=Alt_err2,aes(x=t, y=epsilon),alpha=1,size=1.8)

acf1(as.ts(Alt_err2[,2]), na.action = na.pass)
pacf(as.ts(Alt_err2[,2]), na.action = na.pass)

mean_err_over_miniseas2 = matrix(0, nrow = max(miniseason)+1, ncol = 653)
for(i in 1:(max(miniseason)+1)){
  for(j in 1:653){
    ind = intersect(which(data$miniseason == i), (idx[j]:idx[j+1]))
    if(length(ind) == 0){
      mean_err_over_miniseas2[i,j] = 0
    } else {
      mean_err_over_miniseas2[i,j] = mean(Mean_seas_res_over_it[ind])
    }
  }
}

length(which(mean_err_over_miniseas2 == 0))

errors2=matrix(NA,nrow=76,ncol=2)
errors2[,1]=c(1:76)
errors2[,2]=apply(mean_err_over_miniseas2, 1, function(x){mean(x[x != 0])})
errors2=as.data.frame(errors2)
colnames(errors2)=c("t","epsilon")
col = c("blue","black","black","red")
col = rep(col, 19)
ggplot()+geom_point(data=errors2,aes(x=t, y=epsilon),alpha=1,size=1.8,col = col)

acf1(as.ts(errors2[,2]^2), max.lag = 20, na.action = na.pass)
pacf(as.ts(errors2[,2]^2), na.action = na.pass)


# Non-local basis

tdata = data[, c('tij', 'season_number')]
tdata$tij = tdata$tij / max(tdata$tij)
sqldf("select season_number, min(tij), max(tij), max(tij)-min(tij), count(tij) from tdata
      group by season_number ")


#df=p
BB=bs(data$tij/max(data$tij),degree=3,df=80,intercept = F)
tg=seq(0,1,length.out = 200)
bb=bs(tg,degree=3,df=80,intercept = F)
matplot(bb, type="l")

a = data$tij/max(data$tij)
H = cbind(a,BB)
H = as.data.frame(H)
H = H[order(H$a),]
H = as.matrix(H)
matplot(H[,1], H[,2:31], type = "l")

dim(Lambda)
kk = read.table("k.txt", header = F)
Upperk = median(kk$V1)
# nonlocal basis
Psitilda = array(0, c(Upperk, 200))
for (i in 1:Upperk) {
  Lambdai = as.matrix(Lambda)[,((80*(i-1)+1):(80*i))]
  dim(Lambdai)
  Psitildai = Lambdai%*%t(bb)
  Psitilda[i,] = apply(Psitildai,2,mean)
  windows()
  plot(tg, Psitilda[i,], type = "l", lwd = 2  , main=substitute(paste('Nonlocal basis ', k), list(k = c(1:Upperk)[i])))
  abline(v = seasonchange, col = 'gray')
  abline(h = 0.0, col = 'red')
}

# Expected difference in response to changes in the internal covariates

r1 = length(IntRegressors)

RegFunc = array(0, c(r1, 200))

for (j in 1:r1) {
  BetaEtaj = as.matrix(Int_beta)[, j+r1*(0:(Upperk-1))]
  dim(BetaEtaj)
  RegFuncj = BetaEtaj%*%Psitilda
  RegFunc[j,] = apply(RegFuncj,2,mean)
  windows()
  plot(tg, RegFunc[j,],  main=substitute(paste('Response to Internal Regressor ', k), list(k = c(1:r1)[j])), type = 'l')
  abline(v = seasonchange, col = 'gray')
  abline(h = 0.0, col = 'red')
}


# compute LPML
LPLM=LPLM_comp(as.matrix(data[,c(1,2,3,9)]),as.matrix(data[,c(4,5,7)]),as.matrix(BB26),as.matrix(tt[iddd,]),Si[iddd],
                     as.matrix(psi),as.matrix(Theta),as.matrix(Ex_beta),as.matrix(mu)) 


Fi = matrix(NA, nrow = 1599, ncol = 653)
k = 1
q = 80
for(n in 1:653){
  mean = mean(data$Result[which((data$ID==n))])
  Fi[,k] = matrix(1, nrow = 1599, ncol = 1)
  for(j in ((idx[n]+1):(idx[n+1]))){
    Fi[,k] = Fi[,k] / sqrt(2 * pi * as.matrix(psi$V1)) * exp( - (data$Result[j] - (as.matrix(Beta)[,1:3]%*%(Reg_Matrix[j,]) + 
                                                                                  as.matrix(Theta[,((n-1)*q+1):(n*q)]) %*% as.matrix(BB26[j,]) +
                                                                                  mu[, 19*(n-1) + data$season_number[j]] + mean)  )^2 / (2 * as.matrix(psi$V1)) ) 
  }
  k = k+1
}

FI = 1 / Fi

Inv_CPOi = apply(FI, 2 , mean)

CPOi = 1 / Inv_CPOi

LPLM = sum(log(CPOi))
