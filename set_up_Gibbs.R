library("Rcpp")
library("RcppArmadillo")

# Data inport, consistency checks 
data = read.table("COMPLETE_DATA.txt", sep = " ", header = T)
data = data[which(Result >= 10),]
attach(data)
data = data[order(ID,tij),]
max(data$tij)

# numerate seasons from 1
data$season_number=data$season_number+1

dim(data)
colnames(data)
length(unique(data$ID))

# computation of number of observation per athlete
# nni: total # of shot put in career
ni<-tabulate(data[,1])
nni<-c()
for (i in 1:length(ni))
  if (ni[i]!=0) nni<-c(nni,ni[i]) 

# s[j] = # of observations in season j
si<-tabulate(data[,"season_number"]) 

# number of observation per single season per athlete
tt=table(data[,c("ID","season_number")]) 

# incremental index of observations for each athlete 
idx=cumsum(nni[unique(data$ID)])
idx=c(0,idx)  

# Si: number of seasons per athlete i
Si=c()
for (j in unique(data$ID))
  Si<-c(Si,max(data[data$ID==j,"season_number"])) 

# basis function B-spline
library(splines)
BB=bs(data$tij/max(data$tij),degree=3,df=80,intercept = F) # 80df spline basis evaluated on the time points
tg=seq(0,1,length.out = 200) # 200 pointgrid on [0,1]
bb=bs(tg,degree=3,df=80,intercept = F) # 80df spline basis evaluated on the grid
matplot(bb, type="l")

attach(data)

# Hyperparameter selection
Hyperparameters = list(nit = 20000, burn = 12000, thin = 5, 
                       as = 1, bs = 0.3, df = 9, ad1 = 2.1, bd1 = 1, 
                       ad2 = 2.1, bd2 =  1, aw = 0.5, bw = 0.5,
                       nu0alpha1 = 0.5, sigma0alpha1 = 0.5, nu0alpha2 = 0.2, sigma0alpha2 = 0.5,
                       corr0alpha = 0, nu0beta = 0, sigma0beta = 1, 
                       m0 = -0.2, s0 = 0.0001,	
                       nu0psi = 1.0, sigma0psi = 1.0,
                       prop = 1, b0 = 1.5, b1 = 0.0001, epsilon = 0.0001, initialk = 25)


# prepare covariate matrix
Regressors = c("Gender", "age", "Environment")

InitialAge = matrix(NA, ncol = 1, nrow = 41030)
for(n in 1:653){
  InitialAge[(idx[n]+1):idx[j+1]] = min(data$age[which(data$ID == n)])
}
Reg_Matrix = cbind(as.matrix(data[,"Gender"]), InitialAge ,as.matrix(data[,"Environment"]))


# input GIBBS sampler
sourceCpp(file = "Functional_GARCH_AMH.cpp")

# launch GIBBS samper
gibbs_compl_garch(as.matrix(data[,c("ID","tij","Result","season_number")]),
                  Reg_Matrix,as.matrix(BB), as.matrix(tt), Si[],
                  nit = Hyperparameters$nit, burn = Hyperparameters$burn, thin = Hyperparameters$thin,
                  as = Hyperparameters$as, bs = Hyperparameters$bs, 
                  df = Hyperparameters$df, ad1 = Hyperparameters$ad1, bd1 = Hyperparameters$bd1, 
                  ad2 = Hyperparameters$ad2, bd2 = Hyperparameters$bd2, 
                  aw = Hyperparameters$aw, bw = Hyperparameters$bw,
                  nu0alpha1 = Hyperparameters$nu0alpha1, sigma0alpha1 = Hyperparameters$sigma0alpha1,
                  nu0alpha2 = Hyperparameters$nu0alpha2, sigma0alpha2 = Hyperparameters$sigma0alpha2,
                  corr0alpha = Hyperparameters$corr0alpha, 
                  nu0beta = Hyperparameters$nu0beta, sigma0beta = Hyperparameters$sigma0beta,
                  m0 = Hyperparameters$m0, s0 = Hyperparameters$s0,
                  nu0psi = Hyperparameters$nu0psi, sigma0psi = Hyperparameters$sigma0psi,
                  prop = Hyperparameters$prop, b0 = Hyperparameters$b0,
                  b1 = Hyperparameters$b1, epsilon = Hyperparameters$epsilon,
                  k = Hyperparameters$initialk)

