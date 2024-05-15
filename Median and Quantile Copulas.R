
library(ggplot2)
library(reshape2)
library(latex2exp)
library(copula)
library(viridis)
library(cowplot)
library(MASS)

#################
### Section 1 ###
#################

# No figures in this section


#################
### Section 2 ###
#################


### Frank copula ###

C1=function(u1,u2){
  return(pCopula(c(u1,u2), frankCopula(-10)))
}

C2=function(u1,u2){
  return(pCopula(c(u1,u2), frankCopula(10)))
}

u1=seq(0,1,0.01)
u2=u1
m=length(u1)

Frank_cop_minus10=matrix(nrow=m,ncol=m)
Frank_cop_10=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Frank_cop_minus10[i,j]=C1(u1[i],u2[j])
    Frank_cop_10[i,j]=C2(u1[i],u2[j])
  }
}

Frank_copula=rbind(cbind(melt(Frank_cop_minus10),"Parameter"=rep('-10',nrow(melt(Frank_cop_minus10)))) , cbind(melt(Frank_cop_10),"Parameter"=rep('10',nrow(melt(Frank_cop_10)))))
ggplot(Frank_copula, aes(x=u1[Var1], y=u2[Var2], z=value, linetype=Parameter,  colour=Parameter)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, colour='black') +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("solid", "dotted"), labels=unname(TeX(c('$\\theta = -10', '$\\theta = 10'))))

Frank_cop_av_id=matrix(nrow=m,ncol=m)
Frank_cop_av_sq=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Frank_cop_av_id[i,j]=(1/2)*C1(u1[i],u2[j])+(1/2)*C2(u1[i],u2[j])
    Frank_cop_av_sq[i,j]=(1/3)*C1(u1[i],u2[j])+(2/3)*C2(u1[i],u2[j])
  }
}

Frank_copula_av=rbind(cbind(melt(Frank_cop_av_id),"w"=rep('x',nrow(melt(Frank_cop_av_id)))) , cbind(melt(Frank_cop_av_sq),"w"=rep('x^2',nrow(melt(Frank_cop_av_sq)))))
ggplot(Frank_copula_av, aes(x=u1[Var1], y=u2[Var2], z=value, colour=w)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, linetype='dashed') +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', colour=TeX('w(x)'), linetype=TeX('w(x)')) +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_colour_manual(values = c("#00BFC4", "#F8766D"), labels=unname(TeX(c('x', 'x^2'))))

# Example 2.2 with w(x)=x

p=c(0.1,0.5,0.9)
u=seq(0,1,0.01)
v=u

par=c(1:length(p))
for(k in 1:length(p)){
  par[k]=p[k]
  C=matrix(nrow=length(u), ncol=length(v))
  for (i in 1:length(u)){
    for (j in 1:length(v)){
      C[i,j]=(1-par[k])*C1(u[i],v[j])+(par[k])*C2(u[i],v[j])
    }
  }
  cop <- paste("cop_q_", p[k], sep = "")
  assign(cop, C)
}

Quantile_copula=rbind(cbind(melt(cop_q_0.1),"Quantile"=rep('p=0.1',nrow(melt(cop_q_0.1)))) , cbind(melt(cop_q_0.5),"Quantile"=rep('p=0.5',nrow(melt(cop_q_0.5)))) , cbind(melt(cop_q_0.9),"Quantile"=rep('p=0.9',nrow(melt(cop_q_0.9)))))
ggplot(Quantile_copula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Quantile)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, colour="#00BFC4") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"))

# Example 2.2 with w(x)=x^2

p=c(0.1,0.5,0.9)
u=seq(0,1,0.01)
v=u

par=c(1:length(p))
for(k in 1:length(p)){
  par[k]=p[k]
  C=matrix(nrow=length(u), ncol=length(v))
  for (i in 1:length(u)){
    for (j in 1:length(v)){
      C[i,j]=(1-par[k])^2*C1(u[i],v[j])+par[k]*(2-par[k])*C2(u[i],v[j])
    }
  }
  cop <- paste("cop_q_", p[k], sep = "")
  assign(cop, C)
}

Quantile_copula=rbind(cbind(melt(cop_q_0.1),"Quantile"=rep('p=0.1',nrow(melt(cop_q_0.1)))) , cbind(melt(cop_q_0.5),"Quantile"=rep('p=0.5',nrow(melt(cop_q_0.5)))) , cbind(melt(cop_q_0.9),"Quantile"=rep('p=0.9',nrow(melt(cop_q_0.9)))))
ggplot(Quantile_copula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Quantile)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, colour='#00BFC4') +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"))


### Beta distribution ###

C_A=function(u1,u2, a,b){
  return((a/(a+b))*C1(u1,u2)+(b/(a+b))*C2(u1,u2))
}

C_med=function(u1,u2, a,b){
  if(a==b){
    val=0.5*C1(u1,u2)+0.5*C2(u1,u2)
  }
  if(a==1 & b>0){
    val=(1-2^(-1/b))*C1(u1,u2)+(2^(-1/b))*C2(u1,u2)
  }
  if(a>0 & b==1){
    val=(2^(-1/a))*C1(u1,u2)+(1-2^(-1/a))*C2(u1,u2)
  }
  return(val)
}


u1=seq(0,1,0.01)
u2=u1
m=length(u1)

Av_cop_beta0505=matrix(nrow=m,ncol=m)
Av_cop_beta13=matrix(nrow=m,ncol=m)
Av_cop_beta31=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Av_cop_beta0505[i,j]=C_A(u1[i],u2[j], 0.5,0.5)
    Av_cop_beta13[i,j]=C_A(u1[i],u2[j], 1,3)
    Av_cop_beta31[i,j]=C_A(u1[i],u2[j], 3,1)
  }
}

Av_copula=rbind(cbind(melt(Av_cop_beta0505),"Beta"=rep('0505',nrow(melt(Av_cop_beta0505)))) , cbind(melt(Av_cop_beta13),"Beta"=rep('13',nrow(melt(Av_cop_beta13)))) , cbind(melt(Av_cop_beta31),"Beta"=rep('31',nrow(melt(Av_cop_beta31)))))
ggplot(Av_copula, aes(x=u1[Var1], y=u2[Var2], z=value,  colour=Beta)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1,linetype="dashed") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', colour='Parameters (a,b)') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c('(0.5,0.5)','(1,3)', '(3,1)'))

Med_cop_beta0505=matrix(nrow=m,ncol=m)
Med_cop_beta13=matrix(nrow=m,ncol=m)
Med_cop_beta31=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Med_cop_beta0505[i,j]=C_med(u1[i],u2[j], 0.5,0.5)
    Med_cop_beta13[i,j]=C_med(u1[i],u2[j], 1,3)
    Med_cop_beta31[i,j]=C_med(u1[i],u2[j], 3,1)
  }
}

Med_copula=rbind(cbind(melt(Med_cop_beta0505),"Beta"=rep('0505',nrow(melt(Med_cop_beta0505)))) , cbind(melt(Med_cop_beta13),"Beta"=rep('13',nrow(melt(Med_cop_beta13)))) , cbind(melt(Med_cop_beta31),"Beta"=rep('31',nrow(melt(Med_cop_beta31)))))
ggplot(Med_copula, aes(x=u1[Var1], y=u2[Var2], z=value,  colour=Beta)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, linetype="solid") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='Parameters (a,b)', colour='Parameters (a,b)') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c('(0.5,0.5)','(1,3)', '(3,1)'))



#################
### Section 3 ###
#################

u=0.4
v=0.8
# Let X~U[0,1], then w=F(x)=x
x=seq(0,1,0.01)
w=0.4
t=seq(0,1,0.01)

par=c(-0.5,3)
p=c('neg', 'pos')
for(k in 1:length(par)){
  C=c(1:length(x))
  g=c(1:length(t))
  theta=par[k]
  phi=function (t){
    -log((exp(-theta*t)-1)/(exp(-theta)-1))
  }
  dphi=function(t){
    theta/(1-exp(theta*t))
  }
  ddphi=function(t){
    (theta^2)*exp(theta*t)/((1-exp(theta*t))^2)
  }
  phi_inv=function(s){
    (-1/theta)*log(exp(-s)*(exp(-theta)-1)+1)
  }
  dphi_inv=function(s){
    (1/theta)*log(1-(theta/s))
  }
  C_X_Arch=function(w,u,v){
    dphi(w)/dphi(phi_inv( phi( dphi_inv(dphi(w)/u) ) + phi( dphi_inv(dphi(w)/v)) - phi(w)))
  }
  for (i in 1:length(x)){
    C[i]=C_X_Arch(x[i],u,v)
  }
  cop <- paste("cop_", p[k], sep = "")
  assign(cop, C)
  for (i in 1:length(t)){
    g[i]=ddphi(w)/(ddphi(dphi_inv(dphi(w)/t[i]))*t[i]^2)
  }
  gg <- paste("g_", p[k], sep = "")
  assign(gg, g)
}

Copula=rbind(cbind(melt(cop_neg), x,"Parameter"=rep('-0.5',nrow(melt(cop_neg)))), cbind(melt(cop_pos),x,"Parameter"=rep('3',nrow(melt(cop_pos)))))
ggplot(Copula, aes(x=x, y=value, linetype=Parameter)) +
  geom_line(size=1, colour='black')+
  theme_bw() +
  labs(title="", x=TeX('x'), y=TeX('C_X(0.4,0.8)')) +
  theme(legend.position = "bottom") +
  scale_linetype_manual(values = c("solid", "dotted"), labels=unname(TeX(c('$\\theta = -0.5', '$\\theta = 3'))))

g_fun=rbind(cbind(melt(g_neg), x,"Parameter"=rep('-0.5',nrow(melt(g_neg)))), cbind(melt(g_pos),x,"Parameter"=rep('3',nrow(melt(g_pos)))))
ggplot(g_fun, aes(x=x, y=value, linetype=Parameter)) +
  geom_line(size=1, colour='black')+
  theme_bw() +
  labs(title="", x=TeX('t'), y=TeX('g(t)')) +
  theme(legend.position = "bottom") +
  scale_linetype_manual(values = c("solid", "dotted"), labels=unname(TeX(c('$\\theta = -0.5', '$\\theta = 3'))))



### C is Frank copula with theta=3 ###

C_A=function(u1,u2, theta){
  return((u1*u2/(1-(1-u1)*(1-u2)))*((1/theta)*log((1-u1)*(1-u2)*(exp(-theta)-1)+1)+1))
}

C_med=function(u1,u2, theta){
  return(u1*u2/(1-(1-exp(-theta/2))*(1-u1)*(1-u2)))
}

u1=seq(0,1,0.01)
u2=u1
m=length(u1)

Av_cop_frank=matrix(nrow=m,ncol=m)
Med_cop_frank=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Av_cop_frank[i,j]=C_A(u1[i],u2[j], 3)
    Med_cop_frank[i,j]=C_med(u1[i],u2[j], 3)
  }
}

Av_med_copula=rbind(cbind(melt(Av_cop_frank),"Copula"=rep('average',nrow(melt(Av_cop_frank)))) , cbind(melt(Med_cop_frank),"Copula"=rep('median',nrow(melt(Med_cop_frank)))))
ggplot(Av_med_copula, aes(x=u1[Var1], y=u2[Var2], z=value, linetype=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, colour="#00BFC4") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='Copula', colour='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dashed", "solid"), labels=c('Average','Median'))


p=seq(0,1,0.1)
theta= 3 # >= -log(2)/(1-p) and different from 0
u=seq(0,1,0.01)
v=u

par=c(1:length(p))
for(k in 1:length(p)){
  if (theta > 0){
    par[k]=1-exp(-theta*p[k])
    print(par[k])
  } else{ 
    if ( theta>-log(2)/(1-p[k]) & theta<0 ){
      par[k]=1-exp(-theta*(1-p[k]))
      print(par[k])
    } else{ print('error')}
  }
  
  C=matrix(nrow=length(u), ncol=length(v))
  for (i in 1:length(u)){
    for (j in 1:length(v)){
      C[i,j]=u[i]*v[j]/(1-par[k]*(1-u[i])*(1-v[j]))
    }
  }
  cop <- paste("cop_q_", p[k], sep = "")
  assign(cop, C)
}

Quantile_copula=rbind(cbind(melt(cop_q_0.1),"Quantile"=rep('p=0.1',nrow(melt(cop_q_0.1)))) , cbind(melt(cop_q_0.5),"Quantile"=rep('p=0.5',nrow(melt(cop_q_0.5)))) , cbind(melt(cop_q_0.9),"Quantile"=rep('p=0.9',nrow(melt(cop_q_0.9)))))
ggplot(Quantile_copula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Quantile)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, ,  colour="#00BFC4") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"))



### C is AMH copula with theta=-0.6 ####

p=seq(0,1,0.1)
theta=-0.6

u=seq(0,1,0.01)
v=u

par=c(1:length(p))
for(k in 1:length(p)){
  if (theta > 0 & theta < 1){
    par[k]=p[k]
    print(par[k])
  } else{ 
    if ( theta>=-1 & theta<0 ){
      par[k]=1-p[k]
      print(par[k])
    } else{ print('error')}
  }
  
  C=matrix(nrow=length(u), ncol=length(v))
  for (i in 1:length(u)){
    for (j in 1:length(v)){
      C[i,j]= ( ( 4*sqrt(u[i]*v[j])*(1-theta)*(1-theta +theta*par[k])*(1-theta+sqrt((1-theta)^2+4*theta*u[i]*par[k]*(1-theta+theta*par[k])) ) * (1-theta+sqrt((1-theta)^2+4*theta*v[j]*par[k]*(1-theta+theta*par[k])) ) )/( (1-theta+sqrt((1-theta)^2+4*theta*u[i]*par[k]*(1-theta+theta*par[k])) )^2 *(1-theta+sqrt((1-theta)^2+4*theta*v[j]*par[k]*(1-theta+theta*par[k])) )^2 - 16*u[i]*v[j]*theta*par[k]*(1-theta +theta* par[k])^3) )^2 
      
    }
  }
  cop <- paste("cop_q_", p[k], sep = "")
  assign(cop, C)
}

Quantile_copula=rbind(cbind(melt(cop_q_0.1),"Quantile"=rep('p=0.1',nrow(melt(cop_q_0.1)))) , cbind(melt(cop_q_0.5),"Quantile"=rep('p=0.5',nrow(melt(cop_q_0.5)))) , cbind(melt(cop_q_0.9),"Quantile"=rep('p=0.9',nrow(melt(cop_q_0.9)))))
ggplot(Quantile_copula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Quantile)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, colour="#00BFC4") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"))


### C has generator (theta/t + 1)(1 - t) with theta=1 ###

p=seq(0,1,0.1)
theta= 1 # >= 0
u=seq(0,1,0.01)
v=u

par=c(1:length(p))
for(k in 1:length(p)){
  if (theta > 0){
    par[k]=1-p[k]
    print(par[k])
  } else{print('error')}
  
  C=matrix(nrow=length(u), ncol=length(v))
  for (i in 1:length(u)){
    for (j in 1:length(v)){
      Y= (-theta/sqrt(theta*u[i]/(1+theta*par[k]^(-2)-u[i]))) - (theta/sqrt(theta*v[j]/(1+theta*par[k]^(-2)-v[j]))) + (sqrt(theta*u[i]/(1+theta*par[k]^(-2)-u[i]))) + (sqrt(theta*v[j]/(1+theta*par[k]^(-2)-v[j]))) + theta*par[k]^(-1) - par[k]
      C[i,j]=(theta*par[k]^(-2)+1)/(4*theta*( Y+sqrt(Y^2+4*theta) )^(-2)+1)
    }
  }
  cop <- paste("cop_q_", p[k], sep = "")
  assign(cop, C)
}

Quantile_copula=rbind(cbind(melt(cop_q_0.1),"Quantile"=rep('p=0.1',nrow(melt(cop_q_0.1)))) , cbind(melt(cop_q_0.5),"Quantile"=rep('p=0.5',nrow(melt(cop_q_0.5)))) , cbind(melt(cop_q_0.9),"Quantile"=rep('p=0.9',nrow(melt(cop_q_0.9)))))
ggplot(Quantile_copula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Quantile)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1,  colour="#00BFC4") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"))


#################
### Section 4 ###
#################

### Mixture of normals ###

C1=function(u1,u2){
  return(pCopula(c(u1,u2), frankCopula(-10)))
}

C2=function(u1,u2){
  return(pCopula(c(u1,u2), frankCopula(10)))
}

mu1=0.15
mu2=0.95
sigma1=0.05
sigma2=0.01
a=mu2*(pnorm((1-mu2)/sigma2)-pnorm(-mu2/sigma2))+sigma2*(dnorm(-mu2/sigma2)-dnorm((1-mu2)/sigma2))
b=mu1*(pnorm((1-mu1)/sigma1)-pnorm(-mu1/sigma1))+sigma1*(dnorm(-mu1/sigma1)-dnorm((1-mu1)/sigma1)) - (mu2*(pnorm((1-mu2)/sigma2)-pnorm(-mu2/sigma2))+sigma2*(dnorm(-mu2/sigma2)-dnorm((1-mu2)/sigma2)))

C_A=function(u1,u2, gamma){
  return((a+b*gamma)*C1(u1,u2)+(1-(a+b*gamma))*C2(u1,u2))
}

# Quantiles of w(X)
g=0.6

F = function(y,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  1- gamma*pnorm(1,mu1,sigma1)-(1-gamma)*pnorm(1,mu2,sigma2) + gamma * pnorm(y, mu1, sigma1) + (1-gamma) * pnorm(y, mu2, sigma2)
}
F_inv = function(q,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  minmax <- c(0,1)
  uniroot(function(x) F(x,mu1,sigma1,g,mu2,sigma2)-q,
          interval = minmax,
          tol = 10^{-16})$root  
}

F_inv(0.5, gamma=0.6)

C_med=function(u1,u2, gamma){
  if(gamma==1){
    val=0.15*C1(u1,u2)+(1-0.15)*C2(u1,u2)
  }
  if(gamma==0.8){
    val=0.1659*C1(u1,u2)+(1-0.1659)*C2(u1,u2)
    }
  if(gamma==0.6){
    val=0.1984*C1(u1,u2)+(1-0.1984)*C2(u1,u2)
    }
  return(val)
}

C_q01=function(u1,u2, gamma){
  if(gamma==1){
    val=0.086*C1(u1,u2)+(1-0.086)*C2(u1,u2)
  }
  if(gamma==0.8){
    val=0.092*C1(u1,u2)+(1-0.092)*C2(u1,u2)
  }
  if(gamma==0.6){
    val=0.102*C1(u1,u2)+(1-0.102)*C2(u1,u2)
  }
  return(val)
}

C_q09=function(u1,u2, gamma){
  if(gamma==1){
    val=0.214*C1(u1,u2)+(1-0.214)*C2(u1,u2)
  }
  if(gamma==0.8){
    val=0.950*C1(u1,u2)+(1-0.950)*C2(u1,u2)
  }
  if(gamma==0.6){
    val=0.957*C1(u1,u2)+(1-0.957)*C2(u1,u2)
  }
  return(val)
}

u1=seq(0,1,0.01)
u2=u1
m=length(u1)

Av_cop_gamma1=matrix(nrow=m,ncol=m)
Av_cop_gamma08=matrix(nrow=m,ncol=m)
Av_cop_gamma06=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Av_cop_gamma1[i,j]=C_A(u1[i],u2[j], 1)
    Av_cop_gamma08[i,j]=C_A(u1[i],u2[j], 0.8)
    Av_cop_gamma06[i,j]=C_A(u1[i],u2[j], 0.6)
  }
}

Av_copula=rbind(cbind(melt(Av_cop_gamma1),"Gamma"=rep('1',nrow(melt(Av_cop_gamma1)))) , cbind(melt(Av_cop_gamma08),"Gamma"=rep('0.8',nrow(melt(Av_cop_gamma08)))) , cbind(melt(Av_cop_gamma06),"Gamma"=rep('0.6',nrow(melt(Av_cop_gamma06)))))
ggplot(Av_copula, aes(x=u1[Var1], y=u2[Var2], z=value, linetype=Gamma,  colour=Gamma)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1) +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='', colour='') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6')))) +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))

Med_cop_gamma1=matrix(nrow=m,ncol=m)
Med_cop_gamma08=matrix(nrow=m,ncol=m)
Med_cop_gamma06=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Med_cop_gamma1[i,j]=C_med(u1[i],u2[j], 1)
    Med_cop_gamma08[i,j]=C_med(u1[i],u2[j], 0.8)
    Med_cop_gamma06[i,j]=C_med(u1[i],u2[j], 0.6)
  }
}

Med_copula=rbind(cbind(melt(Med_cop_gamma1),"Gamma"=rep('1',nrow(melt(Med_cop_gamma1)))) , cbind(melt(Med_cop_gamma08),"Gamma"=rep('0.8',nrow(melt(Med_cop_gamma08)))) , cbind(melt(Med_cop_gamma06),"Gamma"=rep('0.6',nrow(melt(Med_cop_gamma06)))))
ggplot(Med_copula, aes(x=u1[Var1], y=u2[Var2], z=value, linetype=Gamma,  colour=Gamma)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1) +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='', colour='') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6')))) +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))

Q01_cop_gamma1=matrix(nrow=m,ncol=m)
Q01_cop_gamma08=matrix(nrow=m,ncol=m)
Q01_cop_gamma06=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Q01_cop_gamma1[i,j]=C_q01(u1[i],u2[j], 1)
    Q01_cop_gamma08[i,j]=C_q01(u1[i],u2[j], 0.8)
    Q01_cop_gamma06[i,j]=C_q01(u1[i],u2[j], 0.6)
  }
}

Q01_copula=rbind(cbind(melt(Q01_cop_gamma1),"Gamma"=rep('1',nrow(melt(Q01_cop_gamma1)))) , cbind(melt(Q01_cop_gamma08),"Gamma"=rep('0.8',nrow(melt(Q01_cop_gamma08)))) , cbind(melt(Q01_cop_gamma06),"Gamma"=rep('0.6',nrow(melt(Q01_cop_gamma06)))))
ggplot(Q01_copula, aes(x=u1[Var1], y=u2[Var2], z=value, linetype=Gamma,  colour=Gamma)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1) +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='', colour='') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6')))) +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))

Q09_cop_gamma1=matrix(nrow=m,ncol=m)
Q09_cop_gamma08=matrix(nrow=m,ncol=m)
Q09_cop_gamma06=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Q09_cop_gamma1[i,j]=C_q09(u1[i],u2[j], 1)
    Q09_cop_gamma08[i,j]=C_q09(u1[i],u2[j], 0.8)
    Q09_cop_gamma06[i,j]=C_q09(u1[i],u2[j], 0.6)
  }
}

Q09_copula=rbind(cbind(melt(Q09_cop_gamma1),"Gamma"=rep('1',nrow(melt(Q09_cop_gamma1)))) , cbind(melt(Q09_cop_gamma08),"Gamma"=rep('0.8',nrow(melt(Q09_cop_gamma08)))) , cbind(melt(Q09_cop_gamma06),"Gamma"=rep('0.6',nrow(melt(Q09_cop_gamma06)))))
ggplot(Q09_copula, aes(x=u1[Var1], y=u2[Var2], z=value, linetype=Gamma,  colour=Gamma)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1) +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='', colour='') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dashed"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6')))) +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))


u1=0.5 ; u2=0.1
round(C_A(u1,u2,1),4); round(C_A(u1,u2,0.8),4); round(C_A(u1,u2,0.6),4)
round(C_q01(u1,u2,1),4); round(C_q01(u1,u2,0.8),4); round(C_q01(u1,u2,0.6),4)
round(C_med(u1,u2,1),4); round(C_med(u1,u2,0.8),4); round(C_med(u1,u2,0.6),4)
round(C_q09(u1,u2,1),4); round(C_q09(u1,u2,0.8),4); round(C_q09(u1,u2,0.6),4)


### Mixture of normals: histogram ###

par(mfrow=c(1,3))
nr=1000
# X= p*X_reg + (1-p)*X_outl
p=0.6
set.seed(12)
X_reg=rnorm(p*nr,0.15,0.05)
X_outl=rnorm(nr-nr*p,0.95,0.01)
X=c(X_reg,X_outl)
hist(X,main="", freq = FALSE, xlim=c(0,1), breaks=15, ylim=c(0,9), xlab="w(X)")


### Mixture of normals: plot quantiles of w(X)
p=c(seq(0,0.01,0.001),seq(0.01001,1,0.00001))
mu1=0.15
mu2=0.95
sigma1=0.05
sigma2=0.01
g=1

F = function(y,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  1- gamma*pnorm(1,mu1,sigma1)-(1-gamma)*pnorm(1,mu2,sigma2) + gamma * pnorm(y, mu1, sigma1) + (1-gamma) * pnorm(y, mu2, sigma2)
}
F_inv = function(q,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  minmax <- c(0,1)
  uniroot(function(x) F(x,mu1,sigma1,g,mu2,sigma2)-q,
          interval = minmax,
          tol = 10^{-20})$root  
}

Quantiles_gamma1=c(1:length(p))
Quantiles_gamma1[1]=0
Quantiles_gamma1[2]=0
Quantiles_gamma1[length(p)]=1
for(i in 3:(length(p)-1)){
  Quantiles_gamma1[i]=F_inv(p[i])
}

g=0.8

F = function(y,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  1- gamma*pnorm(1,mu1,sigma1)-(1-gamma)*pnorm(1,mu2,sigma2) + gamma * pnorm(y, mu1, sigma1) + (1-gamma) * pnorm(y, mu2, sigma2)
}
F_inv = function(q,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  minmax <- c(0,1)
  uniroot(function(x) F(x,mu1,sigma1,g,mu2,sigma2)-q,
          interval = minmax,
          tol = 10^{-16})$root  
}
Quantiles_gamma08=c(1:length(p))
Quantiles_gamma08[1]=0
Quantiles_gamma08[2]=0
Quantiles_gamma08[length(p)]=1
for(i in 3:(length(p)-1)){
  Quantiles_gamma08[i]=F_inv(p[i])
}

g=0.6

F = function(y,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  1- gamma*pnorm(1,mu1,sigma1)-(1-gamma)*pnorm(1,mu2,sigma2) + gamma * pnorm(y, mu1, sigma1) + (1-gamma) * pnorm(y, mu2, sigma2)
}
F_inv = function(q,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  minmax <- c(0,1)
  uniroot(function(x) F(x,mu1,sigma1,g,mu2,sigma2)-q,
          interval = minmax,
          tol = 10^{-16})$root  
}
Quantiles_gamma06=c(1:length(p))
Quantiles_gamma06[1]=0
Quantiles_gamma06[length(p)]=1
for(i in 2:(length(p)-1)){
  Quantiles_gamma06[i]=F_inv(p[i])
}

Quantiles_wx=melt(data.frame( p = p, gamma1 =Quantiles_gamma1, gamma08 =Quantiles_gamma08, gamma06 =Quantiles_gamma06),id.vars = 'p')
ggplot(Quantiles_wx, aes(x=p, y=value, colour=variable)) +
  geom_line(size=1)+
  theme_bw() +
  labs(title="", x='p', y=TeX('F_{w(X)}^{-1}(p)'), colour=" ") +
  theme(axis.title=element_text(size=12), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=12)) +
  scale_colour_manual(name=" ", values = c( "#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))

C1=function(u1,u2){
  return(pCopula(c(u1,u2), frankCopula(-10)))
}

C2=function(u1,u2){
  return(pCopula(c(u1,u2), frankCopula(10)))
}

mu1=0.15
mu2=0.95
sigma1=0.05
sigma2=0.01
a=mu2*(pnorm((1-mu2)/sigma2)-pnorm(-mu2/sigma2))+sigma2*(dnorm(-mu2/sigma2)-dnorm((1-mu2)/sigma2))
b=mu1*(pnorm((1-mu1)/sigma1)-pnorm(-mu1/sigma1))+sigma1*(dnorm(-mu1/sigma1)-dnorm((1-mu1)/sigma1)) - (mu2*(pnorm((1-mu2)/sigma2)-pnorm(-mu2/sigma2))+sigma2*(dnorm(-mu2/sigma2)-dnorm((1-mu2)/sigma2)))

C_A=function(u1,u2, gamma){
  return((a+b*gamma)*C1(u1,u2)+(1-(a+b*gamma))*C2(u1,u2))
}

# Quantiles of w(X)
g=1

F = function(y,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  1- gamma*pnorm(1,mu1,sigma1)-(1-gamma)*pnorm(1,mu2,sigma2) + gamma * pnorm(y, mu1, sigma1) + (1-gamma) * pnorm(y, mu2, sigma2)
}
F_inv = function(q,mu1=0.15,sigma1=0.05,gamma=g,mu2=0.95,sigma2=0.01){
  minmax <- c(0,1)
  uniroot(function(x) F(x,mu1,sigma1,g,mu2,sigma2)-q,
          interval = minmax,
          tol = 10^{-16})$root  
}

C_med=function(u1,u2, gamma){
  if(gamma==1){
    val=0.15*C1(u1,u2)+(1-0.15)*C2(u1,u2)
  }
  if(gamma==0.8){
    val=0.1659*C1(u1,u2)+(1-0.1659)*C2(u1,u2)
  }
  if(gamma==0.6){
    val=0.1984*C1(u1,u2)+(1-0.1984)*C2(u1,u2)
  }
  return(val)
}

C_q09=function(u1,u2, gamma){
  if(gamma==1){
    val=0.0859*C1(u1,u2)+(1-0.0859)*C2(u1,u2)
  }
  if(gamma==0.8){
    val=0.0925*C1(u1,u2)+(1-0.0925)*C2(u1,u2)
  }
  if(gamma==0.6){
    val=0.1016*C1(u1,u2)+(1-0.1016)*C2(u1,u2)
  }
  return(val)
}

C_q01=function(u1,u2, gamma){
  if(gamma==1){
    val=0.2141*C1(u1,u2)+(1-0.2141)*C2(u1,u2)
  }
  if(gamma==0.8){
    val=0.9500*C1(u1,u2)+(1-0.9500)*C2(u1,u2)
  }
  if(gamma==0.6){
    val=0.9567*C1(u1,u2)+(1-0.9567)*C2(u1,u2)
  }
  return(val)
}

u1=seq(0,1,0.01)
u2=u1
m=length(u1)

Av_cop_gamma1=matrix(nrow=m,ncol=m)
Av_cop_gamma08=matrix(nrow=m,ncol=m)
Av_cop_gamma06=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Av_cop_gamma1[i,j]=C_A(u1[i],u2[j], 1)
    Av_cop_gamma08[i,j]=C_A(u1[i],u2[j], 0.8)
    Av_cop_gamma06[i,j]=C_A(u1[i],u2[j], 0.6)
  }
}

Av_copula=rbind(cbind(melt(Av_cop_gamma1),"Gamma"=rep('1',nrow(melt(Av_cop_gamma1)))) , cbind(melt(Av_cop_gamma08),"Gamma"=rep('0.8',nrow(melt(Av_cop_gamma08)))) , cbind(melt(Av_cop_gamma06),"Gamma"=rep('0.6',nrow(melt(Av_cop_gamma06)))))
ggplot(Av_copula, aes(x=u1[Var1], y=u2[Var2], z=value,  colour=Gamma)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, linetype='dashed') +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='', colour='') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))

Med_cop_gamma1=matrix(nrow=m,ncol=m)
Med_cop_gamma08=matrix(nrow=m,ncol=m)
Med_cop_gamma06=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Med_cop_gamma1[i,j]=C_med(u1[i],u2[j], 1)
    Med_cop_gamma08[i,j]=C_med(u1[i],u2[j], 0.8)
    Med_cop_gamma06[i,j]=C_med(u1[i],u2[j], 0.6)
  }
}

Med_copula=rbind(cbind(melt(Med_cop_gamma1),"Gamma"=rep('1',nrow(melt(Med_cop_gamma1)))) , cbind(melt(Med_cop_gamma08),"Gamma"=rep('0.8',nrow(melt(Med_cop_gamma08)))) , cbind(melt(Med_cop_gamma06),"Gamma"=rep('0.6',nrow(melt(Med_cop_gamma06)))))
ggplot(Med_copula, aes(x=u1[Var1], y=u2[Var2], z=value,  colour=Gamma)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1,  linetype='solid') +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='', colour='') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))

Q01_cop_gamma1=matrix(nrow=m,ncol=m)
Q01_cop_gamma08=matrix(nrow=m,ncol=m)
Q01_cop_gamma06=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Q01_cop_gamma1[i,j]=C_q01(u1[i],u2[j], 1)
    Q01_cop_gamma08[i,j]=C_q01(u1[i],u2[j], 0.8)
    Q01_cop_gamma06[i,j]=C_q01(u1[i],u2[j], 0.6)
  }
}

Q01_copula=rbind(cbind(melt(Q01_cop_gamma1),"Gamma"=rep('1',nrow(melt(Q01_cop_gamma1)))) , cbind(melt(Q01_cop_gamma08),"Gamma"=rep('0.8',nrow(melt(Q01_cop_gamma08)))) , cbind(melt(Q01_cop_gamma06),"Gamma"=rep('0.6',nrow(melt(Q01_cop_gamma06)))))
ggplot(Q01_copula, aes(x=u1[Var1], y=u2[Var2], z=value,  colour=Gamma)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, linetype='dotted') +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='', colour='') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))

Q09_cop_gamma1=matrix(nrow=m,ncol=m)
Q09_cop_gamma08=matrix(nrow=m,ncol=m)
Q09_cop_gamma06=matrix(nrow=m,ncol=m)
for (i in 1:m){
  for (j in 1:m){
    Q09_cop_gamma1[i,j]=C_q09(u1[i],u2[j], 1)
    Q09_cop_gamma08[i,j]=C_q09(u1[i],u2[j], 0.8)
    Q09_cop_gamma06[i,j]=C_q09(u1[i],u2[j], 0.6)
  }
}

Q09_copula=rbind(cbind(melt(Q09_cop_gamma1),"Gamma"=rep('1',nrow(melt(Q09_cop_gamma1)))) , cbind(melt(Q09_cop_gamma08),"Gamma"=rep('0.8',nrow(melt(Q09_cop_gamma08)))) , cbind(melt(Q09_cop_gamma06),"Gamma"=rep('0.6',nrow(melt(Q09_cop_gamma06)))))
ggplot(Q09_copula, aes(x=u1[Var1], y=u2[Var2], z=value,  colour=Gamma)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, linetype='dotdash') +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='', colour='') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('$\\gamma = 1', '$\\gamma = 0.8', '$\\gamma = 0.6'))))


u1=0.5 ; u2=0.1
round(C_A(u1,u2,1),4); round(C_A(u1,u2,0.8),4); round(C_A(u1,u2,0.6),4)
round(C_q01(u1,u2,1),4); round(C_q01(u1,u2,0.8),4); round(C_q01(u1,u2,0.6),4)
round(C_med(u1,u2,1),4); round(C_med(u1,u2,0.8),4); round(C_med(u1,u2,0.6),4)
round(C_q09(u1,u2,1),4); round(C_q09(u1,u2,0.8),4); round(C_q09(u1,u2,0.6),4)


#### Uranium data example ####

### original ###

data(uranium)

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57

data=as.data.frame(cbind(Y1,Y2,X))
colnames(data)=c("Co","Sc","Ti")

# with continuous color gradation
# use theme_gray() 
# and use scale_color_viridis(option = "B") 
p1c=ggplot(data, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.5,1.5)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p2c=ggplot(data, aes(x=Ti, y=Co, color=Sc))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(2.8,4.4)+
  ylim(0.5,1.5)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p3c=ggplot(data, aes(x=Ti, y=Sc, color=Co))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(2.8,4.4)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

# with discrete color gradation
# discretization by splitting up the data such that the groups are of equal size or such that the classes are of equal length
n/3;2*n/3
breaks_Co_1=c(min(Y1)-0.0001, Y1[order(Y1)][c(219,437,n)])
breaks_Sc_1=c(min(Y2)-0.0001, Y2[order(Y2)][c(219,437,n)])
breaks_Ti_1=c(min(X)-0.0001, X[order(X)][c(219,437,n)])

inferno(7, alpha=0.5)

p1d=ggplot(data, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.5,1.5)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))
p2d=ggplot(data, aes(x=Ti, y=Co, color=cut(Sc, breaks=breaks_Sc_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.32,0.94]","(0.94,1.09]","(1.09,1.46]"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"),labels=c("(0.32,0.94]","(0.94,1.09]","(1.09,1.46]"))+
  labs(color="Sc")+
  stat_ellipse(aes(fill= cut(Sc, breaks=breaks_Sc_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(2.8,4.4)+
  ylim(0.5,1.5)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))
p3d=ggplot(data, aes(x=Ti, y=Sc, color=cut(Co, breaks=breaks_Co_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.56,0.98]","(0.98,1.08]","(1.08,1.45]"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.56,0.98]","(0.98,1.08]","(1.08,1.45]"))+
  labs(color="Co")+
  stat_ellipse(aes(fill= cut(Co, breaks=breaks_Co_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(2.8,4.4)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

plot_grid(p1c,p2c,p3c,p1d,p2d,p3d,nrow=2)


# Copula's

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57
u=seq(0,1,0.05)
C_X=CondCopula_X(Y1,Y2,X,u,h, 'LL')

v=u
m=length(u)
n=length(X)
med_cop=matrix(ncol=m, nrow=m)
av_cop=med_cop
q25_cop=med_cop  
q75_cop=med_cop
for (i in 1:m){
  for (j in 1:m){
    med_cop[i,j]=median(C_X$CondCopula_X[i,j,])
    av_cop[i,j]=mean(C_X$CondCopula_X[i,j,])
    q25_cop[i,j]=quantile(C_X$CondCopula_X[i,j,],0.25)
    q75_cop[i,j]=quantile(C_X$CondCopula_X[i,j,],0.75)
  }
}

AvMedCopula=rbind(cbind(melt(av_cop),"Copula"=rep('Average',nrow(melt(av_cop)))) , cbind(melt(med_cop),"Copula"=rep('Median',nrow(melt(med_cop)))))
QuartileCopula=rbind(cbind(melt(q25_cop),"Copula"=rep('First quartile',nrow(melt(q25_cop)))), cbind(melt(med_cop),"Copula"=rep('Median',nrow(melt(med_cop)))), cbind(melt(q75_cop),"Copula"=rep('Third quartile',nrow(melt(q75_cop)))))

p6=ggplot(AvMedCopula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, colour="#00BFC4") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='Copula', colour='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dashed","solid"))

p7=ggplot(QuartileCopula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1,  colour="#00BFC4") +
  theme_bw() +
  labs(title=" ", x=TeX('u_1'), y='v', linetype='Copula', colour='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"))

plot_grid(p6,p7,nrow=1)


# Kendall's tau

cor(Y1,Y2,method="kendall")
Tau_X_obs=CondKendallsTau_X_Obs(Y1=Y1,Y2=Y2,X=X,h=h,weights='LL')
AvCondKendallsTau=mean(Tau_X_obs$CondKendallsTau_X_Obs)
MedianCondKendallsTau=median(Tau_X_obs$CondKendallsTau_X_Obs)
PartialKendallsTau=cor(Tau_X_obs$Y1.adj, Tau_X_obs$Y2.adj, method="kendall")

AvCondKendallsTau  # 0.3796397
MedianCondKendallsTau  # 0.4186422
PartialKendallsTau  # 0.3929827

x=seq(3,4.2,0.01)
Tau_X=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Tau1=data.frame(Tau_X, x)  
Tau1_sum <- data.frame( x = x, av =AvCondKendallsTau, med=MedianCondKendallsTau, avl='Average', medl='Median')
ggplot(Tau1, aes(x=x, y=Tau_X)) +
  geom_line(size=1, colour='#00BFC4')+
  theme_bw() +
  ylim(-0.15,0.68)+
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  geom_line(aes( x, av ,  linetype=avl), Tau1_sum, size=0.9, colour='gray70') +
  geom_line(aes( x, med , linetype=medl), Tau1_sum, size=0.9, colour='gray70') +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=13)) +
  scale_linetype_manual(name=" ", values = c("dashed", "solid"), labels=unname(TeX(c('$\\tau_n^A', '$\\tau_n^{Med}'))))


### 10% outliers ###

data(uranium)

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57


### create ouliers in X space and in dependency between Y1 and Y2
### by replacying some datapoints by outliers
# outlier proportion
n_out=floor(n*0.1)
# outliers in dependency between Y1 and Y2
covmat=matrix(c(0.002,-0.002,-0.003 ,-0.002, 0.003,0.003 ,-0.003 ,0.003 , 0.005),3,3)
set.seed(12)
# outliers in X (4)
out=mvrnorm(n = n_out, mu = c(1,0.9,4), Sigma = covmat)  
par(mfrow=c(1,3))
plot(out[,1],out[,2])
plot(out[,3],out[,1])
plot(out[,3],out[,2])
set.seed(12)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=out[,1]
Y2[out_ind]=out[,2]
X[out_ind]=out[,3]
n=length(X)
data=as.data.frame(cbind(Y1,Y2,X))
colnames(data)=c("Co","Sc","Ti")

# with continuous color gradation
# use theme_gray() 
# and use scale_color_viridis(option = "B") 
p1c=ggplot(data, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.5,1.5)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p2c=ggplot(data, aes(x=Ti, y=Co, color=Sc))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(2.8,4.4)+
  ylim(0.5,1.5)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p3c=ggplot(data, aes(x=Ti, y=Sc, color=Co))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(2.8,4.4)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

# with discrete color gradation
# discretization by splitting up the data such that the groups are of equal size or such that the classes are of equal length
n/3;2*n/3
breaks_Co_1=c(min(Y1)-0.0001, Y1[order(Y1)][c(219,437,n)])
breaks_Sc_1=c(min(Y2)-0.0001, Y2[order(Y2)][c(219,437,n)])
breaks_Ti_1=c(min(X)-0.0001, X[order(X)][c(219,437,n)])

inferno(7, alpha=0.5)

p1d=ggplot(data, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.5,1.5)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))
p2d=ggplot(data, aes(x=Ti, y=Co, color=cut(Sc, breaks=breaks_Sc_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.32,0.92]","(0.92,1.07]","(1.07,1.46]"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"),labels=c("(0.32,0.92]","(0.92,1.07]","(1.07,1.46]"))+
  labs(color="Sc")+
  stat_ellipse(aes(fill= cut(Sc, breaks=breaks_Sc_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(2.8,4.4)+
  ylim(0.5,1.5)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))
p3d=ggplot(data, aes(x=Ti, y=Sc, color=cut(Co, breaks=breaks_Co_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.56,0.98]","(0.98,1.07]","(1.07,1.42]"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.56,0.98]","(0.98,1.07]","(1.07,1.42]"))+
  labs(color="Co")+
  stat_ellipse(aes(fill= cut(Co, breaks=breaks_Co_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(2.8,4.4)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

plot_grid(p1c,p2c,p3c,p1d,p2d,p3d,nrow=2)


# Copula's

u=seq(0,1,0.05)
C_X=CondCopula_X(Y1,Y2,X,u,h, 'LL')

v=u
m=length(u)
n=length(X)
med_cop=matrix(ncol=m, nrow=m)
av_cop=med_cop
q25_cop=med_cop  
q75_cop=med_cop
for (i in 1:m){
  for (j in 1:m){
    med_cop[i,j]=median(C_X$CondCopula_X[i,j,])
    av_cop[i,j]=mean(C_X$CondCopula_X[i,j,])
    q25_cop[i,j]=quantile(C_X$CondCopula_X[i,j,],0.25)
    q75_cop[i,j]=quantile(C_X$CondCopula_X[i,j,],0.75)
  }
}

AvMedCopula=rbind(cbind(melt(av_cop),"Copula"=rep('Average',nrow(melt(av_cop)))) , cbind(melt(med_cop),"Copula"=rep('Median',nrow(melt(med_cop)))))
QuartileCopula=rbind(cbind(melt(q25_cop),"Copula"=rep('First quartile',nrow(melt(q25_cop)))), cbind(melt(med_cop),"Copula"=rep('Median',nrow(melt(med_cop)))), cbind(melt(q75_cop),"Copula"=rep('Third quartile',nrow(melt(q75_cop)))))

p6=ggplot(AvMedCopula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Copula,  colour=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1) +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='Copula', colour='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("solid","dotted")) +
  scale_color_manual(values = c("#F8766D", "#00BFC4"))

p7=ggplot(QuartileCopula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Copula,  colour=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1) +
  theme_bw() +
  labs(title=" ", x=TeX('u_1'), y='v', linetype='Copula', colour='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dashed", "dotted", "dotdash")) +
  scale_color_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF"))

plot_grid(p6,p7,nrow=1)


# Kendall's tau

cor(Y1,Y2,method="kendall")
Tau_X_obs=CondKendallsTau_X_Obs(Y1=Y1,Y2=Y2,X=X,h=h,weights='LL')
AvCondKendallsTau=mean(Tau_X_obs$CondKendallsTau_X_Obs)
MedianCondKendallsTau=median(Tau_X_obs$CondKendallsTau_X_Obs)
PartialKendallsTau=cor(Tau_X_obs$Y1.adj, Tau_X_obs$Y2.adj, method="kendall")

AvCondKendallsTau  # 0.4145162
MedianCondKendallsTau  # 0.4338301
PartialKendallsTau  # 0.4294278

x=seq(3,4.2,0.01)
Tau_X_10=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Tau1=data.frame(Tau_X_10, x)  
Tau1_sum <- data.frame( x = x, av =AvCondKendallsTau, med=MedianCondKendallsTau, partial=PartialKendallsTau, avl='Average', medl='Median', partiall='Partial')
ggplot(Tau1, aes(x=x, y=Tau_X_10)) +
  geom_line(size=1, colour='black')+
  theme_bw() +
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  geom_line(aes( x, av ,  linetype=avl), Tau1_sum, size=0.7, colour='gray70') +
  geom_line(aes( x, med , linetype=medl), Tau1_sum, size=0.7, colour='gray70') +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=13)) +
  scale_linetype_manual(name=" ", values = c("dashed", "solid"), labels=unname(TeX(c('$\\tau_n^A', '$\\tau_n^{Med}'))))


### 20% outliers ###

data(uranium)

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57


### create ouliers in X space and in dependency between Y1 and Y2
### by replacying some datapoints by outliers
# outlier proportion
n_out=floor(n*0.2)
# outliers in dependency between Y1 and Y2
covmat=matrix(c(0.002,-0.002,-0.003 ,-0.002, 0.003,0.003 ,-0.003 ,0.003 , 0.005),3,3)
set.seed(12)
# outliers in X (4)
out=mvrnorm(n = n_out, mu = c(1,0.9,4), Sigma = covmat)  
par(mfrow=c(1,3))
plot(out[,1],out[,2])
plot(out[,3],out[,1])
plot(out[,3],out[,2])
set.seed(12)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=out[,1]
Y2[out_ind]=out[,2]
X[out_ind]=out[,3]
n=length(X)
data=as.data.frame(cbind(Y1,Y2,X))
colnames(data)=c("Co","Sc","Ti")

# with continuous color gradation
# use theme_gray() 
# and use scale_color_viridis(option = "B") 
p1c=ggplot(data, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.5,1.5)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p2c=ggplot(data, aes(x=Ti, y=Co, color=Sc))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(2.8,4.4)+
  ylim(0.5,1.5)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p3c=ggplot(data, aes(x=Ti, y=Sc, color=Co))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(2.8,4.4)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

# with discrete color gradation
# discretization by splitting up the data such that the groups are of equal size or such that the classes are of equal length
n/3;2*n/3
breaks_Co_1=c(min(Y1)-0.0001, Y1[order(Y1)][c(219,437,n)])
breaks_Sc_1=c(min(Y2)-0.0001, Y2[order(Y2)][c(219,437,n)])
breaks_Ti_1=c(min(X)-0.0001, X[order(X)][c(219,437,n)])

inferno(7, alpha=0.5)

p1d=ggplot(data, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.5,1.5)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))
p2d=ggplot(data, aes(x=Ti, y=Co, color=cut(Sc, breaks=breaks_Sc_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.32,0.91]","(0.91,1.05]","(1.05,1.46]"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"),labels=c("(0.32,0.91]","(0.91,1.05]","(1.05,1.46]"))+
  labs(color="Sc")+
  stat_ellipse(aes(fill= cut(Sc, breaks=breaks_Sc_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(2.8,4.4)+
  ylim(0.5,1.5)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))
p3d=ggplot(data, aes(x=Ti, y=Sc, color=cut(Co, breaks=breaks_Co_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.56,0.98]","(0.98,1.06]","(1.06,1.42]"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"), labels=c("(0.56,0.98]","(0.98,1.06]","(1.06,1.42]"))+
  labs(color="Co")+
  stat_ellipse(aes(fill= cut(Co, breaks=breaks_Co_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(2.8,4.4)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

plot_grid(p1c,p2c,p3c,p1d,p2d,p3d,nrow=2)


# Copula's

u=seq(0,1,0.05)
C_X=CondCopula_X(Y1,Y2,X,u,h, 'LL')

v=u
m=length(u)
n=length(X)
med_cop=matrix(ncol=m, nrow=m)
av_cop=med_cop
q25_cop=med_cop  
q75_cop=med_cop
for (i in 1:m){
  for (j in 1:m){
    med_cop[i,j]=median(C_X$CondCopula_X[i,j,])
    av_cop[i,j]=mean(C_X$CondCopula_X[i,j,])
    q25_cop[i,j]=quantile(C_X$CondCopula_X[i,j,],0.25)
    q75_cop[i,j]=quantile(C_X$CondCopula_X[i,j,],0.75)
  }
}

AvMedCopula=rbind(cbind(melt(av_cop),"Copula"=rep('Average',nrow(melt(av_cop)))) , cbind(melt(med_cop),"Copula"=rep('Median',nrow(melt(med_cop)))))
QuartileCopula=rbind(cbind(melt(q25_cop),"Copula"=rep('First quartile',nrow(melt(q25_cop)))), cbind(melt(med_cop),"Copula"=rep('Median',nrow(melt(med_cop)))), cbind(melt(q75_cop),"Copula"=rep('Third quartile',nrow(melt(q75_cop)))))

p6=ggplot(AvMedCopula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1,  colour="#00BFC4") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dashed","solid"))

p7=ggplot(QuartileCopula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, colour="#00BFC4") +
  theme_bw() +
  labs(title=" ", x=TeX('u_1'), y='v', linetype='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"))

plot_grid(p6,p7,nrow=1)


# Kendall's tau

cor(Y1,Y2,method="kendall")
Tau_X_obs=CondKendallsTau_X_Obs(Y1=Y1,Y2=Y2,X=X,h=h,weights='LL')
AvCondKendallsTau=mean(Tau_X_obs$CondKendallsTau_X_Obs)
MedianCondKendallsTau=median(Tau_X_obs$CondKendallsTau_X_Obs)
PartialKendallsTau=cor(Tau_X_obs$Y1.adj, Tau_X_obs$Y2.adj, method="kendall")

AvCondKendallsTau  # 0.3453647
MedianCondKendallsTau  # 0.4034579
PartialKendallsTau  # 0.3769031

x=seq(3,4.2,0.01)
Tau_X_20=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Tau1=data.frame(Tau_X_20, x)  
Tau1_sum <- data.frame( x = x, av =AvCondKendallsTau, med=MedianCondKendallsTau, avl='Average', medl='Median')
ggplot(Tau1, aes(x=x, y=Tau_X_20)) +
  geom_line(size=1, colour="#7CAE00")+
  theme_bw() +
  ylim(-0.15,0.68)+
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  geom_line(aes( x, av ,  linetype=avl), Tau1_sum, size=0.9, colour='gray70') +
  geom_line(aes( x, med , linetype=medl), Tau1_sum, size=0.9, colour='gray70') +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=13)) +
  scale_linetype_manual(name=" ", values = c("dashed", "solid"), labels=unname(TeX(c('$\\tau_n^A', '$\\tau_n^{Med}'))))


# Conditional Kendall's tau in three cases on one plot

Tau_X_all=melt(data.frame(x, Tau_X, Tau_X_10, Tau_X_20), id.vars='x')
ggplot(Tau_X_all, aes(x=x, y=value, colour=variable)) +
  geom_line(size=1, linetype='solid')+
  theme_bw() +
  ylim(-0.25,0.68)+
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=12)) +
  scale_color_manual(name=" ", values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c('Original', '10% outliers', '20% outliers'))


#### Take different mu ####

### original ###

data(uranium)

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57

# Kendall's tau

x=seq(3,4.2,0.01)
Tau_X=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

### 10% outliers ###

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57

### create ouliers in X space and in dependency between Y1 and Y2
### by replacying some datapoints by outliers
# outlier proportion
n_out=floor(n*0.1)
# outliers in dependency between Y1 and Y2
covmat=matrix(c(0.002,-0.002,-0.003 ,-0.002, 0.003,0.003 ,-0.003 ,0.003 , 0.005),3,3)
set.seed(12)
# outliers in X (4)
out=mvrnorm(n = n_out, mu = c(1.12,1.15,4), Sigma = covmat)  
par(mfrow=c(1,3))
plot(out[,1],out[,2])
plot(out[,3],out[,1])
plot(out[,3],out[,2])
set.seed(12)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=out[,1]
Y2[out_ind]=out[,2]
X[out_ind]=out[,3]
n=length(X)

# Kendall's tau

x=seq(3,4.2,0.01)
Tau_X_10_2=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

### 20% outliers ###

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57


### create ouliers in X space and in dependency between Y1 and Y2
### by replacying some datapoints by outliers
# outlier proportion
n_out=floor(n*0.2)
# outliers in dependency between Y1 and Y2
covmat=matrix(c(0.002,-0.002,-0.003 ,-0.002, 0.003,0.003 ,-0.003 ,0.003 , 0.005),3,3)
set.seed(12)
# outliers in X (4)
out=mvrnorm(n = n_out, mu = c(1.12,1.15,4), Sigma = covmat)  
par(mfrow=c(1,3))
plot(out[,1],out[,2])
plot(out[,3],out[,1])
plot(out[,3],out[,2])
set.seed(12)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=out[,1]
Y2[out_ind]=out[,2]
X[out_ind]=out[,3]
n=length(X)

# Kendall's tau

x=seq(3,4.2,0.01)
Tau_X_20_2=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Tau_X_all_2=melt(data.frame(x, Tau_X, Tau_X_10_2, Tau_X_20_2), id.vars='x')
ggplot(Tau_X_all_2, aes(x=x, y=value, colour=variable)) +
  geom_line(size=1, linetype='solid')+
  theme_bw() +
  ylim(-0.25,0.68)+
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=12)) +
  scale_color_manual(name=" ", values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=c('Original', '10% outliers', '20% outliers'))


### Same code but including outliers in Y1 ###

data(uranium)

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57

data=as.data.frame(cbind(Y1,Y2,X))  # original data
breaks_Ti=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

n_out=floor(n*0.1)
set.seed(120)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=Y1[out_ind]+0.7

data_Y110=as.data.frame(cbind(Y1,Y2,X)) # original data with 10% outliers in Y1
breaks_Ti_Y110=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau_Y110=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X_Y110=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Y1=uranium$Co
n_out=floor(n*0.2)
set.seed(120)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=Y1[out_ind]+0.7

data_Y120=as.data.frame(cbind(Y1,Y2,X)) # original data with 20% outliers in Y1
breaks_Ti_Y120=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau_Y120=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X_Y120=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57
n_out=floor(n*0.1)
covmat=matrix(c(0.002,-0.002,-0.003 ,-0.002, 0.003,0.003 ,-0.003 ,0.003 , 0.005),3,3)
set.seed(12)
out=mvrnorm(n = n_out, mu = c(1,0.9,4), Sigma = covmat)
set.seed(12)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=out[,1]
Y2[out_ind]=out[,2]
X[out_ind]=out[,3]

data_X10=as.data.frame(cbind(Y1,Y2,X))  # dataset with 10% outliers in X
breaks_Ti_X10=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau_X10=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X_X10=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Y1=data_X10[,1]
Y2=data_X10[,2]
X=data_X10[,3]
n_out=floor(n*0.1)
set.seed(120)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=Y1[out_ind]+0.7

data_X10_Y110=as.data.frame(cbind(Y1,Y2,X)) # dataset with 10% outliers in both X and Y1
breaks_Ti_X10_Y110=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau_X10_Y110=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X_X10_Y110=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Y1=data_X10[,1]
Y2=data_X10[,2]
X=data_X10[,3]
n_out=floor(n*0.2)
set.seed(120)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=Y1[out_ind]+0.7

data_X10_Y120=as.data.frame(cbind(Y1,Y2,X)) # dataset with 10% outliers in X and 20% in Y1
breaks_Ti_X10_Y120=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau_X10_Y120=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X_X10_Y120=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Y1=uranium$Co
Y2=uranium$Sc
X=uranium$Ti
n=length(Y1)
h=0.57
n_out=floor(n*0.2)
covmat=matrix(c(0.002,-0.002,-0.003 ,-0.002, 0.003,0.003 ,-0.003 ,0.003 , 0.005),3,3)
set.seed(12)
out=mvrnorm(n = n_out, mu = c(1,0.9,4), Sigma = covmat)
set.seed(12)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=out[,1]
Y2[out_ind]=out[,2]
X[out_ind]=out[,3]

data_X20=as.data.frame(cbind(Y1,Y2,X))  # dataset with 20% outliers in X
breaks_Ti_X20=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau_X20=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X_X20=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Y1=data_X20[,1]
Y2=data_X20[,2]
X=data_X20[,3]
n_out=floor(n*0.1)
set.seed(120)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=Y1[out_ind]+0.7

data_X20_Y110=as.data.frame(cbind(Y1,Y2,X)) # dataset with 20% outliers in X and 10% in Y1
breaks_Ti_X20_Y110=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau_X20_Y110=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X_X20_Y110=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Y1=data_X20[,1]
Y2=data_X20[,2]
X=data_X20[,3]
n_out=floor(n*0.2)
set.seed(120)
out_ind=sample(x = c(1:655), n_out, replace = FALSE)
Y1[out_ind]=Y1[out_ind]+0.7

data_X20_Y120=as.data.frame(cbind(Y1,Y2,X)) # dataset with 20% outliers in X and 20% in Y1
breaks_Ti_X20_Y120=c(min(X)-0.0001, X[order(X)][c(219,437,n)])
Tau_X20_Y120=cor(Y1,Y2,method="kendall")
x=seq(3,4.2,0.01)
Tau_X_X20_Y120=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

colnames(data)=c("Co","Sc","Ti")
colnames(data_Y110)=c("Co","Sc","Ti")
colnames(data_Y120)=c("Co","Sc","Ti")
colnames(data_X10)=c("Co","Sc","Ti")
colnames(data_X10_Y110)=c("Co","Sc","Ti")
colnames(data_X10_Y120)=c("Co","Sc","Ti")
colnames(data_X20)=c("Co","Sc","Ti")
colnames(data_X20_Y110)=c("Co","Sc","Ti")
colnames(data_X20_Y120)=c("Co","Sc","Ti")

# Plot for original data
p1c=ggplot(data, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p2c=ggplot(data_Y110, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p3c=ggplot(data_Y120, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p1d=ggplot(data, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

p2d=ggplot(data_Y110, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_Y110)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_Y110)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

p3d=ggplot(data_Y120, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_Y120)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_Y120)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

plot_grid(p1c,p2c,p3c,p1d,p2d,p3d,nrow=2)

# Plot for 10% outliers in X

p1c=ggplot(data_X10, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p2c=ggplot(data_X10_Y110, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p3c=ggplot(data_X10_Y120, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p1d=ggplot(data_X10, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_X10)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_X10)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

p2d=ggplot(data_X10_Y110, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_X10_Y110)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_X10_Y110)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

p3d=ggplot(data_X10_Y120, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_X10_Y120)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_X10_Y120)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

plot_grid(p1c,p2c,p3c,p1d,p2d,p3d,nrow=2)

# Plot for 20% outliers in X

p1c=ggplot(data_X20, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p2c=ggplot(data_X20_Y110, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p3c=ggplot(data_X20_Y120, aes(x=Co, y=Sc, color=Ti))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

p1d=ggplot(data_X20, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_X20)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_X20)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

p2d=ggplot(data_X20_Y110, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_X20_Y110)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_X20_Y110)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

p3d=ggplot(data_X20_Y120, aes(x=Co, y=Sc, color=cut(Ti, breaks=breaks_Ti_X20_Y120)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Ti")+
  stat_ellipse(aes(fill= cut(Ti, breaks=breaks_Ti_X20_Y120)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  xlim(0.2,2.3)+
  ylim(0.3,1.6)+
  theme_gray()+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=9))

plot_grid(p1c,p2c,p3c,p1d,p2d,p3d,nrow=2)


###### Kendall's tau ######

Tau_X_all=melt(data.frame(x, Tau_X, Tau_X_Y110, Tau_X_Y120), id.vars='x')
ggplot(Tau_X_all, aes(x=x, y=value, colour= variable)) +
  geom_line(size=1, linetype='solid')+
  geom_hline(yintercept=Tau, linetype="twodash", size=0.7, colour="#00BFC4")+
  geom_hline(yintercept=Tau_Y110, linetype="twodash", size=0.7, colour="#F8766D")+
  geom_hline(yintercept=Tau_Y120, linetype="twodash", size=0.7, colour="#7CAE00")+
  theme_bw() +
  ylim(-0.25,0.68)+
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=12)) +
  scale_color_manual(name=" ", values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('No outliers in Y_1', '10% outliers in Y_1', '20% outliers in Y_1'))))


Tau_X10_all=melt(data.frame(x, Tau_X_X10, Tau_X_X10_Y110, Tau_X_X10_Y120), id.vars='x')
ggplot(Tau_X10_all, aes(x=x, y=value, colour= variable)) +
  geom_line(size=1, linetype='solid')+
  geom_hline(yintercept=Tau_X10, linetype="twodash", size=0.7, colour="#00BFC4")+
  geom_hline(yintercept=Tau_X10_Y110, linetype="twodash", size=0.7, colour="#F8766D")+
  geom_hline(yintercept=Tau_X10_Y120, linetype="twodash", size=0.7, colour="#7CAE00")+
  theme_bw() +
  ylim(-0.25,0.68)+
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=12)) +
  scale_color_manual(name=" ", values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('No outliers in Y_1', '10% outliers in Y_1', '20% outliers in Y_1'))))

Tau_X20_all=melt(data.frame(x, Tau_X_X20, Tau_X_X20_Y110, Tau_X_X20_Y120), id.vars='x')
ggplot(Tau_X20_all, aes(x=x, y=value, colour= variable)) +
  geom_line(size=1, linetype='solid')+
  geom_hline(yintercept=Tau_X20, linetype="twodash", size=0.7, colour="#00BFC4")+
  geom_hline(yintercept=Tau_X20_Y110, linetype="twodash", size=0.7, colour="#F8766D")+
  geom_hline(yintercept=Tau_X20_Y120, linetype="twodash", size=0.7, colour="#7CAE00")+
  theme_bw() +
  ylim(-0.25,0.68)+
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=12)) +
  scale_color_manual(name=" ", values = c("#00BFC4", "#F8766D", "#7CAE00"), labels=unname(TeX(c('No outliers in Y_1', '10% outliers in Y_1', '20% outliers in Y_1'))))


#### Second data example ####


##### Example to illustrate the comparison between two conditional Kendall's taus #####

data=read.csv(file.choose(), sep=';', dec=',', header = TRUE)
Liana_data=data[data$Form=='Liana',]
Tree_data=data[data$Form=='Tree',]

table(Liana_data$Species)
table(Tree_data$Species)

Y1=Liana_data[Liana_data$Species==levels(Liana_data$Species)[7],11]
Y2=Liana_data[Liana_data$Species==levels(Liana_data$Species)[7],14]
X=Liana_data[Liana_data$Species==levels(Liana_data$Species)[7],9]
hist(X)
h=(max(X)-min(X))/4
n=length(X)

data=as.data.frame(cbind(Y1,Y2,X))
colnames(data)=c("Cond","VpdL","Tleaf")

# with continuous color gradation
# use theme_gray() 
# and use scale_color_viridis(option = "B") 
p1c=ggplot(data, aes(x=Cond, y=VpdL, color=Tleaf))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_viridis(option = "B")+
  theme_gray()+
  ylim(0,5.5)+
  xlim(-0.3,0.9)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=11))

# with discrete color gradation
# discretization by splitting up the data such that the groups are of equal size or such that the classes are of equal length
n/3;2*n/3
breaks_X_1=c(min(X)-0.0001, X[order(X)][c(floor(n/3),floor(2*n/3),n)])

p1d=ggplot(data, aes(x=Cond, y=VpdL, color=cut(Tleaf, breaks=breaks_X_1)))+
  geom_point(size=2.5, alpha=0.5)+
  scale_color_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  scale_fill_manual(values=c("#330A5F80", "#BB375480", "#FCB51980"))+
  labs(color="Tleaf")+
  stat_ellipse(aes(fill= cut(Tleaf, breaks=breaks_X_1)), type="t", size=1.5, geom="polygon", level=0.99, alpha=0.1, show.legend = F)+
  theme_gray()+
  ylim(0,5.5)+
  xlim(-0.3,0.9)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.position="bottom", legend.title = element_text(size=13), 
        legend.text = element_text(size=9),legend.box.margin=margin(10,0,7,0))

plot_grid(p1c, p1d,ncol=2)


# Copula's

u=seq(0,1,0.05)
h=(max(X)-min(X))/3

C_X=CondCopula_X(Y1,Y2,X,u,h, 'LL')

v=u
m=length(u)
n=length(X)
med_cop=matrix(ncol=m, nrow=m)
av_cop=med_cop
q25_cop=med_cop  
q75_cop=med_cop
for (i in 1:m){
  for (j in 1:m){
    med_cop[i,j]=median(C_X$CondCopula_X[i,j,])
    av_cop[i,j]=mean(C_X$CondCopula_X[i,j,])
    q25_cop[i,j]=quantile(C_X$CondCopula_X[i,j,],0.25)
    q75_cop[i,j]=quantile(C_X$CondCopula_X[i,j,],0.75)
  }
}

AvMedCopula=rbind(cbind(melt(av_cop),"Copula"=rep('Average',nrow(melt(av_cop)))) , cbind(melt(med_cop),"Copula"=rep('Median',nrow(melt(med_cop)))))
QuartileCopula=rbind(cbind(melt(q25_cop),"Copula"=rep('First quartile',nrow(melt(q25_cop)))), cbind(melt(med_cop),"Copula"=rep('Median',nrow(melt(med_cop)))), cbind(melt(q75_cop),"Copula"=rep('Third quartile',nrow(melt(q75_cop)))))

p6=ggplot(AvMedCopula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1, colour="#00BFC4") +
  theme_bw() +
  labs(title="", x=TeX('u_1'), y='v', linetype='Copula', colour='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dashed","solid"))

p7=ggplot(QuartileCopula, aes(x=u[Var1], y=v[Var2], z=value, linetype=Copula)) +
  geom_contour(breaks=seq(0.1,0.9,length.out = 5), size=1,  colour="#00BFC4") +
  theme_bw() +
  labs(title=" ", x=TeX('u_1'), y='v', linetype='Copula', colour='Copula') +
  scale_y_continuous( TeX('u_2'), sec.axis = sec_axis(~ . * 1, name = " ",  breaks = seq(0.1,0.9,0.2))) +
  theme(axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank(), legend.position = "bottom") +
  scale_linetype_manual(values = c("dotted", "solid", "dotdash"))

plot_grid(p6,p7,nrow=1)


# Kendall's tau

Tau_X_obs=CondKendallsTau_X_Obs(Y1=Y1,Y2=Y2,X=X,h=h,weights='LL')
AvCondKendallsTau=mean(Tau_X_obs$CondKendallsTau_X_Obs)
MedianCondKendallsTau=median(Tau_X_obs$CondKendallsTau_X_Obs)

Tau=data.frame(Tau_X=Tau_X_obs$CondKendallsTau_X_Obs, X)  
ggplot(Tau, aes(x=X, y=Tau_X)) +
  geom_point(size=1, colour='#00BFC4')+
  theme_bw() +
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=13))

AvCondKendallsTau  # 0.1950547
MedianCondKendallsTau  # 0.2933494

min(X); max(X)
x=seq(round(min(X),2)+0.01,round(max(X),2)-0.01,0.01)
Tau_X_nonparam=CondKendallsTau_X(x,Y1,Y2,X,h,'LL')

Tau3=data.frame(Tau_X_nonparam, x)  
Tau3_sum <- data.frame( x = x, av =AvCondKendallsTau, med=MedianCondKendallsTau, avl='Average', medl='Median')
ggplot(Tau3, aes(x=x, y=Tau_X_nonparam)) +
  geom_line(size=1, colour='#00BFC4')+
  theme_bw() +
  labs(title="", x='x', y=TeX('$\\tau_n (x)')) +
  geom_line(aes( x, av ,  linetype=avl), Tau3_sum, size=0.9, colour='gray70') +
  geom_line(aes( x, med , linetype=medl), Tau3_sum, size=0.9, colour='gray70') +
  theme(axis.title=element_text(size=13), legend.position="bottom", legend.title = element_text(size=13), legend.text = element_text(size=13)) +
  scale_linetype_manual(name=" ", values = c("dashed", "solid"), labels=unname(TeX(c('$\\tau_n^A', '$\\tau_n^{Med}'))))


