########################################################
### This function calculates the conditional copula ####
############## in the observations of X. ###############
########################################################

#######################################################
############ Conditional copula function ##############
#######################################################

CondCopula_X=function(Y1,Y2,X,u,h, weights){
  
  n=length(X)
  
  ########## Calculation of the weights ############
  
  # scaled X
  X_scaled=matrix(0, nrow=n, ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      X_scaled[i,j]=(X[i]-X[j])/h
    }
  }
  
  # triweight kernel function
  Kernel=matrix(0,nrow=n,ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      if (abs(X_scaled[i,j]) <= 1){
        Kernel[i,j]=(35/32)*((1-((X_scaled[i,j])^2))^3)
      }
      else{
        Kernel[i,j]=0
      }
    }
  }
  
  
  # calculation of the weights depending on the chosen method (NW or LL)
  # NW weights
  if(weights=='NW'){
    
    weights=matrix(0,nrow=n,ncol=n)
    for (i in 1:n){
      for( j in 1:n){
        weights[j,i]=Kernel[i,j]/sum(Kernel[,j])
      }
    }
  # LL weights  
  } else if(weights=='LL'){
    
    # S function
    S0=c(1:n)
    S1=S0
    S2=S0
    for (i in 1:n){
      S0[i]=sum(Kernel[i,])/(n*h)
      S1[i]=sum(X_scaled[i,]*Kernel[i,])/(n*h)
      S2[i]=sum((X_scaled[i,])^2*Kernel[i,])/(n*h)
    }  
    
    # Local Linear weights
    weights=matrix(0,nrow=n,ncol=n)
    for (i in 1:n){
      for( j in 1:n){
        weights[i,j]=(Kernel[i,j]*(S2[i]-X_scaled[i,j]*S1[i]))/(n*h*(S0[i]*S2[i]-S1[i]^2))
      }
    }
  
    # ##### Via library 'locpol' #####
    # 
    # library(locpol)
    # weights=locLinWeightsC(x=X, xeval=X, bw=h, kernel=TriweigK)$locWeig
    
  } else {
    print('Choose Nadaraya-Watson (NW) or Local linear (LL) weights')
  }
    
  ########## Adjusted random variables Y1 and Y2 ############
  ### (sometimes denoted by \tilde{U}_1 and \tilde{U}_2 ) ###

  I_1=matrix(0,nrow=n, ncol=n);
  I_2=I_1;
  
  for (i in 1:n){
    for (j in 1:n){
      if(Y1[i]<=Y1[j]){
        I_1[i,j]=1
      }
      if(Y2[i]<=Y2[j]){
        I_2[i,j]=1
      }
    }
  }
  
  Y1_adj=c(1:n);
  Y2_adj=c(1:n);
  
  for (i in 1:n){
    Y1_adj[i]=sum(weights[i,]*I_1[,i]);
    Y2_adj[i]=sum(weights[i,]*I_2[,i]);
  }  
  
  
  ############### Calculation of conditional copula ##############

  # calculation of G_{1,X_i,h}(U_{1j}) and G_{2,X_i,h}(U_{2j})
  I_1_adj=matrix(0,nrow=n, ncol=n);
  I_2_adj=I_1_adj;
  
  for (i in 1:n){
    for (j in 1:n){
      if(Y1_adj[i]<=Y1_adj[j]){
        I_1_adj[i,j]=1
      }
      if(Y2_adj[i]<=Y2_adj[j]){
        I_2_adj[i,j]=1
      }
    }
  }
  
  G1xU1=matrix(nrow=n, ncol=n)
  G2xU2=matrix(nrow=n, ncol=n)
  
  for (i in 1:n){
    for (j in 1:n){
      G1xU1[i,j]=sum(weights[i,]*I_1_adj[,j])
      G2xU2[i,j]=sum(weights[i,]*I_2_adj[,j])
    }
    cat('l1: ', i)
  }
  
  # calculation of C_{X_i,h}(u,v)  
  Cond_cop=function(u,v,i){
    Cc=sum(ifelse(G1xU1[i,]<=u & G2xU2[i,]<=v,weights[i,],0))
    return(Cc)
  }
  
  v=u
  m=length(u)
  C_x=array(0, c(m,m,n))
  for(i in 1:m){
    for(j in 1:m){
      for(k in 1:n){
        C_x[i,j,k]=Cond_cop(u[i],v[j],k)
      }
    }
    cat('l2: ', i)
  }
  
  ######## Output ########
  
  out_list = list("CondCopula_X"=C_x, "Y1.adj" = Y1_adj, "Y2.adj" = Y2_adj)
  return(out_list)
  
}
