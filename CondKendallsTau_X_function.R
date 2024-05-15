########################################################
# This script calculates the conditional Kendall's tau #
#### for some value of X, instead of an observation of #
#### X (see CondKendallsTau_X_Obs_function.R) #####
########################################################


#######################################################
######### Conditional Kendall's tau function ##########
#######################################################

CondKendallsTau_X=function(x,Y1,Y2,X,h, weights){
  
  n=length(X)
  m=length(x)
  
  ########## Calculation of the weights ############
  
  # scaled X
  X_scaled=matrix(0, nrow=n, ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      X_scaled[i,j]=(X[i]-X[j])/h
    }
  }
  
  # triweight kernel function
  Kernel_X=matrix(0,nrow=n,ncol=n)
  for (i in 1:n){
    for (j in 1:n){
      if (abs(X_scaled[i,j]) <= 1){
        Kernel_X[i,j]=(35/32)*((1-((X_scaled[i,j])^2))^3)
      }
      else{
        Kernel_X[i,j]=0
      }
    }
  }
  
  # calculation of the weights depending on the chosen method (NW or LL)
  # NW weights
  if(weights=='NW'){
    
    weights_X=matrix(0,nrow=n,ncol=n)
    for (i in 1:n){
      for( j in 1:n){
        weights_X[j,i]=Kernel_X[i,j]/sum(Kernel_X[,j])
      }
    }
    # LL weights  
  } else if(weights=='LL'){
    
    # S function
    S0=c(1:n)
    S1=S0
    S2=S0
    for (i in 1:n){
      S0[i]=sum(Kernel_X[i,])/(n*h)
      S1[i]=sum(X_scaled[i,]*Kernel_X[i,])/(n*h)
      S2[i]=sum((X_scaled[i,])^2*Kernel_X[i,])/(n*h)
    }  
    
    # Local Linear weights ( w_i(X_j,h)=weights[i,j] )
    weights_X=matrix(0,nrow=n,ncol=n)
    for (i in 1:n){
      for( j in 1:n){
        weights_X[i,j]=(Kernel_X[i,j]*(S2[i]-X_scaled[i,j]*S1[i]))/(n*h*(S0[i]*S2[i]-S1[i]^2))
      }
    }
    
    # ##### Via library 'locpol' #####
    # 
    # library(locpol)
    # weights=locLinWeightsC(x=X, xeval=X, bw=h, kernel=TriweigK)$locWeig
    # # weights[i,j]= w_i(X_j,h)
    
  } else {
    print('Choose Nadaraya-Watson (NW) or Local linear (LL) weights')
  }
  
  ########## Adjusted random variables Y1 and Y2 ############
  
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
    Y1_adj[i]=sum(weights_X[i,]*I_1[,i]);
    Y2_adj[i]=sum(weights_X[i,]*I_2[,i]);
  }  
  
  ####### Weights of x #######

    # scaled x with X
    x_scaled=matrix(0, nrow=m, ncol=n)
    for (k in 1:m){
      for (i in 1:n){
        x_scaled[k,i]=(x[k]-X[i])/h
      }
    }
    
    # triweight kernel function
    kernel_x=matrix(0,nrow=m,ncol=n)
    for (k in 1:m){
      for (i in 1:n){
        if (abs(x_scaled[k,i]) <= 1){
          kernel_x[k,i]=(35/32)*((1-((x_scaled[k,i])^2))^3)
        }
        else{
          kernel_x[k,i]=0
        }
      }
    }
    
    
    # calculation of the weights depending on the chosen method (NW or LL)
    # NW weights
    if(weights=='NW'){
      
      weights_x=matrix(nrow=m,ncol=n)
      for (k in 1:m){
        for( i in 1:n){
          weights_x[k,i]=kernel_x[k,i]/sum(kernel_x[k,])
        }
      }
      # LL weights  
    } else if(weights=='LL'){
      
      S0_x=c(1:m)  
      S1_x=c(1:m)
      S2_x=c(1:m)
      
      for(k in 1:m){
        S0_x[k]=sum(kernel_x[k,])/(n*h)
        S1_x[k]=sum(x_scaled[k,]*kernel_x[k,])/(n*h)
        S2_x[k]=sum(x_scaled[k,]^2*kernel_x[k,])/(n*h)
      }
      
      weights_x=matrix(nrow=m, ncol=n)
      for(k in 1:m){
        for(i in 1:n){
          weights_x[k,i]=( kernel_x[k,i]*(S2_x[k]-x_scaled[k,i]*S1_x[k]) )/( n*h*(S0_x[k]*S2_x[k]-S1_x[k]^2) )
        }
      }
      
      # ##### Via library 'locpol' #####
      # 
      # library(locpol)
      # weights=locLinWeightsC(x=X, xeval=x, bw=h, kernel=TriweigK)$locWeig

      
    } else {
      print('Choose Nadaraya-Watson (NW) or Local linear (LL) weights')
    }
  
  
  ####### Calculation of conditional Kendall's tau for x #######
  
  conc=matrix(0,nrow=n,ncol=n) 
  for (i in 1:n){
    for (j in 1:n){
      if((Y1_adj[i]<Y1_adj[j]) & (Y2_adj[i]<Y2_adj[j])){
        conc[i,j]=1
      }
      else{
        conc[i,j]=0
      }
      
    }
  }
  
  # for X_k:  
  SW=c(1:n);
  wconc=matrix(0,nrow=n, ncol=n);
  CondKendallsTau_X=c(1:m);
  
  for (k in 1:m){
    for(i in 1:n){
      SW[i]=weights_x[k,i]^2;
      for(j in 1:n){
        if(conc[i,j]==1){
          wconc[i,j]=weights_x[k,i]*weights_x[k,j]
        }
        else{
          wconc[i,j]=0
        }
      }
    }
    
    Swconc=sum(wconc);
    SSW=sum(SW);
    
    CondKendallsTau_X[k]=(4*Swconc/(1-SSW))-1;
      print(k)
  }
  
  ######## Output ########
  
  return(CondKendallsTau_X)
  
}


