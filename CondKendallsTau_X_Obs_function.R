########################################################
## This function calculates the conditional Kendall's ##
############# tau in the observations of X. ############
### If the interest is in a conditional Kendall's tau ##
########### that is in SOME value of X, use ############  
############# CondKendallsTau_X_function.R #############
########################################################


#######################################################
######### Conditional Kendall's tau function ##########
#######################################################

CondKendallsTau_X_Obs=function(Y1,Y2,X,h, weights){
  
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
    
    # Local Linear weights ( w_i(X_j,h)=weights[i,j] )
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
    Y1_adj[i]=sum(weights[i,]*I_1[,i]);
    Y2_adj[i]=sum(weights[i,]*I_2[,i]);
  }  
  
  ####### Calculation of conditional Kendall's tau for all values of X #######

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
  CondKendallsTau_X_Obs=c(1:n);
  
  for (k in 1:n){
    for(i in 1:n){
      SW[i]=weights[k,i]^2;
      for(j in 1:n){
        if(conc[i,j]==1){
          wconc[i,j]=weights[k,i]*weights[k,j]
        }
        else{
          wconc[i,j]=0
        }
      }
    }
    
    Swconc=sum(wconc);
    SSW=sum(SW);
    
    CondKendallsTau_X_Obs[k]=(4*Swconc/(1-SSW))-1;
    print(k)
  }
  
  ######## Output ########
  
  out_list = list("CondKendallsTau_X_Obs"=CondKendallsTau_X_Obs, "Y1.adj" = Y1_adj, "Y2.adj" = Y2_adj)
  return(out_list)

}

