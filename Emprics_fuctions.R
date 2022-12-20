## ----setup, include=FALSE----------------------------------------------------------
library(selectiveInference)
library(glmnet)
library(lmtest)
library(nlme)
library(zeallot)
library(knockoff)
# library(Rfast)
library(data.table)
library(corpcor)
library(stringr)
library(TruncatedNormal)
library(pracma)
# library(RPtests)
library(nleqslv)

posi_lognorm_pval = function(training_ratio, X, Y, isDebug=F, penalty_omega_inv=NA){
  
  
  #######################################
  # Enumerating thru the test portfolios #
  #######################################
  
  t_orig = nrow(Y)
  t_obs = t_orig
  N = ncol(Y)
  J = ncol(X)
  
  all_portfolio_names=1:ncol(Y)
  
  if(isDebug){
    all_portfolio_names=all_portfolio_names[1:2]
    Y=Y[,1:2]
    
    print('warining!!! Debug Mode')
  }
  
  
 
  
  posi_log_pval_matrix=as.data.frame(matrix(data=Inf,ncol = N,nrow = J))
  print(dim(posi_log_pval_matrix))
  colnames(posi_log_pval_matrix)=all_portfolio_names
  rownames(posi_log_pval_matrix)=colnames(X)
  
  
  naive_t_p_val=posi_log_pval_matrix
  
  
  beta_hat_table=as.data.frame(matrix(data=NA,ncol = length(all_portfolio_names),nrow = length(colnames(X))))
  colnames(beta_hat_table)=all_portfolio_names
  rownames(beta_hat_table)=colnames(X)
  
  residual_matrix=matrix(ncol=ncol(Y),nrow=t_orig)
  for(i_col in 1:length(all_portfolio_names)){
    col = all_portfolio_names[i_col]
    
    
    y = Y[,col]
    k_fold = 10
    # if(is.na(penalty_omega_inv)){
    #   penalty_omega_inv=rep(1,ncol(X))
    # }

    
    # cvob0=cv.glmnet(X, y,nfolds = k_fold,lambda=exp((32:-32)/4)*log(ncol(X))/sqrt(t_obs),penalty.factor = penalty_omega_inv,intercept = T)
    cvob0=cv.glmnet(X, y,nfolds = k_fold,penalty.factor = penalty_omega_inv,intercept = F)
    
    optimal_lambda_index=which(cvob0$lambda==cvob0$lambda.1se)
    optimal_lambda=cvob0$lambda[optimal_lambda_index]
    optimal_lasso=glmnet(X,y,lambda = optimal_lambda,penalty.factor = penalty_omega_inv,intercept = F)
    
    
    optimal_beta=as.vector(optimal_lasso$beta)
    y_hat=predict(optimal_lasso,X)
    residuals_lasso=y-y_hat
    
    
    M_set=optimal_beta!=0
    M_cardi=sum(M_set)
    if(M_cardi==0){
      p_raw_vec=rep(Inf,M_cardi)
      next
    }
    # sigma_squared=var(residuals_lasso)[1]
    sigma_squared=sum(residuals_lasso^2)/(t_obs-M_cardi)
    M_ind=which(M_set)
    beta_hat_LASSO=optimal_beta[M_set]
    X_M=X[,M_ind]
    X_M=as.matrix(X_M)
    s_M=sign(beta_hat_LASSO)
    
    X_min_M=X[,which(optimal_beta==0)]
    X_M_Pinv=pinv(X_M)   # the same as pinv(t(X_M)%*%X_M)%*%t(X_M)
    
    P_M = X_M %*%X_M_Pinv
    
    I_min_P_M=diag(dim(P_M)[1])-P_M
    
    y_shifted=y
    # y_shifted=y
    
    # beta_hat_one_step=X_M_Pinv%*%y
    beta_hat_one_step=X_M_Pinv%*%y_shifted
    beta_hat_table[M_ind,i_col]=beta_hat_one_step
    Fisher_Inv=pinv(t(X_M)%*%X_M)
    covariance_matrix=sigma_squared*Fisher_Inv
    
    dig_var=diag(covariance_matrix)
    dig_sd=sqrt(dig_var)
    
    
    raw_p_t=2*pt(-abs(beta_hat_LASSO)/dig_sd,nrow(X)-length(M_ind))
    
    naive_t_p_val[M_ind,i_col]=raw_p_t
    
    
    # Using Lee, Sun, Sun and Tyalor
    omega_inv_M=penalty_omega_inv[M_ind]
    omega_inv_minus_M=penalty_omega_inv[-M_ind]
    
    component1=t(X_min_M)%*%I_min_P_M
    component2=t(X_min_M)%*%t(X_M_Pinv)%*%(s_M*omega_inv_M)
    
    
    
    
    
    if(M_cardi==1){
      A_mat=rbind(1/optimal_lambda*component1,
                  -1/optimal_lambda*component1,
                  -s_M%*%X_M_Pinv)
      
      
      b_vec=rbind(omega_inv_minus_M-component2,
                  omega_inv_minus_M+component2,
                  -optimal_lambda*s_M%*%Fisher_Inv%*%(s_M*omega_inv_M))
      
    }else{
      A_mat=rbind(1/optimal_lambda*component1,
                  -1/optimal_lambda*component1,
                  -diag(s_M)%*%X_M_Pinv)
      
      
      b_vec=rbind(omega_inv_minus_M-component2,
                  omega_inv_minus_M+component2,
                  -optimal_lambda*diag(s_M)%*%Fisher_Inv%*%(s_M*omega_inv_M))
    }
    
    
    positive_gap = b_vec-A_mat%*%y_shifted
    b_vec[positive_gap<0]=b_vec[positive_gap<0]-min(positive_gap[positive_gap<0])
    p_raw_vec=rep(Inf,M_cardi)
    for(i_of_M in 1:M_cardi){
      this_eta =X_M_Pinv[i_of_M,]
      this_xi = this_eta/sum(this_eta^2)
      xi_eta_outer=this_xi%*%t(this_eta)
      
      this_z=(diag(length(this_eta))-xi_eta_outer)%*%y_shifted
      numerators=b_vec-A_mat%*%this_z
      num_keep_ind=which(numerators>0)
      
      numerators=numerators[num_keep_ind]
      
      denominators=A_mat%*%this_xi
      denominators=denominators[num_keep_ind]
      ratio = numerators/denominators
      
      
      # num_keep_ind2=which(abs(ratio-beta_hat_one_step[i_of_M])>1e-32)
      
      
      # numerators=numerators[num_keep_ind2]
      # denominators=denominators[num_keep_ind2]
      # ratio = numerators/denominators
      negative_ones = which(denominators<0)
      positive_ones = which(denominators>0)
      # print(paste('Ratio chosen from ',length(ratio)))
      V_minus=max(ratio[negative_ones & (beta_hat_one_step[i_of_M]-ratio>1e-32)],na.rm = T)
      
      V_plus=min(ratio[positive_ones& (ratio-beta_hat_one_step[i_of_M]>1e-32)],na.rm = T)
      
      
      # print(paste(dig_sd[i_of_M],V_minus,V_plus))
      
      if(V_minus/dig_sd[i_of_M]>=V_plus/dig_sd[i_of_M]){
        next
        print('bizzare')
        
      }else{
        This_Df=max(1,nrow(X_M)-M_cardi)
        if(beta_hat_one_step[i_of_M]>0){
          right_tail=TruncatedNormal::ptmvnorm(q=-beta_hat_one_step[i_of_M],mu=0,sigma = dig_sd[i_of_M],
                                            lb =-V_plus,ub = -V_minus,type = 'qmc',
                                            log = T,B=20000)
          left_tail=TruncatedNormal::ptmvnorm(q=-beta_hat_one_step[i_of_M],mu=0,sigma = dig_sd[i_of_M],
                                           lb =V_minus,ub = V_plus,type = 'qmc',
                                           log = T,B=20000)
          
        }else{
          right_tail=TruncatedNormal::ptmvnorm(q=beta_hat_one_step[i_of_M],mu=0,sigma = dig_sd[i_of_M],
                                            lb =V_minus,ub = V_plus,type = 'qmc',
                                            log = T,B=20000)
          left_tail=TruncatedNormal::ptmvnorm(q=beta_hat_one_step[i_of_M],mu=0,sigma = dig_sd[i_of_M],
                                           lb =-V_plus,ub = -V_minus,type = 'qmc',
                                           log = T,B=20000)
        }
      }
      if(is.na(right_tail) | is.na(left_tail)){
        p_raw=Inf
        next
      } 
      
      if(right_tail==-Inf & left_tail==-Inf){
        p_raw=Inf
      } else if(abs(right_tail-left_tail)>10){
        p_raw=max(right_tail,left_tail)
      }else{
        p_raw=log(exp(right_tail)+exp(left_tail))
      }
      
      p_raw_vec[i_of_M]=p_raw
    }
    posi_log_pval_matrix[M_ind,i_col]=p_raw_vec
  }
  
  return(list(posi_log_pval_matrix=posi_log_pval_matrix,naive_t_p_val_matrix=naive_t_p_val,ose=beta_hat_table))
}

##########################################################################################################################
panel_unordered_reject = function(posi_log_pval_matrix,significance_level){
  
  
  
  pval_matrix=posi_log_pval_matrix
  pval_matrix[is.na(pval_matrix)]=Inf
  
  ############################################################################################
  
  time_series_Bonf_rejection_table=as.data.frame(matrix(data=F,ncol = 4,nrow = nrow(posi_log_pval_matrix)))
  row.names(time_series_Bonf_rejection_table)=row.names(posi_log_pval_matrix)
  colnames(time_series_Bonf_rejection_table)=c('Reject','bonf_level','n','p_1')
  K_set = rowSums((pval_matrix!=Inf)*1.0)
  M_card = vector(length = nrow(posi_log_pval_matrix),mode = 'numeric')
  for(i_row in 1:nrow(time_series_Bonf_rejection_table)){
    this_pval_row=pval_matrix[i_row,]
    meaningful_ind=!(this_pval_row==Inf)
    meaningful_ind[is.na(meaningful_ind)]=F
    my_df=sum(pval_matrix[,meaningful_ind]!=Inf)
    M_card[i_row]=my_df
  }
  rho_inv=sum(K_set[M_card>0]/M_card[M_card>0])
  rho = 1/rho_inv
  for(i_row in 1:nrow(time_series_Bonf_rejection_table)){
    this_pval_row=pval_matrix[i_row,]
    meaningful_ind=!(this_pval_row==Inf)
    meaningful_ind[is.na(meaningful_ind)]=F
    if(sum(meaningful_ind)>0){
      my_df=M_card[i_row]
      p_1=min(this_pval_row[meaningful_ind])
      
      bonf_level=significance_level/(my_df)
      rejection_char=ifelse((p_1+log(my_df)+log(rho_inv))<log(significance_level) ,'Reject','NotReject')
      time_series_Bonf_rejection_table[i_row,]=c(rejection_char,bonf_level,my_df,p_1)
      
      
    }else{
      rejection_char='NotReject'
      time_series_Bonf_rejection_table[i_row,]=c(rejection_char,NA,my_df,NA)
      
      
    }
    
  }
  
  time_series_Bonf_rejection_table$p_1=exp(as.numeric(time_series_Bonf_rejection_table$p_1))
  time_series_Bonf_rejection_table$Reject=factor(time_series_Bonf_rejection_table$Reject, levels = c("Reject", "NotReject"))
  
  time_series_Bonf_rejection_table =time_series_Bonf_rejection_table[order(time_series_Bonf_rejection_table$Reject,time_series_Bonf_rejection_table$p_1),]
  time_series_Bonf_rejection_table$rho = rho
  return(time_series_Bonf_rejection_table)
}
panel_nested_reject = function(posi_log_pval_matrix,significance_level){
  # Calculate Y_check    
  Y_check_matrix=-posi_log_pval_matrix
  Y_check_matrix[is.na(Y_check_matrix)]=-Inf
  

  Simes_rejection_table=as.data.frame(matrix(data=F,ncol = 2,nrow = nrow(posi_log_pval_matrix)))
  row.names(Simes_rejection_table)=row.names(posi_log_pval_matrix)
  colnames(Simes_rejection_table)=c('q_check','N_check')
  
  # N_check_k iteration
  
  for(i_row in nrow(Simes_rejection_table):1){
    this_pval_row=Y_check_matrix[i_row,]
    meaningful_ind=!(this_pval_row==-Inf)
    meaningful_ind[is.na(meaningful_ind)]=F
    meaningful_ind=which(meaningful_ind) # script_K_i
    
    if(i_row==nrow(Simes_rejection_table)){
      
      script_K_check = meaningful_ind
    }else{
      
      script_K_check = union(meaningful_ind,script_K_check)
    }
    N_check = sum(Y_check_matrix[i_row:nrow(Y_check_matrix),script_K_check]!=-Inf)
    
    Simes_rejection_table[i_row,2]=N_check
    
  }

  # Z_check iteration
  Z_check_cumulative_switch=F
  for(i_row in nrow(Simes_rejection_table):1){
    this_pval_row=Y_check_matrix[i_row,]
    meaningful_ind=!(this_pval_row==-Inf)
    meaningful_ind[is.na(meaningful_ind)]=F
    meaningful_ind=which(meaningful_ind) # script_K_i
    
    if(sum(meaningful_ind)==0){
      Simes_rejection_table[i_row,1]=1
      this_row_contrib_to_Z_check=0
      next
    }
    
    if(!Z_check_cumulative_switch){
      if(i_row==nrow(Simes_rejection_table)){
        denominator=Simes_rejection_table[1,2]
      }else{
        denominator=Simes_rejection_table[1,2]-Simes_rejection_table[i_row+1,2]  
      }
      this_row_contrib_to_Z_check=sum(this_pval_row[meaningful_ind]/denominator)
      Z_check_cumulative=this_row_contrib_to_Z_check
      Z_check_cumulative_switch=T
      
    }else{
      denominator=Simes_rejection_table[1,2]-Simes_rejection_table[i_row+1,2]
      this_row_contrib_to_Z_check=sum(this_pval_row[meaningful_ind]/denominator)
      Z_check_cumulative=this_row_contrib_to_Z_check+Z_check_cumulative
      
      
    }
    q_check = exp(-Z_check_cumulative)
    
    
    Simes_rejection_table[i_row,1]=q_check
    
    
  }
  
  q_check_threshold=Simes_rejection_table$N_check/(nrow(Y_check_matrix)*ncol(Y_check_matrix))*significance_level
  # RHS=(1:nrow(Y_check_matrix))/nrow(Y_check_matrix)*significance_level
  Rejected_indices=which(Simes_rejection_table$q_check<=q_check_threshold)
  Simes_rejection_table$q_check_threshold=q_check_threshold
  Simes_rejection_table$Reject=factor(
    ifelse(Simes_rejection_table$q_check<=q_check_threshold,'Reject','NotReject'), levels = c("Reject", "NotReject"))
  return(Simes_rejection_table)
}
