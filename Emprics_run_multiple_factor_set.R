## ----setup, include=FALSE----------------------------------------------------------
if(length(.libPaths()) == 1){
  # We're in Rscript.exe
  possible_lib_paths <- file.path(Sys.getenv(c('USERPROFILE','R_USER')),
                                  "R","win-library",
                                  paste(R.version$major,
                                        substr(R.version$minor,1,1),
                                        sep='.'))
  indx <- which(file.exists(possible_lib_paths))
  if(length(indx)){
    .libPaths(possible_lib_paths[indx[1]])
  }
  # CLEAN UP
  rm(indx,possible_lib_paths)
}
library(Rfast)
library(data.table)
library(corpcor)
library(stringr)
source(file = '../functions.R')

ptm <- proc.time()


#######################################
# Preparations ########################
#######################################

training_ratio = 1/2
artificial_ratio = 0
isDebug=F
FWER_threshold = 0.05
# T_obs_train_ = as.integer(training_ratio*t_orig)

set.seed(123)
data_source_folder=file.path('/home/jasonzou/data/July2021_Factors_ForUnconditionalSelection')
setwd(data_source_folder)

DS_returns=read.csv('/home/jasonzou/PanelPoSI_Sept2022/data/HouEtAl_processed/DS_resid_df.csv',row.names = 1)
t_orig=nrow(DS_returns)

# Creating PC of HL factors
K=20
pz3 <- prcomp(DS_returns[1:(0.5*t_orig),])
temp=as.matrix(DS_returns)%*%(as.matrix(pz3$rotation[,1:K]))
DS_PC_factors = as.data.frame(sweep(temp, 2, matrix(pz3$sdev[1:K],ncol=1), FUN = '/'))

HL_factors = read.csv('/home/jasonzou/PanelPoSI_Sept2022/data/HouEtAl_processed/HML_no_mkt.csv',row.names = 1)
HL_PC_factors2 = read.csv('/home/jasonzou/PanelPoSI_Sept2022/data/HouEtAl_processed/HML_PC.csv',row.names = 1)
pz3 <- prcomp(HL_factors[1:(0.5*t_orig),])
temp=as.matrix(HL_factors)%*%(as.matrix(pz3$rotation[,1:K]))
HL_PC_factors = as.data.frame(sweep(temp, 2, matrix(pz3$sdev[1:K],ncol=1), FUN = '/'))

adj_r2 = function(y_actual, y_pred, k){
  n = length(y_actual)
  rss <- sum((y_pred - y_actual) ^ 2)  ## residual sum of squares
  tss <- sum((y_pred - mean(y_actual)) ^ 2)  ## total sum of squares
  # ADJ_R2 = 1-rss/tss*(n-1)/(n-k-1)
  ADJ_R2 = 1-rss/tss
  return(ADJ_R2)
}
# selected_factors_index=P_POSI
eval_performance=function(selected_factors_index,all_factors,ff_portfolios){
  n_selection = length(selected_factors_index)
  
  covariates= as.matrix(all_factors[train_index,])
  selected_covariates = covariates[,selected_factors_index]
  response_oos = as.matrix(ff_portfolios[-train_index,])
  response_in_sample = as.matrix(ff_portfolios[train_index,])
  if(n_selection==0){
    predicted = response_oos*0.0
    alpha_vec_ins=colMeans(response_in_sample)
    alpha_vec_oos=colMeans(response_oos)
    
    
    # alpha_vec_ins=sqrt(colMeans(response_in_sample^2))
    # alpha_vec_oos=sqrt(colMeans(response_oos^2))
  }else{
    new_covaraites=as.matrix(all_factors[-train_index,selected_factors_index])
    fitted_beta=matrix(n_selection,ncol(response_in_sample),data=0)
    alpha_vec_ins=rep(0,ncol(response_in_sample))
    alpha_vec_oos=alpha_vec_ins
    for (i_col in 1:ncol(response_in_sample)){
      temp_df=cbind(response_in_sample[,i_col],selected_covariates)
      colnames(temp_df)=c('Y',colnames(selected_covariates))
      curr_ols=lm(Y~.-1,data=as.data.frame(temp_df))
      
      fitted_beta[,i_col]=coef(curr_ols)
      
      alpha_vec_ins[i_col]=mean(residuals(curr_ols))
      alpha_vec_oos[i_col]=mean(predict(curr_ols,as.data.frame(new_covaraites)))
      
      # alpha_vec_ins[i_col]=sqrt(mean((residuals(curr_ols))^2))
      # alpha_vec_oos[i_col]=sqrt(mean((predict(curr_ols,as.data.frame(new_covaraites)))^2))
      
    }
    # fitted_beta=pinv(t(selected_covariates)%*%selected_covariates)%*%t(selected_covariates)%*%response_in_sample
    # alpha_vec_ins=colMeans(response_in_sample-selected_covariates%*%fitted_beta)
    
    predicted = new_covaraites%*%fitted_beta
    
  }
  
  predicted_error = response_oos-predicted
  # alpha_vec_oos = colMeans(predicted_error)
  RMS_alpha=sqrt(sum(alpha_vec_oos^2)/ncol(ff_portfolios))
  
  
  selected_covariates = cbind(rep(1,nrow(covariates)),covariates[,selected_factors_index])
  if(n_selection==0){
    predicted = response_oos*0.0
  }else{
    fitted_beta=pinv(t(selected_covariates)%*%selected_covariates)%*%t(selected_covariates)%*%response_in_sample
    predicted = cbind(rep(1,nrow(new_covaraites)),new_covaraites)%*%fitted_beta
  }
  
  adj_r2_vec=alpha_vec_oos
  
  for (i_unit in 1:ncol(response_oos)){
    adj_r2_vec[i_unit]=adj_r2(response_oos[,i_unit],predicted[,i_unit] ,n_selection)
    
  }
  
  return(list(abs_alpha_INS=abs(alpha_vec_ins),abs_alpha_OOS=abs(alpha_vec_oos),
              adj_r2s=adj_r2_vec,n_selection=n_selection))
}


run_inference=function(ff_portfolios,all_factors,factor_choice){
      
      
  #######################################
  # Enumerating thru the test portfolios #
  #######################################
  
  all_portfolio_names=colnames(ff_portfolios)
  
  
  
  T_obs_train_ = as.integer(training_ratio*t_orig)
  covariates= as.matrix(all_factors[train_index,])
  response = as.matrix(ff_portfolios[train_index,])
  N = ncol(response)
  d = ncol(covariates)
  denom_vec=diag(t(covariates)%*%covariates)
  numerator_ols=t(covariates)%*%response
  if (factor_choice == 'HL wFF3'){
    penalty_omega_inv = rep(1,d)
    # penalty_omega_inv[(d-4):(d-2)]=0
    penalty_omega_inv[(d-3):(d-2)]=0
    penalty_omega_inv=penalty_omega_inv/sum(penalty_omega_inv)*d
  } else if (factor_choice == 'HL wFF5'){
    penalty_omega_inv = rep(1,d)
    penalty_omega_inv[(d-3):(d)]=0
    penalty_omega_inv=penalty_omega_inv/sum(penalty_omega_inv)*d
  } else {
    penalty_omega_inv = rep(1,d)
  }
  
  beta_matrix_ols=numerator_ols/denom_vec
  
  sigmahat_matrix = matrix(data=0,nrow=ncol(covariates),ncol=ncol(response))
  for(i_unit in 1:ncol(beta_matrix_ols)){
    repeated_response_ = t(matrix(rep(response[,i_unit],ncol(covariates)),ncol=T_obs_train_,byrow = T))
    temp=t(t(covariates)*(beta_matrix_ols[,i_unit]))
    
    sigmahat_matrix[,i_unit]=sqrt(colSums((repeated_response_-temp)^2)/(T_obs_train_-1))
  }
  
  
  OLS_t_stat=beta_matrix_ols/sigmahat_matrix* sqrt(denom_vec)
  
  OLS_log_t_stat_vec_ = pt(q=-rowMaxs(abs(OLS_t_stat),value = T),df=T_obs_train_-1,log.p=T)+log(2)
  
  this_PoSI_output = posi_lognorm_pval(training_ratio=training_ratio, 
                                       X=as.matrix(all_factors[train_index,]), 
                                       Y=as.matrix(ff_portfolios[train_index,]),
                                       penalty_omega_inv = penalty_omega_inv)
  
  FWER_thresholds =  c(0.01,0.02,0.05,0.1)
  adjr2_dt_list=list()
  xsec_alpha_INS_dt_list=list()
  xsec_alpha_OOS_dt_list=list()
  N_selection_list=list()
  for (FWER_threshold in FWER_thresholds){
    N_OLS = which(OLS_log_t_stat_vec_<=log(FWER_threshold))
    B_OLS = which(OLS_log_t_stat_vec_+log(d)+log(N)<=log(FWER_threshold))
    
    
    independent_p_val=this_PoSI_output$posi_log_pval_matrix
    naive_t_p_val=this_PoSI_output$naive_t_p_val_matrix
    
    
    LASSO_min_p = rowMins(as.matrix(naive_t_p_val),value = T)
    N_LASSO=which(LASSO_min_p<=FWER_threshold)
    B_LASSO=which(LASSO_min_p*d*N<=FWER_threshold)
    
    
    PoSI_min_log_p = rowMins(as.matrix(independent_p_val),value = T)
    B_POSI = which(PoSI_min_log_p+log(d)+log(N)<=log(FWER_threshold))
    
    P_POSI_df=panel_unordered_reject(posi_log_pval_matrix = independent_p_val,significance_level = FWER_threshold)
    P_POSI_rejected_names = rownames(P_POSI_df)[P_POSI_df$Reject=='Reject']
    P_POSI = match(P_POSI_rejected_names,colnames(covariates))
    
    print(paste('FWER ',FWER_threshold,'; factor set',factor_choice))
    print(P_POSI_df[1:(length(P_POSI_rejected_names)+5),])
    
    
    
    xsec_alpha_INS = data.frame(matrix(nrow = 6,ncol=N))
    xsec_alpha_OOS = data.frame(matrix(nrow = 6,ncol=N))
    N_selection  = data.frame(matrix(nrow = 6,ncol=1))
    fullxsectional_adjr2_result_df = data.frame(matrix(nrow = 6,ncol=N))
    
    print('N_OLS')
    eval_result = eval_performance(N_OLS,all_factors ,ff_portfolios )
    xsec_alpha_INS[1,]=eval_result$abs_alpha_INS
    xsec_alpha_OOS[1,]=eval_result$abs_alpha_OOS
    N_selection[1,]=eval_result$n_selection
    fullxsectional_adjr2_result_df[1,]=eval_result$adj_r2s
    
    print('B_OLS')
    eval_result = eval_performance(B_OLS,all_factors ,ff_portfolios )
    xsec_alpha_INS[2,]=eval_result$abs_alpha_INS
    xsec_alpha_OOS[2,]=eval_result$abs_alpha_OOS
    N_selection[2,]=eval_result$n_selection
    fullxsectional_adjr2_result_df[2,]=eval_result$adj_r2s

    eval_result = eval_performance(N_LASSO,all_factors ,ff_portfolios )
    print('N_LASSO')
    xsec_alpha_INS[3,]=eval_result$abs_alpha_INS
    xsec_alpha_OOS[3,]=eval_result$abs_alpha_OOS
    N_selection[3,]=eval_result$n_selection
    fullxsectional_adjr2_result_df[3,]=eval_result$adj_r2s

    eval_result = eval_performance(B_LASSO,all_factors ,ff_portfolios )
    print('B_LASSO')
    xsec_alpha_INS[4,]=eval_result$abs_alpha_INS
    xsec_alpha_OOS[4,]=eval_result$abs_alpha_OOS
    N_selection[4,]=eval_result$n_selection
    fullxsectional_adjr2_result_df[4,]=eval_result$adj_r2s

    eval_result = eval_performance(B_POSI,all_factors ,ff_portfolios )
    print('B_POSI')
    xsec_alpha_INS[5,]=eval_result$abs_alpha_INS
    xsec_alpha_OOS[5,]=eval_result$abs_alpha_OOS
    N_selection[5,]=eval_result$n_selection
    fullxsectional_adjr2_result_df[5,]=eval_result$adj_r2s

    eval_result = eval_performance(P_POSI,all_factors ,ff_portfolios )
    print('P_POSI')
    xsec_alpha_INS[6,]=eval_result$abs_alpha_INS
    xsec_alpha_OOS[6,]=eval_result$abs_alpha_OOS
    N_selection[6,]=eval_result$n_selection
    fullxsectional_adjr2_result_df[6,]=eval_result$adj_r2s
    
    colnames(xsec_alpha_INS) = colnames(response)
    colnames(xsec_alpha_OOS) = colnames(response)
    xsec_alpha_INS$Method = c('N_OLS','B_OLS','N_LASSO','B_LASSO','B_POSI','P_POSI')
    xsec_alpha_OOS$Method = c('N_OLS','B_OLS','N_LASSO','B_LASSO','B_POSI','P_POSI')
    N_selection$Method = c('N_OLS','B_OLS','N_LASSO','B_LASSO','B_POSI','P_POSI')
    fullxsectional_adjr2_result_df$Method = c('N_OLS','B_OLS','N_LASSO','B_LASSO','B_POSI','P_POSI')
    if(factor_choice%in%c('HL PC','DS PC')){
      O_POSI_df = panel_nested_reject(independent_p_val,FWER_threshold)  
      
      O_POSI_rejected_names = rownames(O_POSI_df)[O_POSI_df$Reject=='Reject']
      O_POSI = match(O_POSI_rejected_names,colnames(covariates))
      eval_result = eval_performance(P_POSI,all_factors ,ff_portfolios )
      print('O_POSI')
      xsec_alpha_INS[7,]=c(eval_result$abs_alpha_INS,'O_POSI')
      xsec_alpha_OOS[7,]=c(eval_result$abs_alpha_OOS,'O_POSI')
      N_selection[7,]=c(eval_result$n_selection,'O_POSI')
      fullxsectional_adjr2_result_df[7,]=c(eval_result$adj_r2s,'O_POSI')
      
      # result_df[7,]=c(eval_result$abs_alpha,'O_POSI')
      # fullxsectional_adjr2_result_df[7,]=eval_result$adj_r2s
      # fullxsectional_adjr2_result_df$Method = c('N_OLS','B_OLS','N_LASSO','B_LASSO','B_POSI','P_POSI','O_POSI')  
    }
    
    xsec_alpha_INS_dt = as.data.table(xsec_alpha_INS)
    xsec_alpha_OOS_dt = as.data.table(xsec_alpha_OOS)
    adjr2_dt = as.data.table(fullxsectional_adjr2_result_df)
    
    xsec_alpha_INS_dt$FWER_threshold=FWER_threshold
    xsec_alpha_OOS_dt$FWER_threshold=FWER_threshold
    adjr2_dt$FWER_threshold=FWER_threshold
    
    xsec_alpha_INS_dt_list[[as.character(FWER_threshold)]]=xsec_alpha_INS_dt
    xsec_alpha_OOS_dt_list[[as.character(FWER_threshold)]]=xsec_alpha_OOS_dt
    adjr2_dt_list[[as.character(FWER_threshold)]]=adjr2_dt
    
    N_selection$FWER_threshold=FWER_threshold
    N_selection_list[[as.character(FWER_threshold)]]=N_selection
  }
  xsec_alpha_INS_joint = rbindlist(xsec_alpha_INS_dt_list)
  xsec_alpha_OOS_joint = rbindlist(xsec_alpha_OOS_dt_list)
  adjr2_dt_join = rbindlist(adjr2_dt_list)
  N_selection_joint = rbindlist(N_selection_list)
  
  return(list(alpha_INS=xsec_alpha_INS_joint,alpha_OOS=xsec_alpha_OOS_joint,adjr2=adjr2_dt_join,N_selections=N_selection_joint))
}

alpha_INS_list= list()
alpha_OOS_list= list()
xsectional_adjr2_list=list()
N_selections_list=list()
mainDir = '/home/jasonzou/PoSI_Nov_2021/q1q2'
factor_choices=c('HL','HL PC','HL + HL PC','DS PC','HL wFF3','HL wFF5')
t_orig = nrow(DS_returns)
train_index = 1:as.integer(training_ratio*t_orig)

for (i_factor_set in 1:length(factor_choices)) {
  factor_choice = factor_choices[i_factor_set]
  
  if(factor_choice %in% c('HL','HL wFF3','HL wFF5')){
    all_factors=HL_factors
  }else if(factor_choice=='HL PC'){
    all_factors=HL_PC_factors
  } else if(factor_choice=='HL + HL PC'){
    all_factors=cbind(HL_factors,HL_PC_factors)
  } else if(factor_choice=='DS PC'){
    all_factors=DS_PC_factors
  }
  
  output=run_inference(DS_returns,all_factors,factor_choice)
  temp=output$alpha_INS
  temp$factor_choice  =factor_choice
  alpha_INS_list[[factor_choice]]=temp
  
  temp=output$alpha_OOS
  temp$factor_choice  =factor_choice
  alpha_OOS_list[[factor_choice]]=temp
  
  temp=output$adjr2
  temp$factor_choice  =factor_choice
  xsectional_adjr2_list[[factor_choice]]=temp
  
  temp=output$N_selections
  temp$factor_choice  =factor_choice
  N_selections_list[[factor_choice]]=temp

}

xsectional_adjr2_bound=rbindlist(xsectional_adjr2_list)
colnames(xsectional_adjr2_bound)=c(colnames(DS_returns),c('Method','FWER_threshold','factor_choice'))
xsectional_adjr2_bound$in.out='out'
alpha_INS_bound = rbindlist(alpha_INS_list)
alpha_INS_bound$in.out='in'
alpha_OOS_bound = rbindlist(alpha_OOS_list)
alpha_OOS_bound$in.out='out'
alpha_bounded_oneFWER=rbind(alpha_INS_bound,alpha_OOS_bound)
proc.time() - ptm

# melted_alpha_abs = melt(alpha_bounded_oneFWER,id.vars = c('Method','FWER_threshold','factor_choice','in.out'))
# fwrite(melted_alpha_abs,'~/PanelPoSI_Sept2022/csv_outputs/melted_alphas.csv')
# fwrite(melted_alpha_abs,'~/PanelPoSI_Sept2022/csv_outputs/melted_rmse.csv')

bounded_selection=rbindlist(N_selections_list)
bounded_selection[FWER_threshold==.05,]
fwrite(bounded_selection,'/home/jasonzou/PanelPoSI_Sept2022/csv_outputs/melted_selection_PCs.csv')
melted_xsectional_adjr2 = melt(xsectional_adjr2_bound,id.vars = c('Method','FWER_threshold','factor_choice','in.out'))

fwrite(melted_xsectional_adjr2,'/home/jasonzou/PanelPoSI_Sept2022/csv_outputs/melted_adjr2.csv')
alpha_INS_bound[FWER_threshold==0.05,]
alpha_OOS_bound[FWER_threshold==0.05,]
proc.time() - ptm

result_=dcast(result_bound, Method +n_selection+`OOS AdjR2`+ `OOS RMSE`~factor_choice)
fwrite(result_,paste('~/PanelPoSI_Sept2022/csv_outputs/empirics_selection_multi_factor_choices.csv'))
proc.time() - ptm
