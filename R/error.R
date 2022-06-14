error <- function(knobi_results,env=NULL) {

  X=knobi_results$data$Average_Biomass
  Y=knobi_results$data$SP
  total_n=length(X)

  r_b<-knobi_results$fit$Parameter_estimates[1]
  K_b<-knobi_results$fit$Parameter_estimates[2]
  if(knobi_results$control$pella){
    p_b<-knobi_results$fit$Parameter_estimates[3]
  } else {p_b=1}

  n_params=length(knobi_results$fit$Parameter_estimates)
  df_b=total_n-n_params

  rss_b = sum((Y-((r_b/p_b)*X*(1-(X/K_b)^(p_b))))^2)
  tss_b = sum((Y-mean(Y))^2)

  ser_b=sqrt(rss_b/(df_b))
  r2_b = 1-rss_b/tss_b
  r2_adj_b = 1-(rss_b/(df_b))/(tss_b/(total_n-1))
  AIC_b = total_n*log(rss_b/total_n)+2*n_params
  RMSE_b<-sqrt(rss_b/total_n)
  MAPE_b<-1/total_n*(sum(abs(((Y-((r_b/p_b)*X*(1-(X/K_b)^(p_b)))))/Y)))

  if(is.null(env)){

    error_table=data.frame(ser_b,r2_b,r2_adj_b,AIC_b, RMSE_b, MAPE_b)
    names=names(error_table)=c("SER","R-squared", "adj-R-squared", "AIC", "RMSE", "MAPE")
    rownames(error_table)="KBPM_error"

  } else {

    base_model=c(ser_b,r2_b,r2_adj_b,AIC_b, RMSE_b,
                 MAPE_b, NA, NA)

    model_env_Additive=env$model_env_Additive
    model_env_Multiplicative=env$model_env_Multiplicative

    if(is.null(env$scaled_environmental_var)){
      env=as.matrix(env$environmental_variables)
    } else {
      env=as.matrix(env$scaled_environmental_var)
    }

    r_a<-model_env_Additive[1]
    K_a<-model_env_Additive[2]
    if(knobi_results$control$pella){
      c_a<-model_env_Additive[3:(length(model_env_Additive)-1)]
      p_a<-model_env_Additive[length(model_env_Additive)]
    } else {
      c_a<-model_env_Additive[3:(length(model_env_Additive))]
      p_a=1}

    n_params_env=length(model_env_Additive)
    df_env=total_n-n_params_env

    rss_a = sum((Y-as.numeric((r_a/p_a)*X*(1-(X/K_a)^(p_a))+as.numeric(env %*% c_a)*X))^2)
    tss_a = sum((Y-mean(Y))^2)

    ser_a=sqrt(rss_a/(df_env))
    r2_a = 1-rss_a/tss_a
    r2_adj_a = 1-(rss_a/(df_env))/(tss_a/(total_n-1))
    AIC_a = total_n*log(rss_a/total_n)+2*n_params_env
    RMSE_a<-sqrt(rss_a/total_n)
    MAPE_a<-1/total_n*(sum(abs((Y-as.numeric((r_a/p_a)*X*(1-(X/K_a)^(p_a))+as.numeric(env %*% c_a)+X))/Y)))
    F_test_a=((rss_b-rss_a)/(df_b-df_env))/(rss_b/(df_b))
    pF_a=stats::pf(F_test_a,df_b-df_env,df_env,lower.tail=F)

    env_a_model=c(ser_a,r2_a,r2_adj_a,AIC_a, RMSE_a,
                  MAPE_a,F_test_a,pF_a)


    r_m<-model_env_Multiplicative[1]
    K_m<-model_env_Multiplicative[2]
    if(knobi_results$control$pella){
      c_m<-model_env_Multiplicative[3:(length(model_env_Additive)-1)]
      p_m<-model_env_Multiplicative[length(model_env_Additive)]
    } else {
      c_m<-model_env_Multiplicative[3:(length(model_env_Additive))]
      p_m=1}

    rss_m = sum((Y-as.numeric(exp(1)^{as.numeric(env %*% c_m)}*((r_m/p_m)*X*(1-(X/K_m)^(p_m)))))^2)
    tss_m = sum((Y-mean(Y))^2)

    ser_m=sqrt(rss_m/(df_env))
    r2_m = 1-rss_m/tss_m
    r2_adj_m = 1-(rss_m/(df_env))/(tss_m/(total_n-1))
    AIC_m = total_n*log(rss_m/total_n)+2*n_params_env
    RMSE_m<-sqrt(rss_m/total_n)
    MAPE_m<-1/total_n*(sum(abs((Y-as.numeric(exp(1)^{as.numeric(env %*% c_m)}*((r_m/p_m)*X*(1-(X/K_m)^(p_m)))))/Y)))
    F_test_m=((rss_b-rss_m)/(df_b-df_env))/(rss_b/(df_b))
    pF_m=stats::pf(F_test_m,df_b-df_env,df_env,lower.tail=F)

    env_m_model=c(ser_m,r2_m,r2_adj_m,AIC_m, RMSE_m,
                  MAPE_m,F_test_m,pF_m)

    error_table=rbind(base_model, env_a_model, env_m_model)
    rownames(error_table)=c("base model", "additive model", "multiplicative model")
    colnames(error_table)=c("SER","R-squared", "adj-R-squared", "AIC", "RMSE",
                            "MAPE"," F-value", "Pr(>F)")
  }


  return(error_table)

}

