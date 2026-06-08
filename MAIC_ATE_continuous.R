library(dplyr)
library(sandwich)
library(limSolve)
library(WeightIt)

set.seed(612)

sim_age<-function(n, p_young=0.7) {
  young<-rbinom(n, 1, prob=p_young)
  ifelse(young==1, runif(n, 45, 60), runif(n, 60, 75))
}

b0<-35.00
b_age<- -0.20
b_female<-2.00
b_A<-8.00
b_B<-5.00
mean_age_ctr<-40
b_ageA<-0.10
b_ageB<- -0.10
sigma<-5.00

n1<-300
n2<-500

sim_arm<-function(n, b_trt, b_age_trt, p_young, p_female) {
  age<-sim_age(n, p_young)
  female<-rbinom(n, 1, prob=p_female)
  trt_effect<-b_trt+b_age_trt*(age-mean_age_ctr)
  Y<-b0+b_age*age+b_female*female+trt_effect+rnorm(n, 0, sigma)
  data.frame(age=age, female=female, Y=Y)
}

t1_A<-sim_arm(n1/2, b_trt=b_A, b_age_trt=b_ageA, p_young=0.80, p_female=0.64)
t1_C<-sim_arm(n1/2, b_trt=0, b_age_trt=0, p_young=0.80, p_female=0.64)
t2_B<-sim_arm(n2/2, b_trt=b_B, b_age_trt=b_ageB, p_young=0.20, p_female=0.80)
t2_C<-sim_arm(n2/2, b_trt=0, b_age_trt=0, p_young=0.20, p_female=0.80)

trial1<-bind_rows(mutate(t1_A, arm="A", trial=1),
                  mutate(t1_C, arm="C", trial=1))
trial2<-bind_rows(mutate(t2_B, arm="B", trial=2),
                  mutate(t2_C, arm="C", trial=2))

write.csv(trial1, "trial1.csv", row.names=FALSE)
write.csv(trial2, "trial2.csv", row.names=FALSE)

trial1<-read.csv("trial1.csv")
trial2<-read.csv("trial2.csv")

MAIC<-function(the_data, id, match_mean, agg_means, agg_sds,
               accuracy=sqrt(.Machine$double.eps)) {
  the_data_o<-the_data
  match_sds<-which(1*(is.na(agg_sds))==0)
  also_match_var<-match_mean[match_sds]
  vars_keep<-c(id, match_mean)
  the_data_m<-dplyr::select(the_data, vars_keep)
  complete_case<-complete.cases(the_data_m)
  the_data_m<-na.omit(the_data_m)
  if(length(also_match_var)>0) {
    the_data_v<-dplyr::select(the_data_m, also_match_var)
    the_data_v<-the_data_v^2
    sq_list<-paste(also_match_var, "sq", sep="_")
    colnames(the_data_v)<-sq_list
    the_data_m<-cbind(the_data_m, the_data_v)
  }
  M1<-matrix(agg_means, nrow=nrow(the_data_m), ncol=length(match_mean), byrow=TRUE)
  if(length(also_match_var)>0) {
    mom2<-agg_means^2+agg_sds^2
    mom2<-mom2[!is.na(mom2)]
    M2<-matrix(mom2, nrow=nrow(the_data_m), ncol=length(also_match_var), byrow=TRUE)
    M1<-cbind(M1, M2)
  }
  the_data_m[, 2:ncol(the_data_m)]<-the_data_m[, 2:ncol(the_data_m)]-M1
  X<-as.matrix(the_data_m[, 2:ncol(the_data_m)])
  start<-rep(0, ncol(X))
  Answers<-weight(start, X, accuracy)
  the_data_o$weight.maic<-0
  the_data_o$weight.maic[which(1*(complete_case)==1)]<-Answers$w
  the_data_o$weight.maic<-the_data_o$weight.maic/sum(the_data_o$weight.maic)
  n_complete_obs<-nrow(X)
  A_mat<-diag(n_complete_obs)
  B_vec<-rep(0, n_complete_obs)
  E_mat<-rbind(rep(1, n_complete_obs), t(X))
  F_vec<-c(1, rep(0, ncol(X)))
  G_mat<-A_mat
  H_vec<-rep(0, n_complete_obs)
  weights_optimal<-lsei(A=A_mat, B=B_vec, E=E_mat, F=F_vec, G=G_mat, H=H_vec)$X
  the_data_o$weight.opt<-0
  the_data_o$weight.opt[which(1*(complete_case)==1)]<-weights_optimal
  return(list(the_data=the_data_o))
}

weight<-function(start, X, accuracy) {
  objective<-function(x, X) return(sum(exp(X%*%x)))
  gradient<-function(x, X) return(t(exp(X%*%x))%*%X)
  max_res<-optim(par=start, fn=objective, gr=gradient, method="BFGS", X=X,
                 control=list(maxit=100000, reltol=accuracy))
  weights<-exp(X%*%max_res$par)
  L<-list(weights, max_res$par, max_res$convergence)
  names(L)<-c("w", "alpha", "converge")
  return(L)
}

ESS<-function(weights) ((sum(weights))^2)/sum(weights^2)

wt_lm_eff<-function(Y_trt, Y_ctrl, w_trt, w_ctrl) {
  n_T<-length(Y_trt)
  n_C<-length(Y_ctrl)
  df<-data.frame(Y=c(Y_trt, Y_ctrl),
                 T=c(rep(1, n_T), rep(0, n_C)),
                 w=c(w_trt, w_ctrl))
  df<-df[df$w>0, ]
  fit<-lm(Y~T, data=df, weights=w)
  vcov_hc<-vcovHC(fit, type="HC3")
  list(delta=as.numeric(coef(fit)["T"]), var=vcov_hc["T", "T"])
}

unwt_lm_eff<-function(Y_trt, Y_ctrl) {
  wt_lm_eff(Y_trt, Y_ctrl,
            w_trt=rep(1, length(Y_trt)),
            w_ctrl=rep(1, length(Y_ctrl)))
}

agg1<-list(mean_age=mean(trial1$age), sd_age=sd(trial1$age))
agg2<-list(mean_age=mean(trial2$age), sd_age=sd(trial2$age))

cat("\n=== Aggregate moments: age ===\n")
cat(sprintf("Trial 1 -- mean: %.2f  SD: %.2f\n", agg1$mean_age, agg1$sd_age))
cat(sprintf("Trial 2 -- mean: %.2f  SD: %.2f\n", agg2$mean_age, agg2$sd_age))

trial1_id<-trial1 %>% mutate(ID=row_number())
trial2_id<-trial2 %>% mutate(ID=row_number())

cat("\n--- MAIC: Trial 2 reweighted to Trial 1 ---\n")
maic_t2_to_t1<-MAIC(trial2_id, "ID", "age", agg1$mean_age, agg1$sd_age)
w_t2_to_t1<-maic_t2_to_t1$the_data$weight.maic

cat("--- MAIC: Trial 1 reweighted to Trial 2 ---\n")
maic_t1_to_t2<-MAIC(trial1_id, "ID", "age", agg2$mean_age, agg2$sd_age)
w_t1_to_t2<-maic_t1_to_t2$the_data$weight.maic

t1_A_df<-trial1 %>% filter(arm=="A")
t1_C_df<-trial1 %>% filter(arm=="C")
t2_B_df<-trial2 %>% filter(arm=="B")
t2_C_df<-trial2 %>% filter(arm=="C")

idx_A1<-which(trial1$arm=="A")
idx_C1<-which(trial1$arm=="C")
idx_B2<-which(trial2$arm=="B")
idx_C2<-which(trial2$arm=="C")

w_A1<-w_t1_to_t2[idx_A1]
w_C1<-w_t1_to_t2[idx_C1]
w_B2<-w_t2_to_t1[idx_B2]
w_C2<-w_t2_to_t1[idx_C2]

res_AC1<-unwt_lm_eff(t1_A_df$Y, t1_C_df$Y)
delta_AC1<-res_AC1$delta
V_AC1<-res_AC1$var

res_BC2<-unwt_lm_eff(t2_B_df$Y, t2_C_df$Y)
delta_BC2<-res_BC2$delta
V_BC2<-res_BC2$var

res_BC1<-wt_lm_eff(t2_B_df$Y, t2_C_df$Y, w_B2, w_C2)
delta_BC1<-res_BC1$delta
V_BC1<-res_BC1$var

res_AC2<-wt_lm_eff(t1_A_df$Y, t1_C_df$Y, w_A1, w_C1)
delta_AC2<-res_AC2$delta
V_AC2<-res_AC2$var

cat("\n=== Four delta components ===\n")
cat(sprintf("delta_AC1 (direct,  T1 pop): %7.4f | SE: %.6f\n", delta_AC1, sqrt(V_AC1)))
cat(sprintf("delta_BC2 (direct,  T2 pop): %7.4f | SE: %.6f\n", delta_BC2, sqrt(V_BC2)))
cat(sprintf("delta_BC1 (MAIC, T2->T1):    %7.4f | SE: %.6f\n", delta_BC1, sqrt(V_BC1)))
cat(sprintf("delta_AC2 (MAIC, T1->T2):    %7.4f | SE: %.6f\n", delta_AC2, sqrt(V_AC2)))

N<-nrow(trial1)+nrow(trial2)
w1<-nrow(trial1)/N
w2<-nrow(trial2)/N

ATE_maic_pooled<-w1*(delta_AC1-delta_BC1)+w2*(delta_AC2-delta_BC2)
Var_maic_pooled<-w1^2*(V_AC1+V_BC1)+w2^2*(V_BC2+V_AC2)+
                 2*w1*w2*(V_AC1+V_BC2)
SE_maic_pooled<-sqrt(Var_maic_pooled)
CI_lo_pooled<-ATE_maic_pooled-1.96*SE_maic_pooled
CI_hi_pooled<-ATE_maic_pooled+1.96*SE_maic_pooled

ATE_maic_t1<-delta_AC1-delta_BC1
Var_maic_t1<-V_AC1+V_BC1
SE_maic_t1<-sqrt(Var_maic_t1)
CI_lo_t1<-ATE_maic_t1-1.96*SE_maic_t1
CI_hi_t1<-ATE_maic_t1+1.96*SE_maic_t1

ATE_maic_t2<-delta_AC2-delta_BC2
Var_maic_t2<-V_AC2+V_BC2
SE_maic_t2<-sqrt(Var_maic_t2)
CI_lo_t2<-ATE_maic_t2-1.96*SE_maic_t2
CI_hi_t2<-ATE_maic_t2+1.96*SE_maic_t2

cat("\n=== MAIC-based ATE (pooled population) ===\n")
cat(sprintf("ATE  : %6.4f\n", ATE_maic_pooled))
cat(sprintf("SE   : %6.4f\n", SE_maic_pooled))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", CI_lo_pooled, CI_hi_pooled))

cat("\n=== MAIC-based ATE (Trial 1 population) ===\n")
cat(sprintf("ATE  : %6.4f\n", ATE_maic_t1))
cat(sprintf("SE   : %6.4f\n", SE_maic_t1))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", CI_lo_t1, CI_hi_t1))

cat("\n=== MAIC-based ATE (Trial 2 population) ===\n")
cat(sprintf("ATE  : %6.4f\n", ATE_maic_t2))
cat(sprintf("SE   : %6.4f\n", SE_maic_t2))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", CI_lo_t2, CI_hi_t2))

combined<-bind_rows(trial1, trial2) %>%
  mutate(S=as.integer(trial==1),
         arm_A=as.integer(arm=="A"),
         arm_B=as.integer(arm=="B"))

W_fit<-weightit(S~age, data=combined, method="glm", estimand="ATE")
combined$ipw<-W_fit$weights

fit_ipw<-lm(Y~arm_A+arm_B+as.factor(trial), data=combined, weights=ipw)
vcov_ipw<-vcovHC(fit_ipw, type="HC0")
coef_A<-coef(fit_ipw)["arm_A"]
coef_B<-coef(fit_ipw)["arm_B"]
ATE_ipw<-as.numeric(coef_A-coef_B)

Var_ipw<-vcov_ipw["arm_A","arm_A"]+vcov_ipw["arm_B","arm_B"]-
         2*vcov_ipw["arm_A","arm_B"]
SE_ipw<-sqrt(Var_ipw)
CI_lo_ipw<-ATE_ipw-1.96*SE_ipw
CI_hi_ipw<-ATE_ipw+1.96*SE_ipw

#unadjusted
ATE_un<-delta_AC1-delta_BC2
SE_un<-sqrt(V_AC1+V_BC2)
CI_lo_un<-ATE_un-1.96*SE_un
CI_hi_un<-ATE_un+1.96*SE_un

cat("\n")
cat("==========================================================\n")
cat("       COMPARISON: MAIC-based vs IPW-based ATE\n")
cat("==========================================================\n")
cat(sprintf("%-26s %8s %8s %22s\n", "Method", "ATE", "SE", "95% CI"))
cat("----------------------------------------------------------\n")
cat(sprintf("%-26s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "MAIC pooled", ATE_maic_pooled, SE_maic_pooled, CI_lo_pooled, CI_hi_pooled))
cat(sprintf("%-26s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "MAIC Trial 1 pop", ATE_maic_t1, SE_maic_t1, CI_lo_t1, CI_hi_t1))
cat(sprintf("%-26s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "MAIC Trial 2 pop", ATE_maic_t2, SE_maic_t2, CI_lo_t2, CI_hi_t2))
cat(sprintf("%-26s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "IPTW pooled IPD", ATE_ipw, SE_ipw, CI_lo_ipw, CI_hi_ipw))
cat(sprintf("%-26s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "Unadjusted", ATE_un, SE_un, CI_lo_un, CI_hi_un))
cat("==========================================================\n")
