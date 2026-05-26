
library(dplyr)
library(sandwich)   # vcovHC
library(limSolve)   # lsei() for MAIC max-ESS weights
library(WeightIt)   # IPW weights in Part 3
options(warn = -1)
set.seed(2)

sim_age <- function(n, p_young = 0.7) {
  young <- rbinom(n, 1, prob = p_young)
  ifelse(young == 1, runif(n, 45, 60), runif(n, 60, 75))
}


sim_age <- function(n, p_young = 0.8) {
  young <- rbinom(n, 1, prob = p_young)
  ifelse(young == 1, runif(n, 45, 60), runif(n, 60, 75))
}
#parameters from Jiang et al. (2024) Section 5
a0           <-  0.85
a_male       <-  0.12
a_age        <-  0.05
beta_A       <- -2.10
beta_B       <- -2.50
b_ageA       <- 0.03
b_ageB       <-  -0.08
mean_age_ctr <-  40

n1 <- 300
n2 <- 500

sim_arm <- function(n, beta_t,a_age_trt, p_young, p_female) {
  age    <- sim_age(n, p_young)
  female <- rbinom(n, 1, prob = p_female)
  male   <- 1 - female
  is_trt     <- as.integer(beta_t != 0)
  trt_effect <- is_trt * (beta_t + a_age_trt * (age - mean_age_ctr))
  lp   <- a0 + a_male * male + a_age * (age - mean_age_ctr) + trt_effect
  prob <- exp(lp) / (1 + exp(lp))
  Y    <- rbinom(n, 1, prob = prob)
  data.frame(age = age, female = female, male = male, Y = Y)
}

t1_A <- sim_arm(n1/2, beta_t = beta_A,a_age_trt=b_ageA, p_young = 0.80, p_female = 0.64)
t1_C <- sim_arm(n1/2, beta_t = 0,   a_age_trt=0,   p_young = 0.80, p_female = 0.64)
t2_B <- sim_arm(n2/2, beta_t = beta_B, a_age_trt=b_ageB,p_young = 0.20, p_female = 0.80)
t2_C <- sim_arm(n2/2, beta_t = 0,    a_age_trt=0,  p_young = 0.20, p_female = 0.80)

trial1 <- bind_rows(mutate(t1_A, arm = "A", trial = 1),
                    mutate(t1_C, arm = "C", trial = 1))
trial2 <- bind_rows(mutate(t2_B, arm = "B", trial = 2),
                    mutate(t2_C, arm = "C", trial = 2))


write.csv(trial1, "trial1_binary.csv", row.names = FALSE)
write.csv(trial2, "trial2_binary.csv", row.names = FALSE)



trial1 <- read.csv("trial1_binary.csv")
trial2 <- read.csv("trial2_binary.csv")

MAIC <- function(the_data, id, match_mean, agg_means, agg_sds,
                 accuracy = sqrt(.Machine$double.eps)) {
  the_data_o <- the_data
  match_sds      <- which(1 * (is.na(agg_sds)) == 0)
  also_match_var <- match_mean[match_sds]
  vars_keep      <- c(id, match_mean)
  the_data_m     <- dplyr::select(the_data, vars_keep)
  complete_case  <- complete.cases(the_data_m)
  the_data_m     <- na.omit(the_data_m)
  if (length(also_match_var) > 0) {
    the_data_v           <- dplyr::select(the_data_m, also_match_var)
    the_data_v           <- the_data_v^2
    sq_list              <- paste(also_match_var, "sq", sep = "_")
    colnames(the_data_v) <- sq_list
    the_data_m           <- cbind(the_data_m, the_data_v)
  }
  M1 <- matrix(agg_means, nrow = nrow(the_data_m),
               ncol = length(match_mean), byrow = TRUE)
  if (length(also_match_var) > 0) {
    mom2 <- agg_means^2 + agg_sds^2
    mom2 <- mom2[!is.na(mom2)]
    M2   <- matrix(mom2, nrow = nrow(the_data_m),
                   ncol = length(also_match_var), byrow = TRUE)
    M1   <- cbind(M1, M2)
  }
  the_data_m[, 2:ncol(the_data_m)] <- the_data_m[, 2:ncol(the_data_m)] - M1
  X       <- as.matrix(the_data_m[, 2:ncol(the_data_m)])
  start   <- rep(0, ncol(X))
  Answers <- weight(start, X, accuracy)
  the_data_o$weight.maic <- 0
  the_data_o$weight.maic[which(1*(complete_case)==1)] <- Answers$w
  the_data_o$weight.maic <- the_data_o$weight.maic / sum(the_data_o$weight.maic)
  n_complete_obs <- nrow(X)
  A_mat <- diag(n_complete_obs);    B_vec <- rep(0, n_complete_obs)
  E_mat <- rbind(rep(1, n_complete_obs), t(X))
  F_vec <- c(1, rep(0, ncol(X)));  G_mat <- A_mat;  H_vec <- rep(0, n_complete_obs)
  weights_optimal <- lsei(A=A_mat, B=B_vec, E=E_mat, F=F_vec, G=G_mat, H=H_vec)$X
  the_data_o$weight.opt <- 0
  the_data_o$weight.opt[which(1*(complete_case)==1)] <- weights_optimal
  return(list(the_data = the_data_o))
}

weight <- function(start, X, accuracy) {
  objective <- function(x, X) return(sum(exp(X %*% x)))
  gradient  <- function(x, X) return(t(exp(X %*% x)) %*% X)
  max_res   <- optim(par=start, fn=objective, gr=gradient, method="BFGS", X=X,
                     control=list(maxit=100000, reltol=accuracy))
  weights   <- exp(X %*% max_res$par)
  L <- list(weights, max_res$par, max_res$convergence)
  names(L) <- c("w","alpha","converge")
  return(L)
}

ESS <- function(weights) ((sum(weights))^2) / sum(weights^2)

# ==============================================================================
# Helper: returns weighted proportions and arm sizes needed for all estimators
# ==============================================================================
glm_logOR <- function(Y_trt, Y_ctrl, w_trt, w_ctrl) {
  n_T <- length(Y_trt);  n_C <- length(Y_ctrl)
  df  <- data.frame(Y = c(Y_trt, Y_ctrl),
                    T = c(rep(1, n_T), rep(0, n_C)),
                    w = c(w_trt, w_ctrl))
  df      <- df[df$w > 0, ]
  fit     <- glm(Y ~ T, family = binomial(link = "logit"), data = df, weights = w)
  vcov_hc <- vcovHC(fit, type = "HC3")
  # Hajek weighted proportions (needed for covariance formula from document)
  p_trt  <- sum(w_trt  * Y_trt)  / sum(w_trt)
  p_ctrl <- sum(w_ctrl * Y_ctrl) / sum(w_ctrl)
  list(lor    = as.numeric(coef(fit)["T"]),
       var    = vcov_hc["T", "T"],
       p_trt  = p_trt,
       p_ctrl = p_ctrl,
       n_trt  = n_T,
       n_ctrl = n_C)
}

unwt_glm_logOR <- function(Y_trt, Y_ctrl) {
  glm_logOR(Y_trt, Y_ctrl,
            w_trt  = rep(1, length(Y_trt)),
            w_ctrl = rep(1, length(Y_ctrl)))
}


glm_logRR <- function(Y_trt, Y_ctrl, w_trt, w_ctrl) {
  n_T <- length(Y_trt);  n_C <- length(Y_ctrl)
  df  <- data.frame(Y = c(Y_trt, Y_ctrl),
                    T = c(rep(1, n_T), rep(0, n_C)),
                    w = c(w_trt, w_ctrl))
  df      <- df[df$w > 0, ]
  fit     <- glm(Y ~ T, family = binomial(link = "log"), data = df, weights = w,
                 start = c(log(mean(Y_ctrl)), 0))
  vcov_hc <- vcovHC(fit, type = "HC3")
  p_trt  <- sum(w_trt  * Y_trt)  / sum(w_trt)
  p_ctrl <- sum(w_ctrl * Y_ctrl) / sum(w_ctrl)
  list(lrr    = as.numeric(coef(fit)["T"]),
       var    = vcov_hc["T", "T"],
       p_trt  = p_trt,
       p_ctrl = p_ctrl,
       n_trt  = n_T,
       n_ctrl = n_C)
}

unwt_glm_logRR <- function(Y_trt, Y_ctrl) {
  glm_logRR(Y_trt, Y_ctrl,
            w_trt  = rep(1, length(Y_trt)),
            w_ctrl = rep(1, length(Y_ctrl)))
}


glm_RD <- function(Y_trt, Y_ctrl, w_trt, w_ctrl) {
  n_T <- length(Y_trt);  n_C <- length(Y_ctrl)
  df  <- data.frame(Y = c(Y_trt, Y_ctrl),
                    T = c(rep(1, n_T), rep(0, n_C)),
                    w = c(w_trt, w_ctrl))
  df      <- df[df$w > 0, ]
  fit     <- glm(Y ~ T, family = binomial(link = "identity"), data = df, weights = w,
                 start = c(mean(Y_ctrl), 0))
  vcov_hc <- vcovHC(fit, type = "HC3")
  p_trt  <- sum(w_trt  * Y_trt)  / sum(w_trt)
  p_ctrl <- sum(w_ctrl * Y_ctrl) / sum(w_ctrl)
  list(rd     = as.numeric(coef(fit)["T"]),
       var    = vcov_hc["T", "T"],
       p_trt  = p_trt,
       p_ctrl = p_ctrl,
       n_trt  = n_T,
       n_ctrl = n_C)
}

unwt_glm_RD <- function(Y_trt, Y_ctrl) {
  glm_RD(Y_trt, Y_ctrl,
         w_trt  = rep(1, length(Y_trt)),
         w_ctrl = rep(1, length(Y_ctrl)))
}


agg1 <- list(mean_age = mean(trial1$age), sd_age = sd(trial1$age))
agg2 <- list(mean_age = mean(trial2$age), sd_age = sd(trial2$age))

cat("\n=== Aggregate moments: age ===\n")
cat(sprintf("Trial 1 -- mean: %.2f  SD: %.2f\n", agg1$mean_age, agg1$sd_age))
cat(sprintf("Trial 2 -- mean: %.2f  SD: %.2f\n", agg2$mean_age, agg2$sd_age))


# MAIC weights
trial1_id <- trial1 %>% mutate(ID = row_number())
trial2_id <- trial2 %>% mutate(ID = row_number())
maic_t2_to_t1 <- MAIC(trial2_id, "ID", "age", agg1$mean_age, agg1$sd_age)
w_t2_to_t1    <- maic_t2_to_t1$the_data$weight.maic

maic_t1_to_t2 <- MAIC(trial1_id, "ID", "age", agg2$mean_age, agg2$sd_age)
w_t1_to_t2    <- maic_t1_to_t2$the_data$weight.maic


# Arm-level subsets and weights
t1_A_df <- trial1 %>% filter(arm == "A")
t1_C_df <- trial1 %>% filter(arm == "C")
t2_B_df <- trial2 %>% filter(arm == "B")
t2_C_df <- trial2 %>% filter(arm == "C")

idx_A1 <- which(trial1$arm == "A");  idx_C1 <- which(trial1$arm == "C")
idx_B2 <- which(trial2$arm == "B");  idx_C2 <- which(trial2$arm == "C")

w_A1 <- w_t1_to_t2[idx_A1];  w_C1 <- w_t1_to_t2[idx_C1]
w_B2 <- w_t2_to_t1[idx_B2];  w_C2 <- w_t2_to_t1[idx_C2]

# Population weights
N  <- nrow(trial1) + nrow(trial2)
w1 <- nrow(trial1) / N
w2 <- nrow(trial2) / N


#log OR

res_AC1_or   <- unwt_glm_logOR(t1_A_df$Y, t1_C_df$Y)
theta_AC1_or <- res_AC1_or$lor;  V_AC1_or <- res_AC1_or$var

res_BC2_or   <- unwt_glm_logOR(t2_B_df$Y, t2_C_df$Y)
theta_BC2_or <- res_BC2_or$lor;  V_BC2_or <- res_BC2_or$var

res_BC1_or   <- glm_logOR(t2_B_df$Y, t2_C_df$Y, w_B2, w_C2)
theta_BC1_or <- res_BC1_or$lor;  V_BC1_or <- res_BC1_or$var

res_AC2_or   <- glm_logOR(t1_A_df$Y, t1_C_df$Y, w_A1, w_C1)
theta_AC2_or <- res_AC2_or$lor;  V_AC2_or <- res_AC2_or$var

cat("\n=== Four delta components (log OR) ===\n")
cat(sprintf("theta_AC1 (direct,  T1 pop): %7.4f | SE: %.6f\n", theta_AC1_or, sqrt(V_AC1_or)))
cat(sprintf("theta_BC2 (direct,  T2 pop): %7.4f | SE: %.6f\n", theta_BC2_or, sqrt(V_BC2_or)))
cat(sprintf("theta_BC1 (MAIC, T2->T1):    %7.4f | SE: %.6f\n", theta_BC1_or, sqrt(V_BC1_or)))
cat(sprintf("theta_AC2 (MAIC, T1->T2):    %7.4f | SE: %.6f\n", theta_AC2_or, sqrt(V_AC2_or)))

# Pooled log OR 
logOR_maic <- w1*(theta_AC1_or - theta_BC1_or) + w2*(theta_AC2_or - theta_BC2_or)

Cov_AC_or <- 1/(res_AC1_or$n_trt  * res_AC2_or$p_trt  * (1 - res_AC2_or$p_trt))  +
             1/(res_AC1_or$n_ctrl * res_AC2_or$p_ctrl * (1 - res_AC2_or$p_ctrl))

Cov_BC_or <- 1/(res_BC2_or$n_trt  * res_BC1_or$p_trt  * (1 - res_BC1_or$p_trt))  +
             1/(res_BC2_or$n_ctrl * res_BC1_or$p_ctrl * (1 - res_BC1_or$p_ctrl))

# Var(theta_AB) = w1^2*(V_AC1 + V_BC1w) + w2^2*(V_AC2w + V_BC2) + 2*w1*w2*(C1+C2)
Var_maic_or <- w1^2*(V_AC1_or + V_BC1_or) + w2^2*(V_AC2_or + V_BC2_or) +
               2*w1*w2*(Cov_AC_or + Cov_BC_or)

SE_maic_or <- sqrt(Var_maic_or)
CI_lo_or   <- logOR_maic - 1.96*SE_maic_or
CI_hi_or   <- logOR_maic + 1.96*SE_maic_or

# (log RR)
res_AC1_rr   <- unwt_glm_logRR(t1_A_df$Y, t1_C_df$Y)
theta_AC1_rr <- res_AC1_rr$lrr;  V_AC1_rr <- res_AC1_rr$var

res_BC2_rr   <- unwt_glm_logRR(t2_B_df$Y, t2_C_df$Y)
theta_BC2_rr <- res_BC2_rr$lrr;  V_BC2_rr <- res_BC2_rr$var

res_BC1_rr   <- glm_logRR(t2_B_df$Y, t2_C_df$Y, w_B2, w_C2)
theta_BC1_rr <- res_BC1_rr$lrr;  V_BC1_rr <- res_BC1_rr$var

res_AC2_rr   <- glm_logRR(t1_A_df$Y, t1_C_df$Y, w_A1, w_C1)
theta_AC2_rr <- res_AC2_rr$lrr;  V_AC2_rr <- res_AC2_rr$var

cat("\n=== Four delta components (log RR) ===\n")
cat(sprintf("theta_AC1 (direct,  T1 pop): %7.4f | SE: %.6f\n", theta_AC1_rr, sqrt(V_AC1_rr)))
cat(sprintf("theta_BC2 (direct,  T2 pop): %7.4f | SE: %.6f\n", theta_BC2_rr, sqrt(V_BC2_rr)))
cat(sprintf("theta_BC1 (MAIC, T2->T1):    %7.4f | SE: %.6f\n", theta_BC1_rr, sqrt(V_BC1_rr)))
cat(sprintf("theta_AC2 (MAIC, T1->T2):    %7.4f | SE: %.6f\n", theta_AC2_rr, sqrt(V_AC2_rr)))

# Pooled log RR
logRR_maic <- w1*(theta_AC1_rr - theta_BC1_rr) + w2*(theta_AC2_rr - theta_BC2_rr)

Cov_AC_rr <- (1 - res_AC1_rr$p_trt)  / (res_AC1_rr$n_trt  * res_AC2_rr$p_trt)  +
             (1 - res_AC1_rr$p_ctrl) / (res_AC1_rr$n_ctrl * res_AC2_rr$p_ctrl)

Cov_BC_rr <- (1 - res_BC2_rr$p_trt)  / (res_BC2_rr$n_trt  * res_BC1_rr$p_trt)  +
             (1 - res_BC2_rr$p_ctrl) / (res_BC2_rr$n_ctrl * res_BC1_rr$p_ctrl)

Var_maic_rr <- w1^2*(V_AC1_rr + V_BC1_rr) + w2^2*(V_AC2_rr + V_BC2_rr) +
               2*w1*w2*(Cov_AC_rr + Cov_BC_rr)

SE_maic_rr <- sqrt(Var_maic_rr)
CI_lo_rr   <- logRR_maic - 1.96*SE_maic_rr
CI_hi_rr   <- logRR_maic + 1.96*SE_maic_rr


# RD

res_AC1_rd   <- unwt_glm_RD(t1_A_df$Y, t1_C_df$Y)
theta_AC1_rd <- res_AC1_rd$rd;  V_AC1_rd <- res_AC1_rd$var

res_BC2_rd   <- unwt_glm_RD(t2_B_df$Y, t2_C_df$Y)
theta_BC2_rd <- res_BC2_rd$rd;  V_BC2_rd <- res_BC2_rd$var

res_BC1_rd   <- glm_RD(t2_B_df$Y, t2_C_df$Y, w_B2, w_C2)
theta_BC1_rd <- res_BC1_rd$rd;  V_BC1_rd <- res_BC1_rd$var

res_AC2_rd   <- glm_RD(t1_A_df$Y, t1_C_df$Y, w_A1, w_C1)
theta_AC2_rd <- res_AC2_rd$rd;  V_AC2_rd <- res_AC2_rd$var


# Pooled RD
RD_maic <- w1*(theta_AC1_rd - theta_BC1_rd) + w2*(theta_AC2_rd - theta_BC2_rd)

Cov_AC_rd <- V_AC1_rd   
Cov_BC_rd <- V_BC2_rd   

Var_maic_rd <- w1^2*(V_AC1_rd + V_BC1_rd) + w2^2*(V_AC2_rd + V_BC2_rd) +
               2*w1*w2*(Cov_AC_rd + Cov_BC_rd)

SE_maic_rd <- sqrt(Var_maic_rd)
CI_lo_rd   <- RD_maic - 1.96*SE_maic_rd
CI_hi_rd   <- RD_maic + 1.96*SE_maic_rd

cat("\n=== Four delta components (RD) ===\n")
cat(sprintf("theta_AC1 (direct,  T1 pop): %7.4f | SE: %.6f\n", theta_AC1_rd, sqrt(V_AC1_rd)))
cat(sprintf("theta_BC2 (direct,  T2 pop): %7.4f | SE: %.6f\n", theta_BC2_rd, sqrt(V_BC2_rd)))
cat(sprintf("theta_BC1 (MAIC, T2->T1):    %7.4f | SE: %.6f\n", theta_BC1_rd, sqrt(V_BC1_rd)))
cat(sprintf("theta_AC2 (MAIC, T1->T2):    %7.4f | SE: %.6f\n", theta_AC2_rd, sqrt(V_AC2_rd)))

# IPW
combined <- bind_rows(trial1, trial2) %>%
  mutate(S     = as.integer(trial == 1),
         arm_A = as.integer(arm == "A"),
         arm_B = as.integer(arm == "B"))

W_fit        <- weightit(S ~ age, data = combined,
                         method = "glm", estimand = "ATE")
combined$ipw <- W_fit$weights


# IPW: Log Odds Ratio

fit_ipw_or <- glm(Y ~ arm_A + arm_B + trial, data = combined,
                  family = binomial(link = "logit"), weights = ipw)

vcov_ipw_or  <- vcovHC(fit_ipw_or, type = "HC0")
coef_A_or    <- coef(fit_ipw_or)["arm_A"]
coef_B_or    <- coef(fit_ipw_or)["arm_B"]
logOR_ipw    <- as.numeric(coef_A_or - coef_B_or)

Var_ipw_or   <- vcov_ipw_or["arm_A","arm_A"] + vcov_ipw_or["arm_B","arm_B"] -
  2*vcov_ipw_or["arm_A","arm_B"]
SE_ipw_or    <- sqrt(Var_ipw_or)
CI_lo_ipw_or <- logOR_ipw - 1.96*SE_ipw_or
CI_hi_ipw_or <- logOR_ipw + 1.96*SE_ipw_or



# IPW: Log Risk Ratio
fit_ipw_rr <- glm(Y ~ arm_A + arm_B + trial, data = combined,
                  family = binomial(link = "log"), weights = ipw,
                  start = c(log(mean(combined$Y[combined$arm_A == 0 & combined$arm_B == 0])),
                            0, 0, 0))

vcov_ipw_rr  <- vcovHC(fit_ipw_rr, type = "HC0")
coef_A_rr    <- coef(fit_ipw_rr)["arm_A"]
coef_B_rr    <- coef(fit_ipw_rr)["arm_B"]
logRR_ipw    <- as.numeric(coef_A_rr - coef_B_rr)

Var_ipw_rr   <- vcov_ipw_rr["arm_A","arm_A"] + vcov_ipw_rr["arm_B","arm_B"] -
  2*vcov_ipw_rr["arm_A","arm_B"]
SE_ipw_rr    <- sqrt(Var_ipw_rr)
CI_lo_ipw_rr <- logRR_ipw - 1.96*SE_ipw_rr
CI_hi_ipw_rr <- logRR_ipw + 1.96*SE_ipw_rr


# IPW: Risk Difference
fit_ipw_rd <- glm(Y ~ arm_A + arm_B + trial, data = combined,
                  family = binomial(link = "identity"), weights = ipw,
                  start = c(mean(combined$Y[combined$arm_A == 0 & combined$arm_B == 0]),
                            0, 0, 0))

vcov_ipw_rd  <- vcovHC(fit_ipw_rd, type = "HC0")
coef_A_rd    <- coef(fit_ipw_rd)["arm_A"]
coef_B_rd    <- coef(fit_ipw_rd)["arm_B"]
RD_ipw       <- as.numeric(coef_A_rd - coef_B_rd)

Var_ipw_rd   <- vcov_ipw_rd["arm_A","arm_A"] + vcov_ipw_rd["arm_B","arm_B"] -
  2*vcov_ipw_rd["arm_A","arm_B"]
SE_ipw_rd    <- sqrt(Var_ipw_rd)
CI_lo_ipw_rd <- RD_ipw - 1.96*SE_ipw_rd
CI_hi_ipw_rd <- RD_ipw + 1.96*SE_ipw_rd


################################################################################
# COMPARISON TABLE
################################################################################
cat("\n")
cat("=============================================================================\n")
cat("        COMPARISON: MAIC-based vs IPW-based — all effect measures\n")
cat("=============================================================================\n")
cat(sprintf("%-22s %-8s %8s %8s %26s\n", "Method", "TE", "Estimate", "SE", "95% CI"))
cat("-----------------------------------------------------------------------------\n")
# MAIC log OR
cat(sprintf("%-22s %-8s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "MAIC", "log OR", logOR_maic, SE_maic_or, CI_lo_or, CI_hi_or))
# MAIC log RR
cat(sprintf("%-22s %-8s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "MAIC", "log RR", logRR_maic, SE_maic_rr, CI_lo_rr, CI_hi_rr))
# MAIC RD
cat(sprintf("%-22s %-8s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "MAIC", "RD", RD_maic, SE_maic_rd, CI_lo_rd, CI_hi_rd))
cat("-----------------------------------------------------------------------------\n")
# IPW log OR
cat(sprintf("%-22s %-8s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "IPW", "log OR", logOR_ipw, SE_ipw_or, CI_lo_ipw_or, CI_hi_ipw_or))
# IPW log RR
cat(sprintf("%-22s %-8s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "IPW", "log RR", logRR_ipw, SE_ipw_rr, CI_lo_ipw_rr, CI_hi_ipw_rr))
# IPW RD
cat(sprintf("%-22s %-8s %8.4f %8.4f  [%7.4f, %7.4f]\n",
            "IPW", "RD", RD_ipw, SE_ipw_rd, CI_lo_ipw_rd, CI_hi_ipw_rd))
cat("=============================================================================\n")
