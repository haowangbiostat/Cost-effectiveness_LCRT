############################################################################################################
#                                                                                                          #
#        program:        utils_SWCRT_incomplete.R                                                          #
#                                                                                                          #
#        usage:          to propose LOD and MMD in stepped-wedge CRTs for cost-effectiveness              #
#                        supports incomplete designs with NA in treatment matrix                           #
#                                                                                                          #
#        assumption:     nested exchangeable correlation structure                                         #
#                        equal cluster size across time periods                                            #
#                        cross-sectional design (different individuals in each cluster-period)             #
#                                                                                                          #
#        required:       base R, nloptr package for optimization                                           #
#                                                                                                          #
############################################################################################################

# load required libraries
if (!require(nloptr)) {
  install.packages("nloptr")
  library(nloptr)
}

# ============================================================================================================
# helper functions
# ============================================================================================================

# compute variance of INMB estimator for stepped-wedge design (numerical approach for incomplete designs)
compute_variance_SWCRT <- function(Z, K, rho0_E, rho1_E, rho0_C, rho1_C, 
                                   rho0_EC, rho1_EC, rho2_EC, lambda, sigma_E, sigma_C) {
  I <- nrow(Z)
  J <- ncol(Z)
  
  # rho0_E >= rho1_E
  if (rho0_E < rho1_E) return(NA_real_)
  
  # rho0_C >= rho1_C
  if (rho0_C < rho1_C) return(NA_real_)
  
  # rho0_EC <= min(rho0_E, rho0_C)
  if (rho0_EC > rho0_E || rho0_EC > rho0_C) return(NA_real_)
  
  # rho1_EC <= min(rho1_E, rho1_C)
  if (rho1_EC > rho1_E || rho1_EC > rho1_C) return(NA_real_)
  
  # rho0_EC >= rho1_EC
  if (rho0_EC < rho1_EC) return(NA_real_)
  
  # rho2_EC >= rho0_EC
  if (rho2_EC < rho0_EC) return(NA_real_)
  
  # compute kappa values
  kappa_E <- 1 + (K - 1) * rho0_E - K * rho1_E
  kappa_C <- 1 + (K - 1) * rho0_C - K * rho1_C
  kappa_EC <- rho2_EC + (K - 1) * rho0_EC - K * rho1_EC
  
  # check all eigenvalues are positive
  xi_1 <- ((K - 1) * (rho0_E - rho0_C) + (J - 1) * K * (rho1_E - rho1_C))^2 + 
    4 * (rho2_EC + (K - 1) * rho0_EC + (J - 1) * K * rho1_EC)^2
  
  lambda_1_plus <- 0.5 * ((2 + (K - 1) * (rho0_E + rho0_C) + (J - 1) * K * (rho1_E + rho1_C)) + sqrt(xi_1))
  lambda_1_minus <- 0.5 * ((2 + (K - 1) * (rho0_E + rho0_C) + (J - 1) * K * (rho1_E + rho1_C)) - sqrt(xi_1))
  lambda_2_plus <- 0.5 * ((kappa_E + kappa_C) + sqrt((kappa_E - kappa_C)^2 + 4 * kappa_EC^2))
  lambda_2_minus <- 0.5 * ((kappa_E + kappa_C) - sqrt((kappa_E - kappa_C)^2 + 4 * kappa_EC^2))
  lambda_3_plus <- 0.5 * ((2 - rho0_E - rho0_C) + sqrt((rho0_E - rho0_C)^2 + 4 * (rho2_EC - rho0_EC)^2))
  lambda_3_minus <- 0.5 * ((2 - rho0_E - rho0_C) - sqrt((rho0_E - rho0_C)^2 + 4 * (rho2_EC - rho0_EC)^2))
  
  if (lambda_1_plus <= 0 || lambda_1_minus <= 0 || 
      lambda_2_plus <= 0 || lambda_2_minus <= 0 || 
      lambda_3_plus <= 0 || lambda_3_minus <= 0) return(NA_real_)
  
  # covariance matrices for random effects
  Sigma_b <- matrix(c(rho1_E * sigma_E^2, rho1_EC * sigma_E * sigma_C,
                      rho1_EC * sigma_E * sigma_C, rho1_C * sigma_C^2), 2, 2)
  Sigma_s <- matrix(c((rho0_E - rho1_E) * sigma_E^2, (rho0_EC - rho1_EC) * sigma_E * sigma_C,
                      (rho0_EC - rho1_EC) * sigma_E * sigma_C, (rho0_C - rho1_C) * sigma_C^2), 2, 2)
  Sigma_e <- matrix(c((1 - rho0_E) * sigma_E^2, (rho2_EC - rho0_EC) * sigma_E * sigma_C,
                      (rho2_EC - rho0_EC) * sigma_E * sigma_C, (1 - rho0_C) * sigma_C^2), 2, 2)
  Sigma_within <- Sigma_s + Sigma_e / K
  
  # accumulate information matrix
  info_mat <- matrix(0, 2 * J + 2, 2 * J + 2)
  
  for (i in 1:I) {
    obs <- which(!is.na(Z[i, ]))
    J_i <- length(obs)
    if (J_i == 0) next
    
    # V_i for observed periods only
    I_Ji <- diag(J_i)
    ones_Ji <- matrix(1, J_i, J_i)
    V_i <- kronecker(I_Ji, Sigma_within) + kronecker(ones_Ji, Sigma_b)
    
    V_inv <- tryCatch(solve(V_i), error = function(e) NULL)
    if (is.null(V_inv)) return(NA_real_)
    
    # D_i with period indicators for observed periods
    period_ind <- matrix(0, J_i, J)
    for (k in seq_along(obs)) period_ind[k, obs[k]] <- 1
    X_i <- Z[i, obs]
    D_i <- kronecker(cbind(period_ind, X_i), diag(2))
    
    info_mat <- info_mat + t(D_i) %*% V_inv %*% D_i
  }
  
  # check invertibility
  cov_mat <- tryCatch(solve(info_mat), error = function(e) NULL)
  if (is.null(cov_mat)) return(NA_real_)
  
  # extract 2x2 block for (alpha_1, gamma_1)
  idx <- (2 * J + 1):(2 * J + 2)
  Sigma_trt <- cov_mat[idx, idx]
  
  # variance of INMB: Var(lambda * alpha_1 - gamma_1)
  var_beta <- lambda^2 * Sigma_trt[1, 1] + Sigma_trt[2, 2] - 2 * lambda * Sigma_trt[1, 2]
  
  if (var_beta <= 0) return(NA_real_)
  
  return(var_beta)
}

# inequality constraints for optimization
create_constraint_function_SWCRT <- function(K, J) {
  function(x) {
    rho0_E <- x[1]
    rho1_E <- x[2]
    rho0_C <- x[3]
    rho1_C <- x[4]
    rho0_EC <- x[5]
    rho1_EC <- x[6]
    rho2_EC <- x[7]
    
    g <- numeric(14)
    
    # rho0_E >= rho1_E
    g[1] <- -(rho0_E - rho1_E)
    
    # rho0_C >= rho1_C
    g[2] <- -(rho0_C - rho1_C)
    
    # rho0_EC <= min(rho0_E, rho0_C)
    g[3] <- rho0_EC - rho0_E
    g[4] <- rho0_EC - rho0_C
    
    # rho1_EC <= min(rho1_E, rho1_C)
    g[5] <- rho1_EC - rho1_E
    g[6] <- rho1_EC - rho1_C
    
    # rho0_EC >= rho1_EC
    g[7] <- -(rho0_EC - rho1_EC)
    
    # rho2_EC >= rho0_EC
    g[8] <- -(rho2_EC - rho0_EC)
    
    # All 6 eigenvalues > 0 (strictly positive)
    kappa_E <- 1 + (K - 1) * rho0_E - K * rho1_E
    kappa_C <- 1 + (K - 1) * rho0_C - K * rho1_C
    kappa_EC <- rho2_EC + (K - 1) * rho0_EC - K * rho1_EC
    
    xi_1 <- ((K - 1) * (rho0_E - rho0_C) + (J - 1) * K * (rho1_E - rho1_C))^2 + 
      4 * (rho2_EC + (K - 1) * rho0_EC + (J - 1) * K * rho1_EC)^2
    
    lambda_1_plus <- 0.5 * ((2 + (K - 1) * (rho0_E + rho0_C) + (J - 1) * K * (rho1_E + rho1_C)) + sqrt(xi_1))
    lambda_1_minus <- 0.5 * ((2 + (K - 1) * (rho0_E + rho0_C) + (J - 1) * K * (rho1_E + rho1_C)) - sqrt(xi_1))
    lambda_2_plus <- 0.5 * ((kappa_E + kappa_C) + sqrt((kappa_E - kappa_C)^2 + 4 * kappa_EC^2))
    lambda_2_minus <- 0.5 * ((kappa_E + kappa_C) - sqrt((kappa_E - kappa_C)^2 + 4 * kappa_EC^2))
    lambda_3_plus <- 0.5 * ((2 - rho0_E - rho0_C) + sqrt((rho0_E - rho0_C)^2 + 4 * (rho2_EC - rho0_EC)^2))
    lambda_3_minus <- 0.5 * ((2 - rho0_E - rho0_C) - sqrt((rho0_E - rho0_C)^2 + 4 * (rho2_EC - rho0_EC)^2))
    
    g[9] <- -lambda_1_plus + 1e-3
    g[10] <- -lambda_1_minus + 1e-3
    g[11] <- -lambda_2_plus + 1e-3
    g[12] <- -lambda_2_minus + 1e-3
    g[13] <- -lambda_3_plus + 1e-3
    g[14] <- -lambda_3_minus + 1e-3
    
    return(g)
  }
}

# ============================================================================================================
# design matrix construction
# ============================================================================================================

# construct Z matrix for incomplete extended SW design
# L: number of sequences, h: clusters per sequence, J: number of periods
# incomplete: first half clusters miss last 1 periods, second half miss first 2
construct_Z_incomplete <- function(L, h, J) {
  I <- L * h
  
  # base sw design
  Z <- matrix(0, I, J)
  for (s in 1:L) {
    rows <- ((s - 1) * h + 1):(s * h)
    if (s + 1 <= J) Z[rows, (s + 1):J] <- 1
  }
  
  # incomplete: first half miss last 1, second half miss first 2
  Z[1:(I/2), J] <- NA
  Z[(I/2 + 1):I, 1:2] <- NA
  
  return(Z)
}

# ============================================================================================================
# locally optimal design (LOD) function
# ============================================================================================================

compute_LOD_integer_SWCRT <- function(rho0_E, rho1_E, rho0_C, rho1_C, rho0_EC, rho1_EC, rho2_EC, 
                                      trteff, zalpha, J, L, lambda, sigma_E, sigma_C, 
                                      c1, c2, B, max_I, max_K) {
  maxpower <- 0
  lod_result <- NULL
  
  for (h_cur in 1:ceiling(max_I / L)) {
    I <- L * h_cur
    if (I > max_I) break
    
    K <- 2
    Z <- construct_Z_incomplete(L, h_cur, J)
    n_obs <- sum(!is.na(Z))
    actcost <- I * c1 + n_obs * c2 * K
    
    while (actcost <= B && K <= max_K) {
      var <- compute_variance_SWCRT(Z, K, rho0_E, rho1_E, rho0_C, rho1_C, 
                                    rho0_EC, rho1_EC, rho2_EC, lambda, sigma_E, sigma_C)
      
      if (!is.na(var) && var > 0) {
        zbeta <- abs(trteff) / sqrt(var) - zalpha
        power <- pnorm(zbeta)
        
        if (power > maxpower) {
          maxpower <- power
          lod_result <- list(
            rho0_E = rho0_E,
            rho1_E = rho1_E,
            rho0_C = rho0_C,
            rho1_C = rho1_C,
            rho0_EC = rho0_EC,
            rho1_EC = rho1_EC,
            rho2_EC = rho2_EC,
            num_periods = J,
            num_steps = L,
            num_clusters = I,
            clusters_per_step = h_cur,
            cluster_size = K,
            power = power,
            variance = var,
            cost = actcost,
            Z = Z
          )
        }
      }
      
      K <- K + 1
      actcost <- I * c1 + n_obs * c2 * K
    }
  }
  
  return(lod_result)
}

# ============================================================================================================
# maximin design (MMD) functions
# ============================================================================================================

create_mmd_objective_SWCRT <- function(Z, temp_K, J, L, lambda, sigma_E, sigma_C, c1, c2, B, max_K, max_I) {
  
  function(x) {
    rho0_E <- x[1]
    rho1_E <- x[2]
    rho0_C <- x[3]
    rho1_C <- x[4]
    rho0_EC <- x[5]
    rho1_EC <- x[6]
    rho2_EC <- x[7]
    
    var_int <- compute_variance_SWCRT(Z, temp_K, rho0_E, rho1_E, rho0_C, rho1_C,
                                      rho0_EC, rho1_EC, rho2_EC, lambda, sigma_E, sigma_C)
    
    if (is.na(var_int) || var_int <= 0) return(Inf)
    
    I <- nrow(Z)
    n_obs <- sum(!is.na(Z))
    
    obj_fn_inner <- function(K_cont) {
      # K_cont is at least 2 (minimum cluster size)
      if (K_cont < 2) return(NA_real_)
      
      # K_cont does not exceed max_K
      if (K_cont > max_K) return(NA_real_)
      
      # compute decimal I from budget constraint
      I_cont <- B / (c1 + c2 * J * K_cont)
      
      # at least L clusters for stepped-wedge design
      if (I_cont < L) return(NA_real_)
      
      # rho0_E >= rho1_E
      if (rho0_E < rho1_E) return(NA_real_)
      
      # rho0_C >= rho1_C
      if (rho0_C < rho1_C) return(NA_real_)
      
      # rho0_EC <= min(rho0_E, rho0_C)
      if (rho0_EC > rho0_E || rho0_EC > rho0_C) return(NA_real_)
      
      # rho1_EC <= min(rho1_E, rho1_C)
      if (rho1_EC > rho1_E || rho1_EC > rho1_C) return(NA_real_)
      
      # rho0_EC >= rho1_EC
      if (rho0_EC < rho1_EC) return(NA_real_)
      
      # rho2_EC >= rho0_EC
      if (rho2_EC < rho0_EC) return(NA_real_)
      
      # compute kappa values
      kappa_E <- 1 + (K_cont - 1) * rho0_E - K_cont * rho1_E
      kappa_C <- 1 + (K_cont - 1) * rho0_C - K_cont * rho1_C
      kappa_EC <- rho2_EC + (K_cont - 1) * rho0_EC - K_cont * rho1_EC
      
      xi_1 <- ((K_cont - 1) * (rho0_E - rho0_C) + (J - 1) * K_cont * (rho1_E - rho1_C))^2 + 
        4 * (rho2_EC + (K_cont - 1) * rho0_EC + (J - 1) * K_cont * rho1_EC)^2
      
      lambda_1_plus <- 0.5 * ((2 + (K_cont - 1) * (rho0_E + rho0_C) + (J - 1) * K_cont * (rho1_E + rho1_C)) + sqrt(xi_1))
      lambda_1_minus <- 0.5 * ((2 + (K_cont - 1) * (rho0_E + rho0_C) + (J - 1) * K_cont * (rho1_E + rho1_C)) - sqrt(xi_1))
      lambda_2_plus <- 0.5 * ((kappa_E + kappa_C) + sqrt((kappa_E - kappa_C)^2 + 4 * kappa_EC^2))
      lambda_2_minus <- 0.5 * ((kappa_E + kappa_C) - sqrt((kappa_E - kappa_C)^2 + 4 * kappa_EC^2))
      lambda_3_plus <- 0.5 * ((2 - rho0_E - rho0_C) + sqrt((rho0_E - rho0_C)^2 + 4 * (rho2_EC - rho0_EC)^2))
      lambda_3_minus <- 0.5 * ((2 - rho0_E - rho0_C) - sqrt((rho0_E - rho0_C)^2 + 4 * (rho2_EC - rho0_EC)^2))
      
      if (lambda_1_plus <= 0 || lambda_1_minus <= 0 || 
          lambda_2_plus <= 0 || lambda_2_minus <= 0 || 
          lambda_3_plus <= 0 || lambda_3_minus <= 0) return(NA_real_)
      
      var_cont <- compute_variance_SWCRT(Z, K_cont, rho0_E, rho1_E, rho0_C, rho1_C,
                                         rho0_EC, rho1_EC, rho2_EC, lambda, sigma_E, sigma_C)
      
      if (is.na(var_cont) || var_cont <= 0) return(Inf)
      return(var_cont)
    }
    
    K_lower <- max(2, (B - I * c1) / (n_obs * c2 * max_K))
    K_upper <- min(max_K, (B - I * c1) / (n_obs * c2))
    
    if (K_lower > K_upper) return(Inf)
    
    opt_result <- tryCatch({
      optimize(obj_fn_inner, interval = c(K_lower, K_upper), tol = 1e-8)
    }, error = function(e) NULL)
    
    if (is.null(opt_result)) return(Inf)
    
    var_cont_opt <- opt_result$objective
    if (is.na(var_cont_opt) || var_cont_opt <= 0) return(Inf)
    
    RE <- var_cont_opt / var_int
    return(RE)
  }
}

compute_MMD_integer_SWCRT <- function(rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
                                      rho0_EC_min, rho1_EC_min, rho2_EC_min,
                                      rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
                                      rho0_EC_max, rho1_EC_max, rho2_EC_max,
                                      rho0_E_init, rho1_E_init, rho0_C_init, rho1_C_init,
                                      rho0_EC_init, rho1_EC_init, rho2_EC_init,
                                      J, L, lambda, sigma_E, sigma_C,
                                      c1, c2, B, max_I, max_K) {
  reminmax <- 0
  mmd_result <- NULL
  
  lb <- c(rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
          rho0_EC_min, rho1_EC_min, rho2_EC_min)
  ub <- c(rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
          rho0_EC_max, rho1_EC_max, rho2_EC_max)
  
  for (h_cur in 1:ceiling(max_I / L)) {
    I <- L * h_cur
    if (I > max_I) break
    
    Z <- construct_Z_incomplete(L, h_cur, J)
    n_obs <- sum(!is.na(Z))
    
    K <- 2
    actcost <- I * c1 + n_obs * c2 * K
    
    while (actcost <= B && K <= max_K) {
      temp_K <- K
      
      obj_fn <- create_mmd_objective_SWCRT(Z, temp_K, J, L, lambda, sigma_E, sigma_C, c1, c2, B, max_K, max_I)
      eval_g_ineq <- create_constraint_function_SWCRT(K, J)
      
      x0 <- c(rho0_E_init, rho1_E_init, rho0_C_init, rho1_C_init,
              rho0_EC_init, rho1_EC_init, rho2_EC_init)
      
      for (i in 1:length(x0)) {
        if (x0[i] < lb[i] || x0[i] > ub[i]) {
          x0[i] <- (lb[i] + ub[i]) / 2
        }
      }
      
      opts <- list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel = 1e-4,
        maxeval = 1e6
      )
      
      res <- tryCatch({
        nloptr(
          x0 = x0,
          eval_f = obj_fn,
          eval_g_ineq = eval_g_ineq,
          lb = lb,
          ub = ub,
          opts = opts
        )
      }, error = function(e) NULL)
      
      if (!is.null(res) && res$status > 0) {
        remin <- res$objective
        
        if (remin > reminmax && remin > 0 && remin <= 1) {
          reminmax <- remin
          
          var_worst <- compute_variance_SWCRT(Z, K, res$solution[1], res$solution[2], 
                                              res$solution[3], res$solution[4],
                                              res$solution[5], res$solution[6], 
                                              res$solution[7], lambda, sigma_E, sigma_C)
          
          mmd_result <- list(
            rho0_E_worst = res$solution[1],
            rho1_E_worst = res$solution[2],
            rho0_C_worst = res$solution[3],
            rho1_C_worst = res$solution[4],
            rho0_EC_worst = res$solution[5],
            rho1_EC_worst = res$solution[6],
            rho2_EC_worst = res$solution[7],
            num_periods = J,
            num_steps = L,
            num_clusters = I,
            clusters_per_step = h_cur,
            cluster_size = K,
            relative_efficiency = remin,
            variance_worst_case = var_worst,
            cost = actcost,
            optimization_status = res$status,
            budget = B,
            budget_used_pct = (actcost / B) * 100,
            Z = Z
          )
        }
      }
      
      K <- K + 1
      actcost <- I * c1 + n_obs * c2 * K
    }
  }
  
  return(mmd_result)
}

# ============================================================================================================
# main function
# ============================================================================================================

OD_CE_SW_CRT <- function(
    alpha = 0.05,
    trteff = 100,
    lambda = 20000,
    sigma_E = 1,
    sigma_C = 3000,
    rho0_E_min = NULL,
    rho0_E_max = NULL,
    rho1_E_min = NULL,
    rho1_E_max = NULL,
    rho0_C_min = NULL,
    rho0_C_max = NULL,
    rho1_C_min = NULL,
    rho1_C_max = NULL,
    rho0_EC_min = NULL,
    rho0_EC_max = NULL,
    rho1_EC_min = NULL,
    rho1_EC_max = NULL,
    rho2_EC_min = NULL,
    rho2_EC_max = NULL,
    rho0_E_init = 0.06,
    rho1_E_init = 0.03,
    rho0_C_init = 0.06,
    rho1_C_init = 0.03,
    rho0_EC_init = 0.02,
    rho1_EC_init = 0.01,
    rho2_EC_init = 0.5,
    J = 4,
    L = 2,
    c1 = 3000,
    c2 = 250,
    B = 300000,
    max_I = 100,
    max_K = 200
) {
  
  if (is.null(rho0_E_max)) rho0_E_max <- rho0_E_min
  if (is.null(rho1_E_max)) rho1_E_max <- rho1_E_min
  if (is.null(rho0_C_max)) rho0_C_max <- rho0_C_min
  if (is.null(rho1_C_max)) rho1_C_max <- rho1_C_min
  if (is.null(rho0_EC_max)) rho0_EC_max <- rho0_EC_min
  if (is.null(rho1_EC_max)) rho1_EC_max <- rho1_EC_min
  if (is.null(rho2_EC_max)) rho2_EC_max <- rho2_EC_min
  
  if (is.null(rho0_E_min) || is.null(rho1_E_min) || is.null(rho0_C_min) || 
      is.null(rho1_C_min) || is.null(rho0_EC_min) || is.null(rho1_EC_min) || 
      is.null(rho2_EC_min)) {
    stop("All minimum ICC parameters must be provided")
  }
  
  zalpha <- qnorm(1 - alpha/2)
  r <- sigma_E / sigma_C
  
  input_params <- list(
    alpha = alpha,
    trteff = trteff,
    lambda = lambda,
    sigma_E = sigma_E,
    sigma_C = sigma_C,
    r = r,
    rho0_E_range = c(rho0_E_min, rho0_E_max),
    rho1_E_range = c(rho1_E_min, rho1_E_max),
    rho0_C_range = c(rho0_C_min, rho0_C_max),
    rho1_C_range = c(rho1_C_min, rho1_C_max),
    rho0_EC_range = c(rho0_EC_min, rho0_EC_max),
    rho1_EC_range = c(rho1_EC_min, rho1_EC_max),
    rho2_EC_range = c(rho2_EC_min, rho2_EC_max),
    num_periods = J,
    num_steps = L,
    cost_per_cluster = c1,
    cost_per_participant = c2,
    total_budget = B
  )
  
  is_LOD <- (rho0_E_min == rho0_E_max && rho1_E_min == rho1_E_max &&
               rho0_C_min == rho0_C_max && rho1_C_min == rho1_C_max &&
               rho0_EC_min == rho0_EC_max && rho1_EC_min == rho1_EC_max &&
               rho2_EC_min == rho2_EC_max)
  
  if (is_LOD) {
    lod_integer <- compute_LOD_integer_SWCRT(
      rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
      rho0_EC_min, rho1_EC_min, rho2_EC_min,
      trteff, zalpha, J, L, lambda, sigma_E, sigma_C,
      c1, c2, B, max_I, max_K
    )
    
    return(list(
      design_type = "LOD",
      input_parameters = input_params,
      integer_estimates = lod_integer
    ))
    
  } else {
    mmd_integer <- compute_MMD_integer_SWCRT(
      rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
      rho0_EC_min, rho1_EC_min, rho2_EC_min,
      rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
      rho0_EC_max, rho1_EC_max, rho2_EC_max,
      rho0_E_init, rho1_E_init, rho0_C_init, rho1_C_init,
      rho0_EC_init, rho1_EC_init, rho2_EC_init,
      J, L, lambda, sigma_E, sigma_C,
      c1, c2, B, max_I, max_K
    )
    
    return(list(
      design_type = "MMD",
      input_parameters = input_params,
      integer_estimates = mmd_integer
    ))
  }
}