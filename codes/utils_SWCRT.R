############################################################################################################
#                                                                                                          #
#        program:        utils_SWCRT.R                                                                     #
#                                                                                                          #
#        usage:          to propose LOD and MMD in stepped-wedge CRTs for cost-effectiveness              #
#                                                                                                          #
#        assumption:     nested exchangeable correlation structure                                         #
#                        equal cluster size across time periods                                            #
#                        cross-sectional design (different individuals in each cluster-period)             #
#                                                                                                          #
#        required:       base R, nloptr package for optimization                                           #
#                                                                                                          #
#        parameters:                                                                                       #
#          rho0_E:       within-period effect ICC                                                          #
#          rho0_C:       within-period cost ICC                                                            #
#          rho1_E:       between-period effect ICC                                                         #
#          rho1_C:       between-period cost ICC                                                           #
#          rho0_EC:      within-period effect-cost ICC                                                     #
#          rho1_EC:      between-period effect-cost ICC                                                    #
#          rho2_EC:      within-individual effect-cost ICC                                                 #
#          J:            number of time periods                                                            #
#          L:            number of steps (crossover points)                                                #
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

# compute design constants for stepped-wedge design
compute_design_constants_SWCRT <- function(I, J, L) {
  # h represents number of clusters per step
  h <- I / L
  
  # for J >= L + 1
  if (J < L + 1) {
    stop("Number of periods J must be >= L + 1")
  }
  
  # U = sum_i sum_j Z_ij
  U <- (L * (L + 1) * h) / 2 + L * h * (J - L - 1)
  
  # V = sum_i (sum_j Z_ij)^2
  V <- (1/3 * L^3 + 1/2 * L^2 + 1/6 * L) * h + 
    L * h * (J - L - 1)^2 + 
    (J - L - 1) * L * (L + 1) * h
  
  # W = sum_j (sum_i Z_ij)^2
  W <- (1/3 * L^3 + 1/2 * L^2 + 1/6 * L) * h^2 + 
    L^2 * h^2 * (J - L - 1)
  
  return(list(U = U, V = V, W = W))
}

# compute variance of INMB estimator for stepped-wedge design
compute_variance_SWCRT <- function(I, K, rho0_E, rho1_E, rho0_C, rho1_C, 
                                   rho0_EC, rho1_EC, rho2_EC, J, L, lambda, sigma_E, sigma_C) {
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
  
  # get design constants
  design_const <- compute_design_constants_SWCRT(I, J, L)
  U <- design_const$U
  V <- design_const$V
  W <- design_const$W
  
  # compute intermediate terms
  Delta <- (1 / K^2) * sigma_E^2 * sigma_C^2 * (kappa_E * kappa_C - kappa_EC^2)
  
  if (Delta <= 0) return(NA_real_)
  
  Delta_star <- Delta + J * sigma_E^2 * sigma_C^2 * 
    ((1 / K) * (kappa_E * rho1_C + kappa_C * rho1_E - 2 * kappa_EC * rho1_EC) + 
       J * (rho1_E * rho1_C - rho1_EC^2))
  
  if (Delta_star <= 0) return(NA_real_)
  
  varphi <- J * (I * U - W) / Delta + 
    (U^2 - I * V) * (1 / Delta - 1 / Delta_star)
  
  eta <- varphi * J * (I * U - W) + 
    J^2 * (1 / Delta_star) * (U^2 - I * V) * sigma_E^2 * sigma_C^2 * 
    (rho1_E * rho1_C - rho1_EC^2) * 
    (varphi + (1 / Delta_star) * (U^2 - I * V))
  
  if (eta <= 0) return(NA_real_)
  
  # compute variance
  term1 <- (varphi / K) * (lambda^2 * kappa_E * sigma_E / sigma_C - 
                             2 * lambda * kappa_EC + 
                             kappa_C * sigma_C / sigma_E)
  
  term2 <- (J / Delta_star) * (U^2 - I * V) * 
    (lambda^2 * rho1_E * sigma_E / sigma_C - 
       2 * lambda * rho1_EC + 
       rho1_C * sigma_C / sigma_E)
  
  var_beta <- (I * J * sigma_E * sigma_C / eta) * (term1 - term2)
  
  # check if variance is positive
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
# locally optimal design (LOD) function
# ============================================================================================================

# find locally optimal design with integer constraints
compute_LOD_integer_SWCRT <- function(rho0_E, rho1_E, rho0_C, rho1_C, rho0_EC, rho1_EC, rho2_EC, 
                                      trteff, zalpha, J, L, lambda, sigma_E, sigma_C, 
                                      c1, c2, B, max_I, max_K) {
  maxpower <- 0
  lod_result <- NULL
  
  # iterate through possible number of clusters (must be divisible by L)
  for (I in seq(L, max_I, by = L)) {
    K <- 2
    actcost <- I * (c1 + c2 * J * K)
    
    # iterate through cluster sizes
    while (actcost <= B && K <= max_K) {
      var <- compute_variance_SWCRT(I, K, rho0_E, rho1_E, rho0_C, rho1_C, 
                                    rho0_EC, rho1_EC, rho2_EC, J, L, lambda, sigma_E, sigma_C)
      
      if (!is.na(var) && var > 0) {
        # compute power for testing H0: beta_1 = 0
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
            cluster_size = K,
            power = power,
            variance = var,
            cost = actcost
          )
        }
      }
      
      K <- K + 1
      actcost <- I * (c1 + c2 * J * K)
    }
  }
  
  return(lod_result)
}

# ============================================================================================================
# maximin design (MMD) functions
# ============================================================================================================

# objective function for MMD - computes relative efficiency
create_mmd_objective_SWCRT <- function(temp_I, temp_K, J, L, lambda, sigma_E, sigma_C, c1, c2, B, max_K, max_I) {
  
  function(x) {
    # extract parameters (7 ICC parameters for SW-CRT)
    rho0_E <- x[1]   # within-period effect ICC
    rho1_E <- x[2]   # between-period effect ICC  
    rho0_C <- x[3]   # within-period cost ICC
    rho1_C <- x[4]   # between-period cost ICC
    rho0_EC <- x[5]  # within-period effect-cost ICC
    rho1_EC <- x[6]  # between-period effect-cost ICC
    rho2_EC <- x[7]  # within-individual effect-cost ICC
    
    # compute variance for integer design (denominator of RE)
    var_int <- compute_variance_SWCRT(temp_I, temp_K, rho0_E, rho1_E, rho0_C, rho1_C,
                                      rho0_EC, rho1_EC, rho2_EC, J, L, lambda, sigma_E, sigma_C)
    
    if (is.na(var_int) || var_int <= 0) {
      return(NA_real_)
    }
    
    # Inner optimization: find decimal cluster size K that minimizes variance
    # under budget constraint: I_decimal = B / (c1 + c2 * J * K)
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
      
      # compute variance for decimal design
      var_cont <- compute_variance_SWCRT(I_cont, K_cont, rho0_E, rho1_E, rho0_C, rho1_C,
                                         rho0_EC, rho1_EC, rho2_EC, J, L, lambda, sigma_E, sigma_C)
      
      if (is.na(var_cont) || var_cont <= 0) {
        return(NA_real_)  # infeasible decimal design
      }
      
      return(var_cont)
    }
    
    # determine search bounds for K_cont
    # lower bound: max of 2 and budget-driven bound (ensures feasibility when I=max_I)
    K_lower <- max(2, (B/max_I - c1)/(c2 * J))
    
    # upper bound: simply the user-specified maximum
    K_upper <- min(max_K, (B / L - c1) / (c2 * J))
    
    # ensure valid interval
    if (K_lower > K_upper) {
      return(NA_real_)  # no feasible range
    }
    
    # Single optimization call (simple approach)
    opt_result <- tryCatch({
      optimize(obj_fn_inner, interval = c(K_lower, K_upper), tol = 1e-8)
    }, error = function(e) {
      NULL
    })
    
    if (is.null(opt_result)) {
      return(NA_real_)  # optimization failed
    }
    
    var_cont_opt <- opt_result$objective
    
    if (is.na(var_cont_opt) || var_cont_opt <= 0) {
      return(NA_real_)  # no feasible decimal design found
    }
    
    # compute relative efficiency: decimal_variance / integer_variance
    RE <- var_cont_opt / var_int
    
    # we minimize RE to find worst-case correlations
    # (lower RE means integer design is less efficient vs optimal decimal)
    return(RE)
  }
}

# find maximin design with integer constraints
compute_MMD_integer_SWCRT <- function(rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
                                      rho0_EC_min, rho1_EC_min, rho2_EC_min,
                                      rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
                                      rho0_EC_max, rho1_EC_max, rho2_EC_max,
                                      rho0_E_init, rho1_E_init, rho0_C_init, rho1_C_init,
                                      rho0_EC_init, rho1_EC_init, rho2_EC_init,
                                      J, L, lambda, sigma_E, sigma_C,
                                      c1, c2, B, max_I, max_K) {
  reminmax <- 0  # track best (maximum) worst-case RE
  mmd_result <- NULL
  
  # bounds for optimization
  lb <- c(rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
          rho0_EC_min, rho1_EC_min, rho2_EC_min)
  ub <- c(rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
          rho0_EC_max, rho1_EC_max, rho2_EC_max)
  
  # Grid search over integer designs
  # I must be divisible by L (stepped-wedge constraint)
  for (I in seq(L, max_I, by = L)) {
    K <- 2
    actcost <- I * (c1 + c2 * J * K)
    
    while (actcost <= B && K <= max_K) {
      temp_I <- I
      temp_K <- K
      
      # create objective function for this (I, K) pair
      obj_fn <- create_mmd_objective_SWCRT(temp_I, temp_K, J, L, lambda, sigma_E, sigma_C, c1, c2, B, max_K, max_I)
      
      # create constraint function
      eval_g_ineq <- create_constraint_function_SWCRT(K, J)
      
      # initial values - use midpoint if init values are outside bounds
      x0 <- c(rho0_E_init, rho1_E_init, rho0_C_init, rho1_C_init,
              rho0_EC_init, rho1_EC_init, rho2_EC_init)
      
      # bounds
      lb <- c(rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
              rho0_EC_min, rho1_EC_min, rho2_EC_min)
      ub <- c(rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
              rho0_EC_max, rho1_EC_max, rho2_EC_max)
      
      # use midpoint if initial values are outside bounds
      for (i in 1:length(x0)) {
        if (x0[i] < lb[i] || x0[i] > ub[i]) {
          x0[i] <- (lb[i] + ub[i]) / 2
        }
      }
      
      # optimization settings
      opts <- list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel = 1e-4,
        maxeval = 1e6
      )
      
      # solve constrained optimization to find worst-case ICC parameters
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
        
        # we want the design with maximum worst-case RE
        # (best design = highest minimum efficiency)
        if (remin > reminmax && remin > 0 && remin <= 1) {
          reminmax <- remin
          
          # compute variance at worst-case parameters
          var_worst <- compute_variance_SWCRT(I, K, res$solution[1], res$solution[2], 
                                              res$solution[3], res$solution[4],
                                              res$solution[5], res$solution[6], 
                                              res$solution[7], J, L, lambda, sigma_E, sigma_C)
          
          mmd_result <- list(
            # worst-case correlation parameters
            rho0_E_worst = res$solution[1],
            rho1_E_worst = res$solution[2],
            rho0_C_worst = res$solution[3],
            rho1_C_worst = res$solution[4],
            rho0_EC_worst = res$solution[5],
            rho1_EC_worst = res$solution[6],
            rho2_EC_worst = res$solution[7],
            # design parameters
            num_periods = J,
            num_steps = L,
            num_clusters = I,
            cluster_size = K,
            # performance metrics
            relative_efficiency = remin,
            variance_worst_case = var_worst,
            cost = actcost,
            optimization_status = res$status,
            # for diagnostics
            budget = B,
            budget_used_pct = (actcost / B) * 100
          )
        }
      }
      
      K <- K + 1
      actcost <- I * (c1 + c2 * J * K)
    }
  }
  
  return(mmd_result)
}

# ============================================================================================================
# main function for optimal design in cost-effectiveness SW-CRT
# ============================================================================================================

OD_CE_SW_CRT <- function(
    alpha = 0.05,           # type I error
    trteff = 3000,          # treatment effect (INMB)
    lambda = 20000,         # ceiling ratio (willingness to pay per unit effectiveness)
    sigma_E = 1,            # standard deviation of effectiveness outcome
    sigma_C = 3000,         # standard deviation of cost outcome
    # effect ICC bounds
    rho0_E_min = NULL,      # within-period effect ICC (min or exact)
    rho0_E_max = NULL,      # within-period effect ICC (max or missing)
    rho1_E_min = NULL,      # between-period effect ICC (min or exact)
    rho1_E_max = NULL,      # between-period effect ICC (max or missing)
    # cost ICC bounds
    rho0_C_min = NULL,      # within-period cost ICC (min or exact)
    rho0_C_max = NULL,      # within-period cost ICC (max or missing)
    rho1_C_min = NULL,      # between-period cost ICC (min or exact)
    rho1_C_max = NULL,      # between-period cost ICC (max or missing)
    # effect-cost ICC bounds
    rho0_EC_min = NULL,     # within-period effect-cost ICC (min or exact)
    rho0_EC_max = NULL,     # within-period effect-cost ICC (max or missing)
    rho1_EC_min = NULL,     # between-period effect-cost ICC (min or exact)
    rho1_EC_max = NULL,     # between-period effect-cost ICC (max or missing)
    rho2_EC_min = NULL,     # within-individual effect-cost ICC (min or exact)
    rho2_EC_max = NULL,     # within-individual effect-cost ICC (max or missing)
    # initial values for MMD
    rho0_E_init = 0.06,
    rho1_E_init = 0.03,
    rho0_C_init = 0.06,
    rho1_C_init = 0.03,
    rho0_EC_init = 0.02,
    rho1_EC_init = 0.01,
    rho2_EC_init = 0.5,
    # design parameters
    J = 4,                  # number of time periods
    L = 2,                  # number of steps (crossover points)
    c1 = 3000,              # cost per cluster
    c2 = 250,               # cost per participant per period
    B = 300000,             # total budget
    max_I = 100,            # maximum number of clusters
    max_K = 200             # maximum cluster size
) {
  
  # ==========================================================================================================
  # parameter processing
  # ==========================================================================================================
  
  # handle missing max values - set to min if not provided
  if (is.null(rho0_E_max)) rho0_E_max <- rho0_E_min
  if (is.null(rho1_E_max)) rho1_E_max <- rho1_E_min
  if (is.null(rho0_C_max)) rho0_C_max <- rho0_C_min
  if (is.null(rho1_C_max)) rho1_C_max <- rho1_C_min
  if (is.null(rho0_EC_max)) rho0_EC_max <- rho0_EC_min
  if (is.null(rho1_EC_max)) rho1_EC_max <- rho1_EC_min
  if (is.null(rho2_EC_max)) rho2_EC_max <- rho2_EC_min
  
  # check if all min parameters are provided
  if (is.null(rho0_E_min) || is.null(rho1_E_min) || is.null(rho0_C_min) || 
      is.null(rho1_C_min) || is.null(rho0_EC_min) || is.null(rho1_EC_min) || 
      is.null(rho2_EC_min)) {
    stop("All minimum ICC parameters must be provided")
  }
  
  # check J >= L + 1
  if (J < L + 1) {
    stop("Number of periods J must be >= L + 1 for stepped-wedge design")
  }
  
  # calculate derived parameters
  zalpha <- qnorm(1 - alpha/2)
  r <- sigma_E / sigma_C
  
  # store input parameters
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
  
  # ==========================================================================================================
  # determine design type (LOD or MMD)
  # ==========================================================================================================
  
  is_LOD <- (rho0_E_min == rho0_E_max && rho1_E_min == rho1_E_max &&
               rho0_C_min == rho0_C_max && rho1_C_min == rho1_C_max &&
               rho0_EC_min == rho0_EC_max && rho1_EC_min == rho1_EC_max &&
               rho2_EC_min == rho2_EC_max)
  
  # ==========================================================================================================
  # compute optimal design
  # ==========================================================================================================
  
  if (is_LOD) {
    # compute locally optimal design
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
    # compute maximin design
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