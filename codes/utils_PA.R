############################################################################################################
#                                                                                                          #
#        program:        utils_PA.R                                                                       #
#                                                                                                          #
#        usage:          to propose LOD and MMD in parallel-arm LCRTs for cost-effectiveness              #
#                                                                                                          #
#        assumption:     nested exchangeable correlation structure                                         #
#                        equal cluster size across time periods                                            #
#                        longitudinal design (same individuals in each cluster across periods)             #
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

# compute theta parameter for optimal design (parallel-arm specific)
compute_vartheta_PA <- function(rho0_E, rho1_E, rho0_C, rho1_C, rho0_EC, rho1_EC, rho2_EC, lambda, r, J) {
  # numerator: related to within-cluster variance for parallel-arm
  num <- (1 + J * rho1_E - rho1_E) + 
    2 * (rho1_EC - J * rho1_EC - rho2_EC) / (lambda * r) + 
    (1 + J * rho1_C - rho1_C) / (lambda^2 * r^2)
  
  # denominator: related to between-cluster variance for parallel-arm
  denom <- (J * rho1_E + rho0_E - rho1_E) + 
    2 * (-J * rho1_EC - rho0_EC + rho1_EC) / (lambda * r) + 
    (J * rho1_C + rho0_C - rho1_C) / (lambda^2 * r^2)
  
  # compute vartheta
  if (denom == 0) {
    return(NA_real_)
  }
  
  vartheta <- num / denom - 1
  
  # vartheta must be positive
  if (vartheta <= 0) {
    return(NA_real_)
  }
  
  return(vartheta)
}

# compute variance of INMB estimator for parallel-arm LCRT
compute_variance_PA <- function(I, K, rho0_E, rho1_E, rho0_C, rho1_C, 
                                rho0_EC, rho1_EC, rho2_EC, J, lambda, sigma_E, sigma_C, prop) {
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
  
  # variance formula for parallel-arm LCRT
  var_beta <- ((kappa_C * sigma_C^2 - 2 * lambda * kappa_EC * sigma_C * sigma_E + 
                  lambda^2 * kappa_E * sigma_E^2) + 
                 J * K * (rho1_C * sigma_C^2 - 2 * lambda * rho1_EC * sigma_C * sigma_E + 
                            lambda^2 * rho1_E * sigma_E^2)) / (I * J * K * prop * (1 - prop))
  
  # check if variance is positive
  if (var_beta <= 0) return(NA_real_)
  
  return(var_beta)
}

# inequality constraints
create_constraint_function_PA <- function(K, J) {
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
# locally optimal design (LOD) functions
# ============================================================================================================

# find locally optimal design with integer constraints
compute_LOD_integer_PA <- function(rho0_E, rho1_E, rho0_C, rho1_C, rho0_EC, rho1_EC, rho2_EC, 
                                   trteff, zalpha, J, lambda, sigma_E, sigma_C, prop, 
                                   c1, c2, B, max_I, max_K, piden) {
  maxpower <- 0
  lod_result <- NULL
  
  # iterate through possible number of clusters
  for (I in seq(2, max_I, by = piden)) {
    K <- 2
    actcost <- I * (c1 + c2 * J * K)
    
    # iterate through cluster sizes
    while (actcost <= B && K <= max_K) {
      var <- compute_variance_PA(I, K, rho0_E, rho1_E, rho0_C, rho1_C, 
                                 rho0_EC, rho1_EC, rho2_EC, J, lambda, sigma_E, sigma_C, prop)
      
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
            proportion = prop,
            num_periods = J,
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

# compute locally optimal design with decimal estimates
compute_LOD_decimal_PA <- function(rho0_E, rho1_E, rho0_C, rho1_C, rho0_EC, rho1_EC, rho2_EC,
                                   trteff, zalpha, J, lambda, sigma_E, sigma_C, prop,
                                   c1, c2, B, r) {
  # compute vartheta
  vartheta <- compute_vartheta_PA(rho0_E, rho1_E, rho0_C, rho1_C, rho0_EC, rho1_EC, rho2_EC, lambda, r, J)
  
  if (is.na(vartheta)) {
    return(NULL)
  }
  
  # decimal estimates
  K_dec <- sqrt(c1 * vartheta / (c2 * J))
  I_dec <- B / (c1 + sqrt(vartheta * c1 * c2 * J))
  cost_dec <- I_dec * (c1 + c2 * J * K_dec)
  var_dec <- compute_variance_PA(I_dec, K_dec, rho0_E, rho1_E, rho0_C, rho1_C,
                                 rho0_EC, rho1_EC, rho2_EC, J, lambda, sigma_E, sigma_C, prop)
  
  if (is.na(var_dec)) {
    return(NULL)
  }
  
  zbeta_dec <- abs(trteff) / sqrt(var_dec) - zalpha
  power_dec <- pnorm(zbeta_dec)
  
  return(list(
    rho0_E = rho0_E,
    rho1_E = rho1_E,
    rho0_C = rho0_C,
    rho1_C = rho1_C,
    rho0_EC = rho0_EC,
    rho1_EC = rho1_EC,
    rho2_EC = rho2_EC,
    proportion = prop,
    num_periods = J,
    num_clusters = I_dec,
    cluster_size = K_dec,
    power = power_dec,
    variance = var_dec,
    cost = cost_dec,
    vartheta = vartheta
  ))
}

# ============================================================================================================
# maximin design (MMD) functions
# ============================================================================================================

# compute g(theta) for parallel-arm LCRT
compute_g_theta_PA <- function(vartheta, c1, c2, J) {
  g_theta <- J * c2 * (1 + sqrt(c1 / (c2 * J * vartheta))) * (vartheta + sqrt(c1 * vartheta / (c2 * J)))
  return(g_theta)
}

# objective function for MMD (cost-effectiveness version)
create_mmd_objective_PA <- function(temp_I, temp_K, J, lambda, sigma_E, sigma_C, prop, c1, c2, B) {
  r <- sigma_E / sigma_C
  
  function(x) {
    # extract parameters
    rho0_E <- x[1]
    rho1_E <- x[2]
    rho0_C <- x[3]
    rho1_C <- x[4]
    rho0_EC <- x[5]
    rho1_EC <- x[6]
    rho2_EC <- x[7]
    
    # compute vartheta
    vartheta <- compute_vartheta_PA(rho0_E, rho1_E, rho0_C, rho1_C, rho0_EC, rho1_EC, rho2_EC, lambda, r, J)
    
    if (is.na(vartheta)) {
      return(1e20)
    }
    
    # compute continuous optimal design
    K_dec <- sqrt(c1 * vartheta / (c2 * J))
    I_dec <- B / (c1 + sqrt(vartheta * c1 * c2 * J))
    
    # compute variances
    var_dec <- compute_variance_PA(I_dec, K_dec, rho0_E, rho1_E, rho0_C, rho1_C,
                                   rho0_EC, rho1_EC, rho2_EC, J, lambda, sigma_E, sigma_C, prop)
    var_int <- compute_variance_PA(temp_I, temp_K, rho0_E, rho1_E, rho0_C, rho1_C,
                                   rho0_EC, rho1_EC, rho2_EC, J, lambda, sigma_E, sigma_C, prop)
    
    if (is.na(var_dec) || is.na(var_int) || var_dec <= 0 || var_int <= 0) {
      return(1e20)
    }
    
    # relative efficiency
    RE <- var_dec / var_int
    return(RE)
  }
}

# compute maximin design with decimal estimates
compute_MMD_decimal_PA <- function(rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min, 
                                   rho0_EC_min, rho1_EC_min, rho2_EC_min,
                                   rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
                                   rho0_EC_max, rho1_EC_max, rho2_EC_max,
                                   J, lambda, sigma_E, sigma_C, prop, c1, c2, B, r) {
  # compute theoretical bounds for theta
  vartheta_min <- compute_vartheta_PA(rho0_E_max, rho1_E_min, rho0_C_max, rho1_C_min,
                                      rho0_EC_max, rho1_EC_min, rho2_EC_max, lambda, r, J)
  
  vartheta_max <- compute_vartheta_PA(rho0_E_min, rho1_E_max, rho0_C_min, rho1_C_max,
                                      rho0_EC_min, rho1_EC_max, rho2_EC_min, lambda, r, J)
  
  if (rho0_E_min > rho1_E_min && rho0_C_min > rho1_C_min && rho0_EC_min > rho1_EC_min &&
      rho0_EC_min <= rho0_E_min && rho0_EC_min <= rho0_C_min && 
      rho1_EC_min <= rho1_E_min && rho1_EC_min <= rho1_C_min &&
      !is.na(vartheta_min) && !is.na(vartheta_max) && vartheta_min > 0 && vartheta_max > 0) {
    
    # decimal estimates for parallel-arm with modified g_theta
    g_theta_min <- compute_g_theta_PA(vartheta_min, c1, c2, J)
    g_theta_max <- compute_g_theta_PA(vartheta_max, c1, c2, J)
    
    K_dec_mmd <- (g_theta_min * vartheta_max - vartheta_min * g_theta_max) / (g_theta_max - g_theta_min)
    I_dec_mmd <- B / (c1 + c2 * J * K_dec_mmd)
    cost_dec_mmd <- I_dec_mmd * (c1 + c2 * J * K_dec_mmd)
    
    # compute relative efficiency at decimal design
    re_dec <- g_theta_min / (vartheta_min + K_dec_mmd) * K_dec_mmd / (c1 + c2 * J * K_dec_mmd)
    
    return(list(
      proportion = prop,
      num_periods = J,
      num_clusters = I_dec_mmd,
      cluster_size = K_dec_mmd,
      relative_efficiency = re_dec,
      cost = cost_dec_mmd,
      vartheta_min = vartheta_min,
      vartheta_max = vartheta_max
    ))
  }
  
  return(NULL)
}

# find maximin design with integer constraints
compute_MMD_integer_PA <- function(rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
                                   rho0_EC_min, rho1_EC_min, rho2_EC_min,
                                   rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
                                   rho0_EC_max, rho1_EC_max, rho2_EC_max,
                                   rho0_E_init, rho1_E_init, rho0_C_init, rho1_C_init,
                                   rho0_EC_init, rho1_EC_init, rho2_EC_init,
                                   J, lambda, sigma_E, sigma_C, prop,
                                   c1, c2, B, max_I, max_K, piden) {
  reminmax <- 0
  mmd_result <- NULL
  
  # iterate through integer designs
  for (I in seq(2, max_I, by = piden)) {
    K <- 2
    actcost <- I * (c1 + c2 * J * K)
    
    while (actcost <= B && K <= max_K) {
      temp_I <- I
      temp_K <- K
      
      # create objective function
      obj_fn <- create_mmd_objective_PA(temp_I, temp_K, J, lambda, sigma_E, sigma_C, prop, c1, c2, B)
      
      # create constraint function
      eval_g_ineq <- create_constraint_function_PA(K, J)
      
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
      
      # solve constrained optimization
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
          mmd_result <- list(
            rho0_E = res$solution[1],
            rho1_E = res$solution[2],
            rho0_C = res$solution[3],
            rho1_C = res$solution[4],
            rho0_EC = res$solution[5],
            rho1_EC = res$solution[6],
            rho2_EC = res$solution[7],
            proportion = prop,
            num_periods = J,
            num_clusters = I,
            cluster_size = K,
            relative_efficiency = remin,
            cost = actcost,
            optimization_status = res$status
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
# main function for optimal design in cost-effectiveness parallel-arm LCRT trials
# ============================================================================================================

OD_CE_PA_LCRT <- function(
    alpha = 0.05,           # type I error
    trteff = 4000,          # treatment effect (INMB)
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
    pinum = 1,              # fraction numerator of prop
    piden = 2,              # fraction denominator of prop
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
  
  # calculate derived parameters
  prop <- pinum / piden
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
    proportion = prop,
    num_periods = J,
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
    lod_decimal <- compute_LOD_decimal_PA(
      rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
      rho0_EC_min, rho1_EC_min, rho2_EC_min,
      trteff, zalpha, J, lambda, sigma_E, sigma_C, prop,
      c1, c2, B, r
    )
    
    lod_integer <- compute_LOD_integer_PA(
      rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
      rho0_EC_min, rho1_EC_min, rho2_EC_min,
      trteff, zalpha, J, lambda, sigma_E, sigma_C, prop,
      c1, c2, B, max_I, max_K, piden
    )
    
    return(list(
      design_type = "LOD",
      input_parameters = input_params,
      decimal_estimates = lod_decimal,
      integer_estimates = lod_integer
    ))
    
  } else {
    # compute maximin design
    mmd_decimal <- compute_MMD_decimal_PA(
      rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
      rho0_EC_min, rho1_EC_min, rho2_EC_min,
      rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
      rho0_EC_max, rho1_EC_max, rho2_EC_max,
      J, lambda, sigma_E, sigma_C, prop, c1, c2, B, r
    )
    
    mmd_integer <- compute_MMD_integer_PA(
      rho0_E_min, rho1_E_min, rho0_C_min, rho1_C_min,
      rho0_EC_min, rho1_EC_min, rho2_EC_min,
      rho0_E_max, rho1_E_max, rho0_C_max, rho1_C_max,
      rho0_EC_max, rho1_EC_max, rho2_EC_max,
      rho0_E_init, rho1_E_init, rho0_C_init, rho1_C_init,
      rho0_EC_init, rho1_EC_init, rho2_EC_init,
      J, lambda, sigma_E, sigma_C, prop,
      c1, c2, B, max_I, max_K, piden
    )
    
    return(list(
      design_type = "MMD",
      input_parameters = input_params,
      decimal_estimates = mmd_decimal,
      integer_estimates = mmd_integer
    ))
  }
}