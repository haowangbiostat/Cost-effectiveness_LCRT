source("../codes/utils_CRXO.R")
library(xtable)

# ============================================================================================================
# fixed parameters
# ============================================================================================================

B <- 300000        # total budget
c1 <- 3000         # cluster-level cost
c2 <- 250          # participant-level cost
lambda <- 20000    # ceiling ratio
sigma_E <- 1       # effectiveness sd
sigma_C <- 3000    # cost sd
trteff <- 4000     # average INMB (0.2 * 20000)
alpha <- 0.05      # significance level
pinum <- 1         # pi numerator
piden <- 2         # pi denominator (pi = 0.5)
max_I <- 100       # max clusters
max_K <- 200       # max cluster size

# ============================================================================================================
# design parameters
# ============================================================================================================

J_values <- c(2, 4, 6)  # periods

# ============================================================================================================
# generate parameter combinations for MMD (ICC ranges)
# ============================================================================================================

param_sets <- list(
  # (rho0_E: 0.05-0.10, rho1_E: 0.025-0.040)
  list(rho0_E_min = 0.05, rho0_E_max = 0.10, 
       rho1_E_min = 0.025, rho1_E_max = 0.040),
  
  # (rho0_E: 0.05-0.10, rho1_E: 0.025-0.045)
  list(rho0_E_min = 0.05, rho0_E_max = 0.10,
       rho1_E_min = 0.025, rho1_E_max = 0.045),
  
  # (rho0_E: 0.05-0.20, rho1_E: 0.025-0.040)
  list(rho0_E_min = 0.05, rho0_E_max = 0.20,
       rho1_E_min = 0.025, rho1_E_max = 0.040),
  
  # (rho0_E: 0.05-0.20, rho1_E: 0.025-0.045)
  list(rho0_E_min = 0.05, rho0_E_max = 0.20,
       rho1_E_min = 0.025, rho1_E_max = 0.045),
  
  # (rho0_E: 0.10-0.20, rho1_E: 0.050-0.080)
  list(rho0_E_min = 0.10, rho0_E_max = 0.20,
       rho1_E_min = 0.050, rho1_E_max = 0.080),
  
  # (rho0_E: 0.10-0.20, rho1_E: 0.050-0.090)
  list(rho0_E_min = 0.10, rho0_E_max = 0.20,
       rho1_E_min = 0.050, rho1_E_max = 0.090)
)

# ============================================================================================================
# run MMD for all combinations
# ============================================================================================================

results <- data.frame()

for (i in 1:length(param_sets)) {
  params <- param_sets[[i]]
  
  # initialize row with parameters
  row_result <- data.frame(
    rho0_E_min = params$rho0_E_min,
    rho0_E_max = params$rho0_E_max,
    rho1_E_min = params$rho1_E_min,
    rho1_E_max = params$rho1_E_max,
    stringsAsFactors = FALSE
  )
  
  # run for each period value
  for (J in J_values) {
    result <- OD_CE_CRXO_MPD(
      alpha = alpha,
      trteff = trteff,
      lambda = lambda,
      sigma_E = sigma_E,
      sigma_C = sigma_C,
      # effect ICC ranges
      rho0_E_min = params$rho0_E_min,
      rho0_E_max = params$rho0_E_max,
      rho1_E_min = params$rho1_E_min,
      rho1_E_max = params$rho1_E_max,
      # cost ICCs - proportional to effect ICCs
      rho0_C_min = params$rho0_E_min * 0.8,
      rho0_C_max = params$rho0_E_max * 0.8,
      rho1_C_min = params$rho1_E_min * 0.8,
      rho1_C_max = params$rho1_E_max * 0.8,
      # effect-cost ICCs - set to small values
      rho0_EC_min = 0.01,
      rho0_EC_max = 0.02,
      rho1_EC_min = 0.005,
      rho1_EC_max = 0.01,
      rho2_EC_min = 0.5,
      rho2_EC_max = 0.8,
      # initial values for optimization
      rho0_E_init = (params$rho0_E_min + params$rho0_E_max) / 2,
      rho1_E_init = (params$rho1_E_min + params$rho1_E_max) / 2,
      rho0_C_init = (params$rho0_E_min + params$rho0_E_max) / 2 * 0.8,
      rho1_C_init = (params$rho1_E_min + params$rho1_E_max) / 2 * 0.8,
      rho0_EC_init = 0.015,
      rho1_EC_init = 0.0075,
      rho2_EC_init = 0.65,
      # design parameters
      J = J,
      pinum = pinum,
      piden = piden,
      c1 = c1,
      c2 = c2,
      B = B,
      max_I = max_I,
      max_K = max_K
    )
    
    # extract integer estimates from MMD results
    if (!is.null(result$integer_estimates)) {
      row_result[[paste0("I_J", J)]] <- sprintf("%.0f", result$integer_estimates$num_clusters)
      row_result[[paste0("K_J", J)]] <- sprintf("%.0f", result$integer_estimates$cluster_size)
      row_result[[paste0("RE_J", J)]] <- sprintf("%.3f", result$integer_estimates$relative_efficiency)
    } else {
      row_result[[paste0("I_J", J)]] <- NA
      row_result[[paste0("K_J", J)]] <- NA
      row_result[[paste0("RE_J", J)]] <- NA
    }
  }
  
  results <- rbind(results, row_result)
}

# ============================================================================================================
# format table for latex output
# ============================================================================================================

display_table <- data.frame(
  "$(\\rho_{0\\min}^E, \\rho_{0\\max}^E)$" = sprintf("(%.2f, %.2f)", results$rho0_E_min, results$rho0_E_max),
  "$(\\rho_{1\\min}^E, \\rho_{1\\max}^E)$" = sprintf("(%.3f, %.3f)", results$rho1_E_min, results$rho1_E_max),
  "$I_{\\text{MMD}}$" = results$I_J2,
  "$K_{\\text{MMD}}$" = results$K_J2,
  "RE" = results$RE_J2,
  "$I_{\\text{MMD}}$.1" = results$I_J4,
  "$K_{\\text{MMD}}$.1" = results$K_J4,
  "RE.1" = results$RE_J4,
  "$I_{\\text{MMD}}$.2" = results$I_J6,
  "$K_{\\text{MMD}}$.2" = results$K_J6,
  "RE.2" = results$RE_J6,
  check.names = FALSE
)

# rename columns properly
names(display_table) <- c(
  "$(\\rho_{0\\min}^E, \\rho_{0\\max}^E)$",
  "$(\\rho_{1\\min}^E, \\rho_{1\\max}^E)$",
  "$I_{\\text{MMD}}$", "$K_{\\text{MMD}}$", "RE",
  "$I_{\\text{MMD}}$", "$K_{\\text{MMD}}$", "RE",
  "$I_{\\text{MMD}}$", "$K_{\\text{MMD}}$", "RE"
)

# ============================================================================================================
# generate latex table
# ============================================================================================================

xt <- xtable(display_table,
             caption = "Maximin optimal designs (MMDs) for PA trials with $\\pi = 0.5$",
             label = "tab:PA_MMD_0.5",
             align = c("l", "l", "l", rep("c", 9)))

# add column groups header
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- "& & \\multicolumn{3}{c}{$J = 2$} & \\multicolumn{3}{c}{$J = 4$} & \\multicolumn{3}{c}{$J = 6$} \\\\\n \\cmidrule(lr){3-5} \\cmidrule(lr){6-8} \\cmidrule(lr){9-11}\n"

# print latex code
print(xt, 
      include.rownames = FALSE,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      hline.after = c(-1, 0, nrow(display_table)),
      add.to.row = addtorow,
      comment = FALSE)