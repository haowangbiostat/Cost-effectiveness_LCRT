source("../codes/utils_SWCRT.R")
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
trteff <- 4000     # average INMB
alpha <- 0.05      # significance level
max_I <- 100       # max clusters
max_K <- 200       # max cluster size

# ============================================================================================================
# design parameters
# ============================================================================================================

L_values <- c(3, 5, 7)                   # number of steps
J_range <- 2:9                           # range of periods to search for optimal J
J_fixed <- 9                             # fixed J for comparison

# parameter combinations based on the example table
param_sets <- list(
  list(rho0_E = 0.05, rho1_E = 0.025, rho0_C = 0.05, rho1_C = 0.025, rho0_EC = 0.020, rho1_EC = 0.010, rho2_EC = 0.5),
  list(rho0_E = 0.05, rho1_E = 0.040, rho0_C = 0.05, rho1_C = 0.040, rho0_EC = 0.020, rho1_EC = 0.016, rho2_EC = 0.5),
  list(rho0_E = 0.10, rho1_E = 0.050, rho0_C = 0.10, rho1_C = 0.050, rho0_EC = 0.040, rho1_EC = 0.020, rho2_EC = 0.5),
  list(rho0_E = 0.10, rho1_E = 0.080, rho0_C = 0.10, rho1_C = 0.080, rho0_EC = 0.040, rho1_EC = 0.032, rho2_EC = 0.5),
  list(rho0_E = 0.20, rho1_E = 0.100, rho0_C = 0.20, rho1_C = 0.100, rho0_EC = 0.080, rho1_EC = 0.040, rho2_EC = 0.5),
  list(rho0_E = 0.20, rho1_E = 0.160, rho0_C = 0.20, rho1_C = 0.160, rho0_EC = 0.080, rho1_EC = 0.064, rho2_EC = 0.5)
)

# ============================================================================================================
# run LOD for all combinations
# ============================================================================================================

results <- data.frame()

for (i in 1:length(param_sets)) {
  params <- param_sets[[i]]
  
  # two rows per parameter set: one for optimal J, one for fixed J = 9
  for (j_type in c("optimal", "fixed")) {
    row_result <- data.frame(
      rho_E = sprintf("(%.2f, %.3f)", params$rho0_E, params$rho1_E),
      J_type = ifelse(j_type == "optimal", "Undecided", as.character(J_fixed)),
      stringsAsFactors = FALSE
    )
    
    # for each L value
    for (L in L_values) {
      if (j_type == "optimal") {
        # find optimal J by searching through J_range
        best_power <- -Inf
        best_J <- NA
        best_I <- NA
        best_K <- NA
        
        for (J in J_range) {
          # skip if J < L + 1
          if (J < L + 1) next
          
          result <- OD_CE_SW_CRT(
            alpha = alpha,
            trteff = trteff,
            lambda = lambda,
            sigma_E = sigma_E,
            sigma_C = sigma_C,
            rho0_E_min = params$rho0_E,
            rho0_E_max = params$rho0_E,
            rho1_E_min = params$rho1_E,
            rho1_E_max = params$rho1_E,
            rho0_C_min = params$rho0_C,
            rho0_C_max = params$rho0_C,
            rho1_C_min = params$rho1_C,
            rho1_C_max = params$rho1_C,
            rho0_EC_min = params$rho0_EC,
            rho0_EC_max = params$rho0_EC,
            rho1_EC_min = params$rho1_EC,
            rho1_EC_max = params$rho1_EC,
            rho2_EC_min = params$rho2_EC,
            rho2_EC_max = params$rho2_EC,
            J = J,
            L = L,
            c1 = c1,
            c2 = c2,
            B = B,
            max_I = max_I,
            max_K = max_K
          )
          
          if (!is.null(result$integer_estimates)) {
            current_power <- result$integer_estimates$power
            if (current_power > best_power) {
              best_power <- current_power
              best_J <- J
              best_I <- result$integer_estimates$num_clusters
              best_K <- result$integer_estimates$cluster_size
            }
          }
        }
        
        row_result[[paste0("J_L", L)]] <- best_J
        row_result[[paste0("I_L", L)]] <- best_I
        row_result[[paste0("K_L", L)]] <- best_K
        row_result[[paste0("Power_L", L)]] <- ifelse(is.finite(best_power), sprintf("%.3f", best_power), NA)
        
      } else {
        # use fixed J = 9
        J <- J_fixed
        
        # skip if J < L + 1
        if (J < L + 1) {
          row_result[[paste0("J_L", L)]] <- NA
          row_result[[paste0("I_L", L)]] <- NA
          row_result[[paste0("K_L", L)]] <- NA
          row_result[[paste0("Power_L", L)]] <- NA
          next
        }
        
        result <- OD_CE_SW_CRT(
          alpha = alpha,
          trteff = trteff,
          lambda = lambda,
          sigma_E = sigma_E,
          sigma_C = sigma_C,
          rho0_E_min = params$rho0_E,
          rho0_E_max = params$rho0_E,
          rho1_E_min = params$rho1_E,
          rho1_E_max = params$rho1_E,
          rho0_C_min = params$rho0_C,
          rho0_C_max = params$rho0_C,
          rho1_C_min = params$rho1_C,
          rho1_C_max = params$rho1_C,
          rho0_EC_min = params$rho0_EC,
          rho0_EC_max = params$rho0_EC,
          rho1_EC_min = params$rho1_EC,
          rho1_EC_max = params$rho1_EC,
          rho2_EC_min = params$rho2_EC,
          rho2_EC_max = params$rho2_EC,
          J = J,
          L = L,
          c1 = c1,
          c2 = c2,
          B = B,
          max_I = max_I,
          max_K = max_K
        )
        
        if (!is.null(result$integer_estimates)) {
          row_result[[paste0("J_L", L)]] <- NA  # J is fixed, so we don't show it
          row_result[[paste0("I_L", L)]] <- sprintf("%.0f", result$integer_estimates$num_clusters)
          row_result[[paste0("K_L", L)]] <- sprintf("%.0f", result$integer_estimates$cluster_size)
          row_result[[paste0("Power_L", L)]] <- sprintf("%.3f", result$integer_estimates$power)
        } else {
          row_result[[paste0("J_L", L)]] <- NA
          row_result[[paste0("I_L", L)]] <- NA
          row_result[[paste0("K_L", L)]] <- NA
          row_result[[paste0("Power_L", L)]] <- NA
        }
      }
    }
    
    results <- rbind(results, row_result)
  }
}

# ============================================================================================================
# format table for latex output
# ============================================================================================================

display_table <- data.frame(
  "$(\\rho_0^E, \\rho_1^E)$" = results$rho_E,
  "$J$" = results$J_type,
  "$J_{\\text{LOD}}$" = results$J_L3,
  "$I_{\\text{LOD}}$" = results$I_L3,
  "$K_{\\text{LOD}}$" = results$K_L3,
  "Power" = results$Power_L3,
  "$J_{\\text{LOD}}$.1" = results$J_L5,
  "$I_{\\text{LOD}}$.1" = results$I_L5,
  "$K_{\\text{LOD}}$.1" = results$K_L5,
  "Power.1" = results$Power_L5,
  "$J_{\\text{LOD}}$.2" = results$J_L7,
  "$I_{\\text{LOD}}$.2" = results$I_L7,
  "$K_{\\text{LOD}}$.2" = results$K_L7,
  "Power.2" = results$Power_L7,
  check.names = FALSE
)

# rename columns properly
names(display_table) <- c(
  "$(\\rho_0^E, \\rho_1^E)$",
  "$J$",
  "$J_{\\text{LOD}}$", "$I_{\\text{LOD}}$", "$K_{\\text{LOD}}$", "Power",
  "$J_{\\text{LOD}}$", "$I_{\\text{LOD}}$", "$K_{\\text{LOD}}$", "Power",
  "$J_{\\text{LOD}}$", "$I_{\\text{LOD}}$", "$K_{\\text{LOD}}$", "Power"
)

# ============================================================================================================
# generate latex table
# ============================================================================================================

xt <- xtable(display_table,
             caption = "Local optimal designs (LODs) for SW-CRTs with optimal and fixed number of periods",
             label = "tab:SWCRT_LOD_optimal_J",
             align = c("l", "l", "c", rep("c", 12)))

# add column groups header
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- "& & \\multicolumn{4}{c}{$L = 3$} & \\multicolumn{4}{c}{$L = 5$} & \\multicolumn{4}{c}{$L = 7$} \\\\\n \\cmidrule(lr){3-6} \\cmidrule(lr){7-10} \\cmidrule(lr){11-14}\n"

# print latex code
print(xt,
      include.rownames = FALSE,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      hline.after = c(-1, 0, nrow(display_table)),
      add.to.row = addtorow,
      comment = FALSE)