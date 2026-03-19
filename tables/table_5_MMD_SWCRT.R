source("utils_SWCRT.R")
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

L_values <- c(3, 5, 7)  # number of steps

# ICC parameter combinations for MMD based on the example table
param_sets <- list(
  # (0.05, 0.10) with different rho1 ranges
  list(rho0_E_min = 0.05, rho0_E_max = 0.10, rho1_E_min = 0.025, rho1_E_max = 0.040,
       rho0_C_min = 0.04, rho0_C_max = 0.08, rho1_C_min = 0.020, rho1_C_max = 0.032,
       rho0_EC_min = 0.01, rho0_EC_max = 0.02, rho1_EC_min = 0.005, rho1_EC_max = 0.01,
       rho2_EC_min = 0.5, rho2_EC_max = 0.8),
  list(rho0_E_min = 0.05, rho0_E_max = 0.10, rho1_E_min = 0.025, rho1_E_max = 0.045,
       rho0_C_min = 0.04, rho0_C_max = 0.08, rho1_C_min = 0.020, rho1_C_max = 0.036,
       rho0_EC_min = 0.01, rho0_EC_max = 0.02, rho1_EC_min = 0.005, rho1_EC_max = 0.01,
       rho2_EC_min = 0.5, rho2_EC_max = 0.8),
  
  # (0.05, 0.20) with different rho1 ranges
  list(rho0_E_min = 0.05, rho0_E_max = 0.20, rho1_E_min = 0.025, rho1_E_max = 0.040,
       rho0_C_min = 0.04, rho0_C_max = 0.16, rho1_C_min = 0.020, rho1_C_max = 0.032,
       rho0_EC_min = 0.01, rho0_EC_max = 0.02, rho1_EC_min = 0.005, rho1_EC_max = 0.01,
       rho2_EC_min = 0.5, rho2_EC_max = 0.8),
  list(rho0_E_min = 0.05, rho0_E_max = 0.20, rho1_E_min = 0.025, rho1_E_max = 0.045,
       rho0_C_min = 0.04, rho0_C_max = 0.16, rho1_C_min = 0.020, rho1_C_max = 0.036,
       rho0_EC_min = 0.01, rho0_EC_max = 0.02, rho1_EC_min = 0.005, rho1_EC_max = 0.01,
       rho2_EC_min = 0.5, rho2_EC_max = 0.8),
  
  # (0.10, 0.20) with different rho1 ranges
  list(rho0_E_min = 0.10, rho0_E_max = 0.20, rho1_E_min = 0.050, rho1_E_max = 0.080,
       rho0_C_min = 0.08, rho0_C_max = 0.16, rho1_C_min = 0.040, rho1_C_max = 0.064,
       rho0_EC_min = 0.01, rho0_EC_max = 0.02, rho1_EC_min = 0.005, rho1_EC_max = 0.01,
       rho2_EC_min = 0.5, rho2_EC_max = 0.8),
  list(rho0_E_min = 0.10, rho0_E_max = 0.20, rho1_E_min = 0.050, rho1_E_max = 0.090,
       rho0_C_min = 0.08, rho0_C_max = 0.16, rho1_C_min = 0.040, rho1_C_max = 0.072,
       rho0_EC_min = 0.01, rho0_EC_max = 0.02, rho1_EC_min = 0.005, rho1_EC_max = 0.01,
       rho2_EC_min = 0.5, rho2_EC_max = 0.8)
)

# ============================================================================================================
# run MMD for all combinations
# ============================================================================================================

results <- data.frame()

for (i in 1:length(param_sets)) {
  params <- param_sets[[i]]
  
  # two rows per parameter set: one for J = L+1, one for J = 9
  for (j_type in c("L+1", "9")) {
    row_result <- data.frame(
      rho0_E_min = params$rho0_E_min,
      rho0_E_max = params$rho0_E_max,
      rho1_E_min = params$rho1_E_min,
      rho1_E_max = params$rho1_E_max,
      J_type = j_type,
      stringsAsFactors = FALSE
    )
    
    # for each L value
    for (L in L_values) {
      if (j_type == "L+1") {
        J <- L + 1
      } else {
        J <- 9
      }
      
      # skip if J < L + 1
      if (J < L + 1) {
        row_result[[paste0("I_L", L)]] <- NA
        row_result[[paste0("K_L", L)]] <- NA
        row_result[[paste0("RE_L", L)]] <- NA
        next
      }
      
      result <- OD_CE_SW_CRT(
        alpha = alpha,
        trteff = trteff,
        lambda = lambda,
        sigma_E = sigma_E,
        sigma_C = sigma_C,
        rho0_E_min = params$rho0_E_min,
        rho0_E_max = params$rho0_E_max,
        rho1_E_min = params$rho1_E_min,
        rho1_E_max = params$rho1_E_max,
        rho0_C_min = params$rho0_C_min,
        rho0_C_max = params$rho0_C_max,
        rho1_C_min = params$rho1_C_min,
        rho1_C_max = params$rho1_C_max,
        rho0_EC_min = params$rho0_EC_min,
        rho0_EC_max = params$rho0_EC_max,
        rho1_EC_min = params$rho1_EC_min,
        rho1_EC_max = params$rho1_EC_max,
        rho2_EC_min = params$rho2_EC_min,
        rho2_EC_max = params$rho2_EC_max,
        rho0_E_init = (params$rho0_E_min + params$rho0_E_max) / 2,
        rho1_E_init = (params$rho1_E_min + params$rho1_E_max) / 2,
        rho0_C_init = (params$rho0_C_min + params$rho0_C_max) / 2,
        rho1_C_init = (params$rho1_C_min + params$rho1_C_max) / 2,
        rho0_EC_init = 0.015,
        rho1_EC_init = 0.0075,
        rho2_EC_init = 0.7,
        J = J,
        L = L,
        c1 = c1,
        c2 = c2,
        B = B,
        max_I = max_I,
        max_K = max_K
      )
      
      if (!is.null(result$integer_estimates)) {
        row_result[[paste0("I_L", L)]] <- result$integer_estimates$num_clusters
        row_result[[paste0("K_L", L)]] <- result$integer_estimates$cluster_size
        row_result[[paste0("RE_L", L)]] <- result$integer_estimates$relative_efficiency
      } else {
        row_result[[paste0("I_L", L)]] <- NA
        row_result[[paste0("K_L", L)]] <- NA
        row_result[[paste0("RE_L", L)]] <- NA
      }
    }
    
    results <- rbind(results, row_result)
  }
}

# ============================================================================================================
# format table for latex output
# ============================================================================================================

# create formatted strings for ICC ranges
results$rho0_E_str <- sprintf("(%.2f, %.2f)", results$rho0_E_min, results$rho0_E_max)
results$rho1_E_str <- sprintf("(%.3f, %.3f)", results$rho1_E_min, results$rho1_E_max)

# format I and K as integers, RE with 3 decimals
format_integer <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.0f", x))
}

format_decimal <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.3f", x))
}

display_table <- data.frame(
  "$(\\rho_{0\\min}^E, \\rho_{0\\max}^E)$" = results$rho0_E_str,
  "$(\\rho_{1\\min}^E, \\rho_{1\\max}^E)$" = results$rho1_E_str,
  "$J$" = results$J_type,
  "$I_{\\text{MMD}}$" = format_integer(results$I_L3),
  "$K_{\\text{MMD}}$" = format_integer(results$K_L3),
  "RE" = format_decimal(results$RE_L3),
  "spacer1" = "",  # empty column for spacing
  "$I_{\\text{MMD}}$.1" = format_integer(results$I_L5),
  "$K_{\\text{MMD}}$.1" = format_integer(results$K_L5),
  "RE.1" = format_decimal(results$RE_L5),
  "spacer2" = "",  # empty column for spacing
  "$I_{\\text{MMD}}$.2" = format_integer(results$I_L7),
  "$K_{\\text{MMD}}$.2" = format_integer(results$K_L7),
  "RE.2" = format_decimal(results$RE_L7),
  check.names = FALSE
)

# rename columns properly
names(display_table) <- c(
  "$(\\rho_{0\\min}^E, \\rho_{0\\max}^E)$",
  "$(\\rho_{1\\min}^E, \\rho_{1\\max}^E)$",
  "$J$",
  "$I_{\\text{MMD}}$", "$K_{\\text{MMD}}$", "RE",
  "",
  "$I_{\\text{MMD}}$", "$K_{\\text{MMD}}$", "RE",
  "",
  "$I_{\\text{MMD}}$", "$K_{\\text{MMD}}$", "RE"
)

# ============================================================================================================
# generate latex table
# ============================================================================================================

xt <- xtable(display_table,
             caption = "Maximin optimal designs (MMDs) for SW-CRTs with varying number of steps",
             label = "tab:SWCRT_MMD",
             align = c("l", "c", "c", "c", rep("c", 3), "c", rep("c", 3), "c", rep("c", 3)))

# add column groups header
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- "& & & \\multicolumn{3}{c}{$L = 3$} & & \\multicolumn{3}{c}{$L = 5$} & & \\multicolumn{3}{c}{$L = 7$} \\\\\n \\cmidrule(lr){4-6} \\cmidrule(lr){8-10} \\cmidrule(lr){12-14}\n"

# print latex code
print(xt,
      include.rownames = FALSE,
      sanitize.colnames.function = identity,
      sanitize.text.function = function(x) x,  # don't sanitize to preserve formatting
      hline.after = c(-1, 0, nrow(display_table)),
      add.to.row = addtorow,
      comment = FALSE)