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
trteff <- 4000     # average INMB
alpha <- 0.05      # significance level
pinum <- 1         # pi numerator
piden <- 2         # pi denominator (pi = 0.5)
max_I <- 100       # max clusters
max_K <- 200       # max cluster size
rho2_EC <- 0.5     # within-individual effect-cost correlation

# ============================================================================================================
# design parameters
# ============================================================================================================

J_values <- c(2, 4, 6)                    # periods
rho0_E_values <- c(0.05, 0.10, 0.20)      # within-period effectiveness ICC
CAC_values <- c(0.5, 0.8)                 # cluster autocorrelation
cost_multipliers <- c(1)                  # cost ICC multipliers

# ============================================================================================================
# generate parameter combinations
# ============================================================================================================

param_sets <- list()
for (rho0_E in rho0_E_values) {
  for (CAC in CAC_values) {
    for (mult in cost_multipliers) {
      # calculate correlations
      rho1_E <- CAC * rho0_E
      rho0_C <- mult * rho0_E
      rho1_C <- CAC * rho0_C
      rho0_EC <- 0.4 * rho0_E
      rho1_EC <- CAC * rho0_EC
      
      param_sets <- append(param_sets, list(list(
        rho0_E = rho0_E, rho1_E = rho1_E,
        rho0_C = rho0_C, rho1_C = rho1_C,
        rho0_EC = rho0_EC, rho1_EC = rho1_EC,
        rho2_EC = rho2_EC
      )))
    }
  }
}

# ============================================================================================================
# run LOD for all combinations
# ============================================================================================================

results <- data.frame()

for (i in 1:length(param_sets)) {
  params <- param_sets[[i]]
  
  # initialize row with parameters
  row_result <- data.frame(
    rho_E = sprintf("(%.2f, %.3f)", params$rho0_E, params$rho1_E),
    rho_C = sprintf("(%.2f, %.3f)", params$rho0_C, params$rho1_C),
    rho_EC = sprintf("(%.2f, %.3f, %.1f)", params$rho0_EC, params$rho1_EC, params$rho2_EC),
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
      pinum = pinum,
      piden = piden,
      c1 = c1,
      c2 = c2,
      B = B,
      max_I = max_I,
      max_K = max_K
    )
    
    # extract integer estimates from LOD results
    if (!is.null(result$integer_estimates)) {
      row_result[[paste0("I_J", J)]] <- sprintf("%.0f", result$integer_estimates$num_clusters)
      row_result[[paste0("K_J", J)]] <- sprintf("%.0f", result$integer_estimates$cluster_size)
      row_result[[paste0("Power_J", J)]] <- sprintf("%.3f", result$integer_estimates$power)
    } else {
      row_result[[paste0("I_J", J)]] <- NA
      row_result[[paste0("K_J", J)]] <- NA
      row_result[[paste0("Power_J", J)]] <- NA
    }
  }
  
  results <- rbind(results, row_result)
}

# ============================================================================================================
# format table for latex output
# ============================================================================================================

display_table <- data.frame(
  "$(\\rho_0^E, \\rho_1^E)$" = results$rho_E,
  "$(\\rho_0^C, \\rho_1^C)$" = results$rho_C,
  "$(\\rho_0^{EC}, \\rho_1^{EC}, \\rho_2^{EC})$" = results$rho_EC,
  "$I_{\\text{LOD}}$" = results$I_J2,
  "$K_{\\text{LOD}}$" = results$K_J2,
  "Power" = results$Power_J2,
  "$I_{\\text{LOD}}$.1" = results$I_J4,
  "$K_{\\text{LOD}}$.1" = results$K_J4,
  "Power.1" = results$Power_J4,
  "$I_{\\text{LOD}}$.2" = results$I_J6,
  "$K_{\\text{LOD}}$.2" = results$K_J6,
  "Power.2" = results$Power_J6,
  check.names = FALSE
)

# rename columns properly
names(display_table) <- c(
  "$(\\rho_0^E, \\rho_1^E)$",
  "$(\\rho_0^C, \\rho_1^C)$",
  "$(\\rho_0^{EC}, \\rho_1^{EC}, \\rho_2^{EC})$",
  "$I_{\\text{LOD}}$", "$K_{\\text{LOD}}$", "Power",
  "$I_{\\text{LOD}}$", "$K_{\\text{LOD}}$", "Power",
  "$I_{\\text{LOD}}$", "$K_{\\text{LOD}}$", "Power"
)

# ============================================================================================================
# generate latex table
# ============================================================================================================

xt <- xtable(display_table,
             caption = "Local optimal designs (LODs) for PA-LCRTs with $\\pi = 0.5$",
             label = "tab:PA_LOD_0.5",
             align = c("l", "l", "l", "l", rep("c", 9)))

# add column groups header
addtorow <- list()
addtorow$pos <- list(0)
addtorow$command <- "& & & \\multicolumn{3}{c}{$J = 2$} & \\multicolumn{3}{c}{$J = 4$} & \\multicolumn{3}{c}{$J = 6$} \\\\\n \\cmidrule(lr){4-6} \\cmidrule(lr){7-9} \\cmidrule(lr){10-12}\n"

# print latex code
print(xt,
      include.rownames = FALSE,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      hline.after = c(-1, 0, nrow(display_table)),
      add.to.row = addtorow,
      comment = FALSE)