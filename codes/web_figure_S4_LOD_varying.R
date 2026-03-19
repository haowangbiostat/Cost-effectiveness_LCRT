library(ggplot2)
library(tidyr)
library(dplyr)
library(ggh4x)

# ============================================================================================================
# fixed parameters (same across all scenarios)
# ============================================================================================================

B <- 300000        # total budget
c1 <- 3000         # cluster-level cost
c2 <- 250          # participant-level cost
lambda <- 20000    # ceiling ratio (fixed)
sigma_E <- 1       # effectiveness sd (fixed)
trteff <- 4000     # average INMB
alpha <- 0.05      # significance level
pinum <- 1         # pi numerator
piden <- 2         # pi denominator (pi = 0.5)
max_I <- 100       # max clusters
max_K <- 200       # max cluster size

# ============================================================================================================
# specific parameters
# ============================================================================================================

J <- 4
L <- J - 1         # number of steps for SW-CRT
rho0_E <- 0.1
rho2_EC <- 0.5
rho0_EC <- 0.04
rho0_C <- rho0_E

pi <- 0.5
zalpha <- qnorm(1 - alpha/2)

# ============================================================================================================
# three scenarios for lambda*r by changing sigma_C
# ============================================================================================================

scenarios <- list(
  list(name = "lr < 1", sigma_C = 200000, label = "lr0.1"),
  list(name = "lr = 1", sigma_C = 20000, label = "lr1"),
  list(name = "lr > 1", sigma_C = 10000, label = "lr2")
)

# ============================================================================================================
# cac range and analysis
# ============================================================================================================

CAC_values <- seq(0.1, 0.8, by = 0.01)
all_results_crxo <- data.frame()
all_results_pa <- data.frame()
all_results_swcrt <- data.frame()

cat("Processing LOD analysis with varying lambda*r...\n")

for (scenario in scenarios) {
  sigma_C <- scenario$sigma_C
  r <- sigma_E / sigma_C
  scenario_name <- scenario$name
  scenario_label <- scenario$label
  
  cat(sprintf("\nScenario: %s (lambda*r = %.1f)\n", scenario_name, lambda * r))
  
  results_crxo <- data.frame()
  results_pa <- data.frame()
  results_swcrt <- data.frame()
  
  # ============================================================================================================
  # CRXO analysis
  # ============================================================================================================
  source("utils_CRXO.R")
  cat("  Processing CRXO...\n")
  for (CAC in CAC_values) {
    # calculate between-period correlations
    rho1_E <- CAC * rho0_E
    rho1_C <- CAC * rho0_C
    rho1_EC <- CAC * rho0_EC
    
    result_crxo <- OD_CE_CRXO_MPD(
      alpha = alpha, trteff = trteff, lambda = lambda,
      sigma_E = sigma_E, sigma_C = sigma_C,
      rho0_E_min = rho0_E, rho0_E_max = rho0_E,
      rho1_E_min = rho1_E, rho1_E_max = rho1_E,
      rho0_C_min = rho0_C, rho0_C_max = rho0_C,
      rho1_C_min = rho1_C, rho1_C_max = rho1_C,
      rho0_EC_min = rho0_EC, rho0_EC_max = rho0_EC,
      rho1_EC_min = rho1_EC, rho1_EC_max = rho1_EC,
      rho2_EC_min = rho2_EC, rho2_EC_max = rho2_EC,
      J = J, pinum = pinum, piden = piden,
      c1 = c1, c2 = c2, B = B,
      max_I = max_I, max_K = max_K
    )
    
    # store CRXO results
    results_crxo <- rbind(results_crxo, data.frame(
      CAC = CAC,
      K_int = ifelse(!is.null(result_crxo$integer_estimates), result_crxo$integer_estimates$cluster_size, NA),
      I_int = ifelse(!is.null(result_crxo$integer_estimates), result_crxo$integer_estimates$num_clusters, NA),
      K_dec = ifelse(!is.null(result_crxo$decimal_estimates), result_crxo$decimal_estimates$cluster_size, NA),
      I_dec = ifelse(!is.null(result_crxo$decimal_estimates), result_crxo$decimal_estimates$num_clusters, NA),
      Scenario = scenario_name,
      ScenarioLabel = scenario_label,
      Design = "CRXO"
    ))
  }
  
  # ============================================================================================================
  # PA-LCRT analysis
  # ============================================================================================================
  source("utils_PA.R")
  cat("  Processing PA-LCRT...\n")
  for (CAC in CAC_values) {
    # calculate between-period correlations
    rho1_E <- CAC * rho0_E
    rho1_C <- CAC * rho0_C
    rho1_EC <- CAC * rho0_EC
    
    result_pa <- OD_CE_PA_LCRT(
      alpha = alpha, trteff = trteff, lambda = lambda,
      sigma_E = sigma_E, sigma_C = sigma_C,
      rho0_E_min = rho0_E, rho0_E_max = rho0_E,
      rho1_E_min = rho1_E, rho1_E_max = rho1_E,
      rho0_C_min = rho0_C, rho0_C_max = rho0_C,
      rho1_C_min = rho1_C, rho1_C_max = rho1_C,
      rho0_EC_min = rho0_EC, rho0_EC_max = rho0_EC,
      rho1_EC_min = rho1_EC, rho1_EC_max = rho1_EC,
      rho2_EC_min = rho2_EC, rho2_EC_max = rho2_EC,
      J = J, pinum = pinum, piden = piden,
      c1 = c1, c2 = c2, B = B,
      max_I = max_I, max_K = max_K
    )
    
    # store PA-LCRT results
    results_pa <- rbind(results_pa, data.frame(
      CAC = CAC,
      K_int = ifelse(!is.null(result_pa$integer_estimates), result_pa$integer_estimates$cluster_size, NA),
      I_int = ifelse(!is.null(result_pa$integer_estimates), result_pa$integer_estimates$num_clusters, NA),
      K_dec = ifelse(!is.null(result_pa$decimal_estimates), result_pa$decimal_estimates$cluster_size, NA),
      I_dec = ifelse(!is.null(result_pa$decimal_estimates), result_pa$decimal_estimates$num_clusters, NA),
      Scenario = scenario_name,
      ScenarioLabel = scenario_label,
      Design = "PA-LCRT"
    ))
  }
  
  # ============================================================================================================
  # SW-CRT analysis
  # ============================================================================================================
  source("utils_SWCRT.R")
  cat("  Processing SW-CRT...\n")
  for (CAC in CAC_values) {
    # calculate between-period correlations
    rho1_E <- CAC * rho0_E
    rho1_C <- CAC * rho0_C
    rho1_EC <- CAC * rho0_EC
    
    result_swcrt <- OD_CE_SW_CRT(
      alpha = alpha, trteff = trteff, lambda = lambda,
      sigma_E = sigma_E, sigma_C = sigma_C,
      rho0_E_min = rho0_E, rho0_E_max = rho0_E,
      rho1_E_min = rho1_E, rho1_E_max = rho1_E,
      rho0_C_min = rho0_C, rho0_C_max = rho0_C,
      rho1_C_min = rho1_C, rho1_C_max = rho1_C,
      rho0_EC_min = rho0_EC, rho0_EC_max = rho0_EC,
      rho1_EC_min = rho1_EC, rho1_EC_max = rho1_EC,
      rho2_EC_min = rho2_EC, rho2_EC_max = rho2_EC,
      J = J, L = L,
      c1 = c1, c2 = c2, B = B,
      max_I = max_I, max_K = max_K
    )
    
    # store SW-CRT results (integer only)
    results_swcrt <- rbind(results_swcrt, data.frame(
      CAC = CAC,
      K_int = ifelse(!is.null(result_swcrt$integer_estimates), result_swcrt$integer_estimates$cluster_size, NA),
      I_int = ifelse(!is.null(result_swcrt$integer_estimates), result_swcrt$integer_estimates$num_clusters, NA),
      K_dec = NA,  # no decimal estimates for SW-CRT
      I_dec = NA,
      Scenario = scenario_name,
      ScenarioLabel = scenario_label,
      Design = "SW-CRT"
    ))
  }
  
  all_results_crxo <- rbind(all_results_crxo, results_crxo)
  all_results_pa <- rbind(all_results_pa, results_pa)
  all_results_swcrt <- rbind(all_results_swcrt, results_swcrt)
}

# combine all results
all_results <- rbind(all_results_crxo, all_results_pa, all_results_swcrt)

# ============================================================================================================
# prepare data for plotting
# ============================================================================================================

df_design <- data.frame(
  CAC = rep(all_results$CAC, 4),
  Value = c(all_results$K_int, all_results$K_dec, all_results$I_int, all_results$I_dec),
  Parameter = rep(c("K", "K", "I", "I"), each = nrow(all_results)),
  Type = rep(c("Integer", "Decimal", "Integer", "Decimal"), each = nrow(all_results)),
  Scenario = rep(all_results$Scenario, 4),
  ScenarioLabel = rep(all_results$ScenarioLabel, 4),
  Design = rep(all_results$Design, 4)
)

# create panel labels for 3x3 layout
df_design$Panel <- paste0(df_design$Design, "_", df_design$ScenarioLabel)

# convert to factor with proper ordering and labels
df_design$Panel <- factor(df_design$Panel,
                          levels = c("CRXO_lr0.1", "CRXO_lr1", "CRXO_lr2",
                                     "PA-LCRT_lr0.1", "PA-LCRT_lr1", "PA-LCRT_lr2",
                                     "SW-CRT_lr0.1", "SW-CRT_lr1", "SW-CRT_lr2"),
                          labels = c("(a)~'CRXO trials:'~lambda*r==0.1",
                                     "(b)~'CRXO trials:'~lambda*r==1",
                                     "(c)~'CRXO trials:'~lambda*r==2",
                                     "(d)~'PA-LCRTs:'~lambda*r==0.1",
                                     "(e)~'PA-LCRTs:'~lambda*r==1",
                                     "(f)~'PA-LCRTs:'~lambda*r==2",
                                     "(g)~'SW-CRTs:'~lambda*r==0.1",
                                     "(h)~'SW-CRTs:'~lambda*r==1",
                                     "(i)~'SW-CRTs:'~lambda*r==2"))

# ============================================================================================================
# create combined 3x3 plot
# ============================================================================================================

p_combined <- ggplot(df_design, aes(x = CAC)) +
  geom_line(data = subset(df_design, Type == "Integer"), 
            aes(y = Value, color = interaction(Parameter, "Integer")), size = 1) +
  geom_line(data = subset(df_design, Type == "Decimal"), 
            aes(y = Value, color = interaction(Parameter, "Decimal")), size = 1, linetype = "dotted") +
  scale_color_manual(name = "",
                     values = c("I.Integer" = "#0072B2", 
                                "K.Integer" = "#D55E00",
                                "I.Decimal" = "#0072B2",
                                "K.Decimal" = "#D55E00"),
                     labels = c("I.Integer" = expression(I[LOD]~(Integer)),
                                "K.Integer" = expression(K[LOD]~(Integer)),
                                "I.Decimal" = expression(I[LOD]~(Decimal)),
                                "K.Decimal" = expression(K[LOD]~(Decimal))),
                     breaks = c("I.Integer", "K.Integer",
                                "I.Decimal", "K.Decimal")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "CAC", y = expression(paste("Optimal sample size combination (", I[LOD], ", ", K[LOD], ")"))) +
  scale_x_continuous(breaks = seq(0.1, 0.8, by = 0.1)) +
  facet_wrap(~Panel, ncol = 3, nrow = 3,
             scales = "free_y",
             labeller = label_parsed) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.text = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 11),
        strip.placement = "outside",
        panel.spacing = unit(0.5, "lines"))

# add y-axis scales for each panel
p_combined <- p_combined + 
  ggh4x::facetted_pos_scales(
    y = list(
      Panel == "(a)~'CRXO trials:'~lambda*r==0.1" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)),
      Panel == "(b)~'CRXO trials:'~lambda*r==1" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)),
      Panel == "(c)~'CRXO trials:'~lambda*r==2" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)),
      Panel == "(d)~'PA-LCRTs:'~lambda*r==0.1" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)),
      Panel == "(e)~'PA-LCRTs:'~lambda*r==1" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)),
      Panel == "(f)~'PA-LCRTs:'~lambda*r==2" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)),
      Panel == "(g)~'SW-CRTs:'~lambda*r==0.1" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)),
      Panel == "(h)~'SW-CRTs:'~lambda*r==1" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)),
      Panel == "(i)~'SW-CRTs:'~lambda*r==2" ~ scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10))
    )
  )

# ============================================================================================================
# save plot
# ============================================================================================================

print(p_combined)
ggsave("../figures/web_figure_s4.pdf", p_combined, width = 12, height = 9)