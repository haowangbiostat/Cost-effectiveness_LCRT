source("utils_CRXO.R")
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggh4x)

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

# ============================================================================================================
# specific parameters
# ============================================================================================================

J <- 4
L <- J - 1         # number of steps for SW-CRT
rho0_E <- 0.1
rho2_EC <- 0.5
rho0_EC <- 0.04
rho0_C <- rho0_E

r <- sigma_E / sigma_C
pi <- 0.5
zalpha <- qnorm(1 - alpha/2)

# ============================================================================================================
# cac range and analysis for CRXO
# ============================================================================================================

CAC_values <- seq(0.1, 0.8, by = 0.01)
results_crxo <- data.frame()

cat("Processing CRXO analysis...\n")

for (CAC in CAC_values) {
  # calculate between-period correlations
  rho1_E <- CAC * rho0_E
  rho1_C <- CAC * rho0_C
  rho1_EC <- CAC * rho0_EC
  
  # run optimization for CRXO
  result <- OD_CE_CRXO_MPD(
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
    K_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$cluster_size, NA),
    I_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$num_clusters, NA),
    Power_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$power, NA),
    K_dec = ifelse(!is.null(result$decimal_estimates), result$decimal_estimates$cluster_size, NA),
    I_dec = ifelse(!is.null(result$decimal_estimates), result$decimal_estimates$num_clusters, NA),
    Power_dec = ifelse(!is.null(result$decimal_estimates), result$decimal_estimates$power, NA)
  ))
}

# ============================================================================================================
# cac range and analysis for PA-LCRT
# ============================================================================================================
source("utils_PA.R")
results_pa <- data.frame()

cat("Processing PA-LCRT analysis...\n")

for (CAC in CAC_values) {
  # calculate between-period correlations
  rho1_E <- CAC * rho0_E
  rho1_C <- CAC * rho0_C
  rho1_EC <- CAC * rho0_EC
  
  # run optimization for PA-LCRT
  result <- OD_CE_PA_LCRT(
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
    K_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$cluster_size, NA),
    I_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$num_clusters, NA),
    Power_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$power, NA),
    K_dec = ifelse(!is.null(result$decimal_estimates), result$decimal_estimates$cluster_size, NA),
    I_dec = ifelse(!is.null(result$decimal_estimates), result$decimal_estimates$num_clusters, NA),
    Power_dec = ifelse(!is.null(result$decimal_estimates), result$decimal_estimates$power, NA)
  ))
}

# ============================================================================================================
# cac range and analysis for SW-CRT
# ============================================================================================================
source("utils_SWCRT.R")
results_swcrt <- data.frame()

cat("Processing SW-CRT analysis...\n")

for (CAC in CAC_values) {
  # calculate between-period correlations
  rho1_E <- CAC * rho0_E
  rho1_C <- CAC * rho0_C
  rho1_EC <- CAC * rho0_EC
  
  # run optimization for SW-CRT
  result <- OD_CE_SW_CRT(
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
  
  # store SW-CRT results
  results_swcrt <- rbind(results_swcrt, data.frame(
    CAC = CAC,
    K_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$cluster_size, NA),
    I_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$num_clusters, NA),
    Power_int = ifelse(!is.null(result$integer_estimates), result$integer_estimates$power, NA),
    K_dec = NA,  # no decimal estimates for SW-CRT
    I_dec = NA,
    Power_dec = NA
  ))
}

# ============================================================================================================
# prepare combined data for plotting
# ============================================================================================================

# CRXO data
df_design_crxo <- data.frame(
  CAC = rep(results_crxo$CAC, 4),
  Value = c(results_crxo$K_int, results_crxo$K_dec, results_crxo$I_int, results_crxo$I_dec),
  Parameter = rep(c("K", "K", "I", "I"), each = nrow(results_crxo)),
  Type = rep(c("Integer", "Decimal", "Integer", "Decimal"), each = nrow(results_crxo)),
  Design = "CRXO"
)

df_power_crxo <- data.frame(
  CAC = rep(results_crxo$CAC, 2),
  Power = c(results_crxo$Power_int, results_crxo$Power_dec),
  Type = rep(c("Integer", "Decimal"), each = nrow(results_crxo)),
  Design = "CRXO"
)

# PA-LCRT data
df_design_pa <- data.frame(
  CAC = rep(results_pa$CAC, 4),
  Value = c(results_pa$K_int, results_pa$K_dec, results_pa$I_int, results_pa$I_dec),
  Parameter = rep(c("K", "K", "I", "I"), each = nrow(results_pa)),
  Type = rep(c("Integer", "Decimal", "Integer", "Decimal"), each = nrow(results_pa)),
  Design = "PA-LCRT"
)

df_power_pa <- data.frame(
  CAC = rep(results_pa$CAC, 2),
  Power = c(results_pa$Power_int, results_pa$Power_dec),
  Type = rep(c("Integer", "Decimal"), each = nrow(results_pa)),
  Design = "PA-LCRT"
)

# SW-CRT data (integer only)
df_design_swcrt <- data.frame(
  CAC = rep(results_swcrt$CAC, 2),
  Value = c(results_swcrt$K_int, results_swcrt$I_int),
  Parameter = rep(c("K", "I"), each = nrow(results_swcrt)),
  Type = rep("Integer", 2 * nrow(results_swcrt)),
  Design = "SW-CRT"
)

df_power_swcrt <- data.frame(
  CAC = results_swcrt$CAC,
  Power = results_swcrt$Power_int,
  Type = rep("Integer", nrow(results_swcrt)),
  Design = "SW-CRT"
)

# combine all data
df_combined <- rbind(
  # CRXO LOD
  data.frame(
    CAC = df_design_crxo$CAC,
    Value = df_design_crxo$Value,
    Parameter = df_design_crxo$Parameter,
    Type = df_design_crxo$Type,
    Panel = "CRXO_LOD",
    Design = df_design_crxo$Design
  ),
  # CRXO Power
  data.frame(
    CAC = df_power_crxo$CAC,
    Value = df_power_crxo$Power,
    Parameter = "Power",
    Type = df_power_crxo$Type,
    Panel = "CRXO_Power",
    Design = df_power_crxo$Design
  ),
  # PA-LCRT LOD
  data.frame(
    CAC = df_design_pa$CAC,
    Value = df_design_pa$Value,
    Parameter = df_design_pa$Parameter,
    Type = df_design_pa$Type,
    Panel = "PA_LOD",
    Design = df_design_pa$Design
  ),
  # PA-LCRT Power
  data.frame(
    CAC = df_power_pa$CAC,
    Value = df_power_pa$Power,
    Parameter = "Power",
    Type = df_power_pa$Type,
    Panel = "PA_Power",
    Design = df_power_pa$Design
  ),
  # SW-CRT LOD
  data.frame(
    CAC = df_design_swcrt$CAC,
    Value = df_design_swcrt$Value,
    Parameter = df_design_swcrt$Parameter,
    Type = df_design_swcrt$Type,
    Panel = "SWCRT_LOD",
    Design = df_design_swcrt$Design
  ),
  # SW-CRT Power
  data.frame(
    CAC = df_power_swcrt$CAC,
    Value = df_power_swcrt$Power,
    Parameter = "Power",
    Type = df_power_swcrt$Type,
    Panel = "SWCRT_Power",
    Design = df_power_swcrt$Design
  )
)

# ============================================================================================================
# create combined 3x2 plot
# ============================================================================================================

p_combined <- ggplot(df_combined, aes(x = CAC)) +
  # CRXO LOD lines
  geom_line(data = subset(df_combined, Panel == "CRXO_LOD" & Type == "Integer"), 
            aes(y = Value, color = interaction(Parameter, "Integer")), size = 1) +
  geom_line(data = subset(df_combined, Panel == "CRXO_LOD" & Type == "Decimal"), 
            aes(y = Value, color = interaction(Parameter, "Decimal")), size = 1, linetype = "dotted") +
  # CRXO Power lines
  geom_line(data = subset(df_combined, Panel == "CRXO_Power" & Type == "Integer"), 
            aes(y = Value, color = "Power.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "CRXO_Power" & Type == "Decimal"), 
            aes(y = Value, color = "Power.Decimal"), size = 1, linetype = "dotted") +
  # PA-LCRT LOD lines
  geom_line(data = subset(df_combined, Panel == "PA_LOD" & Type == "Integer"), 
            aes(y = Value, color = interaction(Parameter, "Integer")), size = 1) +
  geom_line(data = subset(df_combined, Panel == "PA_LOD" & Type == "Decimal"), 
            aes(y = Value, color = interaction(Parameter, "Decimal")), size = 1, linetype = "dotted") +
  # PA-LCRT Power lines
  geom_line(data = subset(df_combined, Panel == "PA_Power" & Type == "Integer"), 
            aes(y = Value, color = "Power.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "PA_Power" & Type == "Decimal"), 
            aes(y = Value, color = "Power.Decimal"), size = 1, linetype = "dotted") +
  # SW-CRT LOD lines (integer only)
  geom_line(data = subset(df_combined, Panel == "SWCRT_LOD" & Type == "Integer"), 
            aes(y = Value, color = interaction(Parameter, "Integer")), size = 1) +
  # SW-CRT Power lines (integer only)
  geom_line(data = subset(df_combined, Panel == "SWCRT_Power" & Type == "Integer"), 
            aes(y = Value, color = "Power.Integer"), size = 1) +
  scale_color_manual(name = "",
                     values = c("I.Integer" = "#0072B2", 
                                "K.Integer" = "#D55E00",
                                "Power.Integer" = "blue",
                                "I.Decimal" = "#0072B2",
                                "K.Decimal" = "#D55E00",
                                "Power.Decimal" = "blue"),
                     labels = c("I.Integer" = expression(I[LOD]~(Integer)),
                                "K.Integer" = expression(K[LOD]~(Integer)),
                                "Power.Integer" = "Power (Integer)",
                                "I.Decimal" = expression(I[LOD]~(Decimal)),
                                "K.Decimal" = expression(K[LOD]~(Decimal)),
                                "Power.Decimal" = "Power (Decimal)"),
                     breaks = c("I.Integer", "K.Integer", "Power.Integer",
                                "I.Decimal", "K.Decimal", "Power.Decimal")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "CAC", y = NULL) +
  scale_x_continuous(breaks = seq(0.1, 0.8, by = 0.1)) +
  facet_wrap(~Panel, ncol = 2, nrow = 3,
             labeller = labeller(Panel = c("CRXO_LOD" = "(a) LODs for CRXO Trials (J = 4)",
                                           "CRXO_Power" = "(b) Power of the LODs for CRXO Trials (J = 4)",
                                           "PA_LOD" = "(c) LODs for PA-LCRTs (J = 4)",
                                           "PA_Power" = "(d) Power of the LODs for PA-LCRTs (J = 4)",
                                           "SWCRT_LOD" = "(e) LODs for SW-CRTs (J = 4, L = 3)",
                                           "SWCRT_Power" = "(f) Power of the LODs for SW-CRTs (J = 4, L = 3)")),
             scales = "free_y") +
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
      Panel == "CRXO_LOD" ~ scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)),
      Panel == "CRXO_Power" ~ scale_y_continuous(limits = c(0.20, 1.0), breaks = seq(0.20, 1.00, by = 0.1), name = "Power"),
      Panel == "PA_LOD" ~ scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)),
      Panel == "PA_Power" ~ scale_y_continuous(limits = c(0.20, 1.0), breaks = seq(0.20, 1.00, by = 0.1), name = "Power"),
      Panel == "SWCRT_LOD" ~ scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)),
      Panel == "SWCRT_Power" ~ scale_y_continuous(limits = c(0.20, 1.0), breaks = seq(0.20, 1.00, by = 0.1), name = "Power")
    )
  )

library(patchwork)

# shared elements
shared_theme <- theme(
  legend.position = "bottom",
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

# left plot: sample size
p_left <- ggplot(subset(df_combined, Panel %in% c("CRXO_LOD", "PA_LOD", "SWCRT_LOD")), aes(x = CAC)) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_LOD", "PA_LOD", "SWCRT_LOD") & Type == "Integer"),
            aes(y = Value, color = interaction(Parameter, "Integer")), size = 1) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_LOD", "PA_LOD", "SWCRT_LOD") & Type == "Decimal"),
            aes(y = Value, color = interaction(Parameter, "Decimal")), size = 1, linetype = "dotted") +
  scale_color_manual(name = "",
                     values = c("I.Integer" = "#0072B2", "K.Integer" = "#D55E00",
                                "I.Decimal" = "#0072B2", "K.Decimal" = "#D55E00"),
                     labels = c("I.Integer" = expression(I[LOD]~(Integer)),
                                "K.Integer" = expression(K[LOD]~(Integer)),
                                "I.Decimal" = expression(I[LOD]~(Decimal)),
                                "K.Decimal" = expression(K[LOD]~(Decimal))),
                     breaks = c("I.Integer", "K.Integer", "I.Decimal", "K.Decimal")) +
  scale_x_continuous(breaks = seq(0.1, 0.8, by = 0.1)) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
  labs(x = "CAC", y = expression(paste("Optimal sample size combination (", I[LOD], ", ", K[LOD], ")"))) +
  facet_wrap(~Panel, ncol = 1,
             labeller = labeller(Panel = c("CRXO_LOD" = "(a) LODs for CRXO trials (J = 4)",
                                           "PA_LOD" = "(c) LODs for PA-LCRTs (J = 4)",
                                           "SWCRT_LOD" = "(e) LODs for SW-CRTs (J = 4, L = 3)"))) +
  shared_theme

# right plot: power
p_right <- ggplot(subset(df_combined, Panel %in% c("CRXO_Power", "PA_Power", "SWCRT_Power")), aes(x = CAC)) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_Power", "PA_Power", "SWCRT_Power") & Type == "Integer"),
            aes(y = Value, color = "Power.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_Power", "PA_Power", "SWCRT_Power") & Type == "Decimal"),
            aes(y = Value, color = "Power.Decimal"), size = 1, linetype = "dotted") +
  scale_color_manual(name = "",
                     values = c("Power.Integer" = "blue", "Power.Decimal" = "blue"),
                     labels = c("Power.Integer" = "Power (Integer)", "Power.Decimal" = "Power (Decimal)"),
                     breaks = c("Power.Integer", "Power.Decimal")) +
  scale_x_continuous(breaks = seq(0.1, 0.8, by = 0.1)) +
  scale_y_continuous(limits = c(0.20, 1.0), breaks = seq(0.20, 1.00, by = 0.1)) +
  labs(x = "CAC", y = "Power") +
  facet_wrap(~Panel, ncol = 1,
             labeller = labeller(Panel = c("CRXO_Power" = "(b) Power of the LODs for CRXO trials (J = 4)",
                                           "PA_Power" = "(d) Power of the LODs for PA-LCRTs (J = 4)",
                                           "SWCRT_Power" = "(f) Power of the LODs for SW-CRTs (J = 4, L = 3)"))) +
  shared_theme

# combine and collect legends
p_combined <- p_left + p_right + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

print(p_combined)
ggsave("../figures/figure_3.pdf", p_combined, width = 12, height = 8)