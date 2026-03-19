library(ggplot2)
library(tidyr)
library(dplyr)
library(ggh4x)

# fixed parameters
J <- 4
L <- J - 1  # number of steps for SW-CRT
sigma_E <- 1
sigma_C <- 3000
trteff <- 4000
lambda <- 20000
B <- 300000
c1 <- 3000
c2 <- 250
alpha <- 0.05
r <- sigma_E / sigma_C
pi <- 0.5

# ICC ranges
rho0_E_min <- 0.05
rho0_E_max <- 0.10
rho1_E_min <- 0.025

rho0_C_min <- 0.04
rho0_C_max <- 0.08
rho1_C_min <- 0.02

rho0_EC_min <- 0.01
rho0_EC_max <- 0.02
rho1_EC_min <- 0.005
rho1_EC_max <- 0.01
rho2_EC_min <- 0.5
rho2_EC_max <- 0.8

# initial values
rho0_E_init <- 0.075
rho0_C_init <- 0.06
rho0_EC_init <- 0.015
rho1_EC_init <- 0.0075
rho2_EC_init <- 0.65

# vary maximum between-period correlation
rho1_max_values <- seq(0.025, 0.045, by = 0.001)

# ============================================================================================================
# MMD analysis for CRXO
# ============================================================================================================

results_crxo <- data.frame()

cat("Starting CRXO MMD analysis...\n")

# Source CRXO utilities
source("utils_CRXO.R")

for (i in seq_along(rho1_max_values)) {
  rho1_max <- rho1_max_values[i]
  
  rho1_E_max <- rho1_max
  rho1_C_max <- rho1_max * 0.8
  
  rho1_E_init <- (rho1_E_min + rho1_E_max) / 2
  rho1_C_init <- (rho1_C_min + rho1_C_max) / 2
  
  cat(sprintf("CRXO: Processing rho1_max = %.3f\n", rho1_max))
  
  # call MMD function for CRXO
  capture.output({
    result <- OD_CE_CRXO_MPD(
      alpha = alpha, trteff = trteff, lambda = lambda,
      sigma_E = sigma_E, sigma_C = sigma_C,
      rho0_E_min = rho0_E_min, rho0_E_max = rho0_E_max,
      rho1_E_min = rho1_E_min, rho1_E_max = rho1_E_max,
      rho0_C_min = rho0_C_min, rho0_C_max = rho0_C_max,
      rho1_C_min = rho1_C_min, rho1_C_max = rho1_C_max,
      rho0_EC_min = rho0_EC_min, rho0_EC_max = rho0_EC_max,
      rho1_EC_min = rho1_EC_min, rho1_EC_max = rho1_EC_max,
      rho2_EC_min = rho2_EC_min, rho2_EC_max = rho2_EC_max,
      rho0_E_init = rho0_E_init, rho1_E_init = rho1_E_init,
      rho0_C_init = rho0_C_init, rho1_C_init = rho1_C_init,
      rho0_EC_init = rho0_EC_init, rho1_EC_init = rho1_EC_init,
      rho2_EC_init = rho2_EC_init,
      J = J, c1 = c1, c2 = c2, B = B,
      max_I = 100, max_K = 200
    )
  })
  
  # extract CRXO results
  if (!is.null(result$integer_estimates) && !is.null(result$decimal_estimates)) {
    results_crxo <- rbind(results_crxo, data.frame(
      rho1_max = rho1_max,
      K_int = result$integer_estimates$cluster_size,
      I_int = result$integer_estimates$num_clusters,
      RE_int = result$integer_estimates$relative_efficiency,
      K_dec = result$decimal_estimates$cluster_size,
      I_dec = result$decimal_estimates$num_clusters,
      RE_dec = result$decimal_estimates$relative_efficiency
    ))
    cat(sprintf("  CRXO MMD found: K_int=%d, I_int=%d, RE_int=%.4f\n", 
                result$integer_estimates$cluster_size, 
                result$integer_estimates$num_clusters, 
                result$integer_estimates$relative_efficiency))
  } else {
    cat("  No CRXO MMD found\n")
  }
}

# ============================================================================================================
# MMD analysis for PA-LCRT
# ============================================================================================================

results_pa <- data.frame()

cat("\nStarting PA-LCRT MMD analysis...\n")

# Source PA-LCRT utilities
source("utils_PA.R")

for (i in seq_along(rho1_max_values)) {
  rho1_max <- rho1_max_values[i]
  
  rho1_E_max <- rho1_max
  rho1_C_max <- rho1_max * 0.8
  
  rho1_E_init <- (rho1_E_min + rho1_E_max) / 2
  rho1_C_init <- (rho1_C_min + rho1_C_max) / 2
  
  cat(sprintf("PA-LCRT: Processing rho1_max = %.3f\n", rho1_max))
  
  # call MMD function for PA-LCRT
  capture.output({
    result <- OD_CE_PA_LCRT(
      alpha = alpha, trteff = trteff, lambda = lambda,
      sigma_E = sigma_E, sigma_C = sigma_C,
      rho0_E_min = rho0_E_min, rho0_E_max = rho0_E_max,
      rho1_E_min = rho1_E_min, rho1_E_max = rho1_E_max,
      rho0_C_min = rho0_C_min, rho0_C_max = rho0_C_max,
      rho1_C_min = rho1_C_min, rho1_C_max = rho1_C_max,
      rho0_EC_min = rho0_EC_min, rho0_EC_max = rho0_EC_max,
      rho1_EC_min = rho1_EC_min, rho1_EC_max = rho1_EC_max,
      rho2_EC_min = rho2_EC_min, rho2_EC_max = rho2_EC_max,
      rho0_E_init = rho0_E_init, rho1_E_init = rho1_E_init,
      rho0_C_init = rho0_C_init, rho1_C_init = rho1_C_init,
      rho0_EC_init = rho0_EC_init, rho1_EC_init = rho1_EC_init,
      rho2_EC_init = rho2_EC_init,
      J = J, c1 = c1, c2 = c2, B = B,
      max_I = 100, max_K = 200
    )
  })
  
  # extract PA-LCRT results
  if (!is.null(result$integer_estimates) && !is.null(result$decimal_estimates)) {
    results_pa <- rbind(results_pa, data.frame(
      rho1_max = rho1_max,
      K_int = result$integer_estimates$cluster_size,
      I_int = result$integer_estimates$num_clusters,
      RE_int = result$integer_estimates$relative_efficiency,
      K_dec = result$decimal_estimates$cluster_size,
      I_dec = result$decimal_estimates$num_clusters,
      RE_dec = result$decimal_estimates$relative_efficiency
    ))
    cat(sprintf("  PA-LCRT MMD found: K_int=%d, I_int=%d, RE_int=%.4f\n", 
                result$integer_estimates$cluster_size, 
                result$integer_estimates$num_clusters, 
                result$integer_estimates$relative_efficiency))
  } else {
    cat("  No PA-LCRT MMD found\n")
  }
}

# ============================================================================================================
# MMD analysis for SW-CRT
# ============================================================================================================

results_swcrt <- data.frame()

cat("\nStarting SW-CRT MMD analysis...\n")

# Source SW-CRT utilities
source("utils_SWCRT.R")

for (i in seq_along(rho1_max_values)) {
  rho1_max <- rho1_max_values[i]
  
  rho1_E_max <- rho1_max
  rho1_C_max <- rho1_max * 0.8
  
  rho1_E_init <- (rho1_E_min + rho1_E_max) / 2
  rho1_C_init <- (rho1_C_min + rho1_C_max) / 2
  
  cat(sprintf("SW-CRT: Processing rho1_max = %.3f\n", rho1_max))
  
  # call MMD function for SW-CRT
  capture.output({
    result <- OD_CE_SW_CRT(
      alpha = alpha, trteff = trteff, lambda = lambda,
      sigma_E = sigma_E, sigma_C = sigma_C,
      rho0_E_min = rho0_E_min, rho0_E_max = rho0_E_max,
      rho1_E_min = rho1_E_min, rho1_E_max = rho1_E_max,
      rho0_C_min = rho0_C_min, rho0_C_max = rho0_C_max,
      rho1_C_min = rho1_C_min, rho1_C_max = rho1_C_max,
      rho0_EC_min = rho0_EC_min, rho0_EC_max = rho0_EC_max,
      rho1_EC_min = rho1_EC_min, rho1_EC_max = rho1_EC_max,
      rho2_EC_min = rho2_EC_min, rho2_EC_max = rho2_EC_max,
      rho0_E_init = rho0_E_init, rho1_E_init = rho1_E_init,
      rho0_C_init = rho0_C_init, rho1_C_init = rho1_C_init,
      rho0_EC_init = rho0_EC_init, rho1_EC_init = rho1_EC_init,
      rho2_EC_init = rho2_EC_init,
      J = J, L = L, c1 = c1, c2 = c2, B = B,
      max_I = 100, max_K = 200
    )
  })
  
  # extract SW-CRT results (integer only)
  if (!is.null(result$integer_estimates)) {
    results_swcrt <- rbind(results_swcrt, data.frame(
      rho1_max = rho1_max,
      K_int = result$integer_estimates$cluster_size,
      I_int = result$integer_estimates$num_clusters,
      RE_int = result$integer_estimates$relative_efficiency
    ))
    cat(sprintf("  SW-CRT MMD found: K_int=%d, I_int=%d, RE_int=%.4f\n", 
                result$integer_estimates$cluster_size, 
                result$integer_estimates$num_clusters, 
                result$integer_estimates$relative_efficiency))
  } else {
    cat("  No SW-CRT MMD found\n")
  }
}

cat("\nMMD analysis completed!\n")

# Check if results are empty
if (nrow(results_crxo) == 0 || nrow(results_pa) == 0 || nrow(results_swcrt) == 0) {
  stop("No MMD results found for one or more designs. Check ICC parameter ranges and initial values.")
}

# ============================================================================================================
# prepare combined data for plotting
# ============================================================================================================

# combine all data
df_combined <- rbind(
  # CRXO MMD
  data.frame(
    rho1_max = rep(results_crxo$rho1_max, 4),
    Value = c(results_crxo$I_int, results_crxo$I_dec, results_crxo$K_int, results_crxo$K_dec),
    Parameter = rep(c("I", "I", "K", "K"), each = nrow(results_crxo)),
    Type = rep(c("Integer", "Decimal", "Integer", "Decimal"), each = nrow(results_crxo)),
    Panel = "CRXO_MMD",
    Design = "CRXO"
  ),
  # CRXO RE
  data.frame(
    rho1_max = rep(results_crxo$rho1_max, 2),
    Value = c(results_crxo$RE_int, results_crxo$RE_dec),
    Parameter = rep("RE", 2 * nrow(results_crxo)),
    Type = rep(c("Integer", "Decimal"), each = nrow(results_crxo)),
    Panel = "CRXO_RE",
    Design = "CRXO"
  ),
  # PA-LCRT MMD
  data.frame(
    rho1_max = rep(results_pa$rho1_max, 4),
    Value = c(results_pa$I_int, results_pa$I_dec, results_pa$K_int, results_pa$K_dec),
    Parameter = rep(c("I", "I", "K", "K"), each = nrow(results_pa)),
    Type = rep(c("Integer", "Decimal", "Integer", "Decimal"), each = nrow(results_pa)),
    Panel = "PA_MMD",
    Design = "PA-LCRT"
  ),
  # PA-LCRT RE
  data.frame(
    rho1_max = rep(results_pa$rho1_max, 2),
    Value = c(results_pa$RE_int, results_pa$RE_dec),
    Parameter = rep("RE", 2 * nrow(results_pa)),
    Type = rep(c("Integer", "Decimal"), each = nrow(results_pa)),
    Panel = "PA_RE",
    Design = "PA-LCRT"
  ),
  # SW-CRT MMD (integer only)
  data.frame(
    rho1_max = rep(results_swcrt$rho1_max, 2),
    Value = c(results_swcrt$I_int, results_swcrt$K_int),
    Parameter = rep(c("I", "K"), each = nrow(results_swcrt)),
    Type = rep("Integer", 2 * nrow(results_swcrt)),
    Panel = "SWCRT_MMD",
    Design = "SW-CRT"
  ),
  # SW-CRT RE (integer only)
  data.frame(
    rho1_max = results_swcrt$rho1_max,
    Value = results_swcrt$RE_int,
    Parameter = rep("RE", nrow(results_swcrt)),
    Type = rep("Integer", nrow(results_swcrt)),
    Panel = "SWCRT_RE",
    Design = "SW-CRT"
  )
)

# ============================================================================================================
# create combined 3x2 plot
# ============================================================================================================

p_combined <- ggplot(df_combined, aes(x = rho1_max)) +
  # CRXO MMD lines
  geom_line(data = subset(df_combined, Panel == "CRXO_MMD" & Type == "Integer" & Parameter == "I"), 
            aes(y = Value, color = "I.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "CRXO_MMD" & Type == "Decimal" & Parameter == "I"), 
            aes(y = Value, color = "I.Decimal"), size = 1, linetype = "dotted") +
  geom_line(data = subset(df_combined, Panel == "CRXO_MMD" & Type == "Integer" & Parameter == "K"), 
            aes(y = Value, color = "K.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "CRXO_MMD" & Type == "Decimal" & Parameter == "K"), 
            aes(y = Value, color = "K.Decimal"), size = 1, linetype = "dotted") +
  # CRXO RE lines
  geom_line(data = subset(df_combined, Panel == "CRXO_RE" & Type == "Integer"), 
            aes(y = Value, color = "RE.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "CRXO_RE" & Type == "Decimal"), 
            aes(y = Value, color = "RE.Decimal"), size = 1, linetype = "dotted") +
  # PA-LCRT MMD lines
  geom_line(data = subset(df_combined, Panel == "PA_MMD" & Type == "Integer" & Parameter == "I"), 
            aes(y = Value, color = "I.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "PA_MMD" & Type == "Decimal" & Parameter == "I"), 
            aes(y = Value, color = "I.Decimal"), size = 1, linetype = "dotted") +
  geom_line(data = subset(df_combined, Panel == "PA_MMD" & Type == "Integer" & Parameter == "K"), 
            aes(y = Value, color = "K.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "PA_MMD" & Type == "Decimal" & Parameter == "K"), 
            aes(y = Value, color = "K.Decimal"), size = 1, linetype = "dotted") +
  # PA-LCRT RE lines
  geom_line(data = subset(df_combined, Panel == "PA_RE" & Type == "Integer"), 
            aes(y = Value, color = "RE.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "PA_RE" & Type == "Decimal"), 
            aes(y = Value, color = "RE.Decimal"), size = 1, linetype = "dotted") +
  # SW-CRT MMD lines (integer only)
  geom_line(data = subset(df_combined, Panel == "SWCRT_MMD" & Type == "Integer" & Parameter == "I"), 
            aes(y = Value, color = "I.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel == "SWCRT_MMD" & Type == "Integer" & Parameter == "K"), 
            aes(y = Value, color = "K.Integer"), size = 1) +
  # SW-CRT RE lines (integer only)
  geom_line(data = subset(df_combined, Panel == "SWCRT_RE" & Type == "Integer"), 
            aes(y = Value, color = "RE.Integer"), size = 1) +
  scale_color_manual(name = "",
                     values = c("I.Integer" = "#0072B2", 
                                "K.Integer" = "#D55E00",
                                "RE.Integer" = "blue",
                                "I.Decimal" = "#0072B2",
                                "K.Decimal" = "#D55E00",
                                "RE.Decimal" = "blue"),
                     labels = c("I.Integer" = expression(I[MMD]~(Integer)),
                                "K.Integer" = expression(K[MMD]~(Integer)),
                                "RE.Integer" = "RE (Integer)",
                                "I.Decimal" = expression(I[MMD]~(Decimal)),
                                "K.Decimal" = expression(K[MMD]~(Decimal)),
                                "RE.Decimal" = "RE (Decimal)"),
                     breaks = c("I.Integer", "K.Integer", "RE.Integer",
                                "I.Decimal", "K.Decimal", "RE.Decimal")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = expression(rho[paste("1,max")]^E), y = "") +
  scale_x_continuous(breaks = seq(0.025, 0.045, by = 0.005)) +
  facet_wrap(~Panel, ncol = 2, nrow = 3,
             labeller = labeller(Panel = c("CRXO_MMD" = "(a) MMDs for CRXO trials (J = 4)",
                                           "CRXO_RE" = "(b) RE of MMDs for CRXO trials (J = 4)",
                                           "PA_MMD" = "(c) MMDs for PA-LCRTs (J = 4)",
                                           "PA_RE" = "(d) RE of MMDs for PA-LCRTs (J = 4)",
                                           "SWCRT_MMD" = "(e) MMDs for SW-CRTs (J = 4, L = 3)",
                                           "SWCRT_RE" = "(f) RE of MMDs for SW-CRTs (J = 4, L = 3)")),
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
      Panel == "CRXO_MMD" ~ scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)),
      Panel == "CRXO_RE" ~ scale_y_continuous(limits = c(0.90, 1.0), breaks = seq(0.90, 1.00, by = 0.02), name = "RE"),
      Panel == "PA_MMD" ~ scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)),
      Panel == "PA_RE" ~ scale_y_continuous(limits = c(0.90, 1.0), breaks = seq(0.90, 1.00, by = 0.02), name = "RE"),
      Panel == "SWCRT_MMD" ~ scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)),
      Panel == "SWCRT_RE" ~ scale_y_continuous(limits = c(0.90, 1.0), breaks = seq(0.90, 1.00, by = 0.02), name = "RE")
    )
  )

library(patchwork)

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
p_left <- ggplot(subset(df_combined, Panel %in% c("CRXO_MMD", "PA_MMD", "SWCRT_MMD")), aes(x = rho1_max)) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_MMD", "PA_MMD", "SWCRT_MMD") & Type == "Integer" & Parameter == "I"),
            aes(y = Value, color = "I.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_MMD", "PA_MMD") & Type == "Decimal" & Parameter == "I"),
            aes(y = Value, color = "I.Decimal"), size = 1, linetype = "dotted") +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_MMD", "PA_MMD", "SWCRT_MMD") & Type == "Integer" & Parameter == "K"),
            aes(y = Value, color = "K.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_MMD", "PA_MMD") & Type == "Decimal" & Parameter == "K"),
            aes(y = Value, color = "K.Decimal"), size = 1, linetype = "dotted") +
  scale_color_manual(name = "",
                     values = c("I.Integer" = "#0072B2", "K.Integer" = "#D55E00",
                                "I.Decimal" = "#0072B2", "K.Decimal" = "#D55E00"),
                     labels = c("I.Integer" = expression(I[MMD]~(Integer)),
                                "K.Integer" = expression(K[MMD]~(Integer)),
                                "I.Decimal" = expression(I[MMD]~(Decimal)),
                                "K.Decimal" = expression(K[MMD]~(Decimal))),
                     breaks = c("I.Integer", "K.Integer", "I.Decimal", "K.Decimal")) +
  scale_x_continuous(breaks = seq(0.025, 0.045, by = 0.005)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  labs(x = expression(rho[paste("1,max")]^E),
       y = expression(paste("Optimal sample size combination (", I[MMD], ", ", K[MMD], ")"))) +
  facet_wrap(~Panel, ncol = 1,
             labeller = labeller(Panel = c("CRXO_MMD" = "(a) MMDs for CRXO trials (J = 4)",
                                           "PA_MMD" = "(c) MMDs for PA-LCRTs (J = 4)",
                                           "SWCRT_MMD" = "(e) MMDs for SW-CRTs (J = 4, L = 3)"))) +
  shared_theme

# right plot: RE
p_right <- ggplot(subset(df_combined, Panel %in% c("CRXO_RE", "PA_RE", "SWCRT_RE")), aes(x = rho1_max)) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_RE", "PA_RE", "SWCRT_RE") & Type == "Integer"),
            aes(y = Value, color = "RE.Integer"), size = 1) +
  geom_line(data = subset(df_combined, Panel %in% c("CRXO_RE", "PA_RE") & Type == "Decimal"),
            aes(y = Value, color = "RE.Decimal"), size = 1, linetype = "dotted") +
  scale_color_manual(name = "",
                     values = c("RE.Integer" = "blue", "RE.Decimal" = "blue"),
                     labels = c("RE.Integer" = "RE (Integer)", "RE.Decimal" = "RE (Decimal)"),
                     breaks = c("RE.Integer", "RE.Decimal")) +
  scale_x_continuous(breaks = seq(0.025, 0.045, by = 0.005)) +
  scale_y_continuous(limits = c(0.90, 1.0), breaks = seq(0.90, 1.00, by = 0.02)) +
  labs(x = expression(rho[paste("1,max")]^E), y = "RE") +
  facet_wrap(~Panel, ncol = 1,
             labeller = labeller(Panel = c("CRXO_RE" = "(b) RE of MMDs for CRXO trials (J = 4)",
                                           "PA_RE" = "(d) RE of MMDs for PA-LCRTs (J = 4)",
                                           "SWCRT_RE" = "(f) RE of MMDs for SW-CRTs (J = 4, L = 3)"))) +
  shared_theme

# combine and collect legends
p_combined <- p_left + p_right + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

print(p_combined)
ggsave("../figures/figure_4.pdf", p_combined, width = 12, height = 8)