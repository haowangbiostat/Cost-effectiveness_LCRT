# sw-crt incomplete designs for australian trial context

source("utils_SWCRT_incomplete.R")

# parameters from trial data
rho0_E <- 0.048
rho1_E <- 0.042
rho0_C <- 0.020
rho1_C <- 0.018
rho0_EC <- 0.007
rho1_EC <- 0.004
rho2_EC <- 0.75
sigma_E <- 6.48
sigma_C <- 11635
lambda <- 216
trteff <- 2089

c1 <- 3000
c2 <- 250
B <- 600000
max_I <- 100
max_K <- 200
L <- 7

# mmd icc ranges
rho0_EC_min <- 0; rho0_EC_max <- 0.01
rho1_EC_min <- 0; rho1_EC_max <- 0.005
rho2_EC_min <- 0.5; rho2_EC_max <- 0.8

# lod for J = 8
swcrt_lod_8 <- OD_CE_SW_CRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho1_E_min = rho1_E,
  rho0_C_min = rho0_C, rho1_C_min = rho1_C,
  rho0_EC_min = rho0_EC, rho1_EC_min = rho1_EC, rho2_EC_min = rho2_EC,
  J = 8, L = L, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(swcrt_lod_8$integer_estimates[names(swcrt_lod_8$integer_estimates) != "Z"])

# lod for J = 9
swcrt_lod_9 <- OD_CE_SW_CRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho1_E_min = rho1_E,
  rho0_C_min = rho0_C, rho1_C_min = rho1_C,
  rho0_EC_min = rho0_EC, rho1_EC_min = rho1_EC, rho2_EC_min = rho2_EC,
  J = 9, L = L, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(swcrt_lod_9$integer_estimates[names(swcrt_lod_9$integer_estimates) != "Z"])

# lod for J = 10
swcrt_lod_10 <- OD_CE_SW_CRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho1_E_min = rho1_E,
  rho0_C_min = rho0_C, rho1_C_min = rho1_C,
  rho0_EC_min = rho0_EC, rho1_EC_min = rho1_EC, rho2_EC_min = rho2_EC,
  J = 10, L = L, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(swcrt_lod_10$integer_estimates[names(swcrt_lod_10$integer_estimates) != "Z"])

# mmd for J = 8
swcrt_mmd_8 <- OD_CE_SW_CRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho0_E_max = rho0_E,
  rho1_E_min = rho1_E, rho1_E_max = rho1_E,
  rho0_C_min = rho0_C, rho0_C_max = rho0_C,
  rho1_C_min = rho1_C, rho1_C_max = rho1_C,
  rho0_EC_min = rho0_EC_min, rho0_EC_max = rho0_EC_max,
  rho1_EC_min = rho1_EC_min, rho1_EC_max = rho1_EC_max,
  rho2_EC_min = rho2_EC_min, rho2_EC_max = rho2_EC_max,
  J = 8, L = L, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(swcrt_mmd_8$integer_estimates[names(swcrt_mmd_8$integer_estimates) != "Z"])