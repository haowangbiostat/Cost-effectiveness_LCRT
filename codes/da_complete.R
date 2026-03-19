source("utils_PA.R")
source("utils_CRXO.R")
source("utils_SWCRT.R")

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
B <-600000
max_I <- 100
max_K <- 200

# mmd icc ranges
rho0_EC_min <- 0; rho0_EC_max <- 0.01
rho1_EC_min <- 0; rho1_EC_max <- 0.005
rho2_EC_min <- 0.5; rho2_EC_max <- 0.8

# lod for crxo
cat("=== crxo lod ===\n")
crxo_lod <- OD_CE_CRXO_MPD(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho1_E_min = rho1_E,
  rho0_C_min = rho0_C, rho1_C_min = rho1_C,
  rho0_EC_min = rho0_EC, rho1_EC_min = rho1_EC, rho2_EC_min = rho2_EC,
  J = 8, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(crxo_lod$integer_estimates)

# mmd for crxo
cat("\n=== crxo mmd ===\n")
crxo_mmd <- OD_CE_CRXO_MPD(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho0_E_max = rho0_E,
  rho1_E_min = rho1_E, rho1_E_max = rho1_E,
  rho0_C_min = rho0_C, rho0_C_max = rho0_C,
  rho1_C_min = rho1_C, rho1_C_max = rho1_C,
  rho0_EC_min = rho0_EC_min, rho0_EC_max = rho0_EC_max,
  rho1_EC_min = rho1_EC_min, rho1_EC_max = rho1_EC_max,
  rho2_EC_min = rho2_EC_min, rho2_EC_max = rho2_EC_max,
  J = 8, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(crxo_mmd$integer_estimates)

# lod for pa-lcrt
cat("\n=== pa-lcrt lod ===\n")
pa_lod <- OD_CE_PA_LCRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho1_E_min = rho1_E,
  rho0_C_min = rho0_C, rho1_C_min = rho1_C,
  rho0_EC_min = rho0_EC, rho1_EC_min = rho1_EC, rho2_EC_min = rho2_EC,
  J = 8, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(pa_lod$integer_estimates)

# mmd for pa-lcrt
cat("\n=== pa-lcrt mmd ===\n")
pa_mmd <- OD_CE_PA_LCRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho0_E_max = rho0_E,
  rho1_E_min = rho1_E, rho1_E_max = rho1_E,
  rho0_C_min = rho0_C, rho0_C_max = rho0_C,
  rho1_C_min = rho1_C, rho1_C_max = rho1_C,
  rho0_EC_min = rho0_EC_min, rho0_EC_max = rho0_EC_max,
  rho1_EC_min = rho1_EC_min, rho1_EC_max = rho1_EC_max,
  rho2_EC_min = rho2_EC_min, rho2_EC_max = rho2_EC_max,
  J = 8, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(pa_mmd$integer_estimates)

# lod for sw-crt (J = L + 1 = 8, L = 7)
swcrt_lod <- OD_CE_SW_CRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho1_E_min = rho1_E,
  rho0_C_min = rho0_C, rho1_C_min = rho1_C,
  rho0_EC_min = rho0_EC, rho1_EC_min = rho1_EC, rho2_EC_min = rho2_EC,
  J = 8, L = 7, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(swcrt_lod$integer_estimates)

# lod for sw-crt (J = L + 2 = 9, L = 7)
swcrt_lod <- OD_CE_SW_CRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho1_E_min = rho1_E,
  rho0_C_min = rho0_C, rho1_C_min = rho1_C,
  rho0_EC_min = rho0_EC, rho1_EC_min = rho1_EC, rho2_EC_min = rho2_EC,
  J = 9, L = 7, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(swcrt_lod$integer_estimates)

# lod for sw-crt (J = L + 3 = 10, L = 7)
swcrt_lod <- OD_CE_SW_CRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho1_E_min = rho1_E,
  rho0_C_min = rho0_C, rho1_C_min = rho1_C,
  rho0_EC_min = rho0_EC, rho1_EC_min = rho1_EC, rho2_EC_min = rho2_EC,
  J = 10, L = 7, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(swcrt_lod$integer_estimates)

# mmd for sw-crt (J = L + 1 = 9, L = 7)
swcrt_mmd <- OD_CE_SW_CRT(
  trteff = trteff, lambda = lambda, sigma_E = sigma_E, sigma_C = sigma_C,
  rho0_E_min = rho0_E, rho0_E_max = rho0_E,
  rho1_E_min = rho1_E, rho1_E_max = rho1_E,
  rho0_C_min = rho0_C, rho0_C_max = rho0_C,
  rho1_C_min = rho1_C, rho1_C_max = rho1_C,
  rho0_EC_min = rho0_EC_min, rho0_EC_max = rho0_EC_max,
  rho1_EC_min = rho1_EC_min, rho1_EC_max = rho1_EC_max,
  rho2_EC_min = rho2_EC_min, rho2_EC_max = rho2_EC_max,
  J = 8, L = 7, c1 = c1, c2 = c2, B = B, max_I = max_I, max_K = max_K
)
print(swcrt_mmd$integer_estimates)