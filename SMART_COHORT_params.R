# -----------------------------------------------------------------#
# SMART COHORT model
# Parameters and storage

# File with all static parameter values and initialization of all storage

# Updated March 2019 - Megan Vahsen in R

# Load required packages
library(pracma)
install.packages("pracma")

# ----------------------------------------------------------------- #

# --------------------- #
# Deposition parameters #
# --------------------- #
tr <- 0.44 # tidal range
amp <- tr/2
numiterations <- 500
P <- 12.5 * 3600 # tidal period in seconds; A tidal period is 12.5 hours long
dt <- P/numiterations # this is just to smooth out the tidal period??
timestep <- 365 * (24/12.5) # tidal cycles per year; number to multiple accretion simulated over a tidal cycle by
# Here, it takes ~ 12.5 hours to get from one high tide (or low tide) to the
# next. There are 700.8 tidal periods ever year because we shift from one period
# to the next every 12.5 hours.

# Storage variables
dep_alloc_slow <- NULL
dep_alloc_fast <- NULL
dep_alloc_min <- NULL


############
# Sediment #
############
Ci <- 5 # concentration (mg/L) of suspended sediment (mineral + organic), all supsended concentration
COrg <- 0.05 # g/g organic content of suspended sediment Ci
depOrg_frac <- c(0.1, 0.9) # fast and slow sediment OM fractions
ws <- 2e-4 # settling velocity (D'Alpaos et al. 2007)
rhos <- 1.99e6 # bulk density of mineral matter (g/m^3) (i.e. sediment); from Morris et al. 2016
rhoo <- 8.5e4 # bulk density of organic matter (g/m^3); from Morris et al. 2016


# ------------------ #
# Biomass parameters #
# ------------------ #

#######
# CO2 #
#######
# Turn on or off CO2 (1 = on; 0 = off)
CO2 <- 1

##################
# Tidal and time #
##################
tr <- 0.44 # tidal range
amp <- tr/2 # tidal amplitude
R0 <- 0.001 # initial rate of sea-level rise in m/yr
RSLRi <- R0 * 10^3 # initial rate of sea level rise in mm/yr 
Z_end <- 0.8 # starting elevation (m); user defined to being marsh building wherever you want in vegetation growth range
tend <- 2600 # end of time (years)
#yr <- 1000:tend # set year vector

###############
# C3 specific #
###############
Z_C3 <- c(0.28, 0.16, 0.05) # elevation profile for C3 (3 entries)
ABG_C3 <- c(0, 600, 0) # aboveground biomass for C3 given Z_C3 elevations (3 entries)
eCO2_C3 <- c(1, 1.3, 1) # multiplicative effect of CO2 on biomass (3 entries)
# eCO_C3 based on Chunwu Zhu's observations at SERC (e.g. ~20% increase in abg
# at optimum elevation)

###############
# C4 specific #
###############
Z_C4 <- c(0.45, 0.38, 0.27, 0.12) # elevation profile for C4 (4 entries)
ABG_C4 <- c(0, 500, 700, 0) # aboveground biomass for C4 given Z_C3 elevations (4 entries)
eCO2_C4 <- rep(1,4) # multiplicative effect of CO2 on biomass (4 entries)

################
# C3, C4 joint #
################
# No CO2
npp <- c(0.175, 0.175) # proportional change in net carbon uptake (gC/g leaf/day)
abg_rh <- c(1,1) # aboveground to rhizome biomass ratio
N_Scalar <- c(0.6, 0.6)
Vmax <- c(3.6e-3, 5.3e-3) # maximum nitrogen uptake velocity
oN <- c(0.02, 0.02, 0.02, 0.02) # oN; optimal [N] in tissue (fraction); (C3, C4, mixed - C3 dom, mixed - C4 dom)

# Belowground parameters
root_ash <- c(0.08, 0.08, 0.08, 0.08) # fraction of root biomass that is mineral; (C3, C4, C3MD, C4MD)
rhizome_ash <- c(0.08, 0.08, 0.08, 0.08) # fraction of rhizome biomass that is mineral; (C3, C4, C3MD, C4MD)
Ro_T <- c(1.1, 1.6) # root turnover rate; number of times root biomass diesback and regrows
Rh_T <- c(0.5, 1) # rhizome turnover rate; number of times rhizome biomass diesback and regrows
bioOrg_frac <- c(0.84, 0.16) # fast and slow sediment OM fractions
max_depth_bgb <- 1 # maximum rooting depth

# Multiplicative CO2 effects
npp_eCO2 <- c(1.2, 1) # proportional change in net carbon uptake (gC/g leaf/day)
abg_rh_eCO2 <- c(2,1) # multiplier for rhizome biomass in relation to aboveground biomass
N_Scalar_eCO2 <- c(1,1)

# Storage
Root <- NULL
Rhizome <- NULL
abg <- NULL
bgb <- NULL
spp_weight <- NULL
species <- NULL

# ----------------------------------------------------------------- #

# ------------------------ #
# Preliminary calculations #
# ------------------------ #

# Preallocate variables for pre-existing stratigraphy

# Calculate mean sea-level (starts increasing at time = 1000)
msl <- rep(0, tend)
for (i in 1000:tend){ # creates msl matric relative to your R0
  msl[i] <- msl[i-1] + R0
}

layer_depth <- rep(1e-3, 1000) # Defines thickness of each layer (0.001 m or 1mm)

# Calculate maximum and minimum depths below tidal amplitude. Negative depths
# mean that the elevation of the marsh is above the height of the tidal
# amplitude.
D_max <- amp - c(min(Z_C3), min(Z_C4))
D_min <- amp - c(max(Z_C3), max(Z_C4))
E_max <- c(max(Z_C3), max(Z_C4))
E_min <- c(min(Z_C3), min(Z_C4))

# Calculate coefficients of polynomials for C3 and C4 by depth below tidal
# amplitude
coef <- c(polyfit((amp - Z_C3), ABG_C3, 2), polyfit((amp - Z_C4), ABG_C4, 3))

# If CO2 is elevated, apply eCO2 multipliers. These will only influence C3 as
# all C4 multipliers are 1.
if (CO2 == 1){
  npp <- npp * npp_eCO2
  abg_rh <- abg_rh * abg_rh_eCO2
  N_Scalar <- N_Scalar * N_Scalar_eCO2
  coef <- c(polyfit((amp - Z_C3), ABG_C3 * eCO2_C3, 2), polyfit((amp - Z_C4), ABG_C4 * eCO2_C4, 3))
}

############################
# Build marsh stratigraphy #
############################
Z_start = Z_end - (1000 * 0.001) # builds stratigraphy; Z increment must be the same as layer_depth
# Create empty columns for soil column, mineral accretion, accretion, and elevation Z
soilcolumn <- rep(0, tend)
minacc <- rep(0, tend)
accretion <- rep(0, tend)
Z <- rep(0, tend)

# Spin-up elevation to time = 1000
for (i in 1:1000){
  Z[i] <- Z_start + i*0.001
  soilcolumn[i] <- Z[i] - Z_start
  if (i > 1){
    minacc[i] <- Z[i] - Z[i-1]
    accretion[i] <- Z[i] - Z[i-1]
    # These are identical... do we need both?
  }
}

# Create vectors to hold depth and elevation values over time
d <- rep(0, tend)
e <- rep(0, tend)

###############################################################
# Evaluate the intersection of the two biomass x depth curves #
############################################################### 
# Create sequences of depths from minimum to maximum. Remember minimum depth ==
# maximum elevation
dd <- seq(min(D_min), max(D_max), 0.0001)
# Create functions to evaluate y values given coefficients and depth vector (f1
# and f2). Also create a function to evaluate differences between f1 and f2.
f1 <- function(dd){
  return(polyval(coef[1:3], dd))
}
f2 <- function(dd){
  return(polyval(coef[4:7], dd))
}
f <- function(dd){
  return(f1(dd) - f2(dd))
}

# Identify where C3 is greater than C4
t <- f(dd) > 0
# Figure out where there is a switch between C3 and C4 being higher. This
# 'switch' is where the intersection is.
i0 <- which(diff(t) != 0)
# Create a matrix of the identified 'switch' and the next step in the sequence
i0 <- matrix(c(i0, i0+1), nrow = 2, ncol = 1)
# Identify the dimensions of the matrix. If there are multiple intersections,
# there could be multiple 'switch' points
n <- dim(i0)[2]
# Find intersection (x-axis; depth (m)) between two curves and store.
intersection <- NULL
for (jj in 1:n){
  intersection[jj] <- fzero(f, dd[i0[,jj]])$x
}

# ----------------------------------------------------------------- #

# ------------------------ #
# Decomposition parameters #
# ------------------------ #
g <- 0.27 # "depth over which root growth is reduced by a factor of e^{-1}" (Mudd et al. 2009)
# NOTE: new documentation by Kirwan says the depth at which root growth is
# reduced by 1/3.. I don't think this makes sense and this doesn't match that
# paper. Although e^-1 = 0.367 so it is close.

layer_depth <- rep(0.001, 2600)
rhos <- 1.99e6
# -----------------------------------------------------------------#


# Initialize
liveroot <- rep(0, tend)
liverhizome <- rep(0, tend)
liveroot_Turnover <- rep(0, tend)
liverhizome_Turnover <- rep(0, tend)
total_Turnover <- rep(0, tend)
root_min_dep <- rep(0, tend)
rhizome_min_dep <- rep(0, tend)
tempbgb_min <- rep(0, tend)
tempbgb_fast <- rep(0, tend)
tempbgb_slow <- rep(0, tend)
bgb_org_dep <- rep(0, tend)
tempal_slow <- rep(0, tend)
tempal_fast <- rep(0, tend)
tempal_min <- rep(0, tend)
mki_fast <- rep(0, tend)
mki_slow <- rep(0, tend)
decomp_al_fast <- rep(0, tend)
al_fast <- rep(0, tend)
decomp_al_slow <- rep(0, tend)
al_slow <- rep(0, tend)
decomp_bgb_fast <- rep(0, tend)
bgb_fast <- rep(0, tend)
decomp_bgb_slow <- rep(0, tend)
bgb_slow <- rep(0, tend)
bgb_min <- rep(0, tend)
organic <- rep(0, tend)
mineral <- rep(0, tend)
loi <- rep(0, tend)
perc_C <- rep(0, tend)
layer_C <- rep(0, tend)

Fout_C <- rep(0, tend)
column_C <- rep(0, tend)
C_accum <- rep(0, tend)
torgacc <- rep(0, tend)
orgacc <- rep(0, tend)
org_in <- rep(0, tend)

# mineral mass in each layer is eqivalent to a layer depth of 1mm/yr
al_min1 <- rep(rhos * 0.001, 1000) 
al_min2 <- rep(0, tend - 1000)
al_min <- c(al_min1, al_min2)

# Calculate decay rate at soil surface based on aboveground biomass
mkiref_fast <- 0.027 # reference/average decomposition constant of labile organic matter
mkiref_slow <- 0.0027 # reference/average decomposition constant of recalcitrant organic matter

mki_decrease <- 0.0001 # defines the lowest decomposition rate (i.e. decomposition rate is never = 0)
x <- seq(127, 900, length.out = 1000) # range of biomass values to fit linear relationship - lowest biomass defined by Mueller et al. 2015
# This means that when biomass is below 127, k = k_decrease

Bref <- 500 # Reference/average aboveground biomass in a given year - i.e. the biomass at which you find your kref value

mki_increase_fast <- mkiref_fast * 2 # The maximum amount by which kref can increase - doubling estimated by Mueller et al. 2015
mki_increase_slow <- mkiref_slow * 2
mkiFast_range <- seq(mki_decrease, mki_increase_fast, length.out = 1000) # Full range of k from k_increase to k_decrease (lowest to highest possible values)

# Fit linear relationship between biomass and decomposition rate (higher biomass means higher decomposition rate)
newx <- mki_increase_fast*x/Bref
# Fit linear model without intercept to get value of a_coef
f_fast <- lm(mkiFast_range ~ -1 + newx) 
a_coef <- as.numeric(coef(f_fast)) # coef of the relationship between biomass and decomp rate; used in calculating k for a given abg

# Calculate coefs for distributin k with depth in the soil column
x2 <- seq(1, 0, length.out = 1000) # defines 1000 points between 0 and 1 - refers to soil column and max rooting depth (i.e. 0 to 1m)
expfunc <- function(a, b, x){mki_increase_fast * a * exp(b*x)} # function for fitting the range of k values to depth - exponential disitribution to max rooting depth
f_depth <- nls(mkiFast_range ~ expfunc(a, b, x = x2), start = list(a = 1, b = 1))
coef_depth <- as.numeric(coef(f_depth)) # coef values used to calculate k at a given depth
a_depth <- coef_depth[1]
b_depth <- coef_depth[2]

