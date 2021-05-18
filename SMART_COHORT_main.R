# -----------------------------------------------------------------#
# SMART COHORT model
# Mother function

# A dynamic point-based soil-cohort model for the moprhological evolution of a
# marsh soil column. The model captures the ecogeomorphic feedbacks between
# flooding, organic matter accumulation, sediment deposition, and marsh surface
# elevation under scenarios of sea level rise and elevated CO2. Changes in
# primary productivity and sedimentation are based elevation above MSL, and the
# decay rate of organic matter scales with aboveground biomass and exponentially
# declines with depth. The model incorporates the differential responses of C3
# and C4 plants, allowing for the switching between parameterizations for each
# plant type, creating weighted mixed communities as defined by the overlapping
# ranges of each species along an elevation continuum, and incorporating the
# empirically derived relationship between C3 productivity and elevated CO2.

# Updated February 2019 - Anthony Rietl in Matlab
# Updated March 2019 - Megan Vahsen in R

# Load required packages
library(pracma)

# ----------------------------------------------------------------- #
rm(list = ls())
# Load parameters and storage vectors
source("code/SMART_COHORT/SMART_COHORT_params.R")

# Load deposition, biomass and decomposition functions
source("code/SMART_COHORT/SMART_COHORT_getdep.R")
source("code/SMART_COHORT/SMART_COHORT_biomass.R")
source("code/SMART_COHORT/SMART_COHORT_decompose.R")

tend <- 1200

# Set timer
old <- Sys.time()

for (yr in 1000:tend){

# temporarily set year
#yr <- 1000

# ---------- #
# Deposition #
# ---------- #

# Execute function
dep_out <- getdep(yr)
# Store outputs
dep_alloc_slow[yr] <- dep_out$dep_slow
dep_alloc_fast[yr] <- dep_out$dep_fast
dep_alloc_min[yr] <- dep_out$dep_min
depth <- dep_out$depth
#plot(depth)
# ------- #
# Biomass #
# ------- #

# Execute function
biomass_out <- biomass(yr)
# Store outputs
species[yr] <- biomass_out$tempspecies
Root[yr] <- biomass_out$temproot
Rhizome[yr] <- biomass_out$temprhizome
abg[yr] <- biomass_out$tempabg
bgb[yr] <- biomass_out$temproot + biomass_out$temprhizome
spp_weight[yr] <- biomass_out$dWeight

# Indices for decomposition function
kk <- biomass_out$kk
Ro_T <- biomass_out$Ro_T
Rh_T <- biomass_out$Rh_T

# ------------- #
# Decomposition #
# ------------- #

# Execute function
decomp_out <- decomp(yr = yr, kk = kk, temproot = Root[yr], temprhizome = Rhizome[yr], tempabg = abg[yr],
                     dep_slow = dep_alloc_slow[yr], dep_fast = dep_alloc_fast[yr], dep_min = dep_alloc_min[yr],
                     Ro_T = Ro_T, Rh_T = Rh_T)
# Store outputs
Fout_C[yr] <- decomp_out$tempFout_C
soilcolumn[yr] <- decomp_out$tempsoilcolumn
column_C[yr] <- decomp_out$tempcolumn_C
torgacc[yr] <- decomp_out$temporgacc
minacc[yr] <- decomp_out$tempminacc
org_in[yr] <- decomp_out$temporgin
orgacc[yr] <- torgacc[yr] - torgacc[yr - 1]

# ----------------------- #
# Accretion and elevation #
# ----------------------- #

accretion[yr] <- soilcolumn[yr] - soilcolumn[yr - 1]
C_accum[yr] <- column_C[yr] - column_C[yr - 1] # Carbon accumulation rate (g/m^2/yr)
Z[yr + 1] <- Z[yr] + accretion[yr]

# Print progress
print(yr)
}

new <- Sys.time() - old

# Summary plots
#png("~/SMART_COHORT/Figs/output_plot_18Nov.png", height = 1800, width = 3000, res = 300)
par(mfrow = c(2,1), mar = c(5,5,2,2))
colors <- RColorBrewer::brewer.pal(4, "Set2")
colors_plot <- colors[species[1000:tend]]
plot(Z[1001:tend], type = "b", xlab = "year", ylab = "elevation (m)", pch = 16, col = colors_plot)
legend("bottomleft",legend = c("C3", "C4", "C3 > C4", "C3 < C4"), pch = 16, col = colors, bty = "n")
plot(C_accum[1001:tend], type = "l", xlab = "year", ylab = "C accumulation")
#dev.off()

Z[tend]
C_accum[tend]
plot(torgacc[1005:1050], pch = 16)

print(new)
