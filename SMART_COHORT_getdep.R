# -----------------------------------------------------------------#
# SMART COHORT model
# Deposition function

# Function to calculate the amount of sediment trapping and settling for a tide
# range and suspended sediment concentration (C, mg/l OR g/m3). D Alpaos A,
# Lanzoni S, Marani M, Rinaldo A (2007) Landscape evolution in tidal embayments:
# Modeling the interplay of erosion, sedimentation, and vegetation dynamics. J
# Geophys Res Earth Surf 112:1-17.

# Updated February 2019 - Anthony Rietl in Matlab
# Updated February 2019 - Megan Vahsen in R

# Load required packages
library(pracma)

# ----------------------------------------------------------------- #

getdep <- function(yr){
# We calculate how much accretion there is for every tidal period. Then we
# multiply up this to the annual timestep to work well with the rest of the
# model.

# ----------------------------------------------------------------- #

# ------------------------- #
# Initiate variable storage #
# ------------------------- #
C <- rep(0, numiterations)
qs <- rep(0, numiterations)
depth_dep <- rep(0, numiterations)
flooding_dur <- rep(0, tend)

# Iterate through single tidal cycle
for (t in 2:numiterations) {
  depth_dep[t] <- amp * sin(2 * pi * (t * dt / P - 0.25)) + (msl[yr] - Z[yr]) # Depth of water over marsh surface, as a fxn of tidal stage
  
  # If tide is out, depth is 0. Basically, if tide is out, we aren't interested
  # in calculating the sediment in the water above the marsh surface and how it
  # is trapped, because the water doesn't exist.
  
  if (depth_dep[t] < 0) {
    depth_dep[t] <- 0
  } else if (depth_dep[t] > 0) {
    flooding_dur[yr] <- flooding_dur[yr] + dt
    
    # Here, ws * Ci * flooding_dur = estimated sediment deposited in a tidal
    # cycle; multiply by the number of tidal cycles in a year to get annual
    # sediment deposition
    
    # Change in water level (above marsh surface) from previous timestep
    dh <- depth_dep[t] - depth_dep[t - 1]
    
    # If there is more water present, we will add sediment suspended in the extra water
    if (dh > 0) {
      mi <- C[t - 1] * depth_dep[t - 1] #initial mass of sediment in the water column above the marsh platform
      dm <- -qs[t - 1] + Ci * dh # change in mass; mass balance between deposition (previous timestep) and sediment coming in with tide (SSC).
      # Sediment coming in from tide always has the same concentration!!
      if (dm + mi < 0) {
        # Cannot remove more sediment from water column than what is there
        dm <- -mi # to keep concentration from becoming negative, set it to zero in this case. This doesn't look like it is setting to 0...?
      }
      C[t] <- (mi + dm) / depth_dep[t] # New suspended sediment concentration is equal to the initial mass of the sediment plus the change
      # in mass of sediment divided by the depth of the water column
      qs[t] <- C[t] * ws * dt # sediment settling (Marani et al. 2007)
    } else{
      mi <- C[t - 1] * depth_dep[t - 1] # initial mass of sediment in water column
      dm <- -qs[t - 1] + C[t - 1] * dh # There's no new water coming in because dh is not greater than 0. When dh < 0, tide is going out.
      # You loose mass given the SSC at the previous timestep.
      # Change in mass: mass balance between deposition and sediment going out with tide
      if (dm + mi < 0) {
        # cannot remove more sediment from the water column than what is there
        dm <- -mi # To keep concentration from becoming negative, set it to zero in this case; this isn't actually setting to 0???
      }
      C[t] <- (mi + dm) / depth_dep[t] # new suspended sediment concentration
      qs[t] <- C[t] * ws * dt # sediment settling (Marani et al. 2007)
    }
  }
}

# Calculate annual deposition rates from trapping and settling
dep <- sum(qs) * timestep # Total deposition: annual settling mass flux (g/m2/yr)
dep_slow <- COrg * dep * depOrg_frac[1] # Mass flux of deposition that is organic (5%) times proportion of organic that is slow (recalcitrant) organic C (10%)
dep_fast <- COrg * dep * depOrg_frac[2] # Mass flux of deposition that is organic (5%) times proportion of organic that is fast (labile) organic C (90%)
dep_min <- dep - (dep_slow + dep_fast) # Mass flux of deposition that is inorganic

# if (abs(dep - (dep_slow + dep_fast + dep_min)) < 1){
#   set <- "Running year"
#   print(paste(set, yr))
# }else{
#   print("check dep calc")
# }

return(list(dep = dep, dep_slow = dep_slow, dep_fast = dep_fast, dep_min = dep_min, flooding_dur = flooding_dur,
            depth_dep = depth_dep))

}
