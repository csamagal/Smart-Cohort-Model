# -----------------------------------------------------------------#
# SMART COHORT model
# Biomass function

# Function to calculate biomass related characteristics -- biomass is a function
# of elevation and depth. Includes vegetation switching and weighted
# parameterizations, and root:shoot calculations.

# Updated February 2019 - Anthony Rietl in Matlab
# Updated March 2019 - Megan Vahsen in R

# Load required packages
library(pracma)

# ----------------------------------------------------------------- #
biomass <- function(yr){

# ----------------- #
# Start annual loop #
# ----------------- #

#for (yr in 1000:tend){

# Calculate current depth and elevation  
d[yr] <- amp + msl[yr] - Z[yr]
e[yr] <- amp - d[yr]

# Identify maximum and minimum elevations
E_max <- c(max(Z_C3), max(Z_C4))
E_min <- c(min(Z_C3), min(Z_C4))

# Calculate shootmass for each species given depths
shootmass <- c(polyval(coef[1:3], d[yr]), polyval(coef[4:7], d[yr]))
# Keep shootmass from going negative. If it is negative for either species,
# force to 0.
if(shootmass[1] < 0){
  shootmass[1] <- 0 
}
if(shootmass[2] < 0){
  shootmass[2] <- 0
}

#####################
# Determine habitat # 
#####################

# There are 6 different scenarios: (1) marsh is too high for either species to
# survive, (2) marsh is too low for either species to survive (3) marsh is in a
# range that only supports the more flood-tolerant species (i.e. C3) (4) marsh
# is in a range that supports both species, but is more flooded, so
# flood-tolerant species dominates (5) marsh is in a range that supports both
# species, but is less flooded, so less flood-tolerant species dominates (6)
# marsh is in a range that only supports the less flood-tolerant species (i.e.
# C4).

#################
# START BIOMASS #
#################

# Scenario 1: marsh is too high for either species to survive
if (e[yr] > max(E_max)){
  kk <- 2 # C4; just a placeholder, because nothing is actually there
  tempabg <- 0
  tempspecies <- 2 # C4; just a placeholder, because nothing is actually there
  temprhizome <- 0
  Root_Shoot <- 0
  dWeight <- 0
  
# Scenario 2: marsh is too low for either species to survive
}else if(e[yr] < min(E_min)){
  kk <- 1 # C3; just a placeholder, because nothing is actually there
  tempabg <- 0
  tempspecies <- 1 # C3; just a placeholder because nothing is actually there
  temprhizome <- 0
  Root_Shoot <- 0
  dWeight <- 0
  
# Scenario 3: flood-tolerant species range (no less flood-tolerant species)
}else if(e[yr] < max(E_min) & e[yr] >= min(E_min)){
  kk <- 1 # C3 species code
  tempabg <- shootmass[kk] # Aboveground biomass
  tempspecies <- 1 # C3 species code
  temprhizome <- abg_rh[kk] * tempabg # Rhizome biomass
  Nup <- Vmax[kk] * N_Scalar[kk] # N-uptake rate
  Root_Shoot <- oN[kk] * npp[kk] / Nup # Root to shoot ratio
  dWeight <- 0 # Species weight; 0 because 100% C3

# Scenario 4: marsh is in a range that supports both species, but is more
# flooded, so flood-tolerant species dominates.
}else if(e[yr] >= max(E_min) & e[yr] <= amp - intersection){
  kk <- 3 # C3MD species code (C3 mixed dominant)
  tempspecies <- 3 # C3MD species code (C3 mixed dominant)

  # New weighting scheme
  c3_biomass <- polyval(coef[1:3], amp - e[yr])
  c4_biomass <- polyval(coef[4:7], amp - e[yr])
  dWeight <- c3_biomass / (c3_biomass + c4_biomass)
  
  # Update biomass values using weighting
  tempabg <- shootmass[1] * dWeight + shootmass[2] * (1-dWeight) # Calculate total aboveground biomass
  temprhizome <- (abg_rh[1]*dWeight + abg_rh[2]*(1 - dWeight))*tempabg # Calculate total rhizome biomass
  # Calculate root:shoot ratio based on N parameters using weighting
  Nup <- (Vmax[1]*dWeight + Vmax[2]*(1-dWeight)) * (N_Scalar[1]*dWeight + N_Scalar[2]*(1 - dWeight))
  Root_Shoot <- oN[kk] * (npp[1]*dWeight + npp[2]*(1-dWeight)) / Nup 
  # Calculate turnover of roots and rhizomes using weighting
  Ro_T[3] <- Ro_T[1] * dWeight + Ro_T[2] * (1 - dWeight)
  Rh_T[3] <- Rh_T[1] * dWeight + Rh_T[2] * (1 - dWeight)
  
# Scenario 5: marsh is in a range that supports both species, but is less
# flooded, so less flood-tolerant species dominates. This is the same as
# Scenario 5 except that weights reverse for species 1 or 2 when updating
# biomass, root:shoot, and turnover.
}else if (e[yr] >= amp - intersection & e[yr] <= min(E_max)){
  kk <- 4 # C4MD species code (C4 mixed dominant)
  tempspecies <- 4 # C4MD species code (C4 mixed dominant)
  
  # New weighting scheme
  c3_biomass <- polyval(coef[1:3], amp - e[yr])
  c4_biomass <- polyval(coef[4:7], amp - e[yr])
  dWeight <- c4_biomass / (c3_biomass + c4_biomass)
  
  # Update biomass using weighting
  tempabg <- shootmass[1]*(1 - dWeight) + shootmass[2]*dWeight
  temprhizome <- (abg_rh[1]*(1 - dWeight) + abg_rh[2]*dWeight)*tempabg
  # Calculate root:shoot ratio based on N parameters using weighting
  Nup <- (Vmax[1]*(1-dWeight) + Vmax[2]*dWeight) * (N_Scalar[1]*(1 - dWeight) + N_Scalar[2]*dWeight)
  Root_Shoot <- oN[kk] * (npp[1]*(1 - dWeight) + npp[2]*dWeight) / Nup
  # Calculate turnover of roots and rhizomes using weighting
  Ro_T[4] <- Ro_T[1] * (1 - dWeight) + Ro_T[2] * dWeight
  Rh_T[4] <- Rh_T[1] * (1 - dWeight) + Rh_T[2] * dWeight

# Scenario 6: less flood tolerant species range (no flood tolerant species)
}else if (e[yr] > min(E_max) & e[yr] <= max(E_max)){
  kk <- 2 # C4 species code
  tempspecies <- 2 # C4 species code
  tempabg <- shootmass[kk] # aboveground biomass
  temprhizome <- abg_rh[kk] * tempabg # rhizome biomass
  Nup <- Vmax[kk] * N_Scalar[kk] # N-uptake rate
  Root_Shoot <- oN[kk] * npp[kk] / Nup # Root to shoot ratio
  dWeight <- 0 # species weight, 0 because 100% C4
}

# Calculate root biomass
temproot <- Root_Shoot * tempabg # Total root biomass in a given year

return(list(tempspecies = tempspecies, temproot = temproot, temprhizome = temprhizome,
            tempabg = tempabg, dWeight = dWeight, kk = kk, D_min = D_min, D_max = D_max,
            Ro_T = Ro_T, Rh_T = Rh_T))

}
