# -----------------------------------------------------------------#
# SMART COHORT model
# Decomposition function

# Function to calculate all decomposition related parameters. Distributes live
# roots exponentially down the soil column, calculates OM in each cohort,
# includes total turnover and depth dependent decay rates  

# Updated February 2019 - Anthony Rietl in Matlab
# Updated February 2019 - Megan Vahsen in R

# Load required packages
library(pracma)

# -----------------------------------------------------------------#

decomp <- function(yr, kk, tempabg, temproot, temprhizome, dep_slow, dep_fast, dep_min, Ro_T, Rh_T){

# Preliminaries
layer_depth[yr] <- layer_depth[yr - 1] # Layer depth must be designated before calculating true layer depth, thus approximated by last year's value
bro <- temproot / g # Root biomass at the surface (Mudd et al. 2009)
brh <- temprhizome / g # Rhizome biomass at the surface (Mudd et al. 2009)
mkitop_fast <- mki_increase_fast*a_coef*(tempabg/Bref)
mkitop_slow <- mki_increase_slow*a_coef*(tempabg/Bref)

if (mkitop_fast < mki_decrease){# Doesn't allow ktop to go below the value set in k_decrease
  mkitop_fast <- mki_decrease
}

if (mkitop_slow < mki_decrease){
  mkitop_slow <- mki_decrease
}

# Add and decompose material in each layer (cohort) of the soil column
depth <- rep(0, length.out = length(seq(yr, 1, -1)))
for (cohort in seq(yr, 1, -1)){
  # if(cohort == yr){
  #   depth[cohort] <- 0
  # }else{
  depth[cohort] <- sum(layer_depth[(cohort + 1):yr])
  # }
  
  if(depth[cohort] < max_depth_bgb){
    topdepth <- depth[cohort]
    bottomdepth <- depth[cohort] + layer_depth[cohort]
    
    fun1 <- function(depth){bro * exp(-depth/g)} # Defines the function that distributes root biomass exponentially with depth in the soil profile
    fun2 <- function(depth){brh * exp(-depth/g)}# Defines the function that distributes rhizome biomass exponentially with depth in the soil profile
    
    liveroot[cohort] <- integral(fun1, topdepth, bottomdepth) # Integrates under the depth distribution curve for each cohort to calculate root biomass
    liverhizome[cohort] <- integral(fun2, topdepth, bottomdepth)# Integrates under the depth distribution curve for each cohort to calculate rhizome biomass
    
    liveroot_Turnover[cohort] <- liveroot[cohort] * Ro_T[kk] # turnover: dead material only
    liverhizome_Turnover[cohort] <- liverhizome[cohort] * Rh_T[kk]
    total_Turnover[cohort] <- liveroot_Turnover[cohort] + liverhizome_Turnover[cohort]
    
    root_min_dep[cohort] <- liveroot_Turnover[cohort] * root_ash[kk]
    rhizome_min_dep[cohort] <- liverhizome_Turnover[cohort] * rhizome_ash[kk]
    tempbgb_min[cohort] <- root_min_dep[cohort] + rhizome_min_dep[cohort] # mineral deposition in the soil column from fraction of biomass that is mineral
    
    bgb_org_dep[cohort] <- total_Turnover[cohort] - tempbgb_min[cohort] # Calculates total organic inputs - turnover minus fraction that is mineral 
    tempbgb_fast[cohort] <- bgb_org_dep[cohort] * bioOrg_frac[1] # Splits total into fast and slow fractions
    tempbgb_slow[cohort] <- bgb_org_dep[cohort] * bioOrg_frac[2]

# If the cohort is on the surface, allocthonous material from suspended sediment is added    
    if (depth[cohort] == 0){
      tempal_slow[cohort] <- dep_slow
      tempal_fast[cohort] <- dep_fast
      tempal_min[cohort] <- dep_min
    }else{
      tempal_slow[cohort] <- 0
      tempal_fast[cohort] <- 0
      tempal_min[cohort] <- 0
    }
    # Determine depth-dependent decomposition parameters
    mki_fast[cohort] <- mkitop_fast * a_depth * exp(b_depth * depth[cohort])
    mki_slow[cohort] <- mkitop_slow * a_depth * exp(b_depth * depth[cohort])
  }else{
    tempbgb_min[cohort] <- 0
    tempbgb_fast[cohort] <- 0
    tempbgb_slow[cohort] <- 0
    bgb_org_dep[cohort] <- 0
    
    mki_fast[cohort] <- 0
    mki_slow[cohort] <- 0
  }
  # Decompose the marsh sediment
  if (depth[cohort] <= max_depth_bgb){
    # Proportion of decomposition in allochthonous fast in each cohort
    decomp_al_fast[cohort] <- (al_fast[cohort] + tempal_fast[cohort]) * mki_fast[cohort]
    # Allochthonous fast in each cohort (total - decomp)
    al_fast[cohort] <- al_fast[cohort] + tempal_fast[cohort] - decomp_al_fast[cohort]
    # Proportion of decomposition in allochthonous slow in each cohort
    decomp_al_slow[cohort] <- (al_slow[cohort] + tempal_slow[cohort]) * mki_slow[cohort]
    # Allochthonous slow in each cohort (total - decomp)
    al_slow[cohort] <- al_slow[cohort] + tempal_slow[cohort] - decomp_al_slow[cohort]
    # Proportion of decomposition in bgb fast in each cohort 
    decomp_bgb_fast[cohort] <- (bgb_fast[cohort] + tempbgb_fast[cohort]) * mki_fast[cohort]
    # bgb fast in each cohort (total - decomp)
    bgb_fast[cohort] <- bgb_fast[cohort] + tempbgb_fast[cohort] - decomp_bgb_fast[cohort]
    # Proportion of decomposition in bgb slow in each cohort
    decomp_bgb_slow[cohort] <- (bgb_slow[cohort] + tempbgb_slow[cohort]) * mki_slow[cohort]
    # bgb slow in each cohort (total - decomp)
    bgb_slow[cohort] <- bgb_slow[cohort] + tempbgb_slow[cohort] - decomp_bgb_slow[cohort]
    
    # update components
    al_min[cohort] <- al_min[cohort] + tempal_min[cohort] # Allochthonous mineral
    bgb_min[cohort] <- bgb_min[cohort] + tempbgb_min[cohort] # bgb mineral
    
    # Calculate organic and inorganic contributions to cohort
    organic[cohort] <- al_fast[cohort] + al_slow[cohort] + bgb_fast[cohort] + bgb_slow[cohort]
    mineral[cohort] <- al_min[cohort] + bgb_min[cohort]
    # Loss-on-ignition: proportion organic out of total (organic + mineral)
    loi[cohort] <- (organic[cohort] + liveroot[cohort] + liverhizome[cohort]) / 
      (organic[cohort] + liveroot[cohort] + liverhizome[cohort] + mineral[cohort])
    
    perc_C[cohort] <- 0.4 * loi[cohort] + 0.0025 * loi[cohort]^2 # Percent caron in each cohort
    layer_C[cohort] <- organic[cohort] * perc_C[cohort] # Total carbon

# It's possible under conditions of little to no mineral inputs for layers to be
# extremely thin - so thin that it falls below matlabs lowest recognized digit
# (realmin) the following keeps this number from being below realmin
  
    # Define lowest floating number
    realmin <- .Machine$double.xmin
    
    if ((organic[cohort]/rhoo) + (mineral[cohort]/rhos) <= realmin * 2){
      layer_depth[cohort] <- realmin
    }else{
      layer_depth[cohort] <- (organic[cohort]/rhoo + mineral[cohort]/rhos)
    }
  }
}

temporgin <- sum(bgb_org_dep[1:yr] + dep_fast + dep_slow) # Total organic in
temporgacc <- (sum(organic[1:yr]) - organic[yr-1])/rhoo # Organic accretion
tempminacc <- (sum(mineral[1:yr]) - sum(mineral[1:yr-1]))/rhos
tempsoilcolumn <- sum(layer_depth[1:yr]) # Soil column depth
tempcolumn_C <- sum(layer_C[1:yr]) # Total Carbon in column
tempFout_C <- sum(decomp_al_fast[1:yr]) + sum(decomp_al_slow[1:yr]) + sum(decomp_bgb_fast[1:yr]) + sum(decomp_bgb_slow[1:yr]) # Total decomposition

return(list(temporgin = temporgin, temporgacc = temporgacc, tempminacc = tempminacc,
            tempsoilcolumn = tempsoilcolumn, tempcolumn_C = tempcolumn_C, tempFout_C = tempFout_C,
            dep_slow = dep_alloc_slow[yr], dep_fast = dep_alloc_fast[yr], dep_min = dep_alloc_min[yr],
            Ro_T = Ro_T, Rh_T = Rh_T))

}
