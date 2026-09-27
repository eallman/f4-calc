include("functions/load_functions.jl")


## ########################################################################### ##
## --------------------------------------------------------------------------- ##
##                        QUARTET CONCORDANCE FACTORS                          ##
## --------------------------------------------------------------------------- ##
## ########################################################################### ##

## PARAMETER NAMES
## - h_1, h_2: for hybrid parameters
## - x_1, .. , x_k: for edges in a cycle but not hybrid
## - l_1, .. , l_y: for cut edges
## - g: hybridization parameter (γ)



## ------------------------- ##
##     3-2 cycle (5 taxa)    ##
## ------------------------- ##

R, (g, h_1, h_2, x, l_1, l_2) = polynomial_ring(QQ, [:g, :h_1, :h_2, :x, :l_1, :l_2])

#### INDEPENDENT INHERITANCE CORRELATION (ρ = 0) ####
cf_5_32_ind = CFs_5_3₂cycle(g, h_1, h_2, x, l_1, l_2, false)
I_ind = vanishing_ideal(cf_5_32_ind)

#### COMMON INHERITANCE CORRELATION (ρ = 1) ####
cf_5_32_com = CFs_5_3₂cycle_com(g, h_1, h_2, x, l_1, l_2, false)
I_com = vanishing_ideal(cf_5_32_com)

### Check if I_ind == I_com
R_ind = base_ring(I_ind); R_com = base_ring(I_com);
inject_Rind = hom(R_ind, R_com, gens(R_com))

inject_Rind(I_ind) == I_com



## ------------------------- ## 
##     3-2 cycle (6 taxa)    ## 
## ------------------------- ##

R, (g, h_1, h_2, x, l_1, l_2, l_3) = polynomial_ring(QQ, [:g, :h_1, :h_2, :x, :l_1, :l_2, :l_3])

#### INDEPENDENT INHERITANCE CORRELATION  (ρ = 0) ####
cf_6_32_ind = CFs_6_3₂cycle(g, h_1, h_2, x, l_1, l_2, l_3, false)
I_ind = vanishing_ideal(cf_6_32_ind)

#### COMMON INHERITANCE CORRELATION  (ρ = 1) ####
cf_6_32_com = CFs_6_3₂cycle_com(g, h_1, h_2, x, l_1, l_2, l_3, false)
I_com = vanishing_ideal(cf_6_32_com)

### Check if I_ind == I_com
R_ind = base_ring(I_ind); R_com = base_ring(I_com);
inject_Rind = hom(R_ind, R_com, gens(R_com))

inject_Rind(I_ind) == I_com



## ------------------------- ## 
##   4-cycle SOUTH (5 taxa)  ## 
## ------------------------- ##

R, (g, h_1, h_2, x_1, x_2, l) = polynomial_ring(QQ, [:g, :h_1, :h_2, :x_1, :x_2, :l])

#### INDEPENDENT INHERITANCE CORRELATION  (ρ = 0) ####
cf_5_south_ind = CFs_5_SOUTH(g, h_1, h_2, x_1, x_2, l, ["a₁", "a₂", "b", "c", "d"], true)
I_ind = vanishing_ideal(cf_5_south_ind)

#### COMMON INHERITANCE CORRELATION  (ρ = 1) ####
cf_5_south_com = CFs_5_SOUTH_com(g, h_1, h_2, x_1, x_2, l, ["a₁", "a₂", "b", "c", "d"], true)
I_com = vanishing_ideal(cf_5_south_com)

### Check if I_ind == I_com ####
R_ind = base_ring(I_ind); R_com = base_ring(I_com);
inject_Rind = hom(R_ind, R_com, gens(R_com))

inject_Rind(I_ind) == I_com


## Identifiability of numerical parameters for ρ = 1
## -------------------------------------------------
elim_ring, elim_gens, I  = elimination_ideal(parametrization(cf_5_south_com))
# elim_gens[1] =  [g, h_1,  h_2,  x_1,  x_2,  l]

eliminate(I, elim_gens[1][[2,3,4,5,6]]) # Find relations for g (γ)
eliminate(I, elim_gens[1][[1,2,3,5,6]]) # Find relations for x_1
eliminate(I, elim_gens[1][[1,2,3,4,6]]) # Find relations for x_2
# -> g, x_1 and x_2 seem identifiable in this case (ρ = 1)

eliminate(I, elim_gens[1][[1,4,5,6]]) # Find relations for h_1 and h_2 
eliminate(I, elim_gens[1][[1,3,4,5]]) # Find relations for h_1 and l
# -> 1 dregree of freedom between h_1, h_2 and l
