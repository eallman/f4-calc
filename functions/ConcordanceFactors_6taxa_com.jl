include("ConcordanceFactors_utilities.jl")

include("ConcordanceFactors_4taxa_com.jl")
include("ConcordanceFactors_5taxa_com.jl")

## Tree 2-2-2
function CFs_6_tree_2_2_2(ℓ₁, ℓ₂, ℓ₃, l::Vector{String}, reduced::Bool)
    CF₁₂ = CF_tree(ℓ₂*ℓ₃, [l[3], l[4], l[5], l[6]], reduced)
    CF₃₄ = CF_tree(ℓ₁*ℓ₃, [l[1], l[2], l[5], l[6]], reduced)
    CF₅₆ = CF_tree(ℓ₁*ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    CF₂₄ = CF_tree(ℓ₃, [l[1], l[3], l[5], l[6]], reduced)
    CF₂₆ = CF_tree(ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    CF₄₆ = CF_tree(ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    merge(CF₁₂, CF₃₄, CF₅₆, CF₂₄, CF₂₆, CF₄₆)
end
CFs_6_tree_2_2_2(ℓ₁, ℓ₂, ℓ₃) = CFs_6_tree_2_2_2(ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b₁", "b₂", "c₁", "c₂"], false);
CFs_6_tree_2_2_2(ℓ₁, ℓ₂, ℓ₃, reduced::Bool) = CFs_6_tree_2_2_2(ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b₁", "b₂", "c₁", "c₂"], reduced);
CFs_6_tree_2_2_2(ℓ₁, ℓ₂, ℓ₃, label::Vector{String}) = CFs_6_tree_2_2_2(ℓ₁, ℓ₂, ℓ₃, label, false);

## 3₂-cycle
function CFs_6_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, l::Vector{String}, reduced::Bool)
    CF₁₂ = CF_tree(x*ℓ₂*ℓ₃, [l[3], l[4], l[5], l[6]], reduced)
    #CF₃₄ = CF_2cycle_com(γ, h₁*x, h₂, ℓ₁, ℓ₃, [l[1], l[2], l[5], l[6]], reduced)
    CF₃₄ = CF_2cycle_com(γ, h₁*x, h₂, ℓ₃, ℓ₁, [l[5], l[6], l[1], l[2]], reduced)
    #CF₅₆ = CF_2cycle_com(γ, h₁, x*h₂, ℓ₁, ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    CF₅₆ = CF_2cycle_com(γ, h₁, x*h₂, ℓ₂, ℓ₁, [l[3], l[4], l[1], l[2]], reduced)
    CF₂₄ = CF_3₁cycle_com(γ, x, ℓ₃, [l[1], l[3], l[5], l[6]], reduced)
    CF₂₆ = CF_3₁cycle_com(1-γ, x, ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    CF₄₆ = CF_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    merge(CF₁₂, CF₃₄, CF₅₆, CF₂₄, CF₂₆, CF₄₆)
end
CFs_6_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃) = CFs_6_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b₁", "b₂", "c₁", "c₂"], false);
CFs_6_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, reduced::Bool) = CFs_6_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b₁", "b₂", "c₁", "c₂"], reduced);
CFs_6_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, label::Vector{String}) = CFs_6_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, label, false);

## 4₂-cycle: SOUTH-WEST
function CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    CF₃ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ₁, [l[1], l[2], l[3], l[5], l[6]], reduced)

    CF₂₆ = CF_3₁cycle_com(1-γ, x₁, ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    CF₂₅ = CF_3₁cycle_com(1-γ, x₁*x₂, ℓ₂, [l[1], l[6], l[3], l[4]], reduced)
    CF₁₂ = CF_tree(ℓ₂*x₁, [l[3], l[4], l[5], l[6]], reduced)

    CF₅₆ = CF_2cycle_com(γ, h₁, h₂*x₁*x₂, ℓ₁, ℓ₂, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF₃, CF₂₆, CF₂₅, CF₁₂, CF₅₆)
end
CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂) = CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c", "d"], false);
CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c", "d"], reduced);
CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label, false);

## 4₂-cycle: SOUTH-NORTH
function CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    CF₁ = CFs_5_NORTH(γ, x₁, x₂, ℓ₂, [l[1], l[3], l[4], l[5], l[6]], reduced)
    CF₃ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ₁, [l[1], l[2], l[3], l[4], l[6]], reduced)
    CF₅₆ = CF_2cycle_com(γ, h₁*x₁, h₂*x₂, ℓ₁, ℓ₂, [l[1], l[2], l[4], l[5]], reduced)
    merge(CF₁, CF₃, CF₅₆)
end
CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂) = CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label, false);

## 5₂-cycle
function CFs_6_SOUTH(γ, h₁, h₂, x₁, x₂, x₃, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂ =  CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃, [l[1], l[3], l[4], l[5], l[6]], reduced)
    CF₃₄ = CF_3₂cycle_com(γ, h₁*x₁*x₂, h₂, x₃, ℓ, [l[1], l[2], l[5], l[6]], reduced)
    CF₃₅ = CF_3₂cycle_com(γ, h₁*x₁, h₂, x₂*x₃, ℓ, [l[1], l[2], l[4], l[6]], reduced)
    CF₃₆ = CF_3₂cycle_com(γ, h₁*x₁, h₂*x₃, x₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF₄₅ = CF_3₂cycle_com(γ, h₁, h₂, x₁*x₂*x₃, ℓ, [l[1], l[2], l[3], l[6]], reduced)
    CF₄₆ = CF_3₂cycle_com(γ, h₁, h₂*x₃, x₁*x₂, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    CF₅₆ = CF_3₂cycle_com(γ, h₁, h₂*x₂*x₃, x₁, ℓ, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF₁₂, CF₃₄, CF₃₅, CF₃₆, CF₄₅, CF₄₆, CF₅₆)
end;
CFs_6_SOUTH(γ, h₁, h₂, x₁, x₂, x₃, ℓ) = CFs_6_SOUTH(γ, h₁, h₂, x₁, x₂, x₃, ℓ, ["a₁", "a₂", "b", "c", "d", "e"], false);
CFs_6_SOUTH(γ, h₁, h₂, x₁, x₂, x₃, ℓ, reduced::Bool) = CFs_6_SOUTH(γ, h₁, h₂, x₁, x₂, x₃, ℓ, ["a₁", "a₂", "b", "c", "d", "e"], reduced);
CFs_6_SOUTH(γ, h₁, h₂, x₁, x₂, x₃, ℓ, label::Vector{String}) = CFs_6_SOUTH(γ, h₁, h₂, x₁, x₂, x₃, ℓ, label, false);


## 4₂: SOUTH-NORTH
function CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄₆ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ₁, [l[1], l[2], l[3], l[4], l[6]], reduced)

    CF₁₃₅₆ = CF_4₁cycle(γ, x₁, x₂,[l[1], l[3], l[5], l[6]], reduced)
    CF₂₃₅₆ = CF_4₁cycle(γ, x₁, x₂,[l[2], l[3], l[5], l[6]], reduced)
    
    CF₁₃₄₅ = CF_3₁cycle_com(γ, x₁, ℓ₂, [l[1], l[3], l[4], l[5]], reduced)
    CF₂₃₄₅ = CF_3₁cycle_com(γ, x₁, ℓ₂, [l[2], l[3], l[4], l[5]], reduced)
    CF₁₄₅₆ = CF_3₁cycle_com(1-γ, x₂, ℓ₂, [l[1], l[6], l[4], l[5]], reduced)
    CF₂₄₅₆ = CF_3₁cycle_com(1-γ, x₂, ℓ₂, [l[2], l[6], l[4], l[5]], reduced)

    CF₁₂₃₅ = CF_3₂cycle_com(γ, h₁, h₂*x₂, x₁, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    CF₁₂₅₆ = CF_3₂cycle_com(γ, h₁*x₁, h₂, x₂, ℓ₁, [l[1], l[2], l[5], l[6]], reduced)
    
    CF₁₂₄₅ = CF_2cycle_com(γ, h₁*x₁, h₂*x₂, ℓ₁, ℓ₂, [l[1], l[2], l[4], l[5]], reduced)
    CF₃₄₅₆ = CF_tree(ℓ₂, [l[3], l[6], l[4], l[5]], reduced)

    merge(CF₁₂₃₄₆, CF₁₂₅₆, CF₁₃₅₆, CF₂₃₅₆, CF₁₂₄₅, CF₁₃₄₅, CF₂₃₄₅, CF₁₄₅₆, CF₂₄₅₆, CF₃₄₅₆, CF₁₂₃₅)
end
CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂) = CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a", "b", "c", "d", "e", "f"], false);
CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a", "b", "c", "d", "e", "f"], reduced);
CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label, false);

# ## 4₂: SOUTH-NORTH
# function CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
#     CF₁₂₃₄₆ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ₁, [l[1], l[2], l[3], l[4], l[6]], reduced)

#     CF₁₃₅₆ = CF_4₁cycle(γ, x₁, x₂,[l[1], l[3], l[5], l[6]], reduced)
#     CF₂₃₅₆ = CF_4₁cycle(γ, x₁, x₂,[l[2], l[3], l[5], l[6]], reduced)
    
#     CF₁₂₃₅ = CF_3₂cycle_com(γ, h₁, h₂*x₂, x₁, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
#     CF₁₂₅₆ = CF_3₂cycle_com(γ, h₁*x₁, h₂, x₂, ℓ₁, [l[1], l[2], l[5], l[6]], reduced)
    
#     CF₃₄₅₆ = CF_tree(ℓ₂, [l[3], l[6], l[4], l[5]], reduced)

#     CF₁₃₄₅ = CF_3₁_3₁_chain(γ, x₁, δ, y, ℓ₂, [l[1], l[3], l[5], l[4]], reduced)
#     CF₂₃₄₅ = CF_3₁_3₁_chain(γ, x₁, δ, y, ℓ₂, [l[2], l[3], l[5], l[4]], reduced)
#     CF₁₄₅₆ = CF_3₁_3₁_chain(1-γ, x₂, δ, y, ℓ₂, [l[1], l[6], l[5], l[4]], reduced)
#     CF₂₄₅₆ = CF_3₁_3₁_chain(1-γ, x₂, δ, y, ℓ₂, [l[2], l[6], l[5], l[4]], reduced)

#     CF₁₂₄₅ = CF_2₁_3₁_chain(γ, h₁*x₁, h₂*x₂, ℓ₁, ℓ₂, δ, y, [l[1], l[2], l[5], l[4]], reduced)
    

#     merge(CF₁₂₃₄₆, CF₁₂₅₆, CF₁₃₅₆, CF₂₃₅₆, CF₁₂₄₅, CF₁₃₄₅, CF₂₃₄₅, CF₁₄₅₆, CF₂₄₅₆, CF₃₄₅₆, CF₁₂₃₅)
# end
# CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂) = CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, ["a", "b", "c", "d", "e", "f"], false);
# CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, ["a", "b", "c", "d", "e", "f"], reduced);
# CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, label, false);




