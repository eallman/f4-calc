include("ConcordanceFactors_utilities.jl")

include("ConcordanceFactors_4taxa.jl")
include("ConcordanceFactors_5taxa.jl")

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
function CFs_6_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, l::Vector{String}, reduced::Bool)
    CF₁₂ = CF_tree(x*ℓ₂*ℓ₃, [l[3], l[4], l[5], l[6]], reduced)
    #CF₃₄ = CF_2cycle(γ, h₁*x, h₂, ℓ₁, ℓ₃, [l[1], l[2], l[5], l[6]], reduced)
    CF₃₄ = CF_2cycle(γ, h₁*x, h₂, ℓ₃, ℓ₁, [l[5], l[6], l[1], l[2]], reduced)
    #CF₅₆ = CF_2cycle(γ, h₁, x*h₂, ℓ₁, ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    CF₅₆ = CF_2cycle(γ, h₁, x*h₂, ℓ₂, ℓ₁, [l[3], l[4], l[1], l[2]], reduced)
    CF₂₄ = CF_3₁cycle(γ, x, ℓ₃, [l[1], l[3], l[5], l[6]], reduced)
    CF₂₆ = CF_3₁cycle(1-γ, x, ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    CF₄₆ = CF_3₂cycle(γ, h₁, h₂, x, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    merge(CF₁₂, CF₃₄, CF₅₆, CF₂₄, CF₂₆, CF₄₆)
end
CFs_6_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃) = CFs_6_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b₁", "b₂", "c₁", "c₂"], false);
CFs_6_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, reduced::Bool) = CFs_6_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b₁", "b₂", "c₁", "c₂"], reduced);
CFs_6_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, label::Vector{String}) = CFs_6_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ℓ₃, label, false);

## 4₂-cycle: SOUTH-WEST
function CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    CF₃ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ₁, [l[1], l[2], l[3], l[5], l[6]], reduced)

    CF₂₆ = CF_3₁cycle(1-γ, x₁, ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    CF₂₅ = CF_3₁cycle(1-γ, x₁*x₂, ℓ₂, [l[1], l[6], l[3], l[4]], reduced)
    CF₁₂ = CF_tree(ℓ₂*x₁, [l[3], l[4], l[5], l[6]], reduced)

    CF₅₆ = CF_2cycle(γ, h₁, h₂*x₁*x₂, ℓ₁, ℓ₂, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF₃, CF₂₆, CF₂₅, CF₁₂, CF₅₆)
end
CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂) = CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c", "d"], false);
CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c", "d"], reduced);
CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_SOUTH_WEST(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label, false);

## 4₂-cycle: SOUTH-NORTH
function CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    CF₁ = CFs_5_NORTH(γ, x₁, x₂, ℓ₂, [l[1], l[3], l[4], l[5], l[6]], reduced)
    CF₃ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ₁, [l[1], l[2], l[3], l[4], l[6]], reduced)
    CF₅₆ = CF_2cycle(γ, h₁*x₁, h₂*x₂, ℓ₁, ℓ₂, [l[1], l[2], l[4], l[5]], reduced)
    merge(CF₁, CF₃, CF₅₆)
end
CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂) = CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_SOUTH_NORTH(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label, false);

## 5₂-cycle
function CFs_6_SOUTH(γ, h₁, h₂, x₁, x₂, x₃, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂ =  CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃, [l[1], l[3], l[4], l[5], l[6]], reduced)
    CF₃₄ = CF_3₂cycle(γ, h₁*x₁*x₂, h₂, x₃, ℓ, [l[1], l[2], l[5], l[6]], reduced)
    CF₃₅ = CF_3₂cycle(γ, h₁*x₁, h₂, x₂*x₃, ℓ, [l[1], l[2], l[4], l[6]], reduced)
    CF₃₆ = CF_3₂cycle(γ, h₁*x₁, h₂*x₃, x₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF₄₅ = CF_3₂cycle(γ, h₁, h₂, x₁*x₂*x₃, ℓ, [l[1], l[2], l[3], l[6]], reduced)
    CF₄₆ = CF_3₂cycle(γ, h₁, h₂*x₃, x₁*x₂, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    CF₅₆ = CF_3₂cycle(γ, h₁, h₂*x₂*x₃, x₁, ℓ, [l[1], l[2], l[3], l[4]], reduced)

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
    
    CF₁₃₄₅ = CF_3₁cycle(γ, x₁, ℓ₂, [l[1], l[3], l[4], l[5]], reduced)
    CF₂₃₄₅ = CF_3₁cycle(γ, x₁, ℓ₂, [l[2], l[3], l[4], l[5]], reduced)
    CF₁₄₅₆ = CF_3₁cycle(1-γ, x₂, ℓ₂, [l[1], l[6], l[4], l[5]], reduced)
    CF₂₄₅₆ = CF_3₁cycle(1-γ, x₂, ℓ₂, [l[2], l[6], l[4], l[5]], reduced)

    CF₁₂₃₅ = CF_3₂cycle(γ, h₁, h₂*x₂, x₁, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    CF₁₂₅₆ = CF_3₂cycle(γ, h₁*x₁, h₂, x₂, ℓ₁, [l[1], l[2], l[5], l[6]], reduced)
    
    CF₁₂₄₅ = CF_2cycle(γ, h₁*x₁, h₂*x₂, ℓ₁, ℓ₂, [l[1], l[2], l[4], l[5]], reduced)
    CF₃₄₅₆ = CF_tree(ℓ₂, [l[3], l[6], l[4], l[5]], reduced)

    merge(CF₁₂₃₄₆, CF₁₂₅₆, CF₁₃₅₆, CF₂₃₅₆, CF₁₂₄₅, CF₁₃₄₅, CF₂₃₄₅, CF₁₄₅₆, CF₂₄₅₆, CF₃₄₅₆, CF₁₂₃₅)
end
CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂) = CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a", "b", "c", "d", "e", "f"], false);
CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, ["a", "b", "c", "d", "e", "f"], reduced);
CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_SOUTH_NORTH_2(γ, h₁, h₂, x₁, x₂, ℓ₁, ℓ₂, label, false);

## 4₂: SOUTH-NORTH
function CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄₆ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ₁, [l[1], l[2], l[3], l[4], l[6]], reduced)

    CF₁₃₅₆ = CF_4₁cycle(γ, x₁, x₂,[l[1], l[3], l[5], l[6]], reduced)
    CF₂₃₅₆ = CF_4₁cycle(γ, x₁, x₂,[l[2], l[3], l[5], l[6]], reduced)
    
    CF₁₂₃₅ = CF_3₂cycle(γ, h₁, h₂*x₂, x₁, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    CF₁₂₅₆ = CF_3₂cycle(γ, h₁*x₁, h₂, x₂, ℓ₁, [l[1], l[2], l[5], l[6]], reduced)
    
    CF₃₄₅₆ = CF_tree(ℓ₂, [l[3], l[6], l[4], l[5]], reduced)

    CF₁₃₄₅ = CF_3₁_3₁_chain(γ, x₁, δ, y, ℓ₂, [l[1], l[3], l[5], l[4]], reduced)
    CF₂₃₄₅ = CF_3₁_3₁_chain(γ, x₁, δ, y, ℓ₂, [l[2], l[3], l[5], l[4]], reduced)
    CF₁₄₅₆ = CF_3₁_3₁_chain(1-γ, x₂, δ, y, ℓ₂, [l[1], l[6], l[5], l[4]], reduced)
    CF₂₄₅₆ = CF_3₁_3₁_chain(1-γ, x₂, δ, y, ℓ₂, [l[2], l[6], l[5], l[4]], reduced)

    CF₁₂₄₅ = CF_2₁_3₁_chain(γ, h₁*x₁, h₂*x₂, ℓ₁, ℓ₂, δ, y, [l[1], l[2], l[5], l[4]], reduced)
    

    merge(CF₁₂₃₄₆, CF₁₂₅₆, CF₁₃₅₆, CF₂₃₅₆, CF₁₂₄₅, CF₁₃₄₅, CF₂₃₄₅, CF₁₄₅₆, CF₂₄₅₆, CF₃₄₅₆, CF₁₂₃₅)
end
CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂) = CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, ["a", "b", "c", "d", "e", "f"], false);
CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, ["a", "b", "c", "d", "e", "f"], reduced);
CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_SOUTH_3₁NORTH(γ, h₁, h₂, x₁, x₂, δ, y, ℓ₁, ℓ₂, label, false);






########################
## ## ## CHAINS ## ## ##
########################

## 4₁-3₁-cherry-cycles: WEST_3₁-cherry
function CFs_6_WEST_3₁cherry(γ, x₁, x₂, δ, y, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    CF₂ = CFs_5_WEST_3₁(γ, x₁, x₂, δ, y, ℓ₁, [l[1], l[3], l[4], l[5], l[6]], reduced)
    CF₃₄ = CF_3₁cycle(1-δ, y, ℓ₂, [l[6], l[5], l[1], l[2]], reduced)
    CF₃₅ = CF_3₁cycle(1-δ, y, ℓ₂, [l[6], l[4], l[1], l[2]], reduced)
    CF₃₆ = CF_3₁cycle(1-γ, x₁*x₂, y*ℓ₁*ℓ₂, [l[5], l[4], l[1], l[2]], reduced)
    CF₄₅ = CF_3₁cycle(1-δ, y, ℓ₂, [l[6], l[3], l[1], l[2]], reduced)
    CF₄₆ = CF_3₁cycle(1-γ, x₁, y*ℓ₁*ℓ₂, [l[5], l[3], l[1], l[2]], reduced)

    CF₅₆ = CF_tree(x₁*y*ℓ₁*ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    merge(CF₂, CF₃₄, CF₃₅, CF₃₆, CF₄₅, CF₄₆, CF₅₆)
end

CFs_6_WEST_3₁cherry(γ, x₁, x₂, δ, y, ℓ₁, ℓ₂) = CFs_6_WEST_3₁cherry(γ, x₁, x₂, δ, y, ℓ₁, ℓ₂, ["a₁", "a₂", "b", "c", "d", "e"], false);
CFs_6_WEST_3₁cherry(γ, x₁, x₂, δ, y, ℓ₁, ℓ₂, reduced::Bool) = CFs_6_WEST_3₁cherry(γ, x₁, x₂, δ, y, ℓ₁, ℓ₂, ["a₁", "a₂", "b", "c", "d", "e"], reduced);
CFs_6_WEST_3₁cherry(γ, x₁, x₂, δ, y, ℓ₁, ℓ₂, label::Vector{String}) = CFs_6_WEST_3₁cherry(γ, x₁, x₂, δ, y, ℓ₁, ℓ₂, label, false);

## Tree
#function CFs_6_tree(ℓ₁, ℓ₂, ℓ₃, l::Vector{String}, reduced::Bool)
#    CF₁₂₄₅₆ = CFs_5_tree(ℓ₃*ℓ₁, ℓ₂, [l[1], l[2], l[4], l[5], l[6]], reduced)
#    CF₁₂₃₄ = CF_tree(ℓ₁, [l[1], l[2], l[3], l[4]], reduced)
#    CF₁₂₃₆ = CF_tree(ℓ₁, [l[1], l[2], l[3], l[6]], reduced)
#    CF₁₃₄₅ = CF_tree(ℓ₃*ℓ₁*ℓ₂, [l[1], l[3], l[4], l[5]], reduced)
#    CF₁₃₄₆ = CF_tree(ℓ₃, [l[1], l[3], l[4], l[6]], reduced)
#    CF₃₄₅₆ = CF_tree(ℓ₂, [l[3], l[6], l[4], l[5]], reduced)
#
#    merge(CF₁₂₄₅₆, CF₁₂₃₄, CF₁₂₃₆, CF₁₃₄₅, CF₁₃₄₆, CF₃₄₅₆)
#end
#CFs_6_tree(ℓ₁, ℓ₂, ℓ) = CFs_6_tree(ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
#CFs_6_tree(ℓ₁, ℓ₂, ℓ, reduced::Bool) = CFs_6_tree(ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
#CFs_6_tree(ℓ₁, ℓ₂, ℓ, label::Vector{String}) = CFs_6_tree(ℓ₁, ℓ₂, ℓ, label, false);

## Tree
function CFs_6_tree(ℓ₁, ℓ₂, ℓ₃, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄₅ = CFs_5_tree(ℓ₁, ℓ₂, [l[1], l[2], l[3], l[4], l[5]], reduced)
    CF₁₂₅₆ = CF_tree(ℓ₁*ℓ₃, [l[1], l[2], l[5], l[6]], reduced)
    CF₃₄₅₆ = CF_tree(ℓ₂*ℓ₃, [l[3], l[4], l[5], l[6]], reduced)
    CF₁₃₅₆ = CF_tree(ℓ₃, [l[1], l[3], l[5], l[6]], reduced)
    merge(CF₁₂₃₄₅, CF₁₂₅₆, CF₃₄₅₆, CF₁₃₅₆)
end
CFs_6_tree(ℓ₁, ℓ₂, ℓ₃) = CFs_6_tree(ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_tree(ℓ₁, ℓ₂, ℓ₃, reduced::Bool) = CFs_6_tree(ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_tree(ℓ₁, ℓ₂, ℓ₃, label::Vector{String}) = CFs_6_tree(ℓ₁, ℓ₂, ℓ₃, label, false);


## Tree
function CFs_6_tree_2cherries(ℓ₁, ℓ₂, ℓ₃, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄₅ = CFs_5_tree(ℓ₁, ℓ₂*ℓ₃, [l[1], l[2], l[3], l[4], l[5]], reduced)
    CF₁₂₅₆ = CF_tree(ℓ₁, [l[1], l[2], l[5], l[6]], reduced)
    CF₃₄₅₆ = CF_tree(ℓ₂, [l[3], l[4], l[5], l[6]], reduced)
    CF₁₃₅₆ = CF_tree(ℓ₃, [l[1], l[5], l[3], l[6]], reduced)
    merge(CF₁₂₃₄₅, CF₁₂₅₆, CF₃₄₅₆, CF₁₃₅₆)
end
CFs_6_tree_2cherries(ℓ₁, ℓ₂, ℓ₃) = CFs_6_tree_2cherries(ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_tree_2cherries(ℓ₁, ℓ₂, ℓ₃, reduced::Bool) = CFs_6_tree_2cherries(ℓ₁, ℓ₂, ℓ₃, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_tree_2cherries(ℓ₁, ℓ₂, ℓ₃, label::Vector{String}) = CFs_6_tree_2cherries(ℓ₁, ℓ₂, ℓ₃, label, false);


## 3₁-3₁-cycles
function CFs_6_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂₄₅₆ = CFs_5_3₁cycle(δ, y, x*ℓ*ℓ₁, ℓ₂, [l[1], l[2], l[4], l[5], l[6]], reduced)
    CF₁₂₃₄ = CF_3₁cycle(γ, x, ℓ₁, [l[3], l[4], l[1], l[2]], reduced)
    CF₁₂₃₆ = CF_3₁cycle(γ, x, ℓ₁, [l[3], l[6], l[1], l[2]], reduced)
    CF₁₃₄₅ = CF_3₁cycle(1-γ, x, ℓ*y*ℓ₂, [l[3], l[1], l[4], l[5]], reduced)
    CF₁₃₄₆ = CF_3₁_3₁_chain(1-γ, x, δ, y, ℓ, [l[3], l[1], l[4], l[6]], reduced)
    CF₃₄₅₆ = CF_3₁cycle(δ, y, ℓ₂, [l[6], l[3], l[4], l[5]], reduced)

    merge(CF₁₂₄₅₆, CF₁₂₃₄, CF₁₂₃₆, CF₁₃₄₅, CF₁₃₄₆, CF₃₄₅₆)
end
CFs_6_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ℓ) = CFs_6_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ℓ, reduced::Bool) = CFs_6_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ℓ, label::Vector{String}) = CFs_6_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ℓ, label, false);

## 3₁-3₂out-cycles
function CFs_6_3₁_3₂o_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂₄₅₆ = CFs_5_3₂cycle(1-δ, k₂, k₁, y, ℓ₂, ℓ*x*ℓ₁, [l[4], l[5], l[1], l[2], l[6]], reduced)
    CF₁₂₃₄ = CF_3₁cycle(γ, x, ℓ₁, [l[3], l[4], l[1], l[2]], reduced)
    CF₁₂₃₆ = CF_3₁cycle(γ, x, ℓ₁, [l[3], l[6], l[1], l[2]], reduced)
    CF₁₃₄₅ = CF_3₁_2₂_chain(1-γ, x, δ, y*k₁, k₂, ℓ, ℓ₂, [l[3], l[1], l[4], l[5]], reduced)
    CF₁₃₄₆ = CF_3₁_3₁_chain(1-γ, x, 1-δ, y, ℓ, [l[3], l[1], l[6], l[4]], reduced)
    CF₃₄₅₆ = CF_3₂cycle(1-δ, k₂, k₁, y, ℓ₂, [l[4], l[5], l[3], l[6]], reduced)

    merge(CF₁₂₄₅₆, CF₁₂₃₄, CF₁₂₃₆, CF₁₃₄₅, CF₁₃₄₆, CF₃₄₅₆)
end
CFs_6_3₁_3₂o_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ) = CFs_6_3₁_3₂o_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_3₁_3₂o_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, reduced::Bool) = CFs_6_3₁_3₂o_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_3₁_3₂o_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, label::Vector{String}) = CFs_6_3₁_3₂o_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, label, false);

## 3₂out-3₂out_cycles
function CFs_6_3₂o_3₂o_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄ = CF_3₂cycle(γ, h₁, h₂, x, ℓ₁, [l[1], l[2], l[3], l[4]], reduced)
    CF₁₂₃₆ = CF_3₂cycle(γ, h₁, h₂, x, ℓ₁, [l[1], l[2], l[3], l[6]], reduced)
    CF₁₂₄₅ = CF_2₂_2₂_chain(γ, h₁*x, h₂, δ, k₁*y, k₂, ℓ₁, ℓ, ℓ₂, [l[1], l[2], l[4], l[5]], reduced)
    CF₁₂₄₆ = CF_3₁_2₂_chain(δ, y, γ, h₁*x, h₂, ℓ, ℓ₁, [l[4], l[6], l[1], l[2]], reduced)
    CF₁₃₄₅ = CF_3₁_2₂_chain(γ, x, δ, k₁*y, k₂, ℓ, ℓ₂, [l[1], l[3], l[4], l[5]], reduced)
    CF₁₃₄₆ = CF_3₁_3₁_chain(γ, x, 1-δ, y, ℓ, [l[1], l[3], l[6], l[4]], reduced)
    CF₁₄₅₆ = CF_3₂cycle(δ, k₁, k₂, y, ℓ₂, [l[4], l[5], l[6], l[1]], reduced)
    CF₃₄₅₆ = CF_3₂cycle(δ, k₁, k₂, y, ℓ₂, [l[4], l[5], l[6], l[3]], reduced)

    merge(CF₁₂₃₄, CF₁₂₃₆, CF₁₂₄₅, CF₁₂₄₆, CF₁₃₄₅, CF₁₃₄₆, CF₁₄₅₆, CF₃₄₅₆)
end
CFs_6_3₂o_3₂o_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ) = CFs_6_3₂o_3₂o_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_3₂o_3₂o_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, reduced::Bool) = CFs_6_3₂o_3₂o_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_3₂o_3₂o_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, label::Vector{String}) = CFs_6_3₂o_3₂o_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, label, false);

## 3₁-3₂_cycles
function CFs_6_3₁_3₂_chain(γ, h₁, h₂, x, δ, y, ℓ₁, ℓ₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄₅ = CFs_5_3₂cycle(γ, h₁, h₂, x, ℓ*y*ℓ₂, ℓ₁, [l[4], l[5], l[1], l[2], l[3]], reduced)

    CF₁₂₃₆ = CF_3₁cycle(1-γ, x, ℓ₁, [l[6], l[3], l[1], l[2]], reduced)
    CF₁₄₅₆ = CF_3₁cycle(δ, y, ℓ₂, [l[6], l[1], l[4], l[5]], reduced)
    CF₃₄₅₆ = CF_3₁cycle(δ, y, ℓ₂, [l[6], l[3], l[4], l[5]], reduced)

    CF₁₂₄₆ = CF_2₁_3₁_chain(δ, y, γ, h₁, x*h₂, ℓ, ℓ₁, [l[6], l[4], l[1], l[2]], reduced)
    CF₁₃₄₆ = CF_3₂_3₁_chain(γ, h₁, h₂, x, δ, y, ℓ, [l[1], l[3], l[4], l[6]], reduced)

    merge(CF₁₂₃₄₅, CF₁₂₃₆, CF₁₄₅₆, CF₃₄₅₆, CF₁₂₄₆, CF₁₃₄₆)
end
CFs_6_3₁_3₂_chain(γ, h₁, h₂, x, δ, y, ℓ₁, ℓ₂, ℓ) = CFs_6_3₁_3₂_chain(γ, h₁, h₂, x, δ, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_3₁_3₂_chain(γ, h₁, h₂, x, δ, y, ℓ₁, ℓ₂, ℓ, reduced::Bool) = CFs_6_3₁_3₂_chain(γ, h₁, h₂, x, δ, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_3₁_3₂_chain(γ, h₁, h₂, x, δ, y, ℓ₁, ℓ₂, ℓ, label::Vector{String}) = CFs_6_3₁_3₂_chain(γ, h₁, h₂, x, δ, y, ℓ₁, ℓ₂, ℓ, label, false);

## 3₂-3₂_cycles
function CFs_6_3₂_3₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, l::Vector{String}, reduced::Bool)

    CF₁₂₃₄ = CF_3₁cycle(1-γ, x, ℓ₁, [l[4], l[3], l[1], l[2]], reduced)
    CF₁₂₃₆ = CF_3₁cycle(1-γ, x, ℓ₁, [l[6], l[3], l[1], l[2]], reduced)

    CF₁₄₅₆ = CF_3₂cycle(δ, k₁, k₂, y, ℓ₂, [l[4], l[5], l[6], l[1]], reduced)
    CF₃₄₅₆ = CF_3₂cycle(δ, k₁, k₂, y, ℓ₂, [l[4], l[5], l[6], l[3]], reduced)

    CF₁₂₄₅ = CF_2₂_2₂_chain(γ, h₁, h₂*x, δ, k₁*y, k₂, ℓ₁, ℓ, ℓ₂, [l[1], l[2], l[4], l[5]], reduced)
    CF₁₂₄₆ = CF_3₁_2₂_chain(δ, y, γ, h₁*x, h₂, ℓ, ℓ₁, [l[4], l[6], l[1], l[2]], reduced)
    CF₁₃₄₆ = CF_3₂_3₁_chain(γ, h₁, h₂, x, 1-δ, y, ℓ, [l[1], l[3], l[6], l[4]], reduced)

    CF₁₃₄₅ = CF_3₂_2₂_chain(γ, h₁, h₂, x, δ, y*k₁, k₂, ℓ, ℓ₂, [l[4], l[5], l[1], l[3]], reduced)

    merge(CF₁₂₃₄, CF₁₂₃₆, CF₁₄₅₆, CF₃₄₅₆, CF₁₂₄₅, CF₁₂₄₆, CF₁₃₄₆, CF₁₃₄₅)
end
CFs_6_3₂_3₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ) = CFs_6_3₂_3₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], false);
CFs_6_3₂_3₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, reduced::Bool) = CFs_6_3₂_3₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, ["a₁", "a₂", "b", "c₁", "c₂", "d"], reduced);
CFs_6_3₂_3₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, label::Vector{String}) = CFs_6_3₂_3₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ℓ, label, false);

## 4₁-4₁: EAST-WEST
function CFs_6_EAST_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[3], l[4]], reduced)
    CF₁₂₃₅ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[3], l[5]], reduced)
    CF₁₂₃₆ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[3], l[6]], reduced)
    CF₁₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[6], l[1], l[4], l[5]], reduced)
    CF₂₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[6], l[2], l[4], l[5]], reduced)
    CF₃₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[6], l[3], l[4], l[5]], reduced)

    CF₁₂₄₅ = CF_3₁cycle(γ, x₁*x₂, ℓ*y₁, [l[1], l[2], l[4], l[5]], reduced)
    CF₁₃₄₅ = CF_3₁cycle(γ, x₂, ℓ*y₁, [l[1], l[3], l[4], l[5]], reduced)
    CF₂₃₄₆ = CF_3₁cycle(1-δ, y₁, ℓ*x₂, [l[6], l[4], l[2], l[3]], reduced)
    CF₂₃₅₆ = CF_3₁cycle(1-δ, y₁*y₂, ℓ*x₂, [l[6], l[5], l[2], l[3]], reduced)

    CF₁₂₄₆ = CF_3₁_3₁_chain(γ, x₁*x₂, δ, y₁, ℓ, [l[1], l[2], l[4], l[6]], reduced)
    CF₁₂₅₆ = CF_3₁_3₁_chain(γ, x₁*x₂, δ, y₁*y₂, ℓ, [l[1], l[2], l[5], l[6]], reduced)
    CF₁₃₄₆ = CF_3₁_3₁_chain(γ, x₂, δ, y₁, ℓ, [l[1], l[3], l[4], l[6]], reduced)
    CF₁₃₅₆ = CF_3₁_3₁_chain(γ, x₂, δ, y₁*y₂, ℓ, [l[1], l[3], l[5], l[6]], reduced)

    CF₂₃₄₅ = CF_tree(x₂*ℓ*y₁, [l[2], l[3], l[4], l[5]], reduced)

    merge(CF₁₂₃₄, CF₁₂₃₅, CF₁₂₃₆, CF₁₂₄₅, CF₁₂₄₆, CF₁₂₅₆, CF₁₃₄₅, CF₁₃₄₆, CF₁₃₅₆, CF₁₄₅₆, CF₂₃₄₅, CF₂₃₄₆, CF₂₃₅₆, CF₂₄₅₆, CF₃₄₅₆)
end
CFs_6_EAST_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ) = CFs_6_EAST_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], false);
CFs_6_EAST_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, reduced::Bool) = CFs_6_EAST_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], reduced);
CFs_6_EAST_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, label::Vector{String}) = CFs_6_EAST_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, label, false);

## 4₁-4₁: SOUTH-NORTH
function CFs_6_SOUTH_NORTH_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄₆ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ, [l[4], l[6], l[1], l[2], l[3]], reduced)

    CF₁₂₅₆ = CF_3₂_3₁_chain(γ, h₁, h₂*x₂, x₁, 1-δ, y₁, ℓ, [l[1], l[2], l[6], l[5]], reduced)
    CF₁₃₅₆ = CF_3₂_3₁_chain(γ, h₁, h₂, x₁*x₂, 1-δ, y₁, ℓ, [l[1], l[3], l[6], l[5]], reduced)
    CF₂₃₅₆ = CF_3₂_3₁_chain(γ, h₁*x₁, h₂, x₂, 1-δ, y₁, ℓ, [l[2], l[3], l[6], l[5]], reduced)
    CF₁₂₄₅ = CF_3₂_3₁_chain(γ, h₁, h₂*x₂, x₁, δ, y₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF₁₃₄₅ = CF_3₂_3₁_chain(γ, h₁, h₂, x₁*x₂, δ, y₂, ℓ, [l[1], l[3], l[4], l[5]], reduced)
    CF₂₃₄₅ = CF_3₂_3₁_chain(γ, h₁*x₁, h₂, x₂, δ, y₂, ℓ, [l[2], l[3], l[4], l[5]], reduced)

    CF₁₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[5], l[6], l[1], l[4]], reduced)
    CF₂₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[5], l[6], l[2], l[4]], reduced)
    CF₃₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[5], l[6], l[3], l[4]], reduced)
    CF₁₂₃₅ = CF_4₁cycle(γ, x₁, x₂, [l[5], l[1], l[2], l[3]], reduced)

    merge(CF₁₂₃₄₆, CF₁₂₅₆, CF₁₃₅₆, CF₂₃₅₆, CF₁₂₄₅, CF₁₃₄₅, CF₂₃₄₅, CF₁₄₅₆, CF₂₄₅₆, CF₃₄₅₆, CF₁₂₃₅)
end
CFs_6_SOUTH_NORTH_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ) = CFs_6_SOUTH_NORTH_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], false);
CFs_6_SOUTH_NORTH_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, reduced::Bool) = CFs_6_SOUTH_NORTH_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], reduced);
CFs_6_SOUTH_NORTH_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, label::Vector{String}) = CFs_6_SOUTH_NORTH_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, label, false);

## 4₁-4₁: SOUTH-WEST
function CFs_6_SOUTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄₅ = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ*y₁, [l[4], l[5], l[1], l[2], l[3]], reduced)

    CF₁₂₅₆ = CF_3₂_3₁_chain(γ, h₁, h₂*x₂, x₁, δ, y₁*y₂, ℓ, [l[1], l[2], l[5], l[6]], reduced)
    CF₁₃₅₆ = CF_3₂_3₁_chain(γ, h₁, h₂, x₁*x₂, δ, y₁*y₂, ℓ, [l[1], l[3], l[5], l[6]], reduced)
    CF₂₃₅₆ = CF_3₂_3₁_chain(γ, h₁*x₁, h₂, x₂, δ, y₁*y₂, ℓ, [l[2], l[3], l[5], l[6]], reduced)
    CF₁₂₄₆ = CF_3₂_3₁_chain(γ, h₁, h₂*x₂, x₁, δ, y₁, ℓ, [l[1], l[2], l[4], l[6]], reduced)
    CF₁₃₄₆ = CF_3₂_3₁_chain(γ, h₁, h₂, x₁*x₂, δ, y₁, ℓ, [l[1], l[3], l[4], l[6]], reduced)
    CF₂₃₄₆ = CF_3₂_3₁_chain(γ, h₁*x₁, h₂, x₂, δ, y₁, ℓ, [l[2], l[3], l[4], l[6]], reduced)

    CF₁₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[6], l[1], l[4], l[5]], reduced)
    CF₂₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[6], l[2], l[4], l[5]], reduced)
    CF₃₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[6], l[3], l[4], l[5]], reduced)
    CF₁₂₃₆ = CF_4₁cycle(γ, x₁, x₂, [l[6], l[1], l[2], l[3]], reduced)

    merge(CF₁₂₃₄₅, CF₁₂₅₆, CF₁₃₅₆, CF₂₃₅₆, CF₁₂₄₆, CF₁₃₄₆, CF₂₃₄₆, CF₁₄₅₆, CF₂₄₅₆, CF₃₄₅₆, CF₁₂₃₆)
end
CFs_6_SOUTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ) = CFs_6_SOUTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], false);
CFs_6_SOUTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, reduced::Bool) = CFs_6_SOUTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], reduced);
CFs_6_SOUTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, label::Vector{String}) = CFs_6_SOUTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, label, false);

## 4₁-4₁: NORTH-WEST
function CFs_6_NORTH_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁₂₃₄₆ = CFs_5_NORTH(γ, x₁, x₂, ℓ*y₁, [l[1], l[2], l[3], l[4], l[6]], reduced)

    CF₁₂₃₅ = CF_3₁_3₁_chain(γ, x₁, δ, y₁, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    CF₁₂₄₅ = CF_3₁_3₁_chain(γ, x₁, δ, y₁*y₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF₁₃₅₆ = CF_3₁_3₁_chain(1-γ, x₂, δ, y₁, ℓ, [l[1], l[6], l[3], l[5]], reduced)
    CF₁₄₅₆ = CF_3₁_3₁_chain(1-γ, x₂, δ, y₁*y₂, ℓ, [l[1], l[6], l[4], l[5]], reduced)

    CF₁₂₅₆ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[5], l[6]], reduced)
    CF₁₃₄₅ = CF_4₁cycle(δ, y₁, y₂, [l[5], l[1], l[3], l[4]], reduced)
    CF₂₃₄₅ = CF_4₁cycle(δ, y₁, y₂, [l[5], l[2], l[3], l[4]], reduced)
    CF₃₄₅₆ = CF_4₁cycle(δ, y₁, y₂, [l[5], l[6], l[3], l[4]], reduced)

    CF₂₃₅₆ = CF_3₁cycle(1-δ, y₁, ℓ, [l[5], l[3], l[2], l[6]], reduced)
    CF₂₄₅₆ = CF_3₁cycle(1-δ, y₁*y₂, ℓ, [l[5], l[4], l[2], l[6]], reduced)

    merge(CF₁₂₃₄₆, CF₁₂₃₅, CF₁₂₄₅, CF₁₃₅₆, CF₁₄₅₆, CF₁₂₅₆, CF₁₃₄₅, CF₂₃₄₅, CF₃₄₅₆, CF₂₃₅₆, CF₂₄₅₆)
end
CFs_6_NORTH_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ) = CFs_6_NORTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], false);
CFs_6_NORTH_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, reduced::Bool) = CFs_6_NORTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], reduced);
CFs_6_NORTH_WEST_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, label::Vector{String}) = CFs_6_NORTH_WEST_chain(γ, h₁, h₂, x₁, x₂, δ, y₁, y₂, ℓ, label, false);

## 4₁-4₁: NORTH-NORTH
function CFs_6_NORTH_NORTH_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₂₃₄₅₆ = CFs_5_NORTH(δ, y₁, y₂, ℓ, [l[4], l[5], l[2], l[6], l[3]], reduced)

    CF₁₂₃₄ = CF_3₁_3₁_chain(γ, x₁, δ, y₂, ℓ,      [l[1], l[2], l[3], l[4]], reduced)
    CF₁₂₄₅ = CF_3₁_3₁_chain(γ, x₁, 1-δ, y₁, ℓ,   [l[1], l[2], l[5], l[4]], reduced)
    CF₁₃₄₆ = CF_3₁_3₁_chain(1-γ, x₂, δ, y₂, ℓ,    [l[1], l[6], l[3], l[4]], reduced)
    CF₁₄₅₆ = CF_3₁_3₁_chain(1-γ, x₂, 1-δ, y₁, ℓ, [l[1], l[6], l[5], l[4]], reduced)

    CF₁₂₃₆ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[3], l[6]], reduced)
    CF₁₂₄₆ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[4], l[6]], reduced)
    CF₁₂₅₆ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[5], l[6]], reduced)
    CF₁₃₄₅ = CF_4₁cycle(δ, y₁, y₂, [l[4], l[5], l[1], l[3]], reduced)

    CF₁₂₃₅ = CF_3₁cycle(γ, x₁, ℓ,  [l[1], l[2], l[3], l[5]], reduced)
    CF₁₃₅₆ = CF_3₁cycle(1-γ, x₂, ℓ,[l[1], l[6], l[3], l[5]], reduced)

    merge(CF₂₃₄₅₆, CF₁₂₃₄, CF₁₂₄₅, CF₁₃₄₆, CF₁₄₅₆, CF₁₂₃₆, CF₁₂₄₆, CF₁₂₅₆, CF₁₃₄₅, CF₁₂₃₅, CF₁₃₅₆)
end
CFs_6_NORTH_NORTH_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ) = CFs_6_NORTH_NORTH_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], false);
CFs_6_NORTH_NORTH_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, reduced::Bool) = CFs_6_NORTH_NORTH_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, ["a", "b", "c", "d", "e", "f"], reduced);
CFs_6_NORTH_NORTH_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, label::Vector{String}) = CFs_6_NORTH_NORTH_chain(γ, x₁, x₂, δ, y₁, y₂, ℓ, label, false);


