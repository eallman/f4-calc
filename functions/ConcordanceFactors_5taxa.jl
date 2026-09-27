include("ConcordanceFactors_utilities.jl")
include("ConcordanceFactors_4taxa.jl")

## Tree
function CFs_5_tree(ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    C₂ = CF_tree(ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    C₄ = CF_tree(ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    C₅ = CF_tree(ℓ₁*ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    merge(C₂, C₄, C₅)
end
CFs_5_tree(ℓ₁, ℓ₂) = CFs_5_tree(ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], false);
CFs_5_tree(ℓ₁, ℓ₂, reduced::Bool) = CFs_5_tree(ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], reduced);
CFs_5_tree(ℓ₁, ℓ₂, label::Vector{String}) = CFs_5_tree(ℓ₁, ℓ₂, label, false);

## 3₁-cycle
function CFs_5_3₁cycle(γ, x, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    C₅ = CF_tree(x*ℓ₁*ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    C₂ = CF_3₁cycle(γ, x, ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    C₄ = CF_3₁cycle(1-γ, x, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    merge(C₂, C₄, C₅)
end
CFs_5_3₁cycle(γ, x, ℓ₁, ℓ₂) = CFs_5_3₁cycle(γ, x, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], false);
CFs_5_3₁cycle(γ, x, ℓ₁, ℓ₂, reduced::Bool) = CFs_5_3₁cycle(γ, x, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], reduced);
CFs_5_3₁cycle(γ, x, ℓ₁, ℓ₂, label::Vector{String}) = CFs_5_3₁cycle(γ, x, ℓ₁, ℓ₂, label, false);

## 3₂-cycle
function CFs_5_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    C₂ = CF_3₁cycle(1-γ, x, ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    C₄ = CF_3₂cycle(γ, h₁, h₂, x, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    C₅ = CF_2cycle(γ, h₁, h₂*x, ℓ₁, ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    merge(C₂, C₄, C₅)
end
CFs_5_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂) = CFs_5_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], false);
CFs_5_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, reduced::Bool) = CFs_5_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], reduced);
CFs_5_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, label::Vector{String}) = CFs_5_3₂cycle(γ, h₁, h₂, x, ℓ₁, ℓ₂, label, false);

## 4₂-cycle: SOUTH
function CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ, l::Vector{String}, reduced::Bool)
    CF_abcd = CF_4₁cycle(γ, x₁, x₂, [l[1], l[3], l[4], l[5]], reduced)
    #CF_a₁a₂cd = CF_3₂cycle(γ, h₁*x₁, x₂, h₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    #CF_a₁a₂bd = CF_3₂cycle(γ, h₁, x₁*x₂, h₂, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    #CF_a₁a₂bc = CF_3₂cycle(γ, h₁, x₁, x₂*h₂, ℓ, [l[1], l[2], l[3], l[4]], reduced)

    CF_a₁a₂cd = CF_3₂cycle(γ, h₁*x₁, h₂, x₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF_a₁a₂bd = CF_3₂cycle(γ, h₁, h₂, x₁*x₂, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    CF_a₁a₂bc = CF_3₂cycle(γ, h₁, x₂*h₂, x₁, ℓ, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF_abcd, CF_a₁a₂cd, CF_a₁a₂bd, CF_a₁a₂bc)
end
CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ) = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ, ["a₁", "a₂", "b", "c", "d"], false);
CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ, reduced::Bool) = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ, ["a₁", "a₂", "b", "c", "d"], reduced);
CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ, label::Vector{String}) = CFs_5_SOUTH(γ, h₁, h₂, x₁, x₂, ℓ, label, false);

## 4₁-cycle: WEST
function CFs_5_WEST(γ, h₁, h₂, x₁, x₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁ = CF_tree(ℓ*x₁, [l[2], l[3], l[4], l[5]], reduced)
    CF₂ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[4], l[5]], reduced)
    CF₃ = CF_3₁cycle(1-γ, x₁*x₂, ℓ, [l[1], l[5], l[2], l[3]], reduced)
    CF₄ = CF_3₁cycle(1-γ, x₁, ℓ, [l[1], l[4], l[2], l[3]], reduced)

    merge(CF₁, CF₂, CF₃, CF₄)
end
CFs_5_WEST(γ, h₁, h₂, x₁, x₂, ℓ) = CFs_5_WEST(γ, h₁, h₂, x₁, x₂, ℓ, ["a", "b₁", "b₂", "c", "d"], false);
CFs_5_WEST(γ, h₁, h₂, x₁, x₂, ℓ, reduced::Bool) = CFs_5_WEST(γ, h₁, h₂, x₁, x₂, ℓ, ["a", "b₁", "b₂", "c", "d"], reduced);
CFs_5_WEST(γ, h₁, h₂, x₁, x₂, ℓ, label::Vector{String}) = CFs_5_WEST(γ, h₁, h₂, x₁, x₂, ℓ, label, false);

## 4₁-cycle: EAST
function CFs_5_EAST(γ, h₁, h₂, x₁, x₂, ℓ, l::Vector{String}, reduced::Bool)
    CFs_4_WEST(1-γ, h₂, h₁, x₂, x₁, ℓ, [l[1], l[4], l[5], l[3], l[2]], reduced)
end
CFs_5_EAST(γ, h₁, h₂, x₁, x₂, ℓ) = CFs_5_EAST(γ, h₁, h₂, x₁, x₂, ℓ, ["a", "b", "c", "d₁", "d₂"], false);
CFs_5_EAST(γ, h₁, h₂, x₁, x₂, ℓ, reduced::Bool) = CFs_5_EAST(γ, h₁, h₂, x₁, x₂, ℓ, ["a", "b", "c", "d₁", "d₂"], reduced);
CFs_5_EAST(γ, h₁, h₂, x₁, x₂, ℓ, label::Vector{String}) = CFs_5_EAST(γ, h₁, h₂, x₁, x₂, ℓ, label, false);

## 4₁-cycle: NORTH
function CFs_5_NORTH(γ, x₁, x₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁ = CF_tree(ℓ, [l[3], l[4], l[2], l[5]], reduced)
    CF₂ = CF_3₁cycle(1-γ, x₂, ℓ, [l[1], l[5], l[3], l[4]], reduced)
    CF₃ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[2], l[3], l[5]], reduced)
    CF₄ = CF_3₁cycle(γ, x₁, ℓ, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF₁, CF₂, CF₃, CF₄)
end
CFs_5_NORTH(γ, x₁, x₂, ℓ) = CFs_5_NORTH(γ, x₁, x₂, ℓ, ["a", "b", "c₁", "c₂", "d"], false);
CFs_5_NORTH(γ, x₁, x₂, ℓ, reduced::Bool) = CFs_5_NORTH(γ, x₁, x₂, ℓ, ["a", "b", "c₁", "c₂", "d"], reduced);
CFs_5_NORTH(γ, x₁, x₂, ℓ, label::Vector{String}) = CFs_5_NORTH(γ, x₁, x₂, ℓ, label, false);

## 5₁-cycle
function CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃, l::Vector{String}, reduced::Bool)
    CF₁ = CF_tree(x₂, [l[2], l[3], l[4], l[5]], reduced)
    CF₂ = CF_4₁cycle(γ, x₂, x₃, [l[1], l[3], l[4], l[5]], reduced)
    CF₃ = CF_4₁cycle(γ, x₁*x₂, x₃, [l[1], l[2], l[4], l[5]], reduced)
    CF₄ = CF_4₁cycle(γ, x₁, x₂*x₃, [l[1], l[2], l[3], l[5]], reduced)
    CF₅ = CF_4₁cycle(γ, x₂, x₃, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF₁, CF₂, CF₃, CF₄, CF₅)
end
CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃) = CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃, ["a", "b", "c", "d", "e"], false);
CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃, reduced::Bool) = CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃, ["a", "b", "c", "d", "e"], reduced);
CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃, label::Vector{String}) = CF_5₁cycle(γ, h₁, h₂, x₁, x₂, x₃, label, false);


## ## ## CHAINS ## ## ##
##

## 3₁-3₁-cycles:
function CFs_5_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    C₅ = CF_3₁cycle(1-δ, y, x*ℓ₁*ℓ₂, [l[4], l[3], l[1], l[2]], reduced)
    C₂ = CF_3₁_3₁_chain(γ, x, δ, y, ℓ₂, [l[5], l[1], l[3], l[4]], reduced)
    C₄ = CF_3₁cycle(1-γ, x, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    merge(C₂, C₄, C₅)
end
CFs_5_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂) = CFs_5_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], false);
CFs_5_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, reduced::Bool) = CFs_5_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], reduced);
CFs_5_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, label::Vector{String}) = CFs_5_3₁_3₁_chain(γ, x, δ, y, ℓ₁, ℓ₂, label, false);

# ## 3₂-3₁-cycles:
function CFs_5_3₁_3₂_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    C₅ = CF_3₂cycle(δ, k₁, k₂, y, ℓ₁, [l[1], l[2], l[3], l[4]], reduced)
    C₂ = CF_3₂_3₁_chain(δ, k₁, k₂, y, 1-γ, x, ℓ₂, [l[3], l[4], l[1], l[5]], reduced)
    C₄ = CF_3₁cycle(1-γ, x, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    merge(C₂, C₄, C₅)
end
CFs_5_3₁_3₂_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂) = CFs_5_3₁_3₂_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], false);
CFs_5_3₁_3₂_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, reduced::Bool) = CFs_5_3₁_3₂_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], reduced);
CFs_5_3₁_3₂_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, label::Vector{String}) = CFs_5_3₁_3₂_chain(γ, x, δ, k₁, k₂, y, ℓ₁, ℓ₂, label, false);



## 4₁-3₁-cycles: WEST_3₁
function CFs_5_WEST_3₁(γ, x₁, x₂, δ, y, ℓ, l::Vector{String}, reduced::Bool)
    CF₂ = CF_4₁cycle(γ, x₁, x₂, [l[5], l[1], l[3], l[4]], reduced)
    CF₃ = CF_3₁_3₁_chain(δ, y, γ, x₁*x₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF₄ = CF_3₁_3₁_chain(δ, y, γ, x₁, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    CF₅ = CF_3₁cycle(δ, y, ℓ*x₁, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF₂, CF₃, CF₄, CF₅)
end
CFs_5_WEST_3₁(γ, x₁, x₂, δ, y, ℓ) = CFs_5_WEST_3₁(γ, x₁, x₂, δ, y, ℓ, ["a₁", "a₂", "b", "c", "d"], false);
CFs_5_WEST_3₁(γ, x₁, x₂, δ, y, ℓ, reduced::Bool) = CFs_5_WEST_3₁(γ, x₁, x₂, δ, y, ℓ, ["a₁", "a₂", "b", "c", "d"], reduced);
CFs_5_WEST_3₁(γ, x₁, x₂, δ, y, ℓ, label::Vector{String}) = CFs_5_WEST_3₁(γ, x₁, x₂, δ, y, ℓ, label, false);

## 4₁-3₂-cycles: WEST_3₂
function CFs_5_WEST_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, l::Vector{String}, reduced::Bool)
    CF₂ = CF_4₁cycle(γ, x₁, x₂, [l[5], l[1], l[3], l[4]], reduced)
    CF₃ = CF_3₂_3₁_chain(δ, k₁, k₂, y, γ, x₁*x₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF₄ = CF_3₂_3₁_chain(δ, k₁, k₂, y, γ, x₁, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    CF₅ = CF_3₂cycle(δ, k₁, k₂, y, x₁*ℓ, [l[3], l[4], l[1], l[2]], reduced)
    merge(CF₂, CF₃, CF₄, CF₅)
end
CFs_5_WEST_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ) = CFs_5_WEST_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, ["a₁", "a₂", "b", "c", "d"], false);
CFs_5_WEST_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, reduced::Bool) = CFs_5_WEST_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, ["a₁", "a₂", "b", "c", "d"], reduced);
CFs_5_WEST_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, label::Vector{String}) = CFs_5_WEST_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, label, false);

## 4₁-3₁-cycles: NORTH_3₁
function CFs_5_NORTH_3₁(γ, x₁, x₂, δ, y, ℓ, l::Vector{String}, reduced::Bool)
    CF₂ = CF_4₁cycle(γ, x₁, x₂, [l[4], l[5], l[1], l[3]], reduced)
    CF₃ = CF_3₁_3₁_chain(δ, y, 1-γ, x₁, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF₅ = CF_3₁_3₁_chain(δ, y, γ, x₂, ℓ, [l[1], l[2], l[3], l[4]], reduced)
    CF₄ = CF_3₁cycle(δ, y, ℓ, [l[1], l[2], l[3], l[5]], reduced)

    merge(CF₂, CF₃, CF₄, CF₅)
end
CFs_5_NORTH_3₁(γ, x₁, x₂, δ, y, ℓ) = CFs_5_NORTH_3₁(γ, x₁, x₂, δ, y, ℓ, ["a₁", "a₂", "b", "c", "d"], false);
CFs_5_NORTH_3₁(γ, x₁, x₂, δ, y, ℓ, reduced::Bool) = CFs_5_NORTH_3₁(γ, x₁, x₂, δ, y, ℓ, ["a₁", "a₂", "b", "c", "d"], reduced);
CFs_5_NORTH_3₁(γ, x₁, x₂, δ, y, ℓ, label::Vector{String}) = CFs_5_NORTH_3₁(γ, x₁, x₂, δ, y, ℓ, label, false);

## 4₁-3₂-cycles: NORTH_3₂
function CFs_5_NORTH_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, l::Vector{String}, reduced::Bool)
    CF₂ = CF_4₁cycle(γ, x₁, x₂, [l[4], l[5], l[1], l[3]], reduced)
    CF₃ = CF_3₂_3₁_chain(δ, k₁, k₂, y, 1-γ, x₁, ℓ, [l[1], l[2], l[5], l[4]], reduced)
    CF₅ = CF_3₂_3₁_chain(δ, k₁, k₂, y, γ, x₂, ℓ, [l[1], l[2], l[3], l[4]], reduced)
    CF₄ = CF_3₂cycle(δ, k₁, k₂, y, ℓ, [l[3], l[5], l[1], l[2]], reduced)

    merge(CF₂, CF₃, CF₄, CF₅)
end
CFs_5_NORTH_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ) = CFs_5_NORTH_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, ["a₁", "a₂", "b", "c", "d"], false);
CFs_5_NORTH_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, reduced::Bool) = CFs_5_NORTH_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, ["a₁", "a₂", "b", "c", "d"], reduced);
CFs_5_NORTH_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, label::Vector{String}) = CFs_5_NORTH_3₂(γ, x₁, x₂, δ, k₁, k₂, y, ℓ, label, false);

## 4₂-3₁-cycles: SOUTH_3₁
function CFs_5_SOUTH_3₁(γ, h₁, h₂, x₁, x₂, δ, y, ℓ, l::Vector{String}, reduced::Bool)
    CF₂ = CF_4₁cycle(γ, x₁, x₂, [l[1], l[3], l[4], l[5]], reduced)
    CF₃ = CF_3₂_3₁_chain(γ, h₁*x₁, h₂, x₂, δ, y, ℓ, [l[4], l[5], l[1], l[2]], reduced)
    CF₅ = CF_3₂_3₁_chain(γ, h₁, h₂, x₁*x₂, δ, y, ℓ, [l[3], l[5], l[1], l[2]], reduced)
    CF₄ = CF_3₂_3₁_chain(γ, h₁, h₂*x₂, x₁, δ, y, ℓ, [l[3], l[4], l[1], l[2]], reduced)

    merge(CF₂, CF₃, CF₄, CF₅)
end
CFs_5_SOUTH_3₁(γ, h₁, h₂, x₁, x₂, δ, y, ℓ) = CFs_5_SOUTH_3₁(γ, h₁, h₂, x₁, x₂, δ, y, ℓ, ["a₁", "a₂", "b", "c", "d"], false);
CFs_5_SOUTH_3₁(γ, h₁, h₂, x₁, x₂, δ, y, ℓ, reduced::Bool) = CFs_5_SOUTH_3₁(γ, h₁, h₂, x₁, x₂, δ, y, ℓ, ["a₁", "a₂", "b", "c", "d"], reduced);
CFs_5_SOUTH_3₁(γ, h₁, h₂, x₁, x₂, δ, y, ℓ, label::Vector{String}) = CFs_5_SOUTH_3₁(γ, h₁, h₂, x₁, x₂, δ, åy, ℓ, label, false);

