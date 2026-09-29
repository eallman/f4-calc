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
function CFs_5_3₁cycle_com(γ, x, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    C₅ = CF_tree(x*ℓ₁*ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    C₂ = CF_3₁cycle_com(γ, x, ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    C₄ = CF_3₁cycle_com(1-γ, x, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    merge(C₂, C₄, C₅)
end
CFs_5_3₁cycle_com(γ, x, ℓ₁, ℓ₂) = CFs_5_3₁cycle_com(γ, x, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], false);
CFs_5_3₁cycle_com(γ, x, ℓ₁, ℓ₂, reduced::Bool) = CFs_5_3₁cycle_com(γ, x, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], reduced);
CFs_5_3₁cycle_com(γ, x, ℓ₁, ℓ₂, label::Vector{String}) = CFs_5_3₁cycle_com(γ, x, ℓ₁, ℓ₂, label, false);

## 3₂-cycle
function CFs_5_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, l::Vector{String}, reduced::Bool)
    C₂ = CF_3₁cycle_com(1-γ, x, ℓ₂, [l[1], l[5], l[3], l[4]], reduced)
    C₄ = CF_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, [l[1], l[2], l[3], l[5]], reduced)
    C₅ = CF_2cycle_com(γ, h₁, h₂*x, ℓ₁, ℓ₂, [l[1], l[2], l[3], l[4]], reduced)
    merge(C₂, C₄, C₅)
end
CFs_5_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂) = CFs_5_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], false);
CFs_5_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, reduced::Bool) = CFs_5_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, ["a₁", "a₂", "b₁", "b₂", "c"], reduced);
CFs_5_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, label::Vector{String}) = CFs_5_3₂cycle_com(γ, h₁, h₂, x, ℓ₁, ℓ₂, label, false);

## 4₂-cycle: SOUTH
function CFs_5_SOUTH_com(γ, h₁, h₂, x₁, x₂, ℓ, l::Vector{String}, reduced::Bool)
    CF_abcd = CF_4₁cycle_com(γ, x₁, x₂, [l[1], l[3], l[4], l[5]], reduced)
    #CF_a₁a₂cd = CF_3₂cycle_com(γ, h₁*x₁, x₂, h₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    #CF_a₁a₂bd = CF_3₂cycle_com(γ, h₁, x₁*x₂, h₂, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    #CF_a₁a₂bc = CF_3₂cycle_com(γ, h₁, x₁, x₂*h₂, ℓ, [l[1], l[2], l[3], l[4]], reduced)

    CF_a₁a₂cd = CF_3₂cycle_com(γ, h₁*x₁, h₂, x₂, ℓ, [l[1], l[2], l[4], l[5]], reduced)
    CF_a₁a₂bd = CF_3₂cycle_com(γ, h₁, h₂, x₁*x₂, ℓ, [l[1], l[2], l[3], l[5]], reduced)
    CF_a₁a₂bc = CF_3₂cycle_com(γ, h₁, x₂*h₂, x₁, ℓ, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF_abcd, CF_a₁a₂cd, CF_a₁a₂bd, CF_a₁a₂bc)
end
CFs_5_SOUTH_com(γ, h₁, h₂, x₁, x₂, ℓ) = CFs_5_SOUTH_com(γ, h₁, h₂, x₁, x₂, ℓ, ["a₁", "a₂", "b", "c", "d"], false);
CFs_5_SOUTH_com(γ, h₁, h₂, x₁, x₂, ℓ, reduced::Bool) = CFs_5_SOUTH_com(γ, h₁, h₂, x₁, x₂, ℓ, ["a₁", "a₂", "b", "c", "d"], reduced);
CFs_5_SOUTH_com(γ, h₁, h₂, x₁, x₂, ℓ, label::Vector{String}) = CFs_5_SOUTH_com(γ, h₁, h₂, x₁, x₂, ℓ, label, false);

## 4₁-cycle: WEST
function CFs_5_WEST_com(γ, h₁, h₂, x₁, x₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁ = CF_tree(ℓ*x₁, [l[2], l[3], l[4], l[5]], reduced)
    CF₂ = CF_4₁cycle_com(γ, x₁, x₂, [l[1], l[2], l[4], l[5]], reduced)
    CF₃ = CF_3₁cycle_com(1-γ, x₁*x₂, ℓ, [l[1], l[5], l[2], l[3]], reduced)
    CF₄ = CF_3₁cycle_com(1-γ, x₁, ℓ, [l[1], l[4], l[2], l[3]], reduced)

    merge(CF₁, CF₂, CF₃, CF₄)
end
CFs_5_WEST_com(γ, h₁, h₂, x₁, x₂, ℓ) = CFs_5_WEST_com(γ, h₁, h₂, x₁, x₂, ℓ, ["a", "b₁", "b₂", "c", "d"], false);
CFs_5_WEST_com(γ, h₁, h₂, x₁, x₂, ℓ, reduced::Bool) = CFs_5_WEST_com(γ, h₁, h₂, x₁, x₂, ℓ, ["a", "b₁", "b₂", "c", "d"], reduced);
CFs_5_WEST_com(γ, h₁, h₂, x₁, x₂, ℓ, label::Vector{String}) = CFs_5_WEST_com(γ, h₁, h₂, x₁, x₂, ℓ, label, false);

## 4₁-cycle: EAST
function CFs_5_EAST_com(γ, h₁, h₂, x₁, x₂, ℓ, l::Vector{String}, reduced::Bool)
    CFs_4_WEST(1-γ, h₂, h₁, x₂, x₁, ℓ, [l[1], l[4], l[5], l[3], l[2]], reduced)
end
CFs_5_EAST_com(γ, h₁, h₂, x₁, x₂, ℓ) = CFs_5_EAST_com(γ, h₁, h₂, x₁, x₂, ℓ, ["a", "b", "c", "d₁", "d₂"], false);
CFs_5_EAST_com(γ, h₁, h₂, x₁, x₂, ℓ, reduced::Bool) = CFs_5_EAST_com(γ, h₁, h₂, x₁, x₂, ℓ, ["a", "b", "c", "d₁", "d₂"], reduced);
CFs_5_EAST_com(γ, h₁, h₂, x₁, x₂, ℓ, label::Vector{String}) = CFs_5_EAST_com(γ, h₁, h₂, x₁, x₂, ℓ, label, false);

## 4₁-cycle: NORTH
function CFs_5_NORTH_com(γ, x₁, x₂, ℓ, l::Vector{String}, reduced::Bool)
    CF₁ = CF_tree(ℓ, [l[3], l[4], l[2], l[5]], reduced)
    CF₂ = CF_3₁cycle_com(1-γ, x₂, ℓ, [l[1], l[5], l[3], l[4]], reduced)
    CF₃ = CF_4₁cycle_com(γ, x₁, x₂, [l[1], l[2], l[3], l[5]], reduced)
    CF₄ = CF_3₁cycle_com(γ, x₁, ℓ, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF₁, CF₂, CF₃, CF₄)
end
CFs_5_NORTH_com(γ, x₁, x₂, ℓ) = CFs_5_NORTH_com(γ, x₁, x₂, ℓ, ["a", "b", "c₁", "c₂", "d"], false);
CFs_5_NORTH_com(γ, x₁, x₂, ℓ, reduced::Bool) = CFs_5_NORTH_com(γ, x₁, x₂, ℓ, ["a", "b", "c₁", "c₂", "d"], reduced);
CFs_5_NORTH_com(γ, x₁, x₂, ℓ, label::Vector{String}) = CFs_5_NORTH_com(γ, x₁, x₂, ℓ, label, false);

## 5₁-cycle
function CF_5₁cycle_com(γ, h₁, h₂, x₁, x₂, x₃, l::Vector{String}, reduced::Bool)
    CF₁ = CF_tree(x₂, [l[2], l[3], l[4], l[5]], reduced)
    CF₂ = CF_4₁cycle_com(γ, x₂, x₃, [l[1], l[3], l[4], l[5]], reduced)
    CF₃ = CF_4₁cycle_com(γ, x₁*x₂, x₃, [l[1], l[2], l[4], l[5]], reduced)
    CF₄ = CF_4₁cycle_com(γ, x₁, x₂*x₃, [l[1], l[2], l[3], l[5]], reduced)
    CF₅ = CF_4₁cycle_com(γ, x₂, x₃, [l[1], l[2], l[3], l[4]], reduced)

    merge(CF₁, CF₂, CF₃, CF₄, CF₅)
end
CF_5₁cycle_com(γ, h₁, h₂, x₁, x₂, x₃) = CF_5₁cycle_com(γ, h₁, h₂, x₁, x₂, x₃, ["a", "b", "c", "d", "e"], false);
CF_5₁cycle_com(γ, h₁, h₂, x₁, x₂, x₃, reduced::Bool) = CF_5₁cycle_com(γ, h₁, h₂, x₁, x₂, x₃, ["a", "b", "c", "d", "e"], reduced);
CF_5₁cycle_com(γ, h₁, h₂, x₁, x₂, x₃, label::Vector{String}) = CF_5₁cycle_com(γ, h₁, h₂, x₁, x₂, x₃, label, false);

