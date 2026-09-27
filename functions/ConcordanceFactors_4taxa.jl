include("ConcordanceFactors_utilities.jl")

## Tree
function CF_tree(x, labels::Vector{String}, reduced::Bool)
    sp = splits(labels)
    CF = Dict(sp[1] => 1-2//3*x,
              sp[2] => 1//3*x,
              sp[3] => 1//3*x)
    returnCF_treeLike(CF, sp, reduced)
end
CF_tree(x) = CF_tree(x, ["a", "b", "c", "d"], false);
CF_tree(x, reduced::Bool) = CF_tree(x, ["a", "b", "c", "d"], reduced);
CF_tree(x, label::Vector{String}) = CF_tree(x, label, false);

## 2-cycle
function CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂, labels::Vector{String}, reduced::Bool)
    sp = splits(labels)
    CF = Dict(sp[1] => (1-γ)^2*(1-2//3*ℓ₁*ℓ₂*h₂) + 2*γ*(1-γ)*(1-2//3*ℓ₁*ℓ₂) + γ^2*(1-2//3*ℓ₁*ℓ₂*h₁),
              sp[2] => (1-γ)^2*(1//3*ℓ₁*ℓ₂*h₂) + 2*γ*(1-γ)*1//3*ℓ₁*ℓ₂ + γ^2*(1//3*ℓ₁*ℓ₂*h₁),
              sp[3] => (1-γ)^2*(1//3*ℓ₁*ℓ₂*h₂) + 2*γ*(1-γ)*1//3*ℓ₁*ℓ₂ + γ^2*(1//3*ℓ₁*ℓ₂*h₁))

    returnCF_treeLike(CF, sp, reduced)
end
CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂) = CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂, ["a", "b", "c", "d"], false);
CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂, reduced::Bool) = CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂, ["a", "b", "c", "d"], reduced);
CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂, label::Vector{String}) = CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂, label, false);

## 3₁-cycle
function CF_3₁cycle(γ, x, ℓ, labels::Vector{String}, reduced::Bool)
    CF_tree_γ = CF_tree(x*ℓ, labels)
    CF_tree_1minusγ = CF_tree(ℓ, labels)

    sp = splits(labels)
    CF = Dict(sp[1] => γ*CF_tree_γ[sp[1]] + (1-γ)*CF_tree_1minusγ[sp[1]],
              sp[2] => γ*CF_tree_γ[sp[2]] + (1-γ)*CF_tree_1minusγ[sp[2]],
              sp[3] => γ*CF_tree_γ[sp[3]] + (1-γ)*CF_tree_1minusγ[sp[3]])

    returnCF_treeLike(CF, sp, reduced)
end
CF_3₁cycle(γ, x, ℓ) = CF_3₁cycle(γ, x, ℓ, ["a", "b", "c", "d"], false);
CF_3₁cycle(γ, x, ℓ, reduced::Bool) = CF_3₁cycle(γ, x, ℓ, ["a", "b", "c", "d"], reduced);
CF_3₁cycle(γ, x, ℓ, label::Vector{String}) = CF_3₁cycle(γ, x, ℓ, label, false);

## 3₂-cycle
function CF_3₂cycle(γ, h₁, h₂, x, ℓ, labels::Vector{String}, reduced::Bool)
    CF_2_γ = CF_tree(h₁, labels)
    CF_2_1minusγ = CF_tree(h₂, labels)
    CF_γ_1minusγ = CF_tree(x, [labels[1],labels[3],labels[2],labels[4]])
    CF_1minusγ_γ = CF_tree(x, [labels[1],labels[4],labels[2],labels[3]])

    sp = splits(labels)

    CF = Dict(sp[1] => (1-ℓ) + ℓ*(γ^2*CF_2_γ[sp[1]] + (1-γ)^2*CF_2_1minusγ[sp[1]] + γ*(1-γ)*(CF_γ_1minusγ[sp[1]] + CF_1minusγ_γ[sp[1]])),
              sp[2] => ℓ*(γ^2*CF_2_γ[sp[2]] + (1-γ)^2*CF_2_1minusγ[sp[2]] + γ*(1-γ)*(CF_γ_1minusγ[sp[2]] + CF_1minusγ_γ[sp[2]])),
              sp[3] => ℓ*(γ^2*CF_2_γ[sp[3]] + (1-γ)^2*CF_2_1minusγ[sp[3]] + γ*(1-γ)*(CF_γ_1minusγ[sp[3]] + CF_1minusγ_γ[sp[3]])))

    returnCF_treeLike(CF, sp, reduced)
end
CF_3₂cycle(γ, h₁, h₂, x, ℓ) = CF_3₂cycle(γ, h₁, h₂, x, ℓ, ["a", "b", "c", "d"], false);
CF_3₂cycle(γ, h₁, h₂, x, ℓ, reduced::Bool) = CF_3₂cycle(γ, h₁, h₂, x, ℓ, ["a", "b", "c", "d"], reduced);
CF_3₂cycle(γ, h₁, h₂, x, ℓ, label::Vector{String}) = CF_3₂cycle(γ, h₁, h₂, x, ℓ, label, false);

## 4₁-cycle
function CF_4₁cycle(γ, x₁, x₂, labels::Vector{String}, reduced::Bool)
    CF_tree_γ = CF_tree(x₁, labels)
    CF_tree_1minusγ = CF_tree(x₂, [labels[1],labels[4],labels[2],labels[3]])

    sp = splits(labels)
    CF = Dict(sp[1] => γ*CF_tree_γ[sp[1]] + (1-γ)*CF_tree_1minusγ[sp[1]],
              sp[2] => γ*CF_tree_γ[sp[2]] + (1-γ)*CF_tree_1minusγ[sp[2]],
              sp[3] => γ*CF_tree_γ[sp[3]] + (1-γ)*CF_tree_1minusγ[sp[3]])

    if reduced
        key = sort([sp[2], sp[3]])[1]
        Dict(sp[1] => CF[sp[1]], key => CF[key])
    else
      CF
    end
end
CF_4₁cycle(γ, x₁, x₂) = CF_4₁cycle(γ, x₁, x₂, ["a", "b", "c", "d"], false);
CF_4₁cycle(γ, x₁, x₂, reduced::Bool) = CF_4₁cycle(γ, x₁, x₂, ["a", "b", "c", "d"], reduced);
CF_4₁cycle(γ, x₁, x₂, label::Vector{String}) = CF_4₁cycle(γ, x₁, x₂, label, false);

########################
## ## ## CHAINS ## ## ##
########################

## 3₁-3₁-cycles
function CF_3₁_3₁_chain(γ, x, δ, y, ℓ, labels::Vector{String}, reduced::Bool)
    CF_tree_γ_δ = CF_tree(x*ℓ, labels)
    CF_tree_1minusγ_δ = CF_tree(ℓ, labels)
    CF_tree_γ_1minusδ = CF_tree(x*ℓ*y, labels)
    CF_tree_1minusγ_1minusδ = CF_tree(ℓ*y, labels)

    sp = splits(labels)
    CF = Dict(sp[1] => γ*δ*CF_tree_γ_δ[sp[1]] + (1-γ)*δ*CF_tree_1minusγ_δ[sp[1]] + γ*(1-δ)*CF_tree_γ_1minusδ[sp[1]] + (1-γ)*(1-δ)*CF_tree_1minusγ_1minusδ[sp[1]],
              sp[2] => γ*δ*CF_tree_γ_δ[sp[2]] + (1-γ)*δ*CF_tree_1minusγ_δ[sp[2]] + γ*(1-δ)*CF_tree_γ_1minusδ[sp[2]] + (1-γ)*(1-δ)*CF_tree_1minusγ_1minusδ[sp[2]],
              sp[3] => γ*δ*CF_tree_γ_δ[sp[3]] + (1-γ)*δ*CF_tree_1minusγ_δ[sp[3]] + γ*(1-δ)*CF_tree_γ_1minusδ[sp[3]] + (1-γ)*(1-δ)*CF_tree_1minusγ_1minusδ[sp[3]])

    returnCF_treeLike(CF, sp, reduced)
end
CF_3₁_3₁_chain(γ, x, δ, y, ℓ) = CF_3₁_3₁_chain(γ, x, δ, y, ℓ, ["a", "b", "c", "d"], false);
CF_3₁_3₁_chain(γ, x, δ, y, ℓ, reduced::Bool) = CF_3₁_3₁_chain(γ, x, δ, y, ℓ, ["a", "b", "c", "d"], reduced);
CF_3₁_3₁_chain(γ, x, δ, y, ℓ, label::Vector{String}) = CF_3₁_3₁_chain(γ, x, δ, y, ℓ, label, false);


## 3₂-3₁-cycles
function CF_3₂_3₁_chain(γ, h₁, h₂, x, δ, y, ℓ, labels::Vector{String}, reduced::Bool)
    CF_3₂_δ = CF_3₂cycle(γ, h₁, h₂, x, ℓ, labels)
    CF_3₂_1minusδ = CF_3₂cycle(γ, h₁, h₂, x, ℓ*y, labels)

    sp = splits(labels)
    CF = Dict(sp[1] => δ*CF_3₂_δ[sp[1]] + (1-δ)*CF_3₂_1minusδ[sp[1]],
              sp[2] => δ*CF_3₂_δ[sp[2]] + (1-δ)*CF_3₂_1minusδ[sp[2]],
              sp[3] => δ*CF_3₂_δ[sp[3]] + (1-δ)*CF_3₂_1minusδ[sp[3]])

    returnCF_treeLike(CF, sp, reduced)
end
CF_3₂_3₁_chain(γ, h₁, h₂, x, δ, y, ℓ) = CF_3₂_3₁_chain(γ, h₁, h₂, x, δ, y, ℓ, ["a", "b", "c", "d"], false);
CF_3₂_3₁_chain(γ, h₁, h₂, x, δ, y, ℓ, reduced::Bool) = CF_3₂_3₁_chain(γ, h₁, h₂, x, δ, y, ℓ, ["a", "b", "c", "d"], reduced);
CF_3₂_3₁_chain(γ, h₁, h₂, x, δ, y, ℓ, label::Vector{String}) = CF_3₂_3₁_chain(γ, h₁, h₂, x, δ, y, ℓ, label, false);

## 2₁-3₁-cycles
function CF_2₁_3₁_chain(γ, h₁, h₂, ℓ₁, ℓ₂, δ, y, labels::Vector{String}, reduced::Bool)
    CF_2cycle_δ = CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂, labels)
    CF_2cycle_1minusδ = CF_2cycle(γ, h₁, h₂, ℓ₁, ℓ₂*y, labels)

    sp = splits(labels)
    CF = Dict(sp[1] => δ*CF_2cycle_δ[sp[1]] + (1-δ)*CF_2cycle_1minusδ[sp[1]],
              sp[2] => δ*CF_2cycle_δ[sp[2]] + (1-δ)*CF_2cycle_1minusδ[sp[2]],
              sp[3] => δ*CF_2cycle_δ[sp[3]] + (1-δ)*CF_2cycle_1minusδ[sp[3]])

    returnCF_treeLike(CF, sp, reduced)
end
CF_2₁_3₁_chain(γ, h₁, h₂, ℓ₁, ℓ₂, δ, y) = CF_2₁_3₁_chain(γ, h₁, h₂, ℓ₁, ℓ₂, δ, y, ["a", "b", "c", "d"], false);
CF_2₁_3₁_chain(γ, h₁, h₂, ℓ₁, ℓ₂, δ, y, reduced::Bool) = CF_2₁_3₁_chain(γ, h₁, h₂, ℓ₁, ℓ₂, δ, y, ["a", "b", "c", "d"], reduced);
CF_2₁_3₁_chain(γ, h₁, h₂, ℓ₁, ℓ₂, δ, y, label::Vector{String}) = CF_2₁_3₁_chain(γ, h₁, h₂, ℓ₁, ℓ₂, δ, y, label, false);

## 3₁-2₂-cycles
function CF_3₁_2₂_chain(γ, x, δ, k₁, k₂, ℓ₁, ℓ₂, labels::Vector{String}, reduced::Bool)
    CF_2cycle_γ = CF_2cycle(δ, k₁, k₂, x*ℓ₁, ℓ₂, labels)
    CF_2cycle_1minusγ = CF_2cycle(δ, k₁, k₂, ℓ₁, ℓ₂, labels)

    sp = splits(labels)
    CF = Dict(sp[1] => γ*CF_2cycle_γ[sp[1]] + (1-γ)*CF_2cycle_1minusγ[sp[1]],
              sp[2] => γ*CF_2cycle_γ[sp[2]] + (1-γ)*CF_2cycle_1minusγ[sp[2]],
              sp[3] => γ*CF_2cycle_γ[sp[3]] + (1-γ)*CF_2cycle_1minusγ[sp[3]])

    returnCF_treeLike(CF, sp, reduced)
end
CF_3₁_2₂_chain(γ, x, δ, k₁, k₂, ℓ₁, ℓ₂) = CF_3₁_2₂_chain(γ, x, δ, k₁, k₂, ℓ₁, ℓ₂, ["a", "b", "c", "d"], false);
CF_3₁_2₂_chain(γ, x, δ, k₁, k₂, ℓ₁, ℓ₂, reduced::Bool) = CF_3₁_2₂_chain(γ, x, δ, k₁, k₂, ℓ₁, ℓ₂, ["a", "b", "c", "d"], reduced);
CF_3₁_2₂_chain(γ, x, δ, k₁, k₂, ℓ₁, ℓ₂, label::Vector{String}) = CF_3₁_2₂_chain(γ, x, δ, k₁, k₂, ℓ₁, ℓ₂, label, false);

# (-4/3)*d*g*l2*l1 + (-4/3)*d*x*l2*l1 + (4/3)*d^2*g*l2*l1 + (4/3)*d^2*x*l2*l1 + (-2/3)*g*k2*l2*l1 + (-2/3)*x*k2*l2*l1 + (4/3)*d*g*k2*l2*l1 + (4/3)*d*g*x*l2*l1 + (4/3)*d*x*k2*l2*l1 + (-2/3)*d^2*g*k1*l2*l1 + (-2/3)*d^2*g*k2*l2*l1 + (-4/3)*d^2*g*x*l2*l1 + (-2/3)*d^2*x*k1*l2*l1 + (-2/3)*d^2*x*k2*l2*l1 + (2/3)*g*x*k2*l2*l1 + (-4/3)*d*g*x*k2*l2*l1 + (2/3)*d^2*g*x*k1*l2*l1 + (2/3)*d^2*g*x*k2*l2*l1

## 3₂-2₂-cycles
function CF_3₂_2₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, ℓ₁, ℓ₂, labels::Vector{String}, reduced::Bool)
    CF_3₂_δ_δ = CF_3₂cycle(γ, h₁, h₂, x, ℓ₁*ℓ₂*k₁, labels)
    CF_3₂_δ_1minusδ = CF_2cycle(γ, h₁, h₂, x, ℓ₁*ℓ₂*k₂, labels)
    CF_3₂_1minusδ_1minusδ = CF_2cycle(γ, h₁, h₂, x, ℓ₁, labels)

    sp = splits(labels)
    CF = Dict(sp[1] => γ^2*CF_3₂_δ_δ[sp[1]] + (1-γ)^2*CF_3₂_δ_1minusδ[sp[1]] + 2*γ*(1-γ)*CF_3₂_1minusδ_1minusδ[sp[1]],
              sp[2] => γ^2*CF_3₂_δ_δ[sp[2]] + (1-γ)^2*CF_3₂_δ_1minusδ[sp[2]] + 2*γ*(1-γ)*CF_3₂_1minusδ_1minusδ[sp[2]],
              sp[3] => γ^2*CF_3₂_δ_δ[sp[3]] + (1-γ)^2*CF_3₂_δ_1minusδ[sp[3]] + 2*γ*(1-γ)*CF_3₂_1minusδ_1minusδ[sp[3]])

    returnCF_treeLike(CF, sp, reduced)
end
CF_3₂_2₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, ℓ₁, ℓ₂) = CF_3₂_2₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, ℓ₁, ℓ₂, ["a", "b", "c", "d"], false);
CF_3₂_2₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, ℓ₁, ℓ₂, reduced::Bool) = CF_3₂_2₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, ℓ₁, ℓ₂, ["a", "b", "c", "d"], reduced);
CF_3₂_2₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, ℓ₁, ℓ₂, label::Vector{String}) = CF_3₂_2₂_chain(γ, h₁, h₂, x, δ, k₁, k₂, ℓ₁, ℓ₂, label, false);


## 2₂-2₂-cycles
function CF_2₂_2₂_chain(γ, h₁, h₂, δ, k₁, k₂, ℓ₁, ℓ₂, ℓ₃, labels::Vector{String}, reduced::Bool)
    CF_γ_γ = CF_2cycle(δ, k₁, k₂, ℓ₁*h₁*ℓ₂, ℓ₃, labels)
    CF_γ_1minusγ = CF_2cycle(δ, k₁, k₂, ℓ₁*ℓ₂, ℓ₃, labels)
    CF_1minusγ_1minusγ = CF_2cycle(δ, k₁, k₂, ℓ₁*h₂*ℓ₂, ℓ₃, labels)

    sp = splits(labels)
    CF = Dict(sp[1] => γ^2*CF_γ_γ[sp[1]] + 2*γ*(1-γ)*CF_γ_1minusγ[sp[1]] + (1-γ)^2*CF_1minusγ_1minusγ[sp[1]],
              sp[2] => γ^2*CF_γ_γ[sp[2]] + 2*γ*(1-γ)*CF_γ_1minusγ[sp[2]] + (1-γ)^2*CF_1minusγ_1minusγ[sp[2]],
              sp[3] => γ^2*CF_γ_γ[sp[3]] + 2*γ*(1-γ)*CF_γ_1minusγ[sp[3]] + (1-γ)^2*CF_1minusγ_1minusγ[sp[3]])

    returnCF_treeLike(CF, sp, reduced)
end
CF_2₂_2₂_chain(γ, h₁, h₂, δ, k₁, k₂, ℓ₁, ℓ₂, ℓ₃) = CF_2₂_2₂_chain(γ, h₁, h₂, δ, k₁, k₂, ℓ₁, ℓ₂, ℓ₃, ["a", "b", "c", "d"], false);
CF_2₂_2₂_chain(γ, h₁, h₂, δ, k₁, k₂, ℓ₁, ℓ₂, ℓ₃, reduced::Bool) = CF_2₂_2₂_chain(γ, h₁, h₂, δ, k₁, k₂, ℓ₁, ℓ₂, ℓ₃, ["a", "b", "c", "d"], reduced);
CF_2₂_2₂_chain(γ, h₁, h₂, δ, k₁, k₂, ℓ₁, ℓ₂, ℓ₃, label::Vector{String}) = CF_2₂_2₂_chain(γ, h₁, h₂, δ, k₁, k₂, ℓ₁, ℓ₂, ℓ₃, label, false);
