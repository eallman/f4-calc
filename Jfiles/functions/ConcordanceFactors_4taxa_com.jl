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
function CF_2cycle_com(γ, h₁, h₂, ℓ₁, ℓ₂, labels::Vector{String}, reduced::Bool)
    sp = splits(labels)
    CF = Dict(sp[1] => (1-γ)*(1-2//3*ℓ₁*ℓ₂*h₂) + γ*(1-2//3*ℓ₁*ℓ₂*h₁),
              sp[2] => (1-γ)*(1//3*ℓ₁*ℓ₂*h₂) + γ*(1//3*ℓ₁*ℓ₂*h₁),
              sp[3] => (1-γ)*(1//3*ℓ₁*ℓ₂*h₂) + γ*(1//3*ℓ₁*ℓ₂*h₁))

    returnCF_treeLike(CF, sp, reduced)
end
CF_2cycle_com(γ, h₁, h₂, ℓ₁, ℓ₂) = CF_2cycle_com(γ, h₁, h₂, ℓ₁, ℓ₂, ["a", "b", "c", "d"], false);
CF_2cycle_com(γ, h₁, h₂, ℓ₁, ℓ₂, reduced::Bool) = CF_2cycle_com(γ, h₁, h₂, ℓ₁, ℓ₂, ["a", "b", "c", "d"], reduced);
CF_2cycle_com(γ, h₁, h₂, ℓ₁, ℓ₂, label::Vector{String}) = CF_2cycle_com(γ, h₁, h₂, ℓ₁, ℓ₂, label, false);

## 3₁-cycle
function CF_3₁cycle_com(γ, x, ℓ, labels::Vector{String}, reduced::Bool)
    CF_tree_γ = CF_tree(x*ℓ, labels)
    CF_tree_1minusγ = CF_tree(ℓ, labels)

    sp = splits(labels)
    CF = Dict(sp[1] => γ*CF_tree_γ[sp[1]] + (1-γ)*CF_tree_1minusγ[sp[1]],
              sp[2] => γ*CF_tree_γ[sp[2]] + (1-γ)*CF_tree_1minusγ[sp[2]],
              sp[3] => γ*CF_tree_γ[sp[3]] + (1-γ)*CF_tree_1minusγ[sp[3]])

    returnCF_treeLike(CF, sp, reduced)
end
CF_3₁cycle_com(γ, x, ℓ) = CF_3₁cycle_com(γ, x, ℓ, ["a", "b", "c", "d"], false);
CF_3₁cycle_com(γ, x, ℓ, reduced::Bool) = CF_3₁cycle_com(γ, x, ℓ, ["a", "b", "c", "d"], reduced);
CF_3₁cycle_com(γ, x, ℓ, label::Vector{String}) = CF_3₁cycle_com(γ, x, ℓ, label, false);

## 3₂-cycle
function CF_3₂cycle_com(γ, h₁, h₂, x, ℓ, labels::Vector{String}, reduced::Bool)
    CF_2_γ = CF_tree(h₁, labels)
    CF_2_1minusγ = CF_tree(h₂, labels)
    CF_γ_1minusγ = CF_tree(x, [labels[1],labels[3],labels[2],labels[4]])
    CF_1minusγ_γ = CF_tree(x, [labels[1],labels[4],labels[2],labels[3]])

    sp = splits(labels)

    CF = Dict(sp[1] => (1-ℓ) + ℓ*(γ*CF_2_γ[sp[1]] + (1-γ)*CF_2_1minusγ[sp[1]]),
              sp[2] => ℓ*(γ*CF_2_γ[sp[2]] + (1-γ)*CF_2_1minusγ[sp[2]]),
              sp[3] => ℓ*(γ*CF_2_γ[sp[3]] + (1-γ)*CF_2_1minusγ[sp[3]]))

    returnCF_treeLike(CF, sp, reduced)
end
CF_3₂cycle_com(γ, h₁, h₂, x, ℓ) = CF_3₂cycle_com(γ, h₁, h₂, x, ℓ, ["a", "b", "c", "d"], false);
CF_3₂cycle_com(γ, h₁, h₂, x, ℓ, reduced::Bool) = CF_3₂cycle_com(γ, h₁, h₂, x, ℓ, ["a", "b", "c", "d"], reduced);
CF_3₂cycle_com(γ, h₁, h₂, x, ℓ, label::Vector{String}) = CF_3₂cycle_com(γ, h₁, h₂, x, ℓ, label, false);

## 4₁-cycle
function CF_4₁cycle_com(γ, x₁, x₂, labels::Vector{String}, reduced::Bool)
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
CF_4₁cycle_com(γ, x₁, x₂) = CF_4₁cycle_com(γ, x₁, x₂, ["a", "b", "c", "d"], false);
CF_4₁cycle_com(γ, x₁, x₂, reduced::Bool) = CF_4₁cycle_com(γ, x₁, x₂, ["a", "b", "c", "d"], reduced);
CF_4₁cycle_com(γ, x₁, x₂, label::Vector{String}) = CF_4₁cycle_com(γ, x₁, x₂, label, false);

