using Oscar

function splits(labels)
    
    #sp1 = join(sort([join(sort(labels[[1,2]])), join(sort(labels[[3,4]]))]))
    #sp2 = join(sort([join(sort(labels[[1,3]])), join(sort(labels[[2,4]]))]))
    #sp3 = join(sort([join(sort(labels[[1,4]])), join(sort(labels[[2,3]]))]))

    sp1 = join([join(sort(labels[[1,2]])), join(sort(labels[[3,4]]))])
    sp2 = join([join(sort(labels[[1,3]])), join(sort(labels[[2,4]]))])
    sp3 = join([join(sort(labels[[1,4]])), join(sort(labels[[2,3]]))])

    return [sp1, sp2, sp3]
end

function CF_coordinates_labels(keys_CF, old_news)
    # Change labels (if any)
    rep_ks = keys_CF
    for k in collect(keys(old_news))
        rep_ks = replace.(rep_ks, k => old_news[k])
    end

    # Specify coordinates names
    Cks = "C" .* rep_ks
end

function translateToM2(CF, nameIdeal="", old_news::Pair...)
    # List of coordinates (Cabcd, ...) and variables (γ, x₁, ...) in CFs
    ks = collect(keys(CF))
    n = length(ks)
    vars = variables(System(collect(values(CF))))

    # Specify coordinates names
    Cks = CF_coordinates_labels(ks, Dict(old_news))

    println("R = QQ[", join(vars, ", "), ", ", join(Cks, ", "), ", MonomialOrder => Lex]", "\n")

    for i in 1:n
        println( "f", i , " = - ", Cks[i], " + ", expand(CF[ks[i]]))
    end

    println("\n", "I", nameIdeal, " = ideal(", join("f".*string.(collect(1:n)), ", "), ")")
    println("J", nameIdeal, " = eliminate({", join(vars, ", "), "}, I", nameIdeal, ")", "\n")
end

function translateToMagma(CF, old_news::Pair...)
    # List of coordinates (Cabcd, ...) and variables (γ, x₁, ...) in CFs
    ks = collect(keys(CF))
    n = length(ks)
    vars = variables(System(collect(values(CF))))

    # Specify coordinates names
    Cks = CF_coordinates_labels(ks, Dict(old_news))

    println("P<", join(vars, ", "), ", ", join(Cks, ", "), "> := PolynomialRing(RationalField(), ", n+length(vars), ");", "\n")

    for i in 1:n
        println( "f", i , " := - ", Cks[i], " + ", CF[ks[i]] , ";")
    end

    println("I := ideal<P | ", join("f".*string.(collect(1:n)), ", "), ">;")
    println("J := EliminationIdeal(I, {", join(Cks, ", ") ,"});")
end

function returnCF_treeLike(CF, sp, reduced)
    if reduced
        key = sort([sp[2], sp[3]])[1]
        Dict(key => CF[key])
    else
        CF
    end
end


function parametrization(CF)
    R = parent(collect(values(CF))[1])

    kCF = collect(keys(CF))
    S, cf = polynomial_ring(QQ,:cf=>kCF)
    phi = hom(S, R, [CF[k] for k in kCF])
    return (phi)
end

function elimination_ideal(phi::Oscar.MPolyAnyMap)
    S = domain(phi)
    R = codomain(phi)

    @req coefficient_ring(R) == coefficient_ring(S) "coefficient fields of parameter and model ring must be the same"

    elim_ring, elim_gens = polynomial_ring(coefficient_ring(R), vcat(symbols(R), symbols(S)))
    inject_R = hom(R, elim_ring, elim_gens[1:ngens(R)])
    inject_S = hom(S, elim_ring, elim_gens[ngens(R) + 1:ngens(elim_ring)])

    I = ideal([inject_S(s) - inject_R(phi(s)) for s in gens(S)])
    
    ## project_S = hom(elim_ring, S, z -> z, [repeat([1], ngens(R)); gens(S)])
    return(elim_ring, [elim_gens[1:ngens(R)], elim_gens[ngens(R) + 1:ngens(elim_ring)]], I)
end


function vanishing_ideal(CF; algorithm::Symbol = :eliminate)
  phi = parametrization(CF)
  
  if algorithm == :kernel
    return kernel(phi)
  else
   
    elim_ring, (params, cf_coord), I  = elimination_ideal(phi)
    S = domain(phi)
    project_S = hom(elim_ring, S, z -> z, [repeat([1], length(params)); gens(S)])

    if algorithm == :eliminate
      # project_S(eliminate(I, params))
      invariants = generating_system(eliminate(I, params))
      
      return ideal(map(project_S, gens(invariants)))
    elseif algorithm == :f4
      invariants = groebner_basis_f4(I, eliminate=length(params))
      emb = hom(base_ring(invariants), S, gens(S))
      
      return ideal(map(emb, invariants))
    else
      error("Didn't recognize  algorithm $algorithm")
    end
  end
end

# function vanishing_ideal(CF_dict)
#     elim_ring, elim_gens, I, project_S  = elimination_ideal(CF_dict)
#     return project_S(eliminate(I, elim_gens[1]))
# end