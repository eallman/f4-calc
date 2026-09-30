# Sept 29, 2026:  f2 computations now using Julia
#
# Aug 26, 2026:   original Singular file

using Oscar

# Timer setup
t_start = time()

# Tree 
#    a         c
#     \__z___/ 
#     /       \ 
#    b         d

# ring declarations
R2, (f2_ac, f2_ad, f2_bc, f2_bd) = 
    polynomial_ring(QQ,[:f2_ac, :f2_ad, :f2_bc, :f2_bd],
    internal_ordering=:degrevlex
    )

R, (f2_ab, f2_ac, f2_ad, f2_bc, f2_bd, f2_cd, a, b, z, c, d) = 
    polynomial_ring(QQ,[:f2_ab, :f2_ac, :f2_ad, :f2_bc, :f2_bd, :f2_cd, :a, :b, :z, :c, :d],
    internal_ordering=:degrevlex
    )

println("Ring R: \n")
println("Variables:    ", symbols(R))
println("Var Count:    ", nvars(R))
println("Base Field:   ", base_ring(R))
println("Term Order:   ", internal_ordering(R))

##
# f2 parameterizations 
p2_ab = a + b
p2_ac = a + z + c
p2_ad = a + z + d
p2_bc = b + z + c
p2_bd = b + z + d
p2_cd = c + d

# f4 parameterizations
p4_abcd = 1//2*(p2_ad + p2_bc - p2_ac - p2_bd)
p4_acdb = 1//2*(p2_ab + p2_cd - p2_ad - p2_bc)
p4_adbc = 1//2*(p2_ac + p2_bd - p2_ab - p2_cd)



#=
S, (f4_abcd, f4_acdb, f4_adbc, a,b,z,c,d) =
    polynomial_ring(QQ,[:f4_abcd, :f4_acdb, :f4_adbc, :a,:b,:z,:c,:d)

phi = hom(R,S,[0,0,0,0,0,0,a,b,z,c,d])

s4_abcd = phi(p4_abcd)
s4_acdb = phi(p4_acdb)
s4_adbc = phi(p4_adbc)

Jbig = ideal(f4_abcd-s4_abcd,f4_acdb-s4_acdb,f4_adbc-s4_adbc)
J = eliminate(Jbig,[a,b,c,d,z])
dim(J)-5
=#


##
# f2 ideal
T = ideal(R, [
    f2_ab - p2_ab, 
    f2_ac - p2_ac, 
    f2_ad - p2_ad, 
    f2_bc - p2_bc, 
    f2_bd - p2_bd, 
    f2_cd - p2_cd
])

# eliminate parameters (BLs)
I = eliminate(T, [a, b, z, c, d])

println("\nIdeal Generators on 4-taxon tree ab|cd:")
println(I)

# Print execution time
println("\nElapsed time: ", time() - t_start, " seconds")

# J = ideal(R, [p2_ab, p2_ad, p2_ad, p2_bc, p2_bd, p2_cd])
# gJ = groebner_basis(J)
# dim(ideal(gJ))

##
# Now do the 4₁-cycle:

# draw picture later
R2, (f2_ab, f2_ac, f2_ad, f2_bc, f2_bd, f2_cd,γ, h1, h2, l1, l2, a, b, c, d) = 
    polynomial_ring(QQ,[:f2_ab, :f2_ac, :f2_ad, :f2_bc, :f2_bd, :f2_cd, 
    :γ, :h1, :h2, :l1, :l2, :a, :b, :c, :d],
    internal_ordering=:degrevlex
    )

p2_ab = a + b + l1 + l2
p2_ac = a + c + γ*(h1+l1) + (1 - γ)*h2+l2
p2_ad = a + d + l1
p2_bc = b + c + γ*(h1+l1) + (1 - γ)*h2+l2
p2_bd = b + d + l2
p2_cd = c + d + γ*(h1+l1) + (1 - γ)*h2+l2

# f2 ideal
T = ideal(R2, [
    f2_ab - p2_ab, 
    f2_ac - p2_ac, 
    f2_ad - p2_ad, 
    f2_bc - p2_bc, 
    f2_bd - p2_bd, 
    f2_cd - p2_cd
])

t_start = time()

# eliminate parameters (BLs)
I = eliminate(T, [a, b, c, d, γ, h1, h2, l1, l2])

println("\nIdeal Generators on 4₁-cycle [c hybrid of a,b, d outgroup]:")
println(I)

# Print execution time
println("\nElapsed time: ", time() - t_start, " seconds")
##
