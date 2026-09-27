# Sept 27, 2026:  looking at CA's notes that admixtools2 is wrong
#
#      CA's file is called simplenet_expectedf2.jl

cd("/Users/eallman/Documents/work/f4-calc/Jfiles")

##
using DataFrames, CSV
using PhyloNetworks, PhyloPlots, RCall

include("interop_admixtools.jl")

const PN = PhyloNetworks
using LinearAlgebra

##
nwk = "(a:0.1,(b:0.05,(((c1:0.001,c2:0.001)i1:0.02)#H1:0.01,#H1:0.02::0.5)i2:0.03)i3:0.05)i4;"
net = readnewick(nwk);

for e in net.edge
  if e.hybrid && e.length > 0
    breakedge!(e, net; lengthratio=1)
  end
end
net_edgelist = net_to_edgelist(net) # convert to an "edge list"
CSV.write("simplenet_edgelist.csv", net_edgelist) # file to read it later
# writenewick(net)


# R"pdf"("simplenet.pdf", height=3, width=5);
# R"par"(mar=[0,0,0,0]);
# plot(net, showedgelength=true, shownodelabel=true, tipoffset=0.2,
# majorhybridedgecolor="deepskyblue", minorlinetype="solid");
# R"dev.off"();

##
f2mat = expectedf2matrix(net)
taxa = tiplabels(net)
DataFrame(f2mat, taxa)

##
f2vec = Float64[]; pop1 =String[]; pop2 = String[]
for i1 in 1:4 
  for i2 in (i1+1):4
    push!(pop1, taxa[i1])
    push!(pop2, taxa[i2])
    push!(f2vec, f2mat[i1,i2])
  end; 
end
df = DataFrame(pop1=pop1, pop2=pop2, f2=round.(f2vec, digits=12))
CSV.write("simplenet_expectedf2_julia.csv", df)

##
m = PN.descendenceweight(net);
m[:all] # 11×11, rows=nodes in preorder, columns=edges in net.edge
m[:tips] # 4×11, rows=tips, ordered as in taxa

Ω = m[:tips] * Diagonal([e.length for e in net.edge]) * transpose(m[:tips]);
DataFrame(Ω, taxa)

##
mc1 = m.V[9:9,:]
mc2 = m.V[8:8,:]

D = Diagonal([e.length for e in net.edge]);

println(mc1 * D * transpose(mc2))

intersect(findall(!iszero,vec(mc1)),findall(!iszero,vec(mc2)))

Oc1_c2 = 0.02 + 0.5^2*0.01 + 0.5^2*0.02 + 0.03 + 0.05