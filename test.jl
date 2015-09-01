using JSON
using DataFrames
using Gadfly

include("utils.jl")

mfile = open("motifs.json")
hashes = JSON.parse(mfile)
close(mfile)

function getIWDB(url::ASCIIString)
  matrix = download(url)
  am = open(matrix)
  A = int(readdlm(am))
  close(am)
  return A
end

function cleanmatrix(A::Array{Int64, 2})
  cumulative_degree = vec(sum(A, 1)) .+ sum(A, 2)
  have_degree = filter((i) -> cumulative_degree[i] > 0, 1:length(cumulative_degree))
  A = A[have_degree, have_degree]
end

function countmotifs(A::Array{Int64, 2}, msize) # TODO motif dict
  nodes_id = combinations(1:size(A, 1), msize)
  motifs = ASCIIString[]
  for nbunch in nodes_id
    sub_A = A[collect(nbunch),collect(nbunch)]
    if nonempty(sub_A)
      V = 3
      E = sum(sub_A)
      m_sub_A = mhash(sub_A)
      m_id = string(V)*"_"*string(E)*"_"*string(hashes[string(V)][string(E)][m_sub_A])
      push!(motifs, m_id)
    end
  end
  mcounts = {x => sum(motifs .== x) for x in unique(motifs)}
  return mcounts
end

function buildIWDBname(name::ASCIIString)
  return "https://www.nceas.ucsb.edu/interactionweb/data/predator_prey/text_matrices/" * name * "txt.txt"
end

function wrapIWDB(x::ASCIIString)
  ok_mat = x |> buildIWDBname |> getIWDB |> cleanmatrix
  return countmotifs(ok_mat, 3)
end

data_names = [
  "AkatoreA",
  "AkatoreB",
  "Catlins",
  "Berwick",
  "Narrowdale",
  "Blackrock",
  "Broad",
  "DempstersSu",
  "DempstersAu",
  "DempstersSp",
  "German",
  "Healy",
  "LilKyeburn",
  "Stony"
]

IWDBmotifs = {x => wrapIWDB(x) for x in data_names}

nrows = sum([length(v) for (k,v) in IWDBmotifs])
iwdb = DataFrame([Symbol, Symbol, Int64], [:fw, :motif, :count], nrows)

cr = 1
for (fw, mc) in IWDBmotifs
  for (mn, cn) in mc
    iwdb[:fw][cr] = symbol(fw)
    iwdb[:motif][cr] = symbol(mn)
    iwdb[:count][cr] = cn
    cr += 1
  end
end

# Functions to deal with proportions
function prop(x::DataArray{Symbol, 1})
  return x
end

function prop(x::DataArray{Int64, 1})
  return x ./ sum(x)
end

# Agregate the dataframe
a_iwdb = aggregate(iwdb, :fw, prop)

draw(
  PDF("motif_iwdb.pdf", 15cm, 10cm),
  plot(a_iwdb, x="fw", color="motif_prop", y="count_prop", Geom.bar)
)
