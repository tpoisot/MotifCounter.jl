using JSON

@everywhere include("utils.jl")

function generatemotifs(V, E)
  t_hashes = pmap(hashmat, from_binarytree(V^2, E))
  t_hashes = filter((x) -> x != nothing, t_hashes)
  hashes = Dict{Any,Any}([x => unhashmat(x) for x in t_hashes])
  return hashes
end

function identifymotifs(motif_set)
  hash_set = Dict()
  known_motifs = 1
  for k in collect(keys(motif_set))
    A = motif_set[k]
    orders = collect(permutations(vec(1:size(A, 1))))
    # We generate a list of hashes
    list = ASCIIString[]
    for o1 in orders
        push!(list, mhash(A[o1, o1]))
    end
    # Is ANY of the hash in the hash_set already?
    is_known = false
    has_id = 0
    for l in list
      if l in keys(hash_set)
        is_known = true
        has_id = hash_set[l]
      end
    end
    if is_known
      for l in list
        hash_set[l] = has_id
      end
    else
      for l in list
        hash_set[l] = known_motifs
      end
      known_motifs += 1
    end
  end
  return hash_set
end

hashes = Dict()

for V in 2:4
  hashes[V] = Dict()
  for E in (V-1):V^2
    println("V: ", V, "\tE: ", E)
    hashes[V][E] = identifymotifs(generatemotifs(V, E));
  end
end

motifs = open("motifs.json", "w")
JSON.print(motifs, hashes)
close(motifs)
