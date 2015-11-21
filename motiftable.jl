using JSON

@everywhere include("utils.jl")

function generatemotifs(V, E)
   
   #=
   The first step is to generate a binary tree of depth V^2
   The from_binarytree function does that, and only returns paths with
      exactly E interactions. This is the entire set of unique permutations
      of E in V^2 elements.
   =#
   t_hashes = pmap(hashmat, from_binarytree(V^2, E))
   
   #=
   TODO maybe this next line will have to be removed
   =#
   t_hashes = filter((x) -> x != nothing, t_hashes)

   #=
   Finally, we return an Dict where every key is a hash (the hash is
      the unfolded matrix), and the value is the matrix itself, at the right
      dimension.
   =#
   hashes = Dict{Any,Any}([x => unhashmat(x) for x in t_hashes])
   return hashes
end

function identifymotifs(motif_set)
   hash_set = Dict()

   #=
   How many motifs do we know at the beginning? 0. Well, 1, because the first
      motif we'll come across will necessarilly be a new one.
   =#
   known_motifs = 1

   #=
   We will then move through all hashes, i.e. all unique matrix conformations.
   =#
   for k in collect(keys(motif_set))

      #=
      The first step is to collect the matrix to be evaluated in A
      =#
      A = motif_set[k]

      #=
      Then, we simply look at all the possible permutations of rows orders.
         The columns HAVE to be re-ordered in the same way since A is the
         adjacency matrix that represents the motif.
      =#
      orders = collect(permutations(vec(1:size(A, 1))))

      #=
      This empty array will store the representations of every possible order.
      =#
      list = ASCIIString[]

      #=
      For every order, we re-order the matrix, hash it, and push it to the list
         of all possible unique conformations. This is an OK overhead because we 
         are often limited to motifs with V <= 5, so that even V! permutations
         are OK.
      =#
      for o1 in orders
         push!(list, mhash(A[o1, o1]))
      end

      #=
      We know want to know whether we already know this motif. By default, we
         assume that no, we don't ...
      =#
      is_known = false

      #=
         ... and it has therefore no known id. But this will change.
      =#
      has_id = 0

      #=
      Then for every possible conformation in the list (the list
         is actually the unique hashes) ...
      =#
      for l in list

         #=
         It is fairly simple. If the motif is know, i.e. if it is in the list
            of all hashes,
         =#
         if l in keys(hash_set)

            #=
            We mark it as known, and give it an ID.
            =#
            is_known = true
            has_id = hash_set[l]
            # TODO we could break of the loop at this point
         end
      end

      #=
      Then, IF we know this motif...
      =#
      if is_known

         #=
         ... we put it's ID in the hash_set ...
         =#
         for l in list
            hash_set[l] = has_id
         end
      else

         #=
         ... but if not, we create a new motif
         =#
         for l in list
            hash_set[l] = known_motifs
         end

         #=
         We implement in the end, because, remember, we started at 1 known motif.
         =#
         known_motifs += 1
      end
   end
   #=
   Done! Good job.
   =#
   return hash_set
end

hashes = Dict()

for V in 2:4
  hashes[V] = Dict()
  for E in (V-1):V^2
    println("V: ", V, "\tE: ", E)
    motifs = generatemotifs(V, E)
    hashes[V][E] = identifymotifs(motifs)
  end
end

motifs = open("motifs.json", "w")
JSON.print(motifs, hashes)
close(motifs)
