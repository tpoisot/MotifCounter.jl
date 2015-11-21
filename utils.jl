function mhash(A::Array{Int64, 2})
  return reduce(*, map(string, A))
end

function nonempty(A::Array{Int64, 2})
  return prod(vec(sum(A, 1)) .+ vec(sum(A, 2))) > 0
end

function hashmat(x)
  V = round(Int64, sqrt(length(x)))
  A = reshape(x, (V, V))
  if nonempty(A)
    return mhash(A)
  end
end

function unhashmat(x)
  c = Int64[]
  for i in 1:length(x)
    push!(c, parse(Int64, x[i]))
  end
  V = round(Int64, sqrt(length(x)))
  return reshape(vec(c), (V, V))
end

function from_binarytree(h, s)
   tree = zeros(Int64, 2^(h+1)-2)
   for i in 1:length(tree)
      if iseven(i)
         tree[i] = 1
      end
   end
   terminal_nodes = (length(tree)-2^h+1):length(tree)
   paths = []
   for n in terminal_nodes
      i = h
      position = n
      sequence = zeros(Int64, h)
      while(i > 0)
         sequence[i] = tree[position]
         position = round(Int,floor((position - 1)/2))
         i = i-1
      end
      if sum(sequence) == s
         push!(paths, sequence)
      end
   end
   return paths
end
