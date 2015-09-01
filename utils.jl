function mhash(A::Array{Int64, 2})
  return reduce(*, map(string, A))
end

function nonempty(A::Array{Int64, 2})
  return prod(vec(sum(A, 1)) .+ vec(sum(A, 2))) > 0
end

function hashmat(x)
  V = int(sqrt(length(x)))
  A = reshape(x, (V, V))
  if nonempty(A)
    return mhash(A)
  end
end

function unhashmat(x)
  c = Int64[]
  for i in 1:length(x)
    push!(c, int(string(x[i])))
  end
  V = int(sqrt(length(x)))
  return reshape(vec(c), (V, V))
end
