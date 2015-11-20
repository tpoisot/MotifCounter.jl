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
  c = vec([parse(Int64, element) for element in x])
  V = round(Int64, sqrt(length(x)))
  return reshape(c, (V, V))
end
