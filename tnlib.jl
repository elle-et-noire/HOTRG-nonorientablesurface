using ITensors, LinearAlgebra

function split(A::ITensor, Linds...; conn = true, kwargs...)
  U, S, V = svd(A, Linds; kwargs...)
  sqrtS = sqrt.(S)
  V *= replaceinds(sqrtS, commonind(U, S) => commonind(S, V)')
  U *= sqrtS
  if conn
    replaceinds!(V, commonind(U, S)' => commonind(U, S))
  end
  U, V
end

function bulk(vweight, hweight = vweight)
  iv = Index(size(vweight)[1])
  ih = Index(size(hweight)[1])

  V1, V2 = split(ITensor(vweight, iv, iv'), iv; conn = false, righttags = "v")
  H1, H2 = split(ITensor(hweight, ih, ih'), ih; conn = false, righttags = "h")
  A = δ(iv, ih, iv', ih') * V1 * V2 * H1 * H2

  A, uniqueind(V1, iv), uniqueind(H1, ih)
end

function isometry(A::ITensor, B::ITensor; kwargs...)
  iAB = commonind(A, B)
  iA = uniqueind(noprime(inds(A)), noprime(iAB), noprime(iAB))
  iB = prime(iA)
  AA = A * prime(A, 4, iA, iAB)
  BB = B * prime(B, 4, iB, iAB)
  U, _, _ = svd(AA * BB, (iA, iB); kwargs...)
  U
end

function hotrg(T::ITensor, iv::Index, ih::Index; maxdim, stepnum, eigvalnum, withsro = true)
  @assert hassameinds([iv, ih, iv', ih'], T)
  norms = zeros(stepnum)
  cftval = Dict(zip(["<C|i>", "<R|i>", "eigval"], [zeros(eigvalnum, stepnum) for _ in 1:3]))

  # init spatial reflection operator
  O = δ(ih, ih')
  o(i, j) = begin
  	@assert hassameinds([ih, ih'], O)
	replaceinds(O, ih => i, ih' => j)
  end

  prime!(T, ih')
  for i in 1:stepnum
    if withsro
      B = T'
    else
      B = prime(T, ih, ih'') * δ(iv, iv'')
    end
    U = isometry(T, B; maxdim, lefttags = "h$i")
    T *= B

    # measure cft data of crosscap and rainbow boundary state
    M = T * δ(iv, iv'')
    U1, S, _ = svd(M, (ih, ih'); maxdim = eigvalnum)
    S = storage(S)
    norms[i] = S[1]; S /= norms[i]; T /= norms[i]
    D = length(S)
    cftval["eigval"][1:D, i] = S

    # if the copy is reflected, correspondence between crosscap/rainbow and contracting δ_ij/O_ij flips
    if withsro
      cftval["<C|i>"][1:D, i] = storage(U1 * δ(ih, ih'))
      cftval["<R|i>"][1:D, i] = storage(U1 * o(ih, ih'))
    else
      cftval["<R|i>"][1:D, i] = storage(U1 * δ(ih, ih'))
      cftval["<C|i>"][1:D, i] = storage(U1 * o(ih, ih'))
    end

    T *= U; T *= prime(U', ih', ih'')
    O = U * o(ih, ih''') * o(ih', ih'') * prime(U', ih', ih'')
    ih = commonind(T, U)
	@assert hassameinds([ih, ih'], O)

    U = isometry(T, T'; maxdim, lefttags = "v$i")
    T *= T'; T *= U; T *= prime(U', iv', iv'')
    iv = commonind(T, U)
  end

  norms, cftval
end


"""
Insert isometry between A and B.
"""
function gilt(A::ITensor, B::ITensor, C::ITensor, D::ITensor; ϵ)
  T = [A, B, C, D]
  iin = []
  for i in 1:4
    @assert length(commoninds(T[i], T[mod1(i + 1, 4)])) == 1
    push!(iin, commonind(T[i], T[mod1(i + 1, 4)]))
  end
  for i in 1:2
    @assert length(commoninds(T[i], T[mod1(i + 2, 4)])) == 0
  end

  iAB = commonind(A, B)
  TT = Vector{ITensor}(undef, 4)
  for i in 1:4
    TT[i] = T[i] * prime(T[i], 2, iin[i], iin[mod1(i - 1, 4)])
  end
  EE = ((TT[3] * TT[4]) * TT[1]) * prime(TT[2], iAB, iAB'')

  U, S, _ = svd(EE, (iAB, iAB'))
  t = storage(U * δ(iAB, iAB'))
  S ./= maximum(S)
  for i in eachindex(t)
    t[i] *= storage(S)[i] |> x -> x^2 / (x^2 + ϵ^2)
  end
  R = ITensor(t, commonind(U, S)) * U
  @assert hassameinds([iAB, iAB'], R)
  U, V = split(R, iAB; righttags = "gilt")
  U, setprime!(V, plev(iAB), iAB')
  # A *= U; B *= noprime(V)
  # U, V
end

function trg(A::ITensor, iv::Index, ih::Index; maxdim, stepnum, eigvalnum)
  @assert hassameinds([iv, iv', ih, ih'], A)
  norms = zeros(stepnum)
  eigval = zeros(eigvalnum, stepnum)

  for i in 1:stepnum
    F1, F3 = split(A, (iv, ih); maxdim, conn = false, righttags = "v$i")
    F2, F4 = split(A, (iv, ih'); maxdim, conn = false, righttags = "h$i")
    iv = uniqueind(F1, A); ih = uniqueind(F2, A)
    A = F1' * F2 * noprime(F3) * F4
    # A = F1 * F2 * F3 * F4
    @assert hassameinds([iv, iv', ih, ih'], A)

    M = A * δ(iv, iv')
    norms[i] = scalar(M * δ(ih, ih'))
    A /= norms[i]
    _, S, _ = svd(M, ih; maxdim = eigvalnum)
    S = storage(S)
    eigval[eachindex(S), i] = S / norms[i]
  end

  norms, eigval
end

function gilttnr(A::ITensor, iv::Index, ih::Index; maxdim, stepnum, eigvalnum, ϵ)
  @assert hassameinds([iv, iv', ih, ih'], A)
  norms = zeros(stepnum)
  eigval = zeros(eigvalnum, stepnum)

  F = Vector{ITensor}(undef, 4)
  for i in 1:stepnum
    F[1], F[3] = split(A, (iv, ih); maxdim, conn = false, righttags = "v$i")
    F[2], F[4] = split(A, (iv, ih'); maxdim, conn = false, righttags = "h$i")
    iv = uniqueind(F[1], A); ih = uniqueind(F[2], A)
    F[1] = F[1]'; F[3] = noprime(F[3])
    for i in 1:4
      U, V = gilt(F[i], F[mod1(i + 1, 4)], F[mod1(i + 2, 4)], F[mod1(i + 3, 4)]; ϵ)
      F[i] *= U; F[mod1(i + 1, 4)] *= V
    end
    # A = F[1]' * F[2] * noprime(F[3]) * F[4]
    A = F[1] * F[2] * F[3] * F[4]
    @assert hassameinds([iv, iv', ih, ih'], A)

    M = A * δ(iv, iv')
    norms[i] = scalar(M * δ(ih, ih'))
    A /= norms[i]
    _, S, _ = svd(M, ih; maxdim = eigvalnum)
    S = storage(S)
    eigval[eachindex(S), i] = S / norms[i]
  end

  norms, eigval
end

function logpartfunc(norms; sitenum_per_step)
  sum = 0.; sitenum = 1; lnz = zero(norms)
  for i in eachindex(norms)
    sitenum *= sitenum_per_step
    lnz[i] = sum += log(norms[i]) / sitenum
  end
  lnz
end