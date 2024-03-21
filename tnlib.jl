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

function gilt_plaq(A1::ITensor, A2::ITensor, iv::Index, ih::Index; ϵ)
  legs = [ih, iv, ih', iv']
  done_legs = Dict(zip(legs, [false for _ in 1:4]))
  for leg in Iterators.cycle(legs)
    A1, A2, done_legs[leg] = apply_gilt(A1, A2, iv, ih, leg; ϵ)
    all(values(done_legs)) && break
  end

  A1, A2
end

function apply_gilt(A1::ITensor, A2::ITensor, iv, ih, leg::Index; ϵ)
  U, S = get_envspec(A1, A2, iv, ih, leg)
  A2 *= δ(leg, leg'')
  S /= sum(storage(S))
  Rp = optimize_Rp(U, S, leg, leg''; ϵ)
  # Rp1, Rp2 = split(Rp, leg; cutoff = ϵ * 1e-3)
  Rp1, Rp2 = split(Rp, leg)
  i12 = commonind(Rp1, Rp2)
  A1 *= Rp1; A2 *= Rp2
  replaceinds!(A1, i12 => leg)
  replaceinds!(A2, i12 => leg)
  A1, A2, true
end

function get_envspec(A1, A2, iv, ih, leg)
  NW = A1 * prime(A1, 3, iv', ih')
  NE = A2 * prime(A2, 3, ih', iv)
  SE = A1 * prime(A1, 3, iv, ih)
  SW = A2 * prime(A2, 3, iv', ih)
  replaceinds!(NE, leg => leg'', leg''' => leg''''')
  replaceinds!(SW, leg => leg'', leg''' => leg''''')

  EE = NW * NE * SE * SW
  D, U = eigen(EE, (leg''', leg'''''), (leg, leg''); ishermitian = true)
  U, sqrt.(abs.(D))
end

function optimize_Rp(U, S, leg1, leg2; ϵ)
  println("optimizing Rp...")
  # leg1, leg2 = uniqueinds(U, S)
  # leg1 = leg; leg2 = leg''
  t = U * δ(leg1, leg2)

  C_err_constterm = norm(t * S)
  function C_err(tp)
    diff = t - tp
    diff *= S
    norm(diff) / C_err_constterm
  end

  weight = zero(collect(storage(S)))
  for i in eachindex(weight)
    ratio = storage(S)[i] / ϵ
    weight[i] = ratio^2 / (1 + ratio^2)
  end
  iUS = commonind(t, S)
  w = diagITensor(weight, iUS, iUS')
  tp = replaceinds(t * w, iUS' => iUS)
  # Rp = U * tp
  Rp = U * t

  # u, s, v = svd(Rp, leg1; cutoff = ϵ * 1e-3)
  # done_recursing = maximum(abs.(storage(s) .- 1)) < 1e-2
  # if !done_recursing
  #   ssqrt = sqrt.(s)
  #   us = u * ssqrt
  #   vs = v * ssqrt
  #   UuvsS = S * U * us * vs
  #   Uinner, Sinner, _ = svd(UuvsS, (commonind(s, v), commonind(u, s)))
  #   Sinner /= sum(storage(Sinner))
  #   Rinner = optimize_Rp(Uinner, Sinner, commonind(s, v), commonind(u, s); ϵ)
  #   Rp = Rinner * us * vs
  # end

  Rp
end

function trgstep(A1, A2, iv, ih, log_fact; maxdim)
  F = Vector{ITensor}(undef, 4)
  F[1], F[3] = split(A1, (iv, ih); maxdim, conn = false, righttags = "i13")
  F[2], F[4] = split(A2, (iv', ih); maxdim, conn = false, righttags = "i24")
  A = F[1] * F[2] * F[3] * F[4]
  iv = commonind(F[1], A); ih = commonind(F[2], A)
  A, iv, ih, log_fact * 2
end

function gilttnr_step(A, iv, ih, log_fact; maxdim, ϵ)
  @assert hassameinds([iv, ih, iv', ih'], A)
  m = norm(A)
  if !iszero(m)
    A /= m
    log_fact += log(m)
  end

  if ϵ > 0
    A1, A2 = gilt_plaq(A, swapprime(A, 0, 1), iv, ih; ϵ)
  else
    A1, A2 = A, swapprime(A, 0, 1)
  end

  display(A1)
  display(A2)
  A, iv, ih, log_fact = trgstep(A1, A2, iv, ih ,log_fact; maxdim)
  A, iv, ih, log_fact = trgstep(A, swapprime(A, 0, 1), iv, ih, log_fact; maxdim)
  display(A)

  replaceinds!(A, iv => ih, ih => iv', iv' => ih', ih' => iv)

  A, iv, ih, log_fact
end

function gilttnr(A::ITensor, iv::Index, ih::Index; maxdim, stepnum, eigvalnum, ϵ)
  @assert hassameinds([iv, iv', ih, ih'], A)
  lnz = zeros(stepnum)
  eigval = zeros(eigvalnum, stepnum)
  log_fact = 0

  for i in 1:stepnum
    println("step $i:")
    A, iv, ih, log_fact = gilttnr_step(A, iv, ih, log_fact; maxdim, ϵ)

    M = A * δ(iv, iv')
    lnz[i] = log(scalar(M * δ(ih, ih'))) + log_fact
    _, S, _ = svd(M, ih; maxdim = eigvalnum)
    S = storage(S)
    eigval[eachindex(S), i] = S / maximum(S)
  end

  lnz, eigval
end

function logpartfunc(norms; sitenum_per_step)
  sum = 0.; sitenum = 1; lnz = zero(norms)
  for i in eachindex(norms)
    sitenum *= sitenum_per_step
    lnz[i] = sum += log(norms[i]) / sitenum
  end
  lnz
end