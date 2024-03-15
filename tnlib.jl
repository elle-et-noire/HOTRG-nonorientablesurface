using ITensors, LinearAlgebra

function isometry(A::ITensor, B::ITensor; maxdim)
  iAB = commonind(A, B)
  iA = noprime(uniqueind(A, iAB, noprime(iAB)))
  iB = noprime(uniqueind(B, iAB, iAB'))
  AA = A * prime(A, 4, iA, iAB)
  BB = B * prime(B, 4, iB, iAB)
  U, _, _ = svd(AA * BB, (iA, iB); maxdim)
  U
end

function split(A::ITensor)
  @assert length(inds(A)) == 2
  U, S, V = svd(A, inds(A)[1])
  sqrtS = sqrt(S)
  V *= replaceinds(sqrtS, commonind(U, S) => commonind(S, V)')
  U *= sqrtS
  U, V
end

function bulk(vweight, hweight = vweight)
  iv = Index(size(vweight)[1])
  ih = Index(size(hweight)[1])

  V1, V2 = split(ITensor(vweight, iv, iv'))
  H1, H2 = split(ITensor(hweight, ih, ih'))
  A = δ(iv, ih, iv', ih') * V1 * V2 * H1 * H2

  A, iv, ih
end

function hotrg(A::ITensor, iv::Index, ih::Index; maxdim, stepnum, eigvalnum, withsro = true)
  @assert hassameinds([iv, ih, iv', ih'], A)
  norms = zeros(stepnum)
  cftval = Dict(zip(["<C|i>", "<R|i>", "eigval"], [zeros(eigvalnum, stepnum) for _ in 1:3]))

  # init spatial reflection operator
  O = δ(ih, ih')
  refl(i, j) = replaceinds(O, ih => i, ih' => j)

  prime!(A, ih')
  for i in 1:stepnum
    if withsro
      B = A'
    else
      B = prime(A, ih, ih'') * δ(iv, iv'')
    end
    AB = A * B

    # measure cft data of crosscap and rainbow boundary state
    M = AB * δ(iv, iv'')
    U1, S, _ = svd(M, (ih, ih'); maxdim = eigvalnum)
    S = storage(S)
    norms[i] = S[1]; S /= norms[i]; AB /= norms[i]
    D = length(S)
    cftval["eigval"][1:D, i] = S

    # if the copy is reflected, correspondence between crosscap/rainbow and contracting δ_ij/O_ij flips
    if withsro
      cftval["<C|i>"][1:D, i] = storage(U1 * δ(ih, ih'))
      cftval["<R|i>"][1:D, i] = storage(U1 * refl(ih, ih'))
    else
      cftval["<R|i>"][1:D, i] = storage(U1 * δ(ih, ih'))
      cftval["<C|i>"][1:D, i] = storage(U1 * refl(ih, ih'))
    end

    U = isometry(A, B; maxdim)
    AB *= U; AB *= U'
    O = U * refl(ih, ih''') * refl(ih', ih'') * U'
    ih = commonind(AB, U)

    AB = A * A'
    U = isometry(A, A'; maxdim)
    AB *= U; AB *= U'
    iv = commonind(AB, U)
  end

  norms, cftval
end