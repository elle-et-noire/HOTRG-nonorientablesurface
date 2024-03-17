include("tnlib.jl")
using QuadGK

abstract type Model end

struct Ising <: Model end

struct Potts <: Model
  q::Int
end

function weight(_::Ising, β)
  [exp(β) exp(-β); exp(-β) exp(β)]
end

function weight(p::Potts, β)
  Diagonal([exp(β) - 1 for _ in 1:p.q]) + ones(p.q, p.q)
end

function Tc(_::Ising)
  2 / log(1 + √2)
end

function Tc(p::Potts)
  1 / log(1 + √p.q)
end

function centralcharge(_::Ising)
  0.5
end

function centralcharge(p::Potts)
  if p.q == 3
    0.8
  elseif p.q == 4
    1.0
  end
end

function quantumdimension(_::Ising)
  1 + 1 / √2
end

function quantumdimension(p::Potts)
  if p.q == 3
    sqrt(3 + 6 / √5)
  elseif p.q == 4
    (3 + 2√2) / 2
  end
end

function freeenergy(_::Ising, β)
  (Main.quadgk(k1 ->
  Main.quadgk(k2 ->
    log(cosh(2β)^2 - sinh(2β) * (cos(k1) + cos(k2))) / 8pi^2,
  0, 2pi)[1],
  0, 2pi)[1] + log(2)) / -β
end