include("tnlib.jl")

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