using ITensors

function gilttest()
  i1, i2, j1, j2, k1, k2 = siteinds(4, 6)
  C1 = combiner(i1, j1); c1 = combinedind(C1)
  C2 = combiner(i2, k1); c2 = combinedind(C2)
  C3 = combiner(j2, k2); c3 = combinedind(C3)
  U = δ(i1, i2) * δ(j1, j2) * δ(k1, k2) * C1 * C2 * C3
  U /= 2
  U, c1, c2, c3
end