using JuMP, Gurobi, Printf
using LinearAlgebra: dot

TL = parse(Int, ARGS[2]) # in secs
nv = 0
edges = Array{Tuple{Int, Int}}(undef, 0)

open(ARGS[1]) do file
  global nv
  nvs, nve = split(readline(file))
  nv = parse(Int, nvs)

  resize!(edges, parse(Int, nve))

  for (i, line) in enumerate(eachline(file))
    u, v = split(line)
    edges[i] = (parse(Int, u), parse(Int, v))
  end
end

MEM = 4 # in GB
m = Model(solver = GurobiSolver(NodefileStart=MEM, TimeLimit=TL))

@variable(m, x[1:nv,1:nv], Bin)
@variable(m, y[1:nv], Bin)

@objective(m, Min, sum(y))

for (u, v) in edges
  for j = 1:nv
    @constraint(m, x[u,j]+x[v,j] <= 1)
  end
end

for i = 1:nv
  @constraint(m, sum(x[i,:]) == 1)
end

for j = 1:nv
  for i = 1:nv
    @constraint(m, y[j] >= x[i,j])
  end
end

println(m)

status = solve(m)

# ----------------------------------------
# Relatório
println("===============================")
if status == :Optimal
  println("Solucão ótima encontrada.")
elseif status == :Unbounded
  println("Problema é ilimitado.")
elseif status == :Infeasible
  println("Problema é inviável.")
elseif status == :UserLimit
  println("Parado por limite de tempo ou iteracões.")
elseif status == :Error
  println("Erro do resolvedor.")
else
  println("Não resolvido.")
end

println("Número de nós explorados: ", getnodecount(m::Model))
D = getobjbound(m::Model)
P = getobjectivevalue(m::Model)
@printf("Melhor limitante dual: %.2f\n", D)
@printf("Melhor limitante primal: %.2f\n", P)
Gap = (abs( D - P )/P)*100
@printf("Gap de otimalidade: %.2f\n", Gap)
@printf("Tempo de execucão: %.2f\n", getsolvetime(m::Model))

# ----------------------------------------
