using JuMP, Gurobi, Printf
using LinearAlgebra: dot

TL = parse(Int, ARGS[2]) # in secs
nu = 0
nv = 0
ne = 0
K = 0
nr = 0

# empty dictionary that will eventually have costs, keys will be used
# as indices for model variables
costs = Dict{Tuple{Int, Int}, Float64}()

open(ARGS[1]) do file
  global nv, nu, ne, K, nr

  sizes = split(readline(file))
  nv = parse(Int, sizes[1])
  nu = parse(Int, sizes[2])
  ne = parse(Int, sizes[3])
  K = parse(Int, sizes[4])
  nr = Int((K >= 2) ? (ceil(nv/2)) : nv) # reduce number of rooms

  for line in eachline(file)
    a, b, c = split(line)
    v = parse(Int, a)
    u = parse(Int, b)
    cost = parse(Float64, c)
    costs[v, u] = cost  
  end
end

MEM = 4 # in GB
mod = Model(solver = GurobiSolver(NodefileStart=MEM, TimeLimit=TL))

@variable(mod, t[keys(costs)], Bin)
@variable(mod, m[1:nr,1:nv], Bin)
@variable(mod, p[1:nr,1:nu], Bin)

@objective(mod, Min, sum(t[k]*costs[k] for k in keys(costs)))

for (v, u) in keys(costs)
  for i = 1:nr
    @constraint(mod, t[(v, u)] >= p[i, u] - m[i, v])
  end
end

for i = 1:nr
  @constraint(mod, sum(m[i, :]) <= K)
end

for v = 1:nv
  @constraint(mod, sum(m[:, v]) == 1)
end

for u = 1:nu
  @constraint(mod, sum(p[:, u]) == 1)
end

status = solve(mod)

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

println("Número de nós explorados: ", getnodecount(mod::Model))
D = getobjbound(mod::Model)
P = getobjectivevalue(mod::Model)
@printf("Melhor limitante dual: %.2f\n", D)
@printf("Melhor limitante primal: %.2f\n", P)
Gap = (abs( D - P )/P)*100
@printf("Gap de otimalidade: %.2f\n", Gap)
@printf("Tempo de execucão: %.2f\n", getsolvetime(mod::Model))

# ----------------------------------------
