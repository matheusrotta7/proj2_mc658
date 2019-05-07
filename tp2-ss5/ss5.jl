using JuMP, Gurobi, Printf
using LinearAlgebra: dot

TL = parse(Int, ARGS[2]) # in secs
n = 0
d = 0
t = Array{Float64}(undef, 0)
d = copy(t)
p = copy(t)

params = [t,d,p]

open(ARGS[1]) do file
  global n
  n = parse(Int, readline(file))

  for prm in params
    resize!(prm, n)
  end

  for (i, line) in enumerate(eachline(file))
    for (j, sn) in enumerate(split(line))
        params[j][i] = parse(Int, sn)
    end
  end
end

M = sum(t)

MEM = 4 # in GB
m = Model(solver = GurobiSolver(NodefileStart=MEM, TimeLimit=TL))

@variable(m, a[1:n] >= 0)
@variable(m, s[1:n] >= 0)
@variable(m, y[1:n, 1:n], Bin)

@objective(m, Min, dot(p,a))

for i=1:n,j=1:n
  @constraint(m, s[i] + t[i] <= s[j] + M*(1-y[i,j]))
  @constraint(m, s[i] + t[i] >= s[j] - M*(1-y[i,j]))
end

for i=1:n,j=1:n
  @constraint(m, y[i,j] + y[j,i] <= 1)
end

for i=1:n
  @constraint(m, sum(y[i,:]) <= 1)
  @constraint(m, sum(y[:,i]) <= 1)
end

@constraint(m, sum(y) == n-1)

for j = 1:n
  @constraint(m, a[j] >= s[j]+t[j]-d[j])
end

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
