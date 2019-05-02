using JuMP, Gurobi, Printf
using LinearAlgebra: dot

TL = parse(Int, ARGS[2]) # in secs
n = 0
d = 0
d = Array{Float64}(undef, 0)
c = copy(d)
b = copy(d)
p = copy(d)
h = copy(d)

params = [d,c,b,p,h]

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

MEM = 4 # in GB
m = Model(solver = GurobiSolver(NodefileStart=MEM, TimeLimit=TL))

@variable(m, e[1:(n+1)] >= 0)
@variable(m, q[1:n] >= 0)
@variable(m, x[1:n], Bin)

@objective(m, Min, dot(b,x)+dot(h,e[1:n])+dot(p,q))

@constraint(m, e[1] == 0)
@constraint(m, e[n+1] == 0)

for i = 1:n
  @constraint(m, q[i] <= c[i]*x[i])
end
for i = 1:n
  @constraint(m, e[i+1]-e[i] == q[i]-d[i])
end

println(m)

status = solve(m)

println(d)
println(getvalue(e))
println(getvalue(q))

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
