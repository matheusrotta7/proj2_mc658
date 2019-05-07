using JuMP, Gurobi, Printf
using LinearAlgebra: dot

TL = parse(Int, ARGS[2]); # in secs
n = 0;
card_S = 0;
S = Array{Tuple{Int64,Int64},1}(undef, 0); #S is an array of tuples
t = Array{Int64,1}(undef, 0); #initializing arrays with size 0 - good practice?
d = Array{Int64,1}(undef, 0);
params = [t, d];
M = 999999

#-------------grab input------------#
open(ARGS[1]) do file
    global n;
    global card_S; #if not global, changes would remain local
    first_line = split(readline(file));
    n = parse(Int, first_line[1]);
    # println(n);
    card_S = parse(Int, first_line[2]);
    # println(card_S);
    resize!(S, card_S);

    for prm in params
        resize!(prm, n); #resize t and d from 0 to n
    end

    for i = 1:n
        cur_line = readline(file);
        for (j, sn) in enumerate(split(cur_line))
            params[j][i] = parse(Int, sn);
        end
    end
    # println(card_S)

    for i = 1:card_S
        # println(i)

        cur_line = split(readline(file));
        # println(cur_line)
        S[i] = (parse(Int, cur_line[1]), parse(Int, cur_line[2]));
    end
end
#-------------grab input(end)------------#

#-------------check input------------#
println("value of n:")
println(n)
println("value of card_S:")
println(card_S)
println("S precedence tuples:")
println(S)
println("values of t:")
println(t)
println("values of d:")
println(d)
#-------------check input(end)------------#


MEM = 4 # in GB
DEADLINE = Model(solver = GurobiSolver(NodefileStart=MEM, TimeLimit=TL)) #model name is DEADLINE

#    Constrói conjunto de pares ordenados de tarefas
m=n*(n-1)  # número de pares ordenados (exclui pares do tipo (i,i))
Pares = Array{Tuple{Int64, Int64}}(undef,m)
k=1
for i in 1:n
    for j in 1:n
        global k
        if (i != j)
            Pares[k] = (i,j)
            k=k+1
        end
    end
end

@variable(DEADLINE, s[1:n] >= 0) #beginning time for task i
@variable(DEADLINE, y[1:n], Bin) #tells whether task is late or not
@variable(DEADLINE, p[e in Pares], Bin)

@objective(DEADLINE, Min, sum(y))

#RESTRIÇÃO 1:
for i in 1:(n-1) #prof. probably wrote like this so as to not have redundancy
   for j in (i+1):n
      f=(i,j)
      g=(j,i)
      @constraint(DEADLINE, p[f]+p[g] == 1)
   end
end

#RESTRIÇÃO 2:
for e in Pares
   i=e[1]
   j=e[2]
   reverse = (e[2], e[1])
    @constraint(DEADLINE, s[j] >= s[i] + t[i] - M*p[reverse])
end

#RESTRIÇÃO 3:
for e in S
   i=e[1]
   j=e[2]
   # reverse = (e[2], e[1])
    @constraint(DEADLINE, p[e] == 1)
end

#RESTRIÇÃO 4:
for j = 1:n
    @constraint(DEADLINE, s[j] + t[j] <= d[j] + M*y[j])
end

status = solve(DEADLINE)

# --------------------------------------------------------------------

# relatório estatístico
println("========================================================================")
if status == :Optimal
  println("Solução ótima encontrada.")
elseif status == :Unbounded
  println("Problema é ilimitado.")
elseif status == :Infeasible
  println("Problema é inviável.")
elseif status == :UserLimit
  println("Parado por limite de tempo ou iterações.")
elseif status == :Error
  println("Erro do resolvedor.")
else
  println("Não resolvido.")
end

println("Número de nós explorados: ", getnodecount(DEADLINE::Model))
D = getobjbound(DEADLINE::Model)
P = getobjectivevalue(DEADLINE::Model)
@printf("Melhor limitante dual: %.2f\n", D)
@printf("Melhor limitante primal: %.2f\n", P)
Gap = (abs( D - P )/P)*100
@printf("Gap de otimalidade: %.2f\n", Gap)
@printf("Tempo de execução: %.2f\n", getsolvetime(DEADLINE::Model))

# end  # => fim do bloco "let"
