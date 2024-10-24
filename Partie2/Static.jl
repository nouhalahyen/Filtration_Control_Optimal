#import Pkg
#Pkg.add("JuMP")
#Pkg.add("Ipopt")
using JuMP
using Ipopt

alpha1 = 1
alpha2 = 1.5
a2 = 1
a1 = 2

# Create a model with the GLPK solver
model = Model(Ipopt.Optimizer)

# Define variables
@variable(model, m1 >= 0)
@variable(model, m2 >= 0)
@variable(model, -1.0 <= u <= 1.0)



# Define the objective function
@objective(model, Max, u/(1+m1 + m2))

# Define the constraints
@constraint(model, alpha1/(1+m1 + m2)*((1+u)/2) - (1-u)/2*a1*m1 == 0.0)
@constraint(model, alpha2/(1+m1 + m2)*((1+u)/2) - (1-u)/2*a2*m2 == 0.0)


# Solve the problem
optimize!(model)

# Get the results
m1_val = value(m1)
m2_val = value(m2)
u_val = value(u)

objective_val = objective_value(model)

println("Optimal solution: m1 = $m1_val, m2 = $m2_val, u = $u_val")
println("Objective value: $objective_val")
