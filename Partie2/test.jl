import Pkg; 
Pkg.add("LSODA")
using LSODA
using JuMP
using Ipopt
using DifferentialEquations
using QuadGK  # For numerical integration

# Define constants
alpha1 = 1.0
alpha2 = 1.5
a1 = 2.0
a2 = 1.0
u_bar = 0.6636975627412119
tf = 40.0

# Define the piecewise control function u(t)
function u_control(t, t1, t2)
    if t <= t1 || t >= t2
        return 1.0
    else
        return u_bar
    end
end

# Define the general dynamics function (common structure)
function my_dynamics!(du, u, p, t, u_t)
    du[1] = alpha1 / (1 + u[1] + u[2]) * ((1 + u_t) / 2) - (1 - u_t) / 2 * a1 * u[1]
    du[2] = alpha2 / (1 + u[1] + u[2]) * ((1 + u_t) / 2) - (1 - u_t) / 2 * a2 * u[2]
end

# Wrapper function for different phases based on control input
function dynamics_1!(du, u, p, t)
    my_dynamics!(du, u, p, t, 1.0)  # u_t is fixed at 1 in this phase
end

function dynamics_2!(du, u, p, t)
    my_dynamics!(du, u, p, t, u_bar)  # u_t is u_bar in this phase
end

# Numerical integration function to integrate the objective
function integrate_objective(t1, t2)
    u0 = [0.0, 0.0]  # Initial conditions for m1, m2
    p = []

    # Define time spans for each phase
    tspan1 = (0.0, t1)
    tspan2 = (t1, t2)
    tspan3 = (t2, tf)

    # Solve ODE for each phase
    prob1 = ODEProblem(dynamics_1!, u0, tspan1, p)
    sol_1 = solve(prob1, AutoVern7(Rodas5()), saveat=0.1)

    prob2 = ODEProblem(dynamics_2!, u0, tspan2, p)
    sol_2 = solve(prob2, AutoVern7(Rodas5()), saveat=0.1)

    prob3 = ODEProblem(dynamics_1!, u0, tspan3, p)
    sol_3 = solve(prob3, AutoVern7(Rodas5()), saveat=0.1)

    # Define the objective integrand (combine into one function)
    function objective_integrand(t, sol)
        u_t = u_control(t, t1, t2)
        return u_t / (1 + sol(t)[1] + sol(t)[2])
    end

    # Numerical integration of the objective across different phases
    integral_value1, _ = QuadGK.quadgk(t -> objective_integrand(t, sol_1), 0.0, t1)
    integral_value2, _ = QuadGK.quadgk(t -> objective_integrand(t, sol_2), t1, t2)
    integral_value3, _ = QuadGK.quadgk(t -> objective_integrand(t, sol_3), t2, tf)

    # Total integral
    return integral_value1 + integral_value2 + integral_value3
end

# Optimization using JuMP
model = JuMP.Model(Ipopt.Optimizer)

# Define variables with initial guesses
@variable(model, 0.0 <= t1 <= tf, start = 2.0)
@variable(model, 0.0 <= t2 <= tf, start = 37.0)

# Define constraints
@constraint(model, t2 >= t1)

# Define the objective function for optimization
@NLobjective(model, Max, integrate_objective(t1, t2))

# Solve the optimization problem
optimize!(model)

# Get and print results
t1_val = value(t1)
t2_val = value(t2)
objective_val = objective_value(model)

println("Optimal solution: t1 = $t1_val, t2 = $t2_val")
println("Objective value: $objective_val")

# Evaluate the objective for specific values of t1 and t2
o = integrate_objective(2.83, 36.48)
println("Objective value of the new method: $o")
