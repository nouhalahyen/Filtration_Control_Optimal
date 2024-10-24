
#import Pkg
#Pkg.activate(".")
#Pkg.update();
#Pkg.add("NLPModelsIpopt")
using NLPModelsIpopt
using Plots
using JuMP
using Ipopt
using DifferentialEquations
using QuadGK  # For numerical integration

t0 = 0     
tf = 40.0

α1 = 1
α2 = 1.5

a1 = 2.0
a2 = 1.0

u_bar = 0.66
m0 = [0,0]

function traiter_controle_masse(m0)
    return m0[1] > 1 || m0[2] >1       
end


# Define the piecewise control function u(t)
function u_control(t, t1, t2,u_bar)
    a = traiter_controle_masse(m0)
    if t <= t1 && !a
        return 1.0
    elseif t <= t1 && a
        return -1
    elseif t1 < t && t <= t2
        return u_bar
    else
        return 1.0
    end
end
 # Define the differential equations
function my_dynamics!(dm, m, p, t)
    t1, t2 = p
    u_t = u_control(t, t1, t2,u_bar)
    dm[1] = α1 / (1 + m[1] + m[2]) * ((1 + u_t) / 2) - (1 - u_t) / 2 * a1 * m[1]  # \dot{m_1}
    dm[2] = α2 / (1 + m[1] + m[2]) * ((1 + u_t) / 2) - (1 - u_t) / 2 * a2 * m[2]  # \dot{m_2}
end
# Numerical integration function to integrate the objective
function integrate_objective(m0,t1, t2)
    #m0 = [0.0, 0.0]  # Initial conditions for m1, m2
    tspan = (0.0, tf)
    p = [t1, t2]  # Parameters to pass t1 and t2

    prob = ODEProblem(my_dynamics!, m0, tspan, p)
    sol = solve(prob, AutoVern7(Rodas5()), saveat=0.01)  # Solve the ODE

    # Define the objective integrand
    function objective_integrand(t)
        u_t = u_control(t, t1, t2,u_bar)
        return u_t / (1 + sol(t)[1] + sol(t)[2])
    end

    # Numerical integration of the objective
    integral_value, _ = QuadGK.quadgk(objective_integrand, 0.0, tf)

    return integral_value
end
function optimum(m0, u_bar)

    # Numerical integration function to integrate the objective
    function integrate_objective2(t1, t2)
        #m0 = [0.0, 0.0]  # Initial conditions for m1, m2
        tspan = (0.0, tf)
        p = [t1, t2]  # Parameters to pass t1 and t2

        prob = ODEProblem(my_dynamics!, m0, tspan, p)
        sol = solve(prob, Tsit5(), saveat=0.01)  # Solve the ODE

        # Define the objective integrand
        function objective_integrand(t)
            u_t = u_control(t, t1, t2,u_bar)
            return u_t / (1 + sol(t)[1] + sol(t)[2])
        end

        # Numerical integration of the objective
        integral_value, _ = QuadGK.quadgk(objective_integrand, 0.0, tf)

        return integral_value
    end

    # Optimization using JuMP
    model = JuMP.Model(Ipopt.Optimizer)

    # Define variables
    @variable(model, 0.0 <= t1 <= tf, start = 2.33 )
    @variable(model, 0.0 <= t2 <= tf, start = 36.69)

    #Define constraints
    @constraint(model, t2 >= t1)
    # Define the objective function that directly uses the optimization variables
    @NLobjective(model, Max, integrate_objective2(t1, t2))

    # Solve the problem
    optimize!(model)

    t1_val = value(t1)
    t2_val = value(t2)
    objective_val = objective_value(model)
    results = [t1_val,t2_val,objective_val]

    return results
    
end
optimum(m0, u_bar)
