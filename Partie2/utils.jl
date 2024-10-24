# ---------------------------------------------------------------------
function plot_traj!(plt, F, z0, tf)
    """ 
    Plot the trajectories in the phase plan of the ode ż(t)=F(z(t)),
    for the initial points given in the z0, until the time tf.
    -------
    
    parameters (input)
    ------------------
    plt : plot object
    F  : function F(z)
    z0 : vector of vector of initial points
         z0[i] = ith initial point
    tf : final time
    returns
    -------
    nothing
    """
    tspan = (0 , tf)
    for i ∈ eachindex(z0)
        prob = ODEProblem((z, _, _) -> F(z), z0[i], tspan)   # défini le problème en Julia
        sol = solve(prob)
        plot!(plt, sol, vars = (1,2), legend = false, color = :blue)
    end
    return nothing
end

# ---------------------------------------------------------------------
function plot_flow!(plt, F, z0, tf; nb_fronts=4)
    """ 
    For times t ∈ range(0, tf, nb_fronts), plot the solution at time t of the Cauchy 
    problems ż(s) = F(z(s)), z(0) = z0[i].
    -------
    
    parameters (input)
    ------------------
    plt : plot object
    F  : function F(z)
    z0 : vector of vector of initial points
         z0[i] = ith initial point
    tf : final time
    returns
    -------
    nothing
    """
    tspan = (0, tf)
    m = nb_fronts
    T = range(0, tf, length=m)
    print("Times = ", [round(t, digits=2) for t ∈ T], "\n")
    N = length(z0)
    V = [ [] for i ∈ 1:m ]
    for i ∈ 1:N
        prob = ODEProblem((z, _, _) -> F(z), z0[i], tspan)   # défini le problème en Julia
        sol = solve(prob)
        for i ∈ 1:m
            t = T[i]
            push!(V[i], sol(t))
        end
    end
    for i ∈ 1:m
        plot!(plt, [V[i][j][1] for j ∈ 1:N], [V[i][j][2] for j ∈ 1:N], legend = false, 
            color = :green, linewidth = 2)
    end
    return nothing
end

# -------------------------------------------------------------------------------
# Fonction qui construit le flot à partir d'un système hamiltonien
# dans le but de résoudre un système hamiltonien 
#
# (dx/dt, dp/dt) = Hv(x, p), (x(t0), p(t0)) = (x0, p0)
#
# où Hv = (dH/dp, -dH/dx).
#
function Flow(Hv)
    """

    parameters (input): Hv

        function Hv(x, p)
            dx = dH/dp
            dp = -dH/dx
            return dx, dp
        end

    returns (output): f

    usage:

        julia> f = Flow(Hv)
        julia> xf, pf = f(t0, x0, p0, tf)

    to get the endpoints. Or

        julia> sol = f((t0, tf), x0, p0)

    to get the solution object, to plot the trajectory for example.
        
        julia> plot(sol, vars=(1,2)

    """
    # 
    function rhs!(dz, z, dummy, t)
        n = size(z, 1)÷2
        dz[1:n], dz[n+1:2n] = Hv(z[1:n], z[n+1:2n])
    end
    
    # cas vectoriel
    function f(tspan::Tuple{Real, Real}, x0::Vector{<:Real}, p0::Vector{<:Real}; abstol=1e-12, reltol=1e-12, saveat=0.01)
        z0 = [ x0 ; p0 ]
        ode = ODEProblem(rhs!, z0, tspan)
        sol = solve(ode, Tsit5(), abstol=abstol, reltol=reltol, saveat=saveat)
        return sol
    end
    
    function f(t0::Real, x0::Vector{<:Real}, p0::Vector{<:Real}, tf::Real; abstol=1e-12, reltol=1e-12, saveat=[])
        sol = f((t0, tf), x0, p0, abstol=abstol, reltol=reltol, saveat=saveat)
        n = size(x0, 1)
        return sol(tf)[1:n], sol(tf)[n+1:2n]
    end
    
    # cas scalaire
    function rhs_scalar!(dz, z, dummy, t)
        dz[1], dz[2] = Hv(z[1], z[2])
    end

    function f(tspan::Tuple{Real, Real}, x0::Real, p0::Real; abstol=1e-12, reltol=1e-12, saveat=0.01)
        z0 = [ x0 ; p0 ]
        ode = ODEProblem(rhs_scalar!, z0, tspan)
        sol = solve(ode, Tsit5(), abstol=abstol, reltol=reltol, saveat=saveat)
        return sol
    end

    function f(t0::Real, x0::Real, p0::Real, tf::Real; abstol=1e-12, reltol=1e-12, saveat=[])
        sol = f((t0, tf), x0, p0, abstol=abstol, reltol=reltol, saveat=saveat)
        return sol(tf)[1], sol(tf)[2]
    end

    return f

end;