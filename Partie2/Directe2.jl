
using OptimalControl
using Plots
t0 = 0     
tf = 40.0
σ1 = 1
σ2 = 1

α1 = 1
α2 = 5

a1 = 1
a2 = 2

# Définition des fonctions dynamiques
O(x1) = 1 / (1 + x1)
#̃m1:
F11(x1) = (α1/σ1 + α2/σ2) * O(x1)
F12(x1,x2) = ((a1^2+a2)*x1 - (a1+a2)*x2)/(1-a1/a2)
#̃m2:
F21(x1) = (α1/σ1 + a1*α2/(σ2*a2)) * O(x1)
F22(x1) =  a1*x1

x0 = [0.0, 0.0]
# Définir le problème de contrôle optimal
Problem(x0) = begin 
    @def ocp begin
        t ∈ [t0, tf], time
        x ∈ R², state
        u ∈ R, control
        m1 = x₁
        m2 = x₂

        x(t0) == x0
        -m1(t) ≤ 0.0
        -m2(t) ≤ 0.0
        -1.0 ≤ u(t) ≤ 1.0
        
        ẋ(t) == [(1 + u(t)) / 2 * F11(m1(t)) - (1 - u(t)) / 2 * F12(m1(t),m2(t)),
                 (1 + u(t)) / 2 * F21(m1(t)) - (1 - u(t)) / 2 * F22(m1(t))]
        
        ∫(u(t) * O(m1(t))) → max
    end
end

# Résoudre le problème de contrôle optimal
prob = Problem(x0)
sol = solve(prob,grid_size=1500)
plotlyjs()
plot(sol, title="σ1 = $σ1, σ2 = $σ2, α1=$α1, α2=$α2, a1 = $a1, a2 = $a2", size=(950, 600))  # Décommenter pour afficher le graphique de la solution