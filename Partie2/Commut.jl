using Plots

α1 = 1
a1 = 1
a2 = 1
α2 = 1

O(x1, x2) = 1 / (1 + x1 + x2 )

function commut(m1,m2,p1,p2)
    return (1 + 1/2*(p1*α1 + p2*α2))*O(m1,m2) + 0.5*(p1*a1*m1 + p2*a2*m2)
end

m_values1 = (0.0:0.1:5.0)
m_values2 = (0.0:0.1:10.0)
p_values1 = (-3.0:0.5:0.0)
p_values2 = (-3.0:0.5:0.0)

plotlyjs()
valeurs_m1, valeurs_m2, valeurs_p1, valeurs_p2 = [],[],[],[]
p2 = -0.1
p1 = -0.1
m1 = -0.0
m2 = 0.0 

for m1 in m_values1
    push!(valeurs_m1, commut(m1,m2,p1,p2))
end

for m2 in m_values2
    push!(valeurs_m2, commut(m1,m2,p1,p2))
 end

 for p1 in p_values1
    push!(valeurs_p1, commut(m1,m2,p1,p2))
 end

 for p2 in p_values2
    push!(valeurs_p2, commut(m1,m2,p1,p2))
 end

 z = [commut(m1, m2, p1, p2) for m1 in m_values1, p1 in p_values1]
 heatmap(m_values1, p_values1, z, xlabel="m1", ylabel="p1", title="Impact of m1 and p1 on commut", colorbar_title="commut")

plt1 = plot()
plt2 = plot()

plot(plt1,m_values1,valeurs_m1)
plot!(m_values2,valeurs_m2)
# so both m1 and m2 have the same behaviour in the commutation function.
plot(plt2,p_values1, valeurs_p1)
plot!(p_values2, valeurs_p2)