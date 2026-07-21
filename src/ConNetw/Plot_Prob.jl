"""
    total_prob(x, distances)

Function to calculate and store all stages of stacked probability. Function is meant for demonstrative purposes.

# Arguments
- `x`: Point where to evaluate function.
- `distances`: Vector of parameter which adjust parameters of RayleighSISO prob. function.
"""
function total_prob(x, distances)
    Prob = RayleighSISO(x, distances[1]/3, eta=2)
    Prob_vec = [Prob]
    for i in 2:length(distances)
        Prob += (1-Prob)*RayleighSISO(x, distances[i], eta=6)
        push!(Prob_vec, Prob)
    end
    return Prob_vec
end

"""
    plot_prob()

Function to visualize stacked probability.
"""
function plot_prob()
    max_distance = 1500
    distances = [max_distance, 300, 800]

    xx = 0:0.1:max_distance
    Werte =zeros(3, length(xx))
    for (i,x) in enumerate(xx) 
        Werte[:,i] = total_prob(x, distances)
    end
    display(Werte)
    p1 = plot(xx, Werte[1,:], label="Dist=$(distances[1])", legend=:outertopright)
    for i in 2:length(distances)
        plot!(p1, xx, Werte[i,:], label="Dist=$(distances[i])")
    end
    display(p1)
end