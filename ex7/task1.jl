include("../utils.jl")

function mc_integrate(n, f, a, b)  
    s = 0
    
    for i in 1:n
        s += f(rand(a:0.0001:b))
    end

    return (b-a) * s / n
end

function mc_integrate_weighted(N, f, a, b, g, G_inv, G)
    s = 0

    for i in 1:N
        x = 0

        while true
            x = G_inv(rand(a:0.0001:b))

            if x < 1
                break
            end
        end

        s += f(x) / g(x)
    end

    return s / N
end

function trap_integrate(n, f, a, b)

    h = (b-a) / n
    
    s = 0
    
    for i in 1:n-1
        s += f(a+i*h)
    end

    return 0.5 * h * f(a) + 0.5 * h * f(b) + h * s
end

function main()  
    Plots.plot()
    f(x) = exp(-x^2)
    I_ex = 0.7468241328124270253994674361318530053544996868126063290276544989
    
    a = 0
    b = 1

    plot_convergence(N=20, sample_size=10, label="uniform mc") do N
        return abs(mc_integrate(N, f, a, b) - I_ex)
    end

    g(x) = exp(-x) / (exp(0) - exp(-1))
    G(x) = -exp(-x)
    G_inv(x) = -log(1 - x)

    plot_convergence(N=20, sample_size=10, label="importance mc") do N
        return abs(mc_integrate_weighted(N, f, a, b, g, G_inv, G) - I_ex)
    end

    plot_convergence(N=10, sample_size=1, label="trapezoidal rule") do N
        return abs(trap_integrate(N, f, a, b) - I_ex)
    end
end