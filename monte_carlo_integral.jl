include("utils.jl")

function mc_integrate(n, f, a, b)
    return (b-a) * sum(map((i) -> f(rand(a:0.0001:b)) , 1:n))/n
end

function main() 
    f(x) = sin(x)
    
    a = -π
    b = π

    error(N) = abs(mc_integrate(N, f, a, b))

    plot_convergence(error, N=20, sample_size=10, expected_convergence=0.5)
end