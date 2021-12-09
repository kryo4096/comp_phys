using Plots

function plot_convergence(error; expected_convergence=-1, N=25, base=2, sample_size=10, label="", xaxis=:log, yaxis=:log)
    errors = zeros(N)
    ns = zeros(N)

    
    for i in 1:N
        ns[i] = base^i
        print("N=$(ns[i])")
        @time errors[i] = sum(error.(ones(sample_size) * ns[i]))
        println("\terr=$(errors[i])")
    end

    plt = Plots.plot!(ns, errors, xaxis=xaxis, yaxis=yaxis, label=label)
    if(expected_convergence > 0)
        Plots.plot!(ns, 1 ./ns.^expected_convergence * errors[1] * ns[1]^expected_convergence, xaxis=:log, yaxis=:log, label=label * " expected")
    end

    display(plt)
end