import Plots

function mc(n)
    count = 0

    for i in 1:n
        if sqrt((rand()-0.5)^2 + (rand()-0.5)^2) < 0.5
            count += 1
        end
    end

    return count / n
end

function main() 

    N = 25
    errors = zeros(N)
    ns = zeros(N)

    for i in 1:N
        ns[i] = 2^(i)
        Threads.@threads for j in 1:10
            errors[i] += abs(mc(ns[i]) - π/4)
        end
        errors[i] /= 10
    end

    display(errors)

    Plots.plot(ns, errors, xaxis=:log, yaxis=:log)
    Plots.plot!(ns, 1 ./ns.^0.5 * errors[1], xaxis=:log, yaxis=:log)

end