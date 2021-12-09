using LinearAlgebra

include("../utils.jl")

function main()

    plot_convergence(xaxis=:lin, yaxis=:lin, N=20) do M
        L = 1
        N = 10
        R = 0.1

        total = Threads.Atomic{Float64}(0.0)

        Threads.@threads for i in 1:M
            xs = zeros(3, N)
            for i in 1:N
                while true
                    xs[:,i] .= rand(3) .* L

                    overlap = false

                    for j in 1:i-1
                        if norm(xs[:, i] - xs[:, j]) <= 2R
                            overlap = true
                            break
                        end
                    end

                    if !overlap
                        break
                    end
                end
            end

            d = 0
            for i in 1:N, j in 1:(i - 1)
                d += norm(xs[:,i]-xs[:,j])
            end

            avg = 2 / (N*(N-1)) * d
            Threads.atomic_add!(total, avg)
        end
        

        return total[] / M
    end

end
