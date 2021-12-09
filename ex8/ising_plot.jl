using GLMakie
using Statistics


function ising(N, beta, t_sweeps, n_samples, J=1)
    sigma = ones(Int8, N, N)

    H = -Float64(2 * N^2 * J)
    M = Float64(N^2)

    function step()
        i = rand(1:N)
        j = rand(1:N)

        deltaH = 2 * (sigma[mod1(i+1,N),j] + sigma[mod1(i-1,N),j] + sigma[i, mod1(j+1, N)] + sigma[i, mod1(j-1, N)]) * sigma[i,j] * J

        if deltaH < 0 || rand() < exp(-deltaH*beta)
            H += deltaH
            M -= 2 * sigma[i,j]
            sigma[i,j] *= -1
        end
    end

    Ms = Vector{Float64}()
    Hs = Vector{Float64}()

    for _ in 1:t_sweeps*N^2
        step()
    end

    for i in 1:n_samples
        for _ in 1:(3*N^2)
            step()
        end

        push!(Ms, M)
        push!(Hs, H)
    end

    return (Ms/N^2, Hs/N^2)
end

function main() 

    N = 64
    
    beta = LinRange(0,1,64)
    magnetization = zeros(size(beta))
    susceptibility = zeros(size(beta))
    energy = zeros(size(beta))

    p = Threads.Atomic{Int64}(0)

    @time begin
        Threads.@threads for i in 1:size(beta)[1]
            Ms, Hs = ising(N, beta[i], 100, 1000, 1)
            magnetization[i] = Statistics.mean(abs.(Ms))
            susceptibility[i] = beta[i] * Statistics.var(abs.(Ms))
            energy[i] = Statistics.mean(Hs)
            
            Threads.atomic_add!(p, 1)
            println("Progress: $(p[]/size(beta)[1] * 100)%")
        end
    end

    fig = Figure()
    Axis(fig[1,1], xlabel="1/T", title="M")
    Axis(fig[1,2], xlabel="1/T", title="χ")
    Axis(fig[2,1], xlabel="1/T", title="H")
    lines!(fig[1,1], beta, magnetization)
    lines!(fig[1,2], beta, susceptibility)
    lines!(fig[2,1], beta, energy)

    display(fig)
end