using GLMakie
using Statistics


function main() 

    N = 128
    sweeps = 100000
    therm_its = sweeps * N^2
    J = 1

    window_time = N^2

    fig = Figure()
    display(fig)

    Axis(fig[1,1], aspect=1)
    magnet_axis = Axis(fig[1,2], aspect=1, xlabel="sweeps", ylabel="magnetization")
    chi_axis = Axis(fig[1,3], aspect=1, xlabel="sweeps", ylabel="susceptibility")
    
    kBT = Slider(fig[2, 1], range=0:0.01:120, startvalue=5)
    
    sigma = ones(N, N)

    B = 0
    H = -Float64(2 * N^2 * J)
    M = Float64(N^2)

    function flip(i, j)
        M -= 2 * sigma[i, j]
        H += 2 * sigma[mod1(i+1,N),j] * sigma[i, j] * J 
        H += 2 * sigma[mod1(i-1,N),j] * sigma[i, j] * J
        H += 2 * sigma[i, mod1(j+1, N)] * sigma[i, j] * J
        H += 2 * sigma[i, mod1(j-1, N)] * sigma[i, j] * J
        sigma[i,j] *= -1
    end

    function step()
        i = rand(1:N)
        j = rand(1:N)

        H0 = H
        flip(i, j)
        deltaH = H - H0

        if deltaH > 0 && rand() > exp(-deltaH/kBT.value[])
            flip(i,j)
        end
    end


    xs = [0.0]
    Ms = zeros(1)

    stat_xs = [0]
    stat_chis = zeros(1)
    

    Hs = zeros(1)

    Ms[1] = M
    Hs[1] = H

    sigma_obs = Observable(sigma)
    Ms_obs = Observable(Ms)
    xs_obs = Observable(xs)
    stat_xs_obs = Observable(stat_xs)
    stat_chis_obs = Observable(stat_chis)

    GLMakie.heatmap!(fig[1,1], sigma_obs)
    lines!(fig[1,2], xs_obs, Ms_obs)
    lines!(fig[1,3], stat_xs_obs, stat_chis_obs)

    for it in 1:therm_its
        step()

        push!(xs, it / N^2)
        push!(Ms, M)
        push!(Hs, H)

        if it % N^2 == 0
            

            if it > N^2
                chi = Statistics.var(Ms[(it-N^2):it])
                
                push!(stat_chis, chi)
                push!(stat_xs, it / N^2)

                stat_chis_obs[] = stat_chis
            end

            sigma_obs[] = sigma
            Ms_obs.val = Ms
            xs_obs[] = xs
           
            autolimits!(magnet_axis)
            autolimits!(chi_axis)
            sleep(0)
        end
    end
end