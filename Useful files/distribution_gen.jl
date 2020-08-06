function distribution_gen(
    n_comparts,
    qt,
    shunt,
    mean_qt_1,
    mean_qt_2,
    log_sd_1,
    log_sd_2,
    qt_ratio,
    dist_min,
    dist_max,
)


    qt_rest = qt - (qt * shunt)

    D = log(dist_max / dist_min) / (n_comparts - 1)
    distribution = []

    for j in range(2, stop = (n_comparts + 1))
        vq = dist_min * exp(D * (j - 2))
        push!(distribution, vq)
    end



    dist_1 =
        exp.(
            -0.5 .* (log(mean_qt_1) .- log.(distribution)) .*
            (log(mean_qt_1) .- log.(distribution)) ./ (log_sd_1 * log_sd_1),
        )
    dist_2 =
        qt_ratio .*
        exp.(
            -0.5 .* (log(mean_qt_2) .- log.(distribution)) .*
            (log(mean_qt_2) .- log.(distribution)) ./ (log_sd_2 * log_sd_2),
        )
    total_dist = dist_1 + dist_2
    volume_raggiunto_qt = (qt_rest .* total_dist) ./ sum(total_dist)
    #volume_raggiunto_va = volume_raggiunto_qt .* distribution

    return distribution, volume_raggiunto_qt

end
