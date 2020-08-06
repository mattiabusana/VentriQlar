
function oxygen_content(hb, sat_100, po2)
    content = 1.39 * hb * (sat_100 / 100) + 0.003 * po2
    return content
end

function severinghaus(po2)
    sat = ((23400 * (po2^3 + 150 * po2)^(-1)) + 1)^(-1)
    return sat
end

function zander_saturation(po2, hb, pco2, ph)
    hco3 = 0.03 * pco2 * 10^(ph - 6.1)
    x = (1 - 0.014 * hb) * (hco3 - 24 + (1.63 * hb + 9.5) * (ph - 7.4))
    c = (po2 / 26.7)^0.184
    po2pp = po2 * (exp((c + 0.003 * x - 2.2) * (7.4 - ph)))
    sat = 100 / (1 + (23400 / (po2pp^3 + 150 * po2pp)))

    return sat
end

function zander_be_blood(ph, pco2, hb, sat_100)
    be =
        (1 - 0.014 * hb) * (
            (0.0304 * pco2 * 10^(ph - 6.1) - 24.26) +
            (9.5 + 1.63 * hb) * (ph - 7.4)
        ) - (0.2 * hb * ((100 - sat_100) / 100))

    return be
end


function douglas_to_content(ph, pco2, hb, sat_100)
    s_frac = sat_100 / 100
    cont_co2 =
        0.03 *
        pco2 *
        (1 + 10^(ph - 6.1)) *
        (1 - ((0.0289 * hb) / ((3.352 - 0.456 * s_frac) * (8.142 - ph))))

    return cont_co2 * 2.24
end


function douglas_to_pressure(ph, content_ml, hb, sat_100)
    s_frac = sat_100 / 100
    cont_co2 = content_ml / 2.24
    pco2 =
        cont_co2 / (
            0.03 *
            (1 + 10^(ph - 6.1)) *
            (1 - ((0.0289 * hb) / ((3.352 - 0.456 * s_frac) * (8.142 - ph))))
        )

    return pco2
end

function ph_from_be_zander(be_set, pco2, hb, sat_100)
    be_iter = -10
    ph = 6
    while be_iter <= be_set
        be_iter =
            (1 - 0.014 * hb) * (
                (0.0304 * pco2 * 10^(ph - 6.1) - 24.26) +
                (9.5 + 1.63 * hb) * (ph - 7.4)
            ) - (0.2 * hb * ((100 - sat_100) / 100))
        ph += 0.0001
    end

    return ph
end


function log_distribution(min, max, n_comparts)

    dist_space = log(max / min) / n_comparts
    distribution = []

    for j in range(1, stop = n_comparts)
        element = min * exp(dist_space * j)
        push!(distribution, element)
    end

    return distribution

end
